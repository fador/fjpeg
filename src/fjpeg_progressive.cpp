/*
FJPEG
BSD 2-Clause License

Copyright (c) 2024, Marko Viitanen

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

#include "fjpeg.h"
#include "fjpeg_bitstream.h"
#include "fjpeg_huffman.h"
#include "fjpeg_transquant.h"

// A progressive scan. "component" selects the coded component(s):
//   -1 = all components interleaved (only valid for DC scans)
//    0 = luma, 1 = Cb, 2 = Cr (single-component scans, required for AC)
typedef struct {
    int component;
    int ss;
    int se;
    int ah;
    int al;
} fjpeg_scan_t;

// Progressive scan script with spectral selection and AC successive approximation (ITU-T T.81 Annex G).
static const fjpeg_scan_t fjpeg_progressive_script[] = {
    {-1,  0,  0, 0, 1}, // 1. DC initial (Ah=0, Al=1)
    {-1,  0,  0, 1, 0}, // 2. DC refinement (Ah=1, Al=0)
    { 0,  1, 63, 0, 1}, // 3. Luma AC initial (Ah=0, Al=1)
    { 0,  1, 63, 1, 0}, // 4. Luma AC refinement (Ah=1, Al=0)
    { 1,  1, 63, 0, 1}, // 5. Cb AC initial (Ah=0, Al=1)
    { 1,  1, 63, 1, 0}, // 6. Cb AC refinement (Ah=1, Al=0)
    { 2,  1, 63, 0, 1}, // 7. Cr AC initial (Ah=0, Al=1)
    { 2,  1, 63, 1, 0}, // 8. Cr AC refinement (Ah=1, Al=0)
};
static const int fjpeg_progressive_script_len =
    (int)(sizeof(fjpeg_progressive_script) / sizeof(fjpeg_progressive_script[0]));

static inline int fjpeg_abs_int(int v) {
    return v < 0 ? -v : v;
}

static inline int fjpeg_bit_length(int v) {
    int a = fjpeg_abs_int(v);
    int s = 0;
    while (a) {
        a >>= 1;
        s++;
    }
    return s;
}

static inline void fjpeg_emit_symbol(fjpeg_bitstream* stream, const fjpeg_huffman_table_t* table, int symbol) {
    stream->writeBits(table[symbol].code, table[symbol].len);
}

// Point transform: divide by 2^al, rounding towards zero.
static inline int fjpeg_point_transform(int v, int al) {
    if (al == 0) return v;
    return v < 0 ? -((-v) >> al) : (v >> al);
}

// Visit every block of one component in raster order.
template <typename F>
static void fjpeg_for_each_block_component(fjpeg_context* context, int channel, F fn) {
    const int w = (channel == 0) ? context->padded_width : context->padded_width / 2;
    const int h = (channel == 0) ? context->padded_height : context->padded_height / 2;
    for (int y = 0; y < h; y += 8) {
        for (int x = 0; x < w; x += 8) {
            fn(x, y, channel);
        }
    }
}

// Visit every block in the frame's interleaved MCU order.
template <typename F>
static void fjpeg_for_each_block_all(fjpeg_context* context, F fn) {
    const int inc_xy = context->channels == 1 ? 8 : 16;
    const int max_uv = context->channels == 1 ? 1 : 2;
    for (int y = 0; y < context->padded_height; y += inc_xy) {
        for (int x = 0; x < context->padded_width; x += inc_xy) {
            for (int v = 0; v < max_uv; v++) {
                for (int u = 0; u < max_uv; u++) {
                    fn(x + u * 8, y + v * 8, 0);
                }
            }
            if (context->channels == 3) {
                fn(x >> 1, y >> 1, 1);
                fn(x >> 1, y >> 1, 2);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Statistics passes
// ---------------------------------------------------------------------------

static void fjpeg_progressive_gather_dc_stats(fjpeg_context* context, const fjpeg_scan_t* scan,
                                              fjpeg_huffman_statistics_t* stat) {
    memset(stat, 0, sizeof(fjpeg_huffman_statistics_t));
    int last_dc[3] = {0, 0, 0};
    fjpeg_coeff_t block[64];
    fjpeg_for_each_block_all(context, [&](int x, int y, int channel) {
        fjpeg_extract_coeff_8x8(context, block, x, y, channel);
        // DC uses an arithmetic right shift so that the decoder's refinement
        // bits (OR-ed in) reconstruct the exact value.
        int val = fjpeg_round_dc(block[0]) >> scan->al;
        int diff = val - last_dc[channel];
        last_dc[channel] = val;
        int size = fjpeg_bit_length(diff);
        if (size > 11) {
            fprintf(stderr, "Error: DC coefficient size overflow\n");
            exit(1);
        }
        if (channel == 0) stat->luma_dc[size]++;
        else stat->chroma_dc[size]++;
    });
}

static void fjpeg_progressive_gather_ac_stats(fjpeg_context* context, const fjpeg_scan_t* scan,
                                              fjpeg_huffman_statistics_t* stat) {
    memset(stat, 0, sizeof(fjpeg_huffman_statistics_t));
    fjpeg_coeff_t block[64];
    uint32_t* ac = (scan->component == 0) ? stat->luma_ac : stat->chroma_ac;
    int eob_run = 0;

    auto flush_eob_run = [&]() {
        if (eob_run == 0) return;
        int r = 0;
        int temp = eob_run >> 1;
        while (temp > 0) {
            r++;
            temp >>= 1;
        }
        ac[r << 4]++;
        eob_run = 0;
    };

    fjpeg_for_each_block_component(context, scan->component, [&](int x, int y, int channel) {
        fjpeg_extract_coeff_8x8(context, block, x, y, channel);

        // Find last non-zero coefficient in scan->ss .. scan->se
        int last_k = scan->ss - 1;
        for (int k = scan->se; k >= scan->ss; k--) {
            if (fjpeg_point_transform(fjpeg_round(block[k]), scan->al) != 0) {
                last_k = k;
                break;
            }
        }

        if (last_k < scan->ss) {
            eob_run++;
            if (eob_run == 32767) {
                flush_eob_run();
            }
            return;
        }

        flush_eob_run();

        int run = 0;
        for (int k = scan->ss; k <= last_k; k++) {
            int shifted = fjpeg_point_transform(fjpeg_round(block[k]), scan->al);
            if (shifted == 0) {
                run++;
                continue;
            }
            while (run > 15) {
                ac[0xF0]++;
                run -= 16;
            }
            int size = fjpeg_bit_length(shifted);
            if (size > 10) {
                fprintf(stderr, "Error: AC coefficient size overflow\n");
                exit(1);
            }
            ac[(run << 4) | size]++;
            run = 0;
        }

        if (last_k < scan->se) {
            eob_run++;
            if (eob_run == 32767) {
                flush_eob_run();
            }
        }
    });

    flush_eob_run();
}

static void fjpeg_progressive_gather_ac_refine_stats(fjpeg_context* context, const fjpeg_scan_t* scan,
                                                     fjpeg_huffman_statistics_t* stat) {
    memset(stat, 0, sizeof(fjpeg_huffman_statistics_t));
    fjpeg_coeff_t block[64];
    uint32_t* ac = (scan->component == 0) ? stat->luma_ac : stat->chroma_ac;
    int eob_run = 0;

    auto flush_eob_run = [&]() {
        if (eob_run == 0) return;
        int r = 0;
        int temp = eob_run >> 1;
        while (temp > 0) {
            r++;
            temp >>= 1;
        }
        ac[r << 4]++;
        eob_run = 0;
    };

    fjpeg_for_each_block_component(context, scan->component, [&](int x, int y, int channel) {
        fjpeg_extract_coeff_8x8(context, block, x, y, channel);

        int last_newly_nonzero_k = scan->ss - 1;
        for (int k = scan->se; k >= scan->ss; k--) {
            int abs_val = fjpeg_abs_int(fjpeg_round(block[k]));
            if (abs_val >= (1 << scan->al) && abs_val < (1 << scan->ah)) {
                last_newly_nonzero_k = k;
                break;
            }
        }

        if (last_newly_nonzero_k < scan->ss) {
            eob_run++;
            if (eob_run == 32767) {
                flush_eob_run();
            }
            return;
        }

        flush_eob_run();

        int run = 0;
        for (int k = scan->ss; k <= last_newly_nonzero_k; k++) {
            int abs_val = fjpeg_abs_int(fjpeg_round(block[k]));
            if (abs_val < (1 << scan->al)) {
                run++;
                if (run == 16) {
                    ac[0xF0]++;
                    run = 0;
                }
                continue;
            }
            if (abs_val >= (1 << scan->ah)) {
                continue;
            }
            ac[(run << 4) | 1]++;
            run = 0;
        }

        if (last_newly_nonzero_k < scan->se) {
            eob_run++;
            if (eob_run == 32767) {
                flush_eob_run();
            }
        }
    });

    flush_eob_run();
}

// ---------------------------------------------------------------------------
// Entropy coding
// ---------------------------------------------------------------------------

static void fjpeg_progressive_encode_dc_first(fjpeg_bitstream* stream, const fjpeg_coeff_t* block,
                                              int al, const fjpeg_huffman_table_t* dc, int* last_dc) {
    int val = fjpeg_round_dc(block[0]) >> al;
    int diff = val - *last_dc;
    *last_dc = val;

    int temp = diff;
    int negative = 0;
    if (temp < 0) {
        temp = -temp;
        negative = 1;
    }
    int size = fjpeg_bit_length(temp);
    fjpeg_emit_symbol(stream, dc, size);
    if (size > 0) {
        int value = negative ? (~temp & ((1 << size) - 1)) : temp;
        stream->writeBits(value, size);
    }
}

static void fjpeg_progressive_encode_dc_refine(fjpeg_bitstream* stream, const fjpeg_coeff_t* block, int al) {
    int bit = (fjpeg_round_dc(block[0]) >> al) & 1;
    stream->writeBits(bit, 1);
}

static void fjpeg_progressive_encode_ac_first_scan(fjpeg_bitstream* stream, fjpeg_context* context,
                                                  const fjpeg_scan_t* scan, const fjpeg_huffman_table_t* ac) {
    fjpeg_coeff_t block[64];
    int eob_run = 0;

    auto flush_eob_run = [&]() {
        if (eob_run == 0) return;
        int r = 0;
        int temp = eob_run >> 1;
        while (temp > 0) {
            r++;
            temp >>= 1;
        }
        fjpeg_emit_symbol(stream, ac, r << 4);
        if (r > 0) {
            int extra = eob_run - (1 << r);
            stream->writeBits(extra, r);
        }
        eob_run = 0;
    };

    fjpeg_for_each_block_component(context, scan->component, [&](int x, int y, int channel) {
        fjpeg_extract_coeff_8x8(context, block, x, y, channel);

        int last_k = scan->ss - 1;
        for (int k = scan->se; k >= scan->ss; k--) {
            if (fjpeg_point_transform(fjpeg_round(block[k]), scan->al) != 0) {
                last_k = k;
                break;
            }
        }

        if (last_k < scan->ss) {
            eob_run++;
            if (eob_run == 32767) {
                flush_eob_run();
            }
            return;
        }

        flush_eob_run();

        int run = 0;
        for (int k = scan->ss; k <= last_k; k++) {
            int shifted = fjpeg_point_transform(fjpeg_round(block[k]), scan->al);
            if (shifted == 0) {
                run++;
                continue;
            }
            while (run > 15) {
                fjpeg_emit_symbol(stream, ac, 0xF0);
                run -= 16;
            }
            int temp = shifted;
            int negative = 0;
            if (temp < 0) {
                temp = -temp;
                negative = 1;
            }
            int size = fjpeg_bit_length(temp);
            fjpeg_emit_symbol(stream, ac, (run << 4) | size);
            int value = negative ? (~temp & ((1 << size) - 1)) : temp;
            stream->writeBits(value, size);
            run = 0;
        }

        if (last_k < scan->se) {
            eob_run++;
            if (eob_run == 32767) {
                flush_eob_run();
            }
        }
    });

    flush_eob_run();
}

static void fjpeg_progressive_encode_ac_refine_scan(fjpeg_bitstream* stream, fjpeg_context* context,
                                                   const fjpeg_scan_t* scan, const fjpeg_huffman_table_t* ac) {
    fjpeg_coeff_t block[64];
    int eob_run = 0;
    std::vector<int> eob_refine_bits;
    std::vector<int> pending_refine_bits;

    auto flush_eob_run = [&]() {
        if (eob_run == 0) return;
        int r = 0;
        int temp = eob_run >> 1;
        while (temp > 0) {
            r++;
            temp >>= 1;
        }
        fjpeg_emit_symbol(stream, ac, r << 4);
        if (r > 0) {
            int extra = eob_run - (1 << r);
            stream->writeBits(extra, r);
        }
        for (int bit : eob_refine_bits) {
            stream->writeBits(bit, 1);
        }
        eob_refine_bits.clear();
        eob_run = 0;
    };

    fjpeg_for_each_block_component(context, scan->component, [&](int x, int y, int channel) {
        fjpeg_extract_coeff_8x8(context, block, x, y, channel);

        int last_newly_nonzero_k = scan->ss - 1;
        for (int k = scan->se; k >= scan->ss; k--) {
            int abs_val = fjpeg_abs_int(fjpeg_round(block[k]));
            if (abs_val >= (1 << scan->al) && abs_val < (1 << scan->ah)) {
                last_newly_nonzero_k = k;
                break;
            }
        }

        if (last_newly_nonzero_k < scan->ss) {
            for (int k = scan->ss; k <= scan->se; k++) {
                int abs_val = fjpeg_abs_int(fjpeg_round(block[k]));
                if (abs_val >= (1 << scan->ah)) {
                    int bit = (abs_val >> scan->al) & 1;
                    eob_refine_bits.push_back(bit);
                }
            }
            eob_run++;
            if (eob_run == 32767) {
                flush_eob_run();
            }
            return;
        }

        flush_eob_run();

        int run = 0;
        pending_refine_bits.clear();

        for (int k = scan->ss; k <= last_newly_nonzero_k; k++) {
            int rounded = fjpeg_round(block[k]);
            int abs_val = fjpeg_abs_int(rounded);
            if (abs_val < (1 << scan->al)) {
                run++;
                if (run == 16) {
                    fjpeg_emit_symbol(stream, ac, 0xF0); // ZRL
                    for (int b : pending_refine_bits) {
                        stream->writeBits(b, 1);
                    }
                    pending_refine_bits.clear();
                    run = 0;
                }
                continue;
            }
            if (abs_val >= (1 << scan->ah)) {
                int bit = (abs_val >> scan->al) & 1;
                pending_refine_bits.push_back(bit);
                continue;
            }
            // Newly non-zero
            fjpeg_emit_symbol(stream, ac, (run << 4) | 1);
            int sign_bit = (rounded >= 0) ? 1 : 0;
            stream->writeBits(sign_bit, 1);
            for (int b : pending_refine_bits) {
                stream->writeBits(b, 1);
            }
            pending_refine_bits.clear();
            run = 0;
        }

        if (last_newly_nonzero_k < scan->se) {
            for (int k = last_newly_nonzero_k + 1; k <= scan->se; k++) {
                int abs_val = fjpeg_abs_int(fjpeg_round(block[k]));
                if (abs_val >= (1 << scan->ah)) {
                    int bit = (abs_val >> scan->al) & 1;
                    eob_refine_bits.push_back(bit);
                }
            }
            eob_run++;
            if (eob_run == 32767) {
                flush_eob_run();
            }
        }
    });

    flush_eob_run();
}

// ---------------------------------------------------------------------------
// Header helpers
// ---------------------------------------------------------------------------

static void fjpeg_progressive_write_dht(fjpeg_bitstream* stream, const fjpeg_short_huffman_table_t* table,
                                        int tc, int th) {
    int count = 0;
    for (int i = 0; i < 16; i++) {
        count += table->bits[i];
    }
    stream->writeBits(0xFFC4, 16);
    stream->writeBits(2 + 1 + 16 + count, 16);
    stream->writeBits(tc, 4);
    stream->writeBits(th, 4);
    for (int i = 0; i < 16; i++) {
        stream->writeBits(table->bits[i], 8);
    }
    for (int i = 0; i < count; i++) {
        stream->writeBits(table->val[i], 8);
    }
}

static void fjpeg_progressive_write_frame_header(fjpeg_bitstream* stream, fjpeg_context* context) {
    uint8_t tmp[64];

    // SOI
    stream->writeBits(0xFFD8, 16);

    // APP0 (JFIF)
    stream->writeBits(0xFFE0, 16);
    stream->writeBits(16, 16);
    stream->writeBits(0x4A4649, 24); // "JFI"
    stream->writeBits(0x4600, 16);   // "F\0"
    stream->writeBits(0x0102, 16);   // Version
    stream->writeBits(0x00, 8);      // Units
    stream->writeBits(0x0001, 16);   // X density
    stream->writeBits(0x0001, 16);   // Y density
    stream->writeBits(0x00, 8);      // X thumbnail
    stream->writeBits(0x00, 8);      // Y thumbnail

    // DQT - luminance
    stream->writeBits(0xFFDB, 16);
    stream->writeBits(67, 16);
    stream->writeBits(0, 4);
    stream->writeBits(0, 4);
    for (int i = 0; i < 64; i++) {
        tmp[fjpeg_zigzag_8x8[i]] = context->fjpeg_luminance_quantization_table[i];
    }
    for (int i = 0; i < 64; i++) {
        stream->writeBits(tmp[i], 8);
    }

    if (context->channels > 1) {
        // DQT - chrominance
        stream->writeBits(0xFFDB, 16);
        stream->writeBits(67, 16);
        stream->writeBits(0, 4);
        stream->writeBits(1, 4);
        for (int i = 0; i < 64; i++) {
            tmp[fjpeg_zigzag_8x8[i]] = context->fjpeg_chrominance_quantization_table[i];
        }
        for (int i = 0; i < 64; i++) {
            stream->writeBits(tmp[i], 8);
        }
    }

    // SOF2 (progressive DCT)
    stream->writeBits(0xFFC2, 16);
    stream->writeBits(context->channels == 1 ? 11 : 17, 16);
    stream->writeBits(8, 8); // 8 bits per sample
    stream->writeBits(context->height, 16);
    stream->writeBits(context->width, 16);
    stream->writeBits(context->channels, 8);
    for (int i = 0; i < context->channels; i++) {
        stream->writeBits(i + 1, 8);
        stream->writeBits(context->channels == 1 ? 0x11 : (i == 0 ? 0x22 : 0x11), 8);
        stream->writeBits(i == 0 ? 0 : 1, 8); // Quant table
    }
}

static void fjpeg_progressive_write_sos(fjpeg_bitstream* stream, fjpeg_context* context, const fjpeg_scan_t* scan) {
    const int interleaved = (scan->component < 0);
    const int ncomp = interleaved ? context->channels : 1;

    stream->writeBits(0xFFDA, 16);
    stream->writeBits(6 + 2 * ncomp, 16);
    stream->writeBits(ncomp, 8);
    if (interleaved) {
        for (int i = 0; i < context->channels; i++) {
            stream->writeBits(i + 1, 8);
            stream->writeBits(i == 0 ? 0x00 : 0x11, 8);
        }
    } else {
        stream->writeBits(scan->component + 1, 8);
        stream->writeBits(scan->component == 0 ? 0x00 : 0x11, 8);
    }
    stream->writeBits(scan->ss, 8);
    stream->writeBits(scan->se, 8);
    stream->writeBits((scan->ah << 4) | scan->al, 8);
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

bool fjpeg_generate_progressive(fjpeg_bitstream* stream, fjpeg_context* context) {
    // Apply optional trellis re-quantization before the scans are planned.
    fjpeg_trellis_optimize(context);

    fjpeg_progressive_write_frame_header(stream, context);

    fjpeg_huffman_statistics_t stats;

    for (int s = 0; s < fjpeg_progressive_script_len; s++) {
        const fjpeg_scan_t* scan = &fjpeg_progressive_script[s];
        if (scan->component >= context->channels) continue;
        const bool is_dc = (scan->component < 0);

        if (is_dc) {
            if (scan->ah == 0) {
                // First DC scan: gather statistics and define the DC tables.
                fjpeg_progressive_gather_dc_stats(context, scan, &stats);
                memset(context->fjpeg_huffman_luma_dc, 0, sizeof(context->fjpeg_huffman_luma_dc));
                memset(context->fjpeg_huffman_chroma_dc, 0, sizeof(context->fjpeg_huffman_chroma_dc));
                fjpeg_short_huffman_table_t luma =
                    fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_luma_dc, stats.luma_dc, 12);
                fjpeg_short_huffman_table_t chroma =
                    fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_chroma_dc, stats.chroma_dc, 12);
                fjpeg_progressive_write_dht(stream, &luma, 0, 0);
                if (context->channels > 1) {
                    fjpeg_progressive_write_dht(stream, &chroma, 0, 1);
                }
            }
            // DC refinement scans carry raw bits and need no Huffman tables.
        } else {
            if (scan->ah == 0) {
                fjpeg_progressive_gather_ac_stats(context, scan, &stats);
            } else {
                fjpeg_progressive_gather_ac_refine_stats(context, scan, &stats);
            }
            if (scan->component == 0) {
                memset(context->fjpeg_huffman_luma_ac, 0, sizeof(context->fjpeg_huffman_luma_ac));
                fjpeg_short_huffman_table_t luma =
                    fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_luma_ac, stats.luma_ac, 256);
                fjpeg_progressive_write_dht(stream, &luma, 1, 0);
            } else {
                memset(context->fjpeg_huffman_chroma_ac, 0, sizeof(context->fjpeg_huffman_chroma_ac));
                fjpeg_short_huffman_table_t chroma =
                    fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_chroma_ac, stats.chroma_ac, 256);
                fjpeg_progressive_write_dht(stream, &chroma, 1, 1);
            }
        }

        fjpeg_progressive_write_sos(stream, context, scan);

        stream->avoidFF = true;

        if (is_dc) {
            int last_dc[3] = {0, 0, 0};
            fjpeg_coeff_t block[64];
            if (scan->ah == 0) {
                fjpeg_for_each_block_all(context, [&](int x, int y, int channel) {
                    fjpeg_extract_coeff_8x8(context, block, x, y, channel);
                    const fjpeg_huffman_table_t* table = (channel == 0)
                        ? context->fjpeg_huffman_luma_dc : context->fjpeg_huffman_chroma_dc;
                    fjpeg_progressive_encode_dc_first(stream, block, scan->al, table, &last_dc[channel]);
                });
            } else {
                fjpeg_for_each_block_all(context, [&](int x, int y, int channel) {
                    fjpeg_extract_coeff_8x8(context, block, x, y, channel);
                    fjpeg_progressive_encode_dc_refine(stream, block, scan->al);
                });
            }
        } else {
            const fjpeg_huffman_table_t* table = (scan->component == 0)
                ? context->fjpeg_huffman_luma_ac : context->fjpeg_huffman_chroma_ac;
            if (scan->ah == 0) {
                fjpeg_progressive_encode_ac_first_scan(stream, context, scan, table);
            } else {
                fjpeg_progressive_encode_ac_refine_scan(stream, context, scan, table);
            }
        }

        stream->padToByte();
        stream->avoidFF = false;
    }

    // EOI
    stream->writeBits(0xFFD9, 16);
    stream->flushToFile();

    return true;
}
