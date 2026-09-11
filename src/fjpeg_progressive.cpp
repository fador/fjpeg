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

// The DC coefficients are sent with one successive-approximation bit, then each
// component's AC coefficients are sent in two spectral bands. (AC successive
// approximation is not implemented; every AC scan carries full precision.)
static const fjpeg_scan_t fjpeg_progressive_script[] = {
    {-1,  0,  0, 0, 1}, // DC high bits
    {-1,  0,  0, 1, 0}, // DC low bit
    { 0,  1,  5, 0, 0}, // luma AC 1..5
    { 0,  6, 63, 0, 0}, // luma AC 6..63
    { 1,  1,  5, 0, 0}, // Cb AC 1..5
    { 1,  6, 63, 0, 0}, // Cb AC 6..63
    { 2,  1,  5, 0, 0}, // Cr AC 1..5
    { 2,  6, 63, 0, 0}, // Cr AC 6..63
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
        int val = fjpeg_round(block[0]) >> scan->al;
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
    fjpeg_for_each_block_component(context, scan->component, [&](int x, int y, int channel) {
        fjpeg_extract_coeff_8x8(context, block, x, y, channel);
        uint32_t* ac = (channel == 0) ? stat->luma_ac : stat->chroma_ac;
        int run = 0;
        for (int k = scan->ss; k <= scan->se; k++) {
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
        if (run > 0) {
            ac[0x00]++; // EOB
        }
    });
}

// ---------------------------------------------------------------------------
// Entropy coding
// ---------------------------------------------------------------------------

static void fjpeg_progressive_encode_dc_first(fjpeg_bitstream* stream, const fjpeg_coeff_t* block,
                                              int al, const fjpeg_huffman_table_t* dc, int* last_dc) {
    int val = fjpeg_round(block[0]) >> al;
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
    int bit = (fjpeg_round(block[0]) >> al) & 1;
    stream->writeBits(bit, 1);
}

static void fjpeg_progressive_encode_ac_first(fjpeg_bitstream* stream, const fjpeg_coeff_t* block,
                                              const fjpeg_scan_t* scan, const fjpeg_huffman_table_t* ac) {
    int run = 0;
    for (int k = scan->ss; k <= scan->se; k++) {
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
    if (run > 0) {
        fjpeg_emit_symbol(stream, ac, 0x00); // EOB
    }
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
            if (scan->ah != 0) {
                fprintf(stderr, "Error: AC successive approximation is not supported\n");
                return false;
            }
            fjpeg_progressive_gather_ac_stats(context, scan, &stats);
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
            fjpeg_coeff_t block[64];
            fjpeg_for_each_block_component(context, scan->component, [&](int x, int y, int channel) {
                fjpeg_extract_coeff_8x8(context, block, x, y, channel);
                const fjpeg_huffman_table_t* table = (channel == 0)
                    ? context->fjpeg_huffman_luma_ac : context->fjpeg_huffman_chroma_ac;
                fjpeg_progressive_encode_ac_first(stream, block, scan, table);
            });
        }

        stream->padToByte();
        stream->avoidFF = false;
    }

    // EOI
    stream->writeBits(0xFFD9, 16);
    stream->flushToFile();

    return true;
}
