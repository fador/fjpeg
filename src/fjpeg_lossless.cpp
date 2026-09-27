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
#include <algorithm>

#include "fjpeg.h"
#include "fjpeg_global.h"
#include "fjpeg_bitstream.h"
#include "fjpeg_huffman.h"
#include "fjpeg_lossless.h"

// Category and VLI bit calculation for a difference value
static inline void fjpeg_diff_to_vli(int diff, int* out_category, int* out_bits) {
    if (diff == 0) {
        *out_category = 0;
        *out_bits = 0;
        return;
    }
    int abs_diff = diff < 0 ? -diff : diff;
    int cat = 0;
    while (abs_diff > 0) {
        cat++;
        abs_diff >>= 1;
    }
    *out_category = cat;
    if (diff < 0) {
        *out_bits = diff + (1 << cat) - 1;
    } else {
        *out_bits = diff;
    }
}

// Compute difference modulo 2^P for a single sample at (px, py) in a plane
template <typename T>
static inline int fjpeg_calc_diff(const T* plane, int width, int height, int px, int py, int predictor, int precision) {
    (void)height;
    int pred;
    if (py == 0 && px == 0) {
        pred = 1 << (precision - 1); // 128 for 8-bit, 2048 for 12-bit
    } else if (py == 0) {
        pred = plane[py * width + (px - 1)];
    } else if (px == 0) {
        pred = plane[(py - 1) * width + px];
    } else {
        int l  = plane[py * width + (px - 1)];
        int t  = plane[(py - 1) * width + px];
        int tl = plane[(py - 1) * width + (px - 1)];
        pred = fjpeg_lossless_predict(tl, t, l, predictor);
    }
    int diff = (int)plane[py * width + px] - pred;
    int mask = (1 << precision) - 1;
    diff = diff & mask;
    if (diff >= (1 << (precision - 1))) {
        diff -= (1 << precision);
    }
    return diff;
}

template <typename T>
static void fjpeg_lossless_gather_freq(const T* y_plane, const T* cb_plane, const T* cr_plane,
                                      int w, int h, int cw, int ch, int channels,
                                      int predictor, int precision,
                                      uint32_t freq_y[17], uint32_t freq_c[17]) {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int diff = fjpeg_calc_diff(y_plane, w, h, x, y, predictor, precision);
            int cat, bits;
            fjpeg_diff_to_vli(diff, &cat, &bits);
            if (cat <= 16) freq_y[cat]++;
        }
    }

    if (channels == 3) {
        for (int y = 0; y < ch; y++) {
            for (int x = 0; x < cw; x++) {
                int diff_cb = fjpeg_calc_diff(cb_plane, cw, ch, x, y, predictor, precision);
                int cat, bits;
                fjpeg_diff_to_vli(diff_cb, &cat, &bits);
                if (cat <= 16) freq_c[cat]++;

                int diff_cr = fjpeg_calc_diff(cr_plane, cw, ch, x, y, predictor, precision);
                fjpeg_diff_to_vli(diff_cr, &cat, &bits);
                if (cat <= 16) freq_c[cat]++;
            }
        }
    }
}

template <typename T>
static void fjpeg_lossless_encode_scan(fjpeg_bitstream* stream, const T* y_plane, const T* cb_plane, const T* cr_plane,
                                       int w, int h, int cw, int ch, int ncomp, int predictor, int precision,
                                       const fjpeg_huffman_table_t* huff_table_y, const fjpeg_huffman_table_t* huff_table_c) {
    int mb_w = (w + 1) / 2;
    int mb_h = (h + 1) / 2;

    for (int mb_y = 0; mb_y < mb_h; mb_y++) {
        for (int mb_x = 0; mb_x < mb_w; mb_x++) {
            // Y component: 2x2 samples per MCU
            for (int dy = 0; dy < 2; dy++) {
                for (int dx = 0; dx < 2; dx++) {
                    int px = mb_x * 2 + dx;
                    int py = mb_y * 2 + dy;
                    if (px < w && py < h) {
                        int diff = fjpeg_calc_diff(y_plane, w, h, px, py, predictor, precision);
                        int cat, bits;
                        fjpeg_diff_to_vli(diff, &cat, &bits);
                        stream->writeBits(huff_table_y[cat].code, huff_table_y[cat].len);
                        if (cat > 0) {
                            stream->writeBits(bits, cat);
                        }
                    }
                }
            }

            // Cb & Cr components: 1 sample each per MCU
            if (ncomp == 3) {
                int px = mb_x;
                int py = mb_y;
                if (px < cw && py < ch) {
                    // Cb
                    int diff_cb = fjpeg_calc_diff(cb_plane, cw, ch, px, py, predictor, precision);
                    int cat_cb, bits_cb;
                    fjpeg_diff_to_vli(diff_cb, &cat_cb, &bits_cb);
                    stream->writeBits(huff_table_c[cat_cb].code, huff_table_c[cat_cb].len);
                    if (cat_cb > 0) {
                        stream->writeBits(bits_cb, cat_cb);
                    }

                    // Cr
                    int diff_cr = fjpeg_calc_diff(cr_plane, cw, ch, px, py, predictor, precision);
                    int cat_cr, bits_cr;
                    fjpeg_diff_to_vli(diff_cr, &cat_cr, &bits_cr);
                    stream->writeBits(huff_table_c[cat_cr].code, huff_table_c[cat_cr].len);
                    if (cat_cr > 0) {
                        stream->writeBits(bits_cr, cat_cr);
                    }
                }
            }
        }
    }
}

// Evaluate total entropy bits for a candidate predictor
static uint64_t fjpeg_eval_predictor_cost(fjpeg_context* context, int predictor) {
    uint32_t freq_y[17] = {0};
    uint32_t freq_c[17] = {0};

    int w = context->width;
    int h = context->height;
    int cw = w / 2;
    int ch = h / 2;
    int prec = context->bit_depth;

    if (prec == 12) {
        fjpeg_lossless_gather_freq(context->fjpeg_y16, context->fjpeg_cb16, context->fjpeg_cr16,
                                  w, h, cw, ch, context->channels, predictor, prec, freq_y, freq_c);
    } else {
        fjpeg_lossless_gather_freq(context->fjpeg_y, context->fjpeg_cb, context->fjpeg_cr,
                                  w, h, cw, ch, context->channels, predictor, prec, freq_y, freq_c);
    }

    // Estimate total bits: Shannon entropy estimate + payload bits
    uint64_t total_bits = 0;
    auto add_cost = [&](const uint32_t* freq, int count) {
        uint64_t sum = 0;
        for (int i = 0; i < count; i++) sum += freq[i];
        if (sum == 0) return;
        double total_syms = (double)sum;
        for (int i = 0; i < count; i++) {
            if (freq[i] > 0) {
                double p = (double)freq[i] / total_syms;
                double huff_len = -log2(p);
                if (huff_len < 1.0) huff_len = 1.0;
                if (huff_len > 16.0) huff_len = 16.0;
                total_bits += (uint64_t)(freq[i] * (round(huff_len) + i));
            }
        }
    };

    add_cost(freq_y, 17);
    if (context->channels == 3) {
        add_cost(freq_c, 17);
    }
    return total_bits;
}

bool fjpeg_generate_lossless(fjpeg_bitstream* stream, fjpeg_context* context, int predictor) {
    if (!stream || !context || (!context->fjpeg_y && !context->fjpeg_y16)) return false;

    // Auto-select predictor if not specified or invalid
    if (predictor < 1 || predictor > 7) {
        static const int candidate_predictors[] = {1, 2, 4, 5, 6, 7};
        uint64_t best_cost = UINT64_MAX;
        int best_p = 1;
        for (int p : candidate_predictors) {
            uint64_t cost = fjpeg_eval_predictor_cost(context, p);
            if (cost < best_cost) {
                best_cost = cost;
                best_p = p;
            }
        }
        predictor = best_p;
    }

    int w = context->width;
    int h = context->height;
    int cw = w / 2;
    int ch = h / 2;
    int prec = context->bit_depth;

    // Gather frequency statistics for the chosen predictor
    uint32_t freq_y[17] = {0};
    uint32_t freq_c[17] = {0};

    if (prec == 12) {
        fjpeg_lossless_gather_freq(context->fjpeg_y16, context->fjpeg_cb16, context->fjpeg_cr16,
                                  w, h, cw, ch, context->channels, predictor, prec, freq_y, freq_c);
    } else {
        fjpeg_lossless_gather_freq(context->fjpeg_y, context->fjpeg_cb, context->fjpeg_cr,
                                  w, h, cw, ch, context->channels, predictor, prec, freq_y, freq_c);
    }

    // Generate optimal Huffman tables
    fjpeg_huffman_table_t huff_table_y[256];
    fjpeg_huffman_table_t huff_table_c[256];
    memset(huff_table_y, 0, sizeof(huff_table_y));
    memset(huff_table_c, 0, sizeof(huff_table_c));

    fjpeg_short_huffman_table_t short_y = fjpeg_generate_huffman_from_stats(huff_table_y, freq_y, 17);
    fjpeg_short_huffman_table_t short_c = short_y;
    if (context->channels == 3) {
        short_c = fjpeg_generate_huffman_from_stats(huff_table_c, freq_c, 17);
    }

    // 1. SOI
    stream->writeBits(0xFFD8, 16);

    // 2. APP0 (JFIF)
    stream->writeBits(0xFFE0, 16);
    stream->writeBits(16, 16);
    stream->writeBits(0x4A4649, 24); // "JFI"
    stream->writeBits(0x4600, 16);   // "F\0"
    stream->writeBits(0x0101, 16);   // Version 1.01
    stream->writeBits(0x00, 8);      // Units
    stream->writeBits(0x0001, 16);   // X density
    stream->writeBits(0x0001, 16);   // Y density
    stream->writeBits(0x00, 8);      // X thumbnail
    stream->writeBits(0x00, 8);      // Y thumbnail

    // 3. SOF3 (Lossless Huffman, marker 0xFFC3)
    int ncomp = context->channels;
    stream->writeBits(0xFFC3, 16);
    stream->writeBits(8 + 3 * ncomp, 16); // Length
    stream->writeBits(prec, 8);           // Sample precision (8 or 12)
    stream->writeBits(h, 16);             // Height
    stream->writeBits(w, 16);             // Width
    stream->writeBits(ncomp, 8);          // Number of components

    for (int i = 0; i < ncomp; i++) {
        stream->writeBits(i + 1, 8); // Component ID (1, 2, 3)
        // Sampling factors: 4:2:0 -> Y=0x22, Cb/Cr=0x11 (or 0x11 if 1 component)
        uint8_t hv = (ncomp == 1) ? 0x11 : (i == 0 ? 0x22 : 0x11);
        stream->writeBits(hv, 8);
        stream->writeBits(0, 8); // Quantization table selector (0 for lossless)
    }

    // 4. DHT (Huffman tables)
    // Table 0: Luma DC / difference
    int count_y = 0;
    for (int i = 0; i < 16; i++) count_y += short_y.bits[i];

    int count_c = 0;
    if (ncomp == 3) {
        for (int i = 0; i < 16; i++) count_c += short_c.bits[i];
    }

    int dht_len = 2 + (1 + 16 + count_y);
    if (ncomp == 3) {
        dht_len += (1 + 16 + count_c);
    }

    stream->writeBits(0xFFC4, 16);
    stream->writeBits(dht_len, 16);

    // Write table 0 (DC Luma)
    stream->writeBits(0, 4); // Tc = 0 (DC)
    stream->writeBits(0, 4); // Th = 0
    for (int i = 0; i < 16; i++) stream->writeBits(short_y.bits[i], 8);
    for (int i = 0; i < count_y; i++) stream->writeBits(short_y.val[i], 8);

    // Write table 1 (DC Chroma)
    if (ncomp == 3) {
        stream->writeBits(0, 4); // Tc = 0 (DC)
        stream->writeBits(1, 4); // Th = 1
        for (int i = 0; i < 16; i++) stream->writeBits(short_c.bits[i], 8);
        for (int i = 0; i < count_c; i++) stream->writeBits(short_c.val[i], 8);
    }

    // 5. SOS (Start of Scan)
    stream->writeBits(0xFFDA, 16);
    stream->writeBits(6 + 2 * ncomp, 16);
    stream->writeBits(ncomp, 8);

    for (int i = 0; i < ncomp; i++) {
        stream->writeBits(i + 1, 8); // Component selector
        uint8_t td = (i == 0) ? 0 : 1;
        stream->writeBits((td << 4), 8); // Td (DC table), Ta = 0
    }

    stream->writeBits(predictor, 8); // Ss = Predictor selection (1..7)
    stream->writeBits(0, 8);         // Se = 0
    stream->writeBits(0, 8);         // Ah = 0, Al = 0 (Pt = 0)

    // 6. Entropy-coded bitstream
    stream->avoidFF = true;

    if (prec == 12) {
        fjpeg_lossless_encode_scan(stream, context->fjpeg_y16, context->fjpeg_cb16, context->fjpeg_cr16,
                                  w, h, cw, ch, ncomp, predictor, prec, huff_table_y, huff_table_c);
    } else {
        fjpeg_lossless_encode_scan(stream, context->fjpeg_y, context->fjpeg_cb, context->fjpeg_cr,
                                  w, h, cw, ch, ncomp, predictor, prec, huff_table_y, huff_table_c);
    }

    stream->padToByte();
    stream->avoidFF = false;

    // 7. EOI
    stream->writeBits(0xFFD9, 16);
    stream->flushToFile();

    return true;
}
