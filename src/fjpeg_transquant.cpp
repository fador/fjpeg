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
#include <cassert>
#include <vector>

#include "fjpeg.h"



fjpeg_coeff_t* fjpeg_zigzag8x8(fjpeg_coeff_t* block, fjpeg_coeff_t* out) {

    for (int i = 0; i < 64; i++) {
        out[fjpeg_zigzag_8x8[i]] = block[i];
    }

    return out;
}

fjpeg_coeff_t* fjpeg_izigzag8x8(fjpeg_coeff_t* block, fjpeg_coeff_t* out) {

    for (int i = 0; i < 64; i++) {
        out[i] = block[fjpeg_zigzag_8x8[i]];
    }

    return out;
}

fjpeg_pixel_t* fjpeg_extract_8x8(fjpeg_context* context, fjpeg_pixel_t* output, int x, int y, int channel) {    

    fjpeg_pixel_t* image = channel==0?context->fjpeg_y:channel==1?context->fjpeg_cb:context->fjpeg_cr;
    const int input_width = channel==0?context->width:context->width/2;
    const int input_height = channel==0?context->height:context->height/2;

    for (int j = 0; j < 8; j++) {
        int sy = y + j;
        if (sy >= input_height) sy = input_height - 1; // replicate bottom edge
        for (int i = 0; i < 8; i++) {
            int sx = x + i;
            if (sx >= input_width) sx = input_width - 1; // replicate right edge
            output[j * 8 + i] = image[sy * input_width + sx];
        }
    }
    return output;
}

fjpeg_coeff_t* fjpeg_extract_coeff_8x8(fjpeg_context* context, fjpeg_coeff_t* output, int x, int y, int channel) {    

    fjpeg_coeff_t* image = channel==0?context->fjpeg_ydct:channel==1?context->fjpeg_cbdct:context->fjpeg_crdct;
    const int input_width = channel==0?context->padded_width:context->padded_width/2;

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            output[j * 8 + i] = image[(y + j) * input_width  + (x + i)];
        }
    }
    return output;
}

bool fjpeg_store_coeff_8x8(fjpeg_context* context, fjpeg_coeff_t* input, int x, int y, int channel) {

    fjpeg_coeff_t* image = channel==0?context->fjpeg_ydct:channel==1?context->fjpeg_cbdct:context->fjpeg_crdct;
    const int input_width = channel==0?context->padded_width:context->padded_width/2;

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            image[(y + j) * input_width + (x + i)] = input[j * 8 + i];
        }
    }
    return true;
}


fjpeg_coeff_t* fjpeg_dct8x8(fjpeg_context* context, fjpeg_pixel_t* block, fjpeg_coeff_t* out) {
    // Separable 2D DCT: transform rows first, then columns. This reduces the
    // work from O(N^4) to O(2 * N^3) per block (4096 -> 1024 multiplies).
    float tmp[FJPEG_BLOCK_SIZE * FJPEG_BLOCK_SIZE];
    const float c = 1.0f / sqrtf(2.0f);

    for (int y = 0; y < FJPEG_BLOCK_SIZE; y++) {
        for (int u = 0; u < FJPEG_BLOCK_SIZE; u++) {
            float sum = 0.0f;
            for (int x = 0; x < FJPEG_BLOCK_SIZE; x++) {
                sum += ((float)block[y * FJPEG_BLOCK_SIZE + x] - 128.0f) * context->precalc_cos[x][u];
            }
            tmp[y * FJPEG_BLOCK_SIZE + u] = sum;
        }
    }

    for (int u = 0; u < FJPEG_BLOCK_SIZE; u++) {
        for (int v = 0; v < FJPEG_BLOCK_SIZE; v++) {
            float sum = 0.0f;
            for (int y = 0; y < FJPEG_BLOCK_SIZE; y++) {
                sum += tmp[y * FJPEG_BLOCK_SIZE + u] * context->precalc_cos[y][v];
            }
            float cu = (u == 0) ? c : 1.0f;
            float cv = (v == 0) ? c : 1.0f;
            out[v * FJPEG_BLOCK_SIZE + u] = 0.25f * sum * cu * cv;
        }
    }
    return out;
}

fjpeg_pixel_t* fjpeg_idct8x8(fjpeg_context* context, fjpeg_coeff_t* block, fjpeg_pixel_t* out) {
    for (int y = 0; y < FJPEG_BLOCK_SIZE; y++) {
        for (int x = 0; x < FJPEG_BLOCK_SIZE; x++) {
            float sum = 0.0f;
            for (int v = 0; v < FJPEG_BLOCK_SIZE; v++) {
                for (int u = 0; u < FJPEG_BLOCK_SIZE; u++) {
                    float cu = (u == 0) ? 1.0f / sqrtf(2.f) : 1.0f;  // Scaling factors
                    float cv = (v == 0) ? 1.0f / sqrtf(2.f) : 1.0f;
                    sum += cu * cv * block[v*FJPEG_BLOCK_SIZE+u] * context->precalc_cos[x][u]  * context->precalc_cos[y][v] ;
                }
            }
            out[y*FJPEG_BLOCK_SIZE+x] = (fjpeg_pixel_t)((0.25f * sum)+128.0f);  // Apply constant factor
        }
    }
    return out;
}

fjpeg_coeff_t* fjpeg_quant8x8(fjpeg_context* context, fjpeg_coeff_t* input, fjpeg_coeff_t *output, int table) {

    uint8_t *quant_table = table == 0 ? context->fjpeg_luminance_quantization_table : context->fjpeg_chrominance_quantization_table;

    for (int i = 0; i < 64; i++) {
        output[i] = input[i] / (float)quant_table[i];
    }

    return output;
}

fjpeg_coeff_t* fjpeg_dequant8x8(fjpeg_context* context, fjpeg_coeff_t* input, fjpeg_coeff_t *output, int table) {

    uint8_t *quant_table = table == 0 ? context->fjpeg_luminance_quantization_table : context->fjpeg_chrominance_quantization_table;

    for (int i = 0; i < 64; i++) {
        output[i] = input[i] * quant_table[i];
    }

    return output;
}


// Number of bits needed for the magnitude of v (JPEG "size" category).
static inline int fjpeg_coeff_size(int v) {
    int a = v < 0 ? -v : v;
    int s = 0;
    while (a) {
        a >>= 1;
        s++;
    }
    return s;
}

// Code length of a symbol, with a safe fallback for symbols that are absent
// from the current (pre-trellis) Huffman table.
static inline int fjpeg_huff_symbol_len(const fjpeg_huffman_table_t* table, int symbol) {
    return table[symbol].len ? table[symbol].len : 16;
}

void fjpeg_trellis_quant_block(fjpeg_context* context, const fjpeg_coeff_t* block,
                               const fjpeg_huffman_table_t* huff_ac,
                               int channel, fjpeg_coeff_t* out) {
    const int ZIGZAG_COUNT = FJPEG_BLOCK_SIZE * FJPEG_BLOCK_SIZE; // 64
    const float INF = 1.0e30f;
    const float eob_rate = (float)fjpeg_huff_symbol_len(huff_ac, 0x00);
    const float zrl_rate = (float)fjpeg_huff_symbol_len(huff_ac, 0xF0);

    // DC is coded independently of the AC run lengths, so round it as usual.
    out[0] = (fjpeg_coeff_t)fjpeg_round(block[0]);

    // The reconstruction error is Q^2 * (c - n)^2, so weight the distortion by
    // the squared quantization step. The coefficients are in zigzag order.
    const uint8_t* qtab = (channel == 0) ? context->fjpeg_luminance_quantization_table
                                         : context->fjpeg_chrominance_quantization_table;
    int inv_zigzag[ZIGZAG_COUNT];
    for (int n = 0; n < ZIGZAG_COUNT; n++) {
        inv_zigzag[fjpeg_zigzag_8x8[n]] = n;
    }
    float weight[ZIGZAG_COUNT];
    float mean_weight = 0.0f;
    for (int p = 0; p < ZIGZAG_COUNT; p++) {
        float q = (float)qtab[inv_zigzag[p]];
        weight[p] = q * q;
        mean_weight += weight[p];
    }
    mean_weight /= (float)ZIGZAG_COUNT;
    // Scale the multiplier with the quantizer so the tradeoff is independent of
    // the quality setting.
    const float lambda = context->trellis_lambda * mean_weight;

    // Suffix distortion of zeroing coefficients p..63.
    float suffix[ZIGZAG_COUNT + 1];
    suffix[ZIGZAG_COUNT] = 0.0f;
    for (int p = ZIGZAG_COUNT - 1; p >= 1; p--) {
        suffix[p] = suffix[p + 1] + weight[p] * block[p] * block[p];
    }

    // cost[r] = best Lagrangian cost ending with r pending zeros (0..15).
    float prev[16];
    float cur[16];
    for (int r = 0; r < 16; r++) {
        prev[r] = INF;
    }
    prev[0] = 0.0f;

    // Back-pointers for path reconstruction.
    int16_t choice[ZIGZAG_COUNT][16];
    int8_t from[ZIGZAG_COUNT][16];
    memset(choice, 0, sizeof(choice));
    memset(from, 0, sizeof(from));

    // Terminating before any AC coefficient (all-zero block).
    float best = suffix[1] + eob_rate;
    int best_p = 0;
    int best_r = 0;

    for (int p = 1; p < ZIGZAG_COUNT; p++) {
        for (int r = 0; r < 16; r++) {
            cur[r] = INF;
        }

        const float cv = block[p];
        const float w = weight[p];
        const float dist_zero = w * cv * cv;

        // Option 1: quantize this coefficient to zero.
        for (int r = 0; r < 16; r++) {
            if (prev[r] >= INF) continue;
            int nr;
            float nc = prev[r] + dist_zero;
            if (r < 15) {
                nr = r + 1;
            } else {
                nr = 0;
                nc += zrl_rate; // a full run of 16 zeros emits ZRL
            }
            if (nc < cur[nr]) {
                cur[nr] = nc;
                choice[p][nr] = 0;
                from[p][nr] = (int8_t)r;
            }
        }

        // Option 2: keep a nonzero coefficient. Candidates bracket the
        // nearest integer and the neighbouring size categories, since the
        // rate jumps when the magnitude size changes.
        float av = cv < 0.0f ? -cv : cv;
        int nearest = (int)(av + 0.5f);
        if (nearest > 0) {
            int sign = cv < 0.0f ? -1 : 1;
            int nearest_size = fjpeg_coeff_size(nearest);
            int cand[3];
            int ncand = 0;
            cand[ncand++] = sign * nearest;
            if (nearest_size > 1) {
                cand[ncand++] = sign * ((1 << (nearest_size - 1)) - 1);
            }
            if (nearest_size < 10) {
                cand[ncand++] = sign * (1 << nearest_size);
            }
            for (int ci = 0; ci < ncand; ci++) {
                int n = cand[ci];
                int sz = fjpeg_coeff_size(n);
                float d = cv - (float)n;
                d = w * d * d;
                for (int r = 0; r < 16; r++) {
                    if (prev[r] >= INF) continue;
                    float rate = (float)fjpeg_huff_symbol_len(huff_ac, (r << 4) | sz) + (float)sz;
                    float nc = prev[r] + d + lambda * rate;
                    if (nc < cur[0]) {
                        cur[0] = nc;
                        choice[p][0] = (int16_t)n;
                        from[p][0] = (int8_t)r;
                    }
                }
            }
        }

        // Consider ending the block here (EOB) with the remaining AC zeros.
        for (int r = 0; r < 16; r++) {
            if (cur[r] >= INF) continue;
            float term = cur[r] + suffix[p + 1];
            if (p < ZIGZAG_COUNT - 1 || r > 0) {
                term += eob_rate;
            }
            if (term < best) {
                best = term;
                best_p = p;
                best_r = r;
            }
        }

        for (int r = 0; r < 16; r++) {
            prev[r] = cur[r];
        }
    }

    // Reconstruct the chosen coefficients.
    for (int p = 1; p < ZIGZAG_COUNT; p++) {
        out[p] = 0.0f;
    }
    int p = best_p;
    int r = best_r;
    while (p >= 1) {
        out[p] = (fjpeg_coeff_t)choice[p][r];
        int pr = from[p][r];
        p--;
        r = pr;
    }
}


bool fjpeg_transquant_input(fjpeg_context* context) {
    
    fjpeg_pixel_t cur_block[64];
    fjpeg_coeff_t dct_block[64];
    fjpeg_coeff_t dct_block2[64];
    
    for(int y = 0; y < context->padded_height; y+=8) {
        for(int x = 0; x < context->padded_width; x+=8) {
            fjpeg_extract_8x8(context, cur_block, x, y, 0);
            fjpeg_dct8x8(context, cur_block, dct_block);
            fjpeg_quant8x8(context, dct_block,dct_block2, 0);
            fjpeg_zigzag8x8(dct_block2, dct_block);
            fjpeg_store_coeff_8x8(context, dct_block, x, y, 0);
            #ifdef FJPEG_DEBUG_DCT_BLOCK
            fjpeg_izigzag8x8(dct_block, dct_block2);
            fjpeg_dequant8x8(context, dct_block2, dct_block, 0);
            fjpeg_idct8x8(context, dct_block, cur_block);
            for (int j = 0; j < 8; j++) {
                for (int i = 0; i < 8; i++) {
                    image[(y + j) * context->width + (x + i)] = cur_block[j * 8 + i];
                }
            }
            #endif
        }
    }

    if(context->channels == 3) {
        for(int y = 0; y < context->padded_height/2; y+=8) {
            for(int x = 0; x < context->padded_width/2; x+=8) {
                fjpeg_extract_8x8(context, cur_block, x, y, 1);
                fjpeg_dct8x8(context, cur_block, dct_block);
                fjpeg_quant8x8(context, dct_block,dct_block2, 1);
                fjpeg_zigzag8x8(dct_block2, dct_block);
                fjpeg_store_coeff_8x8(context, dct_block, x, y, 1);
            }
        }

        for(int y = 0; y < context->padded_height/2; y+=8) {
            for(int x = 0; x < context->padded_width/2; x+=8) {
                fjpeg_extract_8x8(context, cur_block, x, y, 2);
                fjpeg_dct8x8(context, cur_block, dct_block);
                fjpeg_quant8x8(context, dct_block,dct_block2, 2);
                fjpeg_zigzag8x8(dct_block2, dct_block);
                fjpeg_store_coeff_8x8(context, dct_block, x, y, 2);
            }
        }
    }

    return true;
}