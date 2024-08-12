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

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            output[j * 8 + i] = image[(y + j) * input_width  + (x + i)];
        }
    }
    return output;
}

fjpeg_coeff_t* fjpeg_extract_coeff_8x8(fjpeg_context* context, fjpeg_coeff_t* output, int x, int y, int channel) {    

    fjpeg_coeff_t* image = channel==0?context->fjpeg_ydct:channel==1?context->fjpeg_cbdct:context->fjpeg_crdct;
    const int input_width = channel==0?context->width:context->width/2;

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            output[j * 8 + i] = image[(y + j) * input_width  + (x + i)];
        }
    }
    return output;
}

bool fjpeg_store_coeff_8x8(fjpeg_context* context, fjpeg_coeff_t* input, int x, int y, int channel) {

    fjpeg_coeff_t* image = channel==0?context->fjpeg_ydct:channel==1?context->fjpeg_cbdct:context->fjpeg_crdct;
    const int input_width = channel==0?context->width:context->width/2;

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            image[(y + j) * input_width + (x + i)] = input[j * 8 + i];
        }
    }
    return true;
}


fjpeg_coeff_t* fjpeg_dct8x8(fjpeg_context* context, fjpeg_pixel_t* block, fjpeg_coeff_t* out) {
    for (int v = 0; v < FJPEG_BLOCK_SIZE; v++) {
        for (int u = 0; u < FJPEG_BLOCK_SIZE; u++) {
            float sum = 0.0f;
            float cu = (u == 0) ? 1.0f / sqrtf(2.f) : 1.0f;  // Scaling factors
            float cv = (v == 0) ? 1.0f / sqrtf(2.f) : 1.0f;
            for (int x = 0; x < FJPEG_BLOCK_SIZE; x++) {
                for (int y = 0; y < FJPEG_BLOCK_SIZE; y++) {                    
                    sum += (((float)block[y*FJPEG_BLOCK_SIZE+x])-128.f) *context->precalc_cos[x][u]  * context->precalc_cos[y][v] ;
                }
            }
            out[v*FJPEG_BLOCK_SIZE+u] = 0.25f * sum  * cu * cv;  // Apply constant factor
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


bool fjpeg_transquant_input(fjpeg_context* context) {
    
    fjpeg_pixel_t cur_block[64];
    fjpeg_coeff_t dct_block[64];
    fjpeg_coeff_t dct_block2[64];
    
    for(int y = 0; y < context->height; y+=8) {
        for(int x = 0; x < context->width; x+=8) {
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
        for(int y = 0; y < context->height/2; y+=8) {
            for(int x = 0; x < context->width/2; x+=8) {
                fjpeg_extract_8x8(context, cur_block, x, y, 1);
                fjpeg_dct8x8(context, cur_block, dct_block);
                fjpeg_quant8x8(context, dct_block,dct_block2, 1);
                fjpeg_zigzag8x8(dct_block2, dct_block);
                fjpeg_store_coeff_8x8(context, dct_block, x, y, 1);
            }
        }

        for(int y = 0; y < context->height/2; y+=8) {
            for(int x = 0; x < context->width/2; x+=8) {
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