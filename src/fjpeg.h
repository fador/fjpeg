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
#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>

#include "fjpeg_global.h"
#include "fjpeg_huffman.h"



class fjpeg_context {
    public:
    FILE* input;
    FILE* output;
    int width;
    int height;
    int quality;
    int channels;

    uint8_t fjpeg_luminance_quantization_table[64];
    uint8_t fjpeg_chrominance_quantization_table[64];
    fjpeg_huffman_table_t fjpeg_huffman_luma_dc[16];
    fjpeg_huffman_table_t fjpeg_huffman_luma_ac[256];
    fjpeg_huffman_table_t fjpeg_huffman_chroma_dc[16];
    fjpeg_huffman_table_t fjpeg_huffman_chroma_ac[256];

    fjpeg_short_huffman_table_t fjpeg_short_huffman_chroma_dc;
    fjpeg_short_huffman_table_t fjpeg_short_huffman_chroma_ac;
    
    fjpeg_short_huffman_table_t fjpeg_short_huffman_luma_dc;
    fjpeg_short_huffman_table_t fjpeg_short_huffman_luma_ac;

    fjpeg_pixel_t* fjpeg_y;
    fjpeg_pixel_t* fjpeg_cb;
    fjpeg_pixel_t* fjpeg_cr;

    fjpeg_coeff_t* fjpeg_ydct;
    fjpeg_coeff_t* fjpeg_cbdct;
    fjpeg_coeff_t* fjpeg_crdct;

    float precalc_cos[8][8];

    void fjpeg_precals_cos() {
        for (int i = 0; i < 8; i++) {
            for(int j = 0; j < 8; j++) {
                precalc_cos[i][j] = cosf((2.f * i + 1.f) * j * M_PI / 16.f);
            }
        }
        
    }

    fjpeg_context() {
        input = nullptr;
        output = nullptr;
        width = 0;
        height = 0;
        quality = 0;
        channels = 3;
        memset(fjpeg_luminance_quantization_table, 0, 64);
        memset(fjpeg_chrominance_quantization_table, 0, 64);
        memset(fjpeg_huffman_luma_dc, 0, 16 * sizeof(fjpeg_huffman_table_t));
        memset(fjpeg_huffman_luma_ac, 0, 256 * sizeof(fjpeg_huffman_table_t));
        memset(fjpeg_huffman_chroma_dc, 0, 16 * sizeof(fjpeg_huffman_table_t));
        memset(fjpeg_huffman_chroma_ac, 0, 256 * sizeof(fjpeg_huffman_table_t));
        fjpeg_y = nullptr;
        fjpeg_cb = nullptr;
        fjpeg_cr = nullptr;
        fjpeg_ydct = nullptr;
        fjpeg_cbdct = nullptr;
        fjpeg_crdct = nullptr;

        memcpy(fjpeg_luminance_quantization_table, fjpeg_default_luma_quant_table, 64);
        memcpy(fjpeg_chrominance_quantization_table, fjpeg_default_chroma_quant_table, 64);

        fjpeg_generate_tables(fjpeg_huffman_luma_dc, &fjpeg_default_huffman_luma_dc);
        #ifdef FJPEG_DEBUG_HUFFMAN
        for(int i = 0; i < 256; i++) {
            if(fjpeg_huffman_luma_dc[i].len == 0) continue;
            printf("%i\t%i\t", i, fjpeg_huffman_luma_dc[i].len);
            for(int ii = fjpeg_huffman_luma_dc[i].len-1; ii >= 0; ii--) {
                printf("%i",(fjpeg_huffman_luma_dc[i].code>>ii)&1);
            }
            printf("\r\n");
        }
        #endif

        fjpeg_generate_tables(fjpeg_huffman_luma_ac, &fjpeg_default_huffman_luma_ac);
        fjpeg_generate_tables(fjpeg_huffman_chroma_dc, &fjpeg_default_huffman_chroma_dc);
        fjpeg_generate_tables(fjpeg_huffman_chroma_ac, &fjpeg_default_huffman_chroma_ac);

        memcpy(&fjpeg_short_huffman_chroma_dc, &fjpeg_default_huffman_chroma_dc, sizeof(fjpeg_short_huffman_table_t));
        memcpy(&fjpeg_short_huffman_chroma_ac, &fjpeg_default_huffman_chroma_ac, sizeof(fjpeg_short_huffman_table_t));
        memcpy(&fjpeg_short_huffman_luma_dc, &fjpeg_default_huffman_luma_dc, sizeof(fjpeg_short_huffman_table_t));
        memcpy(&fjpeg_short_huffman_luma_ac, &fjpeg_default_huffman_luma_ac, sizeof(fjpeg_short_huffman_table_t));

        fjpeg_precals_cos();
    }

    bool setQuality(int quality) {
        if (quality < 1 || quality > 100) {
            return false;
        }

        this->quality = quality;

        for (int i = 0; i < 64; i++) {
            int scaled_value_y = ((int32_t)fjpeg_default_luma_quant_table[i] * (100-quality) + FJPEG_Q_FACTOR_SCALE / 2) / FJPEG_Q_FACTOR_SCALE;
            int scaled_value_u = ((int32_t)fjpeg_default_chroma_quant_table[i] * (100-quality) + FJPEG_Q_FACTOR_SCALE / 2) / FJPEG_Q_FACTOR_SCALE;
            
            fjpeg_luminance_quantization_table[i] = FJPEG_CLAMP(scaled_value_y, 1, 255);
            fjpeg_chrominance_quantization_table[i] = FJPEG_CLAMP(scaled_value_u, 1, 255);
        }

        return true;
    }

    bool readInput(const char* filename, int width, int height) {
        input = fopen(filename, "rb");
        if (!input) {
            return false;
        }

        this->width = width;
        this->height = height;

        fjpeg_y = (fjpeg_pixel_t*)malloc(width * height * sizeof(fjpeg_pixel_t));
        fjpeg_cb = (fjpeg_pixel_t*)malloc(width * height * sizeof(fjpeg_pixel_t));
        fjpeg_cr = (fjpeg_pixel_t*)malloc(width * height * sizeof(fjpeg_pixel_t));

        fjpeg_ydct = (fjpeg_coeff_t*)malloc(width * height * sizeof(fjpeg_coeff_t));
        fjpeg_cbdct = (fjpeg_coeff_t*)malloc(width * height * sizeof(fjpeg_coeff_t));
        fjpeg_crdct = (fjpeg_coeff_t*)malloc(width * height * sizeof(fjpeg_coeff_t));

        fread(fjpeg_y, 1, width * height, input);
        fread(fjpeg_cb, 1, (width * height) >> 2, input);
        fread(fjpeg_cr, 1, (width * height) >> 2, input);

        return true;
    }

    ~fjpeg_context() {
        if (input) {
            fclose(input);
        }

        if (output) {
            fclose(output);
        }

        if (fjpeg_y) {
            free(fjpeg_y);
        }

        if (fjpeg_cb) {
            free(fjpeg_cb);
        }

        if (fjpeg_cr) {
            free(fjpeg_cr);
        }

        if (fjpeg_ydct) {
            free(fjpeg_ydct);
        }

        if (fjpeg_cbdct) {
            free(fjpeg_cbdct);
        }

        if (fjpeg_crdct) {
            free(fjpeg_crdct);
        }
    }
 
};


static const char *fjpeg_version() {
    return FJPEG_VERSION;
}

void fjpeg_print_usage();
bool fjpeg_generate_header(fjpeg_bitstream* stream, fjpeg_context* context);
