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

#define _CRT_SECURE_NO_WARNINGS

#define FJPEG_VERSION "0.1.0"

#define FJPEG_MAX(a, b) ((a) > (b) ? (a) : (b))
#define FJPEG_MIN(a, b) ((a) < (b) ? (a) : (b))

#define FJPEG_SWAP(a, b) do { \
    typeof(a) _tmp = (a); \
    (a) = (b); \
    (b) = _tmp; \
} while (0)

#define FJPEG_CLAMP(x, min, max) FJPEG_MIN(FJPEG_MAX((x), (min)), (max))

typedef uint8_t fjpeg_pixel_t;
typedef float fjpeg_coeff_t;

#define FJPEG_UINT32_MAX 0xFFFFFFFF
#define FJPEG_BLOCK_SIZE 8

typedef struct {
    uint8_t len;
    uint16_t code;
} fjpeg_huffman_table_t;

#define FJPEG_Q_FACTOR_SCALE 50
// Set default quant
const uint8_t fjpeg_default_luma_quant_table[64] = { 
  16, 11, 10, 16, 124, 140, 151, 161,
  12, 12, 14, 19, 126, 158, 160, 155,
  14, 13, 16, 24, 140, 157, 169, 156,
  14, 17, 22, 29, 151, 187, 180, 162,
  18, 22, 37, 56, 168, 109, 103, 177,
  24, 35, 55, 64, 181, 104, 113, 192,
  49, 64, 78, 87, 103, 121, 120, 101,
  72, 92, 95, 98, 112, 100, 103, 199,
};

const uint8_t fjpeg_default_chroma_quant_table[64] = {
  17, 18, 24, 47, 99, 99, 99, 99,
  18, 21, 26, 66, 99, 99, 99, 99,
  24, 26, 56, 99, 99, 99, 99, 99,
  47, 66, 99, 99, 99, 99, 99, 99,
  99, 99, 99, 99, 99, 99, 99, 99,
  99, 99, 99, 99, 99, 99, 99, 99,
  99, 99, 99, 99, 99, 99, 99, 99,
  99, 99, 99, 99, 99, 99, 99, 99,
};

const uint8_t fjpeg_zigzag_8x8[64] = {
  // 8x8 zigzag
   0,   1,  5,  6, 14, 15, 27, 28,
   2,   4,  7, 13, 16, 26, 29, 42,
   3,   8, 12, 17, 25, 30, 41, 43,
   9,  11, 18, 24, 31, 40, 44, 53,
   10, 19, 23, 32, 39, 45, 52, 54,
   20, 22, 33, 38, 46, 51, 55, 60,
   21, 34, 37, 47, 50, 56, 59, 61,
   35, 36, 48, 49, 57, 58, 62, 63
};

typedef struct {
    uint8_t bits[16]; // BITS
    uint8_t val[163]; // HUFFVAL
} fjpeg_short_huffman_table_t;

const fjpeg_short_huffman_table_t fjpeg_default_huffman_luma_dc = { 
    { 0x00, 0x01, 0x05, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 },
    { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0xFF }
};

const fjpeg_short_huffman_table_t fjpeg_default_huffman_chroma_dc = { 
  { 0x00, 0x03, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00 },
  { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0xFF }
};

const fjpeg_short_huffman_table_t fjpeg_default_huffman_luma_ac = {
  { 0x00, 0x02, 0x01, 0x03, 0x03, 0x02, 0x04, 0x03, 0x05, 0x05, 0x04, 0x04, 0x00, 0x00, 0x01, 0x7D },
  { 0x01, 0x02, 0x03, 0x00, 0x04, 0x11, 0x05, 0x12, 0x21, 0x31, 0x41, 0x06, 0x13, 0x51, 0x61, 0x07,
    0x22, 0x71, 0x14, 0x32, 0x81, 0x91, 0xA1, 0x08, 0x23, 0x42, 0xB1, 0xC1, 0x15, 0x52, 0xD1, 0xF0,
    0x24, 0x33, 0x62, 0x72, 0x82, 0x09, 0x0A, 0x16, 0x17, 0x18, 0x19, 0x1A, 0x25, 0x26, 0x27, 0x28,
    0x29, 0x2A, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3A, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49,
    0x4A, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5A, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69,
    0x6A, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7A, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89,
    0x8A, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9A, 0xA2, 0xA3, 0xA4, 0xA5, 0xA6, 0xA7,
    0xA8, 0xA9, 0xAA, 0xB2, 0xB3, 0xB4, 0xB5, 0xB6, 0xB7, 0xB8, 0xB9, 0xBA, 0xC2, 0xC3, 0xC4, 0xC5,
    0xC6, 0xC7, 0xC8, 0xC9, 0xCA, 0xD2, 0xD3, 0xD4, 0xD5, 0xD6, 0xD7, 0xD8, 0xD9, 0xDA, 0xE1, 0xE2,
    0xE3, 0xE4, 0xE5, 0xE6, 0xE7, 0xE8, 0xE9, 0xEA, 0xF1, 0xF2, 0xF3, 0xF4, 0xF5, 0xF6, 0xF7, 0xF8,
    0xF9, 0xFA, 0xFF}
};

const fjpeg_short_huffman_table_t fjpeg_default_huffman_chroma_ac = {
 {0x00, 0x02, 0x01, 0x02, 0x04, 0x04, 0x03, 0x04, 0x07, 0x05, 0x04, 0x04, 0x00, 0x01, 0x02, 0x77},
 {0x00, 0x01, 0x02, 0x03, 0x11, 0x04, 0x05, 0x21, 0x31, 0x06, 0x12, 0x41, 0x51, 0x07, 0x61, 0x71,
  0x13, 0x22, 0x32, 0x81, 0x08, 0x14, 0x42, 0x91, 0xA1, 0xB1, 0xC1, 0x09, 0x23, 0x33, 0x52, 0xF0,
  0x15, 0x62, 0x72, 0xD1, 0x0A, 0x16, 0x24, 0x34, 0xE1, 0x25, 0xF1, 0x17, 0x18, 0x19, 0x1A, 0x26,
  0x27, 0x28, 0x29, 0x2A, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3A, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48,
  0x49, 0x4A, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5A, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68,
  0x69, 0x6A, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7A, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
  0x88, 0x89, 0x8A, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9A, 0xA2, 0xA3, 0xA4, 0xA5,
  0xA6, 0xA7, 0xA8, 0xA9, 0xAA, 0xB2, 0xB3, 0xB4, 0xB5, 0xB6, 0xB7, 0xB8, 0xB9, 0xBA, 0xC2, 0xC3,
  0xC4, 0xC5, 0xC6, 0xC7, 0xC8, 0xC9, 0xCA, 0xD2, 0xD3, 0xD4, 0xD5, 0xD6, 0xD7, 0xD8, 0xD9, 0xDA,
  0xE2, 0xE3, 0xE4, 0xE5, 0xE6, 0xE7, 0xE8, 0xE9, 0xEA, 0xF2, 0xF3, 0xF4, 0xF5, 0xF6, 0xF7, 0xF8,
  0xF9, 0xFA, 0xFF}
};

static uint8_t fjpeg_generate_tables(fjpeg_huffman_table_t* output_table, const fjpeg_short_huffman_table_t* data) {
    // Size table
    int k = 0;
    int j = 1;
    int i = 1;
    int lastk = 0;
    int table_len = 0;
    fjpeg_huffman_table_t table[256];
    memset(table, 0, 256 * sizeof(fjpeg_huffman_table_t));
    memset(output_table, 0, 256 * sizeof(fjpeg_huffman_table_t));
    
    while(i <= 16) {
        while(j <= data->bits[i-1]) {
            table[k].len = i;
            j++;
            k++;
        }
        i++;
        j = 1;
    }
    table[k].len = 0;

    lastk = k;

    // Code table
    bool done = false;
    int code = 0;
    int si = table[0].len;
    k = 0;
    while(!done) {
        while(table[k].len == si) {
            table[k].code = code;
            code++;
            k++;
        }
        if(table[k].len == 0) {
            done = true;
        } else {
          while(table[k].len != si) {
            code <<= 1;
            si++;
          }
        }
    }

    // Order tables
    for(int i = 0; i < 256; i++) {
      if(data->val[i] == 0xFF) break;
      output_table[data->val[i]] = table[i];
      table_len++;
    }

    return table_len;
}

typedef struct fjpeg_context {
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
 
} fjpeg_context_t;


// Bitstream handling
typedef struct fjpeg_bitstream {
    uint32_t current;
    int offset;
    FILE *fp;
    std::vector<uint8_t> buffer;
    bool avoidFF;

    fjpeg_bitstream(FILE *fp) : current(0), offset(0), fp(fp), avoidFF(false) {}

    void writeBits(uint32_t input, int bits) {
        assert(bits > 0 && bits <= 24);

        // Is there space in the current buffer
        if (offset + bits >= 32) {
            // Flush the buffer
            while(offset >= 8) {
                uint8_t val = (current >> (offset-8)) & 0xff;
                buffer.push_back(val);
                if(avoidFF && val == 0xff) {
                    buffer.push_back(0);
                }
                offset -= 8;
            }
            current &= (FJPEG_UINT32_MAX >> (32-offset));
        }
        // Write the remaining bits
        current <<= bits;
        current |= input;
        offset += bits;
    }

    void flushToFile() {
        if (offset > 0) {
            if(offset&7) current <<= (8-(offset&7));
            while(offset >= 8) {
                uint8_t val = (current >> (offset-8)) & 0xff;
                buffer.push_back(val);
                if(avoidFF && val == 0xff) {
                    buffer.push_back(0);
                }
                offset -= 8;
            }
            current &= ((~0) >> (32-offset));
        }
        fwrite(buffer.data(), 1, buffer.size(), fp);
        buffer.clear();
    }

    ~fjpeg_bitstream() {
        flushToFile();
    }
} fjpeg_bitstream_t;


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
