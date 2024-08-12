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

#include <cstdint>

#define FJPEG_VERSION "0.1.0"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

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


#define FJPEG_Q_FACTOR_SCALE 50


typedef struct {
    uint8_t bits[16]; // BITS
    uint8_t val[163]; // HUFFVAL
} fjpeg_short_huffman_table_t;

typedef struct {
    uint8_t len;
    uint16_t code;
} fjpeg_huffman_table_t;


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

const uint8_t fjpeg_zigzag_reverse_8x8[64] = {
    0,   1, 15, 27, 29, 42, 24, 31, 
    5,  14, 28, 26,  3, 18, 40, 33, 
    6,   2, 16,  8, 11, 44, 22, 38, 
    4,  13, 12,  9, 53, 20, 46, 56, 
    7,  17, 43, 10, 54, 51, 50, 59, 
    25, 41, 19, 52, 55, 47, 61, 57, 
    30, 23, 45, 60, 37, 35, 49, 58, 
    32, 39, 21, 34, 36, 48, 62, 63,
};
