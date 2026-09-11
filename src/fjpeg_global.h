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

// Round to nearest integer (half away from zero). A plain (int)(x + 0.5f)
// truncates towards zero for negative values and biases the coefficients.
// FJPEG_QUANT_DEADZONE widens the zero bin; values below it quantize to zero,
// which trades a little distortion for fewer transmitted coefficients.
#ifndef FJPEG_QUANT_DEADZONE
#define FJPEG_QUANT_DEADZONE 0.65f
#endif

static inline int fjpeg_round(float x) {
    float ax = (x < 0.0f) ? -x : x;
    if (ax < FJPEG_QUANT_DEADZONE) {
        return 0;
    }
    return (x >= 0.0f) ? (int)(x + 0.5f) : (int)(x - 0.5f);
}

#define FJPEG_UINT32_MAX 0xFFFFFFFF
#define FJPEG_BLOCK_SIZE 8


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
  16, 11, 10, 16, 24, 40, 51, 61,
  12, 12, 14, 19, 26, 58, 60, 55,
  14, 13, 16, 24, 40, 57, 69, 56,
  14, 17, 22, 29, 51, 87, 80, 62,
  18, 22, 37, 56, 68, 109, 103, 77,
  24, 35, 55, 64, 81, 104, 113, 92,
  49, 64, 78, 87, 103, 121, 120, 101,
  72, 92, 95, 98, 112, 100, 103, 99,
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
