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
#include "fjpeg_global.h"
#include "fjpeg_bitstream.h"
#include "fjpeg_transquant.h"
#include "fjpeg_arith.h"

#define V(i, a, b, c, d) \
  (((uint32_t)(a) << 16) | ((uint32_t)(c) << 8) | ((uint32_t)(d) << 7) | (uint32_t)(b))

// ITU-T T.81 Table D.2 probability estimation state machine
const uint32_t fjpeg_aritab[114] = {
  V(   0, 0x5a1d,   1,   1, 1 ),
  V(   1, 0x2586,  14,   2, 0 ),
  V(   2, 0x1114,  16,   3, 0 ),
  V(   3, 0x080b,  18,   4, 0 ),
  V(   4, 0x03d8,  20,   5, 0 ),
  V(   5, 0x01da,  23,   6, 0 ),
  V(   6, 0x00e5,  25,   7, 0 ),
  V(   7, 0x006f,  28,   8, 0 ),
  V(   8, 0x0036,  30,   9, 0 ),
  V(   9, 0x001a,  33,  10, 0 ),
  V(  10, 0x000d,  35,  11, 0 ),
  V(  11, 0x0006,   9,  12, 0 ),
  V(  12, 0x0003,  10,  13, 0 ),
  V(  13, 0x0001,  12,  13, 0 ),
  V(  14, 0x5a7f,  15,  15, 1 ),
  V(  15, 0x3f25,  36,  16, 0 ),
  V(  16, 0x2cf2,  38,  17, 0 ),
  V(  17, 0x207c,  39,  18, 0 ),
  V(  18, 0x17b9,  40,  19, 0 ),
  V(  19, 0x1182,  42,  20, 0 ),
  V(  20, 0x0cef,  43,  21, 0 ),
  V(  21, 0x09a1,  45,  22, 0 ),
  V(  22, 0x072f,  46,  23, 0 ),
  V(  23, 0x055c,  48,  24, 0 ),
  V(  24, 0x0406,  49,  25, 0 ),
  V(  25, 0x0303,  51,  26, 0 ),
  V(  26, 0x0240,  52,  27, 0 ),
  V(  27, 0x01b1,  54,  28, 0 ),
  V(  28, 0x0144,  56,  29, 0 ),
  V(  29, 0x00f5,  57,  30, 0 ),
  V(  30, 0x00b7,  59,  31, 0 ),
  V(  31, 0x008a,  60,  32, 0 ),
  V(  32, 0x0068,  62,  33, 0 ),
  V(  33, 0x004e,  63,  34, 0 ),
  V(  34, 0x003b,  32,  35, 0 ),
  V(  35, 0x002c,  33,   9, 0 ),
  V(  36, 0x5ae1,  37,  37, 1 ),
  V(  37, 0x484c,  64,  38, 0 ),
  V(  38, 0x3a0d,  65,  39, 0 ),
  V(  39, 0x2ef1,  67,  40, 0 ),
  V(  40, 0x261f,  68,  41, 0 ),
  V(  41, 0x1f33,  69,  42, 0 ),
  V(  42, 0x19a8,  70,  43, 0 ),
  V(  43, 0x1518,  72,  44, 0 ),
  V(  44, 0x1177,  73,  45, 0 ),
  V(  45, 0x0e74,  74,  46, 0 ),
  V(  46, 0x0bfb,  75,  47, 0 ),
  V(  47, 0x09f8,  77,  48, 0 ),
  V(  48, 0x0861,  78,  49, 0 ),
  V(  49, 0x0706,  79,  50, 0 ),
  V(  50, 0x05cd,  48,  51, 0 ),
  V(  51, 0x04de,  50,  52, 0 ),
  V(  52, 0x040f,  50,  53, 0 ),
  V(  53, 0x0363,  51,  54, 0 ),
  V(  54, 0x02d4,  52,  55, 0 ),
  V(  55, 0x025c,  53,  56, 0 ),
  V(  56, 0x01f8,  54,  57, 0 ),
  V(  57, 0x01a4,  55,  58, 0 ),
  V(  58, 0x0160,  56,  59, 0 ),
  V(  59, 0x0125,  57,  60, 0 ),
  V(  60, 0x00f6,  58,  61, 0 ),
  V(  61, 0x00cb,  59,  62, 0 ),
  V(  62, 0x00ab,  61,  63, 0 ),
  V(  63, 0x008f,  61,  32, 0 ),
  V(  64, 0x5b12,  65,  65, 1 ),
  V(  65, 0x4d04,  80,  66, 0 ),
  V(  66, 0x412c,  81,  67, 0 ),
  V(  67, 0x37d8,  82,  68, 0 ),
  V(  68, 0x2fe8,  83,  69, 0 ),
  V(  69, 0x293c,  84,  70, 0 ),
  V(  70, 0x2379,  86,  71, 0 ),
  V(  71, 0x1edf,  87,  72, 0 ),
  V(  72, 0x1aa9,  87,  73, 0 ),
  V(  73, 0x174e,  72,  74, 0 ),
  V(  74, 0x1424,  72,  75, 0 ),
  V(  75, 0x119c,  74,  76, 0 ),
  V(  76, 0x0f6b,  74,  77, 0 ),
  V(  77, 0x0d51,  75,  78, 0 ),
  V(  78, 0x0bb6,  77,  79, 0 ),
  V(  79, 0x0a40,  77,  48, 0 ),
  V(  80, 0x5832,  80,  81, 1 ),
  V(  81, 0x4d1c,  88,  82, 0 ),
  V(  82, 0x438e,  89,  83, 0 ),
  V(  83, 0x3bdd,  90,  84, 0 ),
  V(  84, 0x34ee,  91,  85, 0 ),
  V(  85, 0x2eae,  92,  86, 0 ),
  V(  86, 0x299a,  93,  87, 0 ),
  V(  87, 0x2516,  86,  71, 0 ),
  V(  88, 0x5570,  88,  89, 1 ),
  V(  89, 0x4ca9,  95,  90, 0 ),
  V(  90, 0x44d9,  96,  91, 0 ),
  V(  91, 0x3e22,  97,  92, 0 ),
  V(  92, 0x3824,  99,  93, 0 ),
  V(  93, 0x32b4,  99,  94, 0 ),
  V(  94, 0x2e17,  93,  86, 0 ),
  V(  95, 0x56a8,  95,  96, 1 ),
  V(  96, 0x4f46, 101,  97, 0 ),
  V(  97, 0x47e5, 102,  98, 0 ),
  V(  98, 0x41cf, 103,  99, 0 ),
  V(  99, 0x3c3d, 104, 100, 0 ),
  V( 100, 0x375e,  99,  93, 0 ),
  V( 101, 0x5231, 105, 102, 0 ),
  V( 102, 0x4c0f, 106, 103, 0 ),
  V( 103, 0x4639, 107, 104, 0 ),
  V( 104, 0x415e, 103,  99, 0 ),
  V( 105, 0x5627, 105, 106, 1 ),
  V( 106, 0x50e7, 108, 107, 0 ),
  V( 107, 0x4b85, 109, 103, 0 ),
  V( 108, 0x5597, 110, 109, 0 ),
  V( 109, 0x504f, 111, 107, 0 ),
  V( 110, 0x5a10, 110, 111, 1 ),
  V( 111, 0x5522, 112, 109, 0 ),
  V( 112, 0x59eb, 112, 111, 1 ),
  V( 113, 0x5a1d, 113, 113, 0 )
};

void fjpeg_arith_encoder::init(fjpeg_bitstream* bs) {
    stream = bs;
    c = 0;
    a = 0x10000;
    sc = 0;
    zc = 0;
    ct = 11;
    buffer = -1;
}

void fjpeg_arith_encoder::emit_byte(int val) {
    stream->writeBits(val & 0xFF, 8);
    if ((val & 0xFF) == 0xFF) {
        stream->writeBits(0x00, 8);
    }
}

void fjpeg_arith_encoder::encode(unsigned char* st, int val) {
    int sv = *st;
    uint32_t qe = fjpeg_aritab[sv & 0x7F];
    uint8_t nl = (uint8_t)(qe & 0xFF); qe >>= 8;
    uint8_t nm = (uint8_t)(qe & 0xFF); qe >>= 8;

    a -= qe;
    if (val != (sv >> 7)) {
        if (a >= qe) {
            c += a;
            a = qe;
        }
        *st = (sv & 0x80) ^ nl;
    } else {
        if (a >= 0x8000) {
            return;
        }
        if (a < qe) {
            c += a;
            a = qe;
        }
        *st = (sv & 0x80) ^ nm;
    }

    do {
        a <<= 1;
        c <<= 1;
        if (--ct == 0) {
            uint32_t temp = c >> 19;
            if (temp > 0xFF) {
                if (buffer >= 0) {
                    while (zc > 0) {
                        emit_byte(0x00);
                        zc--;
                    }
                    emit_byte(buffer + 1);
                }
                zc += sc;
                sc = 0;
                buffer = (int)(temp & 0xFF);
            } else if (temp == 0xFF) {
                sc++;
            } else {
                if (buffer == 0) {
                    zc++;
                } else if (buffer >= 0) {
                    while (zc > 0) {
                        emit_byte(0x00);
                        zc--;
                    }
                    emit_byte(buffer);
                }
                if (sc > 0) {
                    while (zc > 0) {
                        emit_byte(0x00);
                        zc--;
                    }
                    while (sc > 0) {
                        emit_byte(0xFF);
                        sc--;
                    }
                }
                buffer = (int)(temp & 0xFF);
            }
            c &= 0x7FFFF;
            ct += 8;
        }
    } while (a < 0x8000);
}

void fjpeg_arith_encoder::flush() {
    uint32_t temp = (a - 1 + c) & 0xFFFF0000;
    if (temp < c) {
        c = temp + 0x8000;
    } else {
        c = temp;
    }
    c <<= ct;
    if (c & 0xF8000000) {
        if (buffer >= 0) {
            while (zc > 0) { emit_byte(0x00); zc--; }
            emit_byte(buffer + 1);
        }
        zc += sc;
        sc = 0;
    } else {
        if (buffer == 0) {
            zc++;
        } else if (buffer >= 0) {
            while (zc > 0) { emit_byte(0x00); zc--; }
            emit_byte(buffer);
        }
        if (sc > 0) {
            while (zc > 0) { emit_byte(0x00); zc--; }
            while (sc > 0) { emit_byte(0xFF); sc--; }
        }
    }
    if (c & 0x7FFF800) {
        while (zc > 0) { emit_byte(0x00); zc--; }
        emit_byte((c >> 19) & 0xFF);
        if (c & 0x7F800) {
            emit_byte((c >> 11) & 0xFF);
        }
    }
}

// ---------------------------------------------------------------------------
// Arithmetic Sequential DCT Encoding (SOF9)
// ---------------------------------------------------------------------------

bool fjpeg_generate_arithmetic(fjpeg_bitstream* stream, fjpeg_context* context) {
    if (!stream || !context) return false;

    // Optional trellis quantization
    fjpeg_trellis_optimize(context);

    // 1. SOI
    stream->writeBits(0xFFD8, 16);

    // 2. APP0 (JFIF)
    stream->writeBits(0xFFE0, 16);
    stream->writeBits(16, 16);
    stream->writeBits(0x4A4649, 24); // "JFI"
    stream->writeBits(0x4600, 16);   // "F\0"
    stream->writeBits(0x0102, 16);   // Version 1.02
    stream->writeBits(0x00, 8);      // Units
    stream->writeBits(0x0001, 16);   // X density
    stream->writeBits(0x0001, 16);   // Y density
    stream->writeBits(0x00, 8);      // X thumbnail
    stream->writeBits(0x00, 8);      // Y thumbnail

    // 3. DQT (Quantization tables)
    uint8_t tmp[64];
    stream->writeBits(0xFFDB, 16);
    stream->writeBits(67, 16);
    stream->writeBits(0, 4); // Pq = 0 (8-bit)
    stream->writeBits(0, 4); // Tq = 0
    for (int i = 0; i < 64; i++) {
        tmp[fjpeg_zigzag_8x8[i]] = context->fjpeg_luminance_quantization_table[i];
    }
    for (int i = 0; i < 64; i++) {
        stream->writeBits(tmp[i], 8);
    }

    if (context->channels > 1) {
        stream->writeBits(0xFFDB, 16);
        stream->writeBits(67, 16);
        stream->writeBits(0, 4); // Pq = 0
        stream->writeBits(1, 4); // Tq = 1
        for (int i = 0; i < 64; i++) {
            tmp[fjpeg_zigzag_8x8[i]] = context->fjpeg_chrominance_quantization_table[i];
        }
        for (int i = 0; i < 64; i++) {
            stream->writeBits(tmp[i], 8);
        }
    }

    // 4. SOF9 (Extended Sequential DCT, Arithmetic Coding, marker 0xFFC9)
    int ncomp = context->channels;
    stream->writeBits(0xFFC9, 16);
    stream->writeBits(ncomp == 1 ? 11 : 17, 16);
    stream->writeBits(context->bit_depth, 8); // Precision (8 or 12)
    stream->writeBits(context->height, 16);
    stream->writeBits(context->width, 16);
    stream->writeBits(ncomp, 8);

    for (int i = 0; i < ncomp; i++) {
        stream->writeBits(i + 1, 8); // Component ID
        uint8_t hv = (ncomp == 1) ? 0x11 : (i == 0 ? 0x22 : 0x11);
        stream->writeBits(hv, 8);
        stream->writeBits(i == 0 ? 0 : 1, 8); // Quant table ID
    }

    // 5. DAC (Define Arithmetic Conditioning, marker 0xFFCC)
    // Standard default conditioning parameters:
    // DC tables: L=0, U=1 -> byte 0x01
    // AC tables: K=5 -> byte 0x05
    int dac_tables = (ncomp == 1) ? 2 : 4;
    stream->writeBits(0xFFCC, 16);
    stream->writeBits(2 + dac_tables * 2, 16);

    // DC table 0
    stream->writeBits(0, 4); // Class 0 (DC)
    stream->writeBits(0, 4); // Table 0
    stream->writeBits(0x10, 8); // U=1, L=0

    // AC table 0
    stream->writeBits(1, 4); // Class 1 (AC)
    stream->writeBits(0, 4); // Table 0
    stream->writeBits(0x05, 8); // K=5

    if (ncomp > 1) {
        // DC table 1
        stream->writeBits(0, 4); // Class 0 (DC)
        stream->writeBits(1, 4); // Table 1
        stream->writeBits(0x10, 8); // U=1, L=0

        // AC table 1
        stream->writeBits(1, 4); // Class 1 (AC)
        stream->writeBits(1, 4); // Table 1
        stream->writeBits(0x05, 8); // K=5
    }

    // 6. SOS (Start of Scan)
    stream->writeBits(0xFFDA, 16);
    stream->writeBits(6 + 2 * ncomp, 16);
    stream->writeBits(ncomp, 8);

    for (int i = 0; i < ncomp; i++) {
        stream->writeBits(i + 1, 8);
        uint8_t td = (i == 0) ? 0 : 1;
        uint8_t ta = (i == 0) ? 0 : 1;
        stream->writeBits((td << 4) | ta, 8);
    }

    stream->writeBits(0, 8);  // Ss = 0
    stream->writeBits(63, 8); // Se = 63
    stream->writeBits(0, 8);  // Ah = 0, Al = 0

    // 7. Entropy-coded scan
    fjpeg_arith_encoder encoder(stream);

    uint8_t dc_stats[2][FJPEG_ARITH_DC_STAT_BINS];
    uint8_t ac_stats[2][FJPEG_ARITH_AC_STAT_BINS];
    uint8_t fixed_bin[4];

    memset(dc_stats, 0, sizeof(dc_stats));
    memset(ac_stats, 0, sizeof(ac_stats));
    memset(fixed_bin, 0, sizeof(fixed_bin));
    fixed_bin[0] = 113; // Fixed 0.5 probability bin in state table

    int last_dc_val[3] = {0, 0, 0};
    int dc_context[3] = {0, 0, 0};

    const int dc_L[2] = {0, 0};
    const int dc_U[2] = {1, 1};
    const int ac_K[2] = {5, 5};

    auto encode_block = [&](fjpeg_coeff_t* block, int ci) {
        int tbl = (ci == 0) ? 0 : 1;

        // --- DC Coefficient ---
        uint8_t* st = dc_stats[tbl] + dc_context[ci];
        int dc_val = fjpeg_round_dc(block[0]);
        int v = dc_val - last_dc_val[ci];

        if (v == 0) {
            encoder.encode(st, 0);
            dc_context[ci] = 0;
        } else {
            last_dc_val[ci] = dc_val;
            encoder.encode(st, 1);
            if (v > 0) {
                encoder.encode(st + 1, 0);
                st += 2;
                dc_context[ci] = 4;
            } else {
                v = -v;
                encoder.encode(st + 1, 1);
                st += 3;
                dc_context[ci] = 8;
            }

            int m = 0;
            if (v -= 1) {
                encoder.encode(st, 1);
                m = 1;
                int v2 = v;
                st = dc_stats[tbl] + 20;
                while (v2 >>= 1) {
                    encoder.encode(st, 1);
                    m <<= 1;
                    st += 1;
                }
            }
            encoder.encode(st, 0);

            if (m < ((1 << dc_L[tbl]) >> 1)) {
                dc_context[ci] = 0;
            } else if (m > ((1 << dc_U[tbl]) >> 1)) {
                dc_context[ci] += 8;
            }

            st += 14;
            while (m >>= 1) {
                encoder.encode(st, (m & v) ? 1 : 0);
            }
        }

        // --- AC Coefficients ---
        // Find last nonzero coefficient in zigzag order
        int ke = 63;
        while (ke > 0 && fjpeg_round(block[ke]) == 0) {
            ke--;
        }

        for (int k = 1; k <= ke; k++) {
            st = ac_stats[tbl] + 3 * (k - 1);
            encoder.encode(st, 0); // Not EOB

            int ac_val = fjpeg_round(block[k]);
            while (ac_val == 0) {
                encoder.encode(st + 1, 0);
                st += 3;
                k++;
                ac_val = fjpeg_round(block[k]);
            }
            encoder.encode(st + 1, 1);

            if (ac_val > 0) {
                encoder.encode(fixed_bin, 0);
            } else {
                ac_val = -ac_val;
                encoder.encode(fixed_bin, 1);
            }

            st += 2;
            int m = 0;
            if (ac_val -= 1) {
                encoder.encode(st, 1);
                m = 1;
                int v2 = ac_val;
                if (v2 >>= 1) {
                    encoder.encode(st, 1);
                    m <<= 1;
                    st = ac_stats[tbl] + (k <= ac_K[tbl] ? 189 : 217);
                    while (v2 >>= 1) {
                        encoder.encode(st, 1);
                        m <<= 1;
                        st += 1;
                    }
                }
            }
            encoder.encode(st, 0);

            st += 14;
            while (m >>= 1) {
                encoder.encode(st, (m & ac_val) ? 1 : 0);
            }
        }

        if (ke < 63) {
            st = ac_stats[tbl] + 3 * ke;
            encoder.encode(st, 1); // EOB
        }
    };

    fjpeg_coeff_t in_block[64];

    if (ncomp == 1) {
        for (int y = 0; y < context->padded_height; y += 8) {
            for (int x = 0; x < context->padded_width; x += 8) {
                fjpeg_extract_coeff_8x8(context, in_block, x, y, 0);
                encode_block(in_block, 0);
            }
        }
    } else {
        // 4:2:0 interleaved: 4 Y blocks, 1 Cb block, 1 Cr block per MCU
        for (int y = 0; y < context->padded_height; y += 16) {
            for (int x = 0; x < context->padded_width; x += 16) {
                for (int v = 0; v < 2; v++) {
                    for (int u = 0; u < 2; u++) {
                        fjpeg_extract_coeff_8x8(context, in_block, x + u * 8, y + v * 8, 0);
                        encode_block(in_block, 0);
                    }
                }
                fjpeg_extract_coeff_8x8(context, in_block, x / 2, y / 2, 1);
                encode_block(in_block, 1);
                fjpeg_extract_coeff_8x8(context, in_block, x / 2, y / 2, 2);
                encode_block(in_block, 2);
            }
        }
    }

    encoder.flush();

    // 8. EOI
    stream->writeBits(0xFFD9, 16);
    stream->flushToFile();

    return true;
}