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
#include <vector>
#include "fjpeg.h"
#include "fjpeg_bitstream.h"

// Table D.2 probability estimation table size
#define FJPEG_ARITH_DC_STAT_BINS 64
#define FJPEG_ARITH_AC_STAT_BINS 256

extern const uint32_t fjpeg_aritab[114];

class fjpeg_arith_encoder {
public:
    fjpeg_arith_encoder() : stream(nullptr), c(0), a(0x10000), sc(0), zc(0), ct(11), buffer(-1) {}
    explicit fjpeg_arith_encoder(fjpeg_bitstream* bs) : stream(bs), c(0), a(0x10000), sc(0), zc(0), ct(11), buffer(-1) {}

    void init(fjpeg_bitstream* bs);
    void encode(unsigned char* st, int val);
    void flush();

private:
    void emit_byte(int val);

    fjpeg_bitstream* stream;
    uint32_t c;
    uint32_t a;
    int sc;
    int zc;
    int ct;
    int buffer;
};

// Generate ITU-T T.81 SOF9 Sequential DCT JPEG with Arithmetic Coding
bool fjpeg_generate_arithmetic(fjpeg_bitstream* stream, fjpeg_context* context);
