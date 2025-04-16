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
#include <cassert>
#include <vector>

// Bitstream handling
class fjpeg_bitstream {
    public:

    uint32_t current;
    int offset;
    FILE *fp;
    std::vector<uint8_t> buffer;
    bool avoidFF;

    fjpeg_bitstream(FILE *fp) : current(0), offset(0), fp(fp), avoidFF(false) {}

    uint32_t readBits(int bits) {
        if(bits <= 0) return 0;
        assert(bits <= 32);
        uint32_t result = 0;
        while (bits > 0) {
            if (offset == 0) {
                current = fgetc(fp);
                if (current == EOF) {
                    fprintf(stderr, "Error: Unexpected end of file\n");
                    exit(1);
                }
                offset = 8;
            }
            int shift = bits - offset;
            if (shift > 0) {
                result |= current << shift;
                bits -= offset;
                offset = 0;
            } else {
                result |= current >> -shift;
                offset -= bits;
                current &= (1 << offset) - 1;
                bits = 0;
            }
        }
        return result;
    }

    std::vector<uint8_t> readBytes(int count) {
        std::vector<uint8_t> result(count);
        if(fread(result.data(), 1, count, fp) != count) {
            fprintf(stderr, "Error: Unexpected end of file\n");
            exit(1);
        }        
        return result;
    }

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
};
