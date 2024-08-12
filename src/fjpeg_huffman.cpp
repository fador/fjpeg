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
#include "fjpeg_bitstream.h"

uint8_t fjpeg_generate_tables(fjpeg_huffman_table_t* output_table, const fjpeg_short_huffman_table_t* data) {
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

// Function to encode a single block of quantized DCT coefficients
int fjpeg_entropy_encode_block(fjpeg_bitstream* stream, fjpeg_context* context, fjpeg_coeff_t* block, int channel, int last_dc) {
    int last_dc_coeff = last_dc;  // For DC coefficient prediction

    fjpeg_huffman_table_t* huff_dc = channel==0?context->fjpeg_huffman_luma_dc:context->fjpeg_huffman_chroma_dc;
    fjpeg_huffman_table_t* huff_ac = channel==0?context->fjpeg_huffman_luma_ac:context->fjpeg_huffman_chroma_ac;

    // Check for last coeff
    int last_coeff = 0;
    for(int i = FJPEG_BLOCK_SIZE*FJPEG_BLOCK_SIZE-1; i >= 0; i--) {
        if((int)(block[i]+0.5f) != 0) {
            last_coeff = i;
            break;
        }        
    }

    // Code DC coefficient
    int coeff = (int)(block[0]+0.5f); // Quantized DCT coefficients
    #ifdef FJPEG_DEBUG_COEFF
    printf("Coeff %d\r\n", coeff);
    #endif
    int diff = coeff - last_dc_coeff;
    last_dc_coeff = coeff;
    int orig_diff = diff;
    int sign = (diff < 0) ? 1 : 0;
    int size = 0;
    // VLI encoding for DC coefficients
    if (sign) diff = -diff;
    while (diff != 0) {
        diff >>= 1;
        size++;
    }
    // Check for overflow
    if (size > 11) {
        // Handle error or clamp the size
        fprintf(stderr, "Error: DC coefficient size overflow 1\n");
        exit(1);
    }
    #ifdef FJPEG_DEBUG_COEFF
    printf("Writing DC coeff %d size %d huff len %d %d\n", orig_diff, size, huff_dc[size].len, huff_dc[size].code);
    #endif
    stream->writeBits(huff_dc[size].code, huff_dc[size].len); // Write size code
    if(size != 0) {
        diff = orig_diff;
        if(sign) {
            diff = (1 << size) + diff - 1;
        }
        #ifdef FJPEG_DEBUG_COEFF
        printf("Writing DC coeff %d size %d\n", diff, size);
        #endif
        stream->writeBits(diff, size); // Write diff value
    }

    if(last_coeff == 0) {
        // All zero block
        stream->writeBits(huff_ac[0x00].code, huff_ac[0x00].len); // EOB
        return last_dc_coeff;
    }
    int i;
    for (i = 1; i < FJPEG_BLOCK_SIZE*FJPEG_BLOCK_SIZE; i++) {
        
        int run_length = 0;
        coeff = (int)(block[i]+0.5f);

        // Run-length coding for AC coefficients
        while (coeff == 0 && i < FJPEG_BLOCK_SIZE*FJPEG_BLOCK_SIZE - 1) {
            run_length++;
            i++;
            coeff = (int)(block[i]+0.5f);
            if(i > last_coeff) {
                break;
            }
            if(run_length == 16) {
                // ZRL
                #ifdef FJPEG_DEBUG_COEFF
                printf("ZRL\r\n");
                #endif
                stream->writeBits(huff_ac[0xF0].code, huff_ac[0xF0].len); // ZRL
                run_length = 0;
            }
        }
        #ifdef FJPEG_DEBUG_COEFF
        printf("Run length %d\r\n", run_length);
        #endif
        // Check for EOB after encoding each coefficient
        if (i > last_coeff || coeff == 0) {
            break;
        }

         // Encode the AC coefficient
        int sign = (coeff < 0) ? 1 : 0;
        int size = 0;
        // VLI encoding for AC coefficients
        int orig_coeff = coeff;
        if(sign) coeff = -coeff;
        while (coeff != 0) {
            coeff >>= 1;
            size++;                
        }
            // Check for overflow
        if (size > 10) {
            // Handle error or clamp the size
            fprintf(stderr, "Error: DC coefficient size overflow 2\n");
            exit(1);
        }
        #ifdef FJPEG_DEBUG_COEFF
        printf("Writing AC coeff %d size %d huff len %d %d\n", orig_coeff, size, huff_ac[(run_length << 4) + size].len, huff_ac[(run_length << 4) + size].code);
        #endif
        stream->writeBits(huff_ac[(run_length << 4) + size].code,
                            huff_ac[(run_length << 4) + size].len); // Write run-length/size code
        if (size > 0) {
            coeff = orig_coeff;
            if(sign) {
                coeff = (1 << size) + coeff - 1;
            }
            stream->writeBits(coeff, size); // Write the remaining bits
        }
    }
    // Don't write EOB if we are at the end of the block already
    if(i < FJPEG_BLOCK_SIZE*FJPEG_BLOCK_SIZE) {
        // EOB
        #ifdef FJPEG_DEBUG_COEFF
        printf("EOB\r\n");
        #endif
        stream->writeBits(huff_ac[0x00].code, huff_ac[0x00].len); // EOB
    }

    return last_dc_coeff;
}