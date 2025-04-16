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

    for(int i = 0; i < 256; i++) {
      if(data->val[i] == 0xFF) break;
      table_len++;
    }

    memset(output_table, 0, table_len * sizeof(fjpeg_huffman_table_t));
    
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
    }

    return table_len;
}

int fjpeg_entropy_stats(fjpeg_huffman_statistics_t* stat, fjpeg_context* context, fjpeg_coeff_t* block, int channel, int last_dc)
{
  int last_dc_coeff = last_dc;  // For DC coefficient prediction

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
    if(channel==0) stat->luma_dc[size]++;
    else stat->chroma_dc[size]++;

    if(last_coeff == 0) {
        // EOB
        if(channel==0) stat->luma_ac[0x00]++;
        else stat->chroma_ac[0x00]++;
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
                if(channel==0) stat->luma_ac[0xF0]++;
                else stat->chroma_ac[0xF0]++;
                run_length = 0;
            }
        }
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
        if(channel==0) stat->luma_ac[(run_length << 4) + size]++;
        else stat->chroma_ac[(run_length << 4) + size]++;
    }
    // Don't write EOB if we are at the end of the block already
    if(i < FJPEG_BLOCK_SIZE*FJPEG_BLOCK_SIZE) {
        // EOB
        if(channel==0) stat->luma_ac[0x00]++;
        else stat->chroma_ac[0x00]++;
    }
    return last_dc_coeff;
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

fjpeg_short_huffman_table_t fjpeg_generate_huffman_from_stats(fjpeg_huffman_table_t* huff_table, uint32_t* freq_array, int size) {
    // Generate Huffman tree from statistics
    
    fjpeg_short_huffman_table_t huff_short;
    memset(&huff_short, 0, sizeof(fjpeg_short_huffman_table_t));
    uint32_t freq[256] = {0};
    int codesize[256] = {0};
    int32_t others[256];
    memset(others, -1, 256 * sizeof(int32_t));

    bool done = false;

    uint32_t current_least = 0xffffffff;
    int current_V1 = 0;
    int current_V2 = 0;

    memcpy(freq, freq_array, size * sizeof(int));

    for(int i = 0; i < size; i++) {
        printf("DC %d freq %d\r\n", i, freq[i]);
        freq[i] += 1;
    }

    // Find code sizes for all symbols
    while(!done) {
        // Find two smallest frequencies
        bool found = false;
        for(int i = 0; i < size; i++) {
            if(freq[i] != 0 && freq[i] < current_least) {
                current_least = freq[i];
                current_V1 = i;
                found = true;
            }
        }
        if(!found) {
            done = true;
            break;
        }

        found = false;
        current_least = 0xffffffff;
        for(int i = 0; i < size; i++) {
            if(freq[i] != 0 && freq[i] < current_least && i != current_V1) {
                current_least = freq[i];
                current_V2 = i;
                found = true;
            }
        }
        if(!found) {
            done = true;
            break;
        }

        // Combine the two smallest frequencies
        freq[current_V1] += freq[current_V2];
        freq[current_V2] = 0;

        // Increase code size for all symbols in the group

        do {
          codesize[current_V1]++;
          if(others[current_V1] == -1) {
            others[current_V1] = current_V2;
            break;
          }  
          current_V1 = others[current_V1];
        } while(current_V1 != -1);

        do {
          codesize[current_V2]++;
          current_V2 = others[current_V2];
        } while(current_V2 != -1);

        // Reset
        current_least = 0xffffffff;
    }

    // Count the number of codes for each code size
    for(int i = 0; i < size; i++) {
        if(codesize[i] != 0) {
            huff_short.bits[codesize[i]-1]++;
        }
    }
    int sum[17] = {0};
    for(int i = 1; i < 17; i++) {
        sum[i] = (sum[i-1] + huff_short.bits[i-1]);
        printf("%d ",huff_short.bits[i-1]);
    }
    printf("\r\n");

    // Print sums
    for(int i = 0; i < 17; i++) {
        printf("%d ",sum[i]);
    }
    printf("\r\n");

    // Print codesizes
    for(int i = 0; i < size; i++) {
        printf("%d ",codesize[i]);
    }
    printf("\r\n");

    // Leftover values get assigned 16-bit code
    huff_short.bits[15] = size-sum[16];
    
    // Assign values to codes
    int next_code[17] = {0};
    for(int i = 0; i < size; i++) {
        if(codesize[i] != 0) {
            huff_short.val[sum[codesize[i]-1]+next_code[codesize[i]]] = i;
            next_code[codesize[i]]++;
        } else {
            huff_short.val[sum[16]+next_code[16]] = i;
            next_code[16]++;
        }
    }

    // Print hex of all values
    for(int i = 0; i < size; i++) {
        printf("0x%02X ", huff_short.val[i]);
    }
    printf("\r\n");

    huff_short.val[size] = 0xFF; // End of table

    fjpeg_generate_tables(huff_table, &huff_short);

    // Assign codes
    for(int i = 0; i < size; i++) {
        printf("DC %d ", i);
        for(int ii = 0; ii < huff_table[i].len; ii++) printf("%d",(huff_table[i].code>>(huff_table[i].len-ii-1))&1);
        printf("\r\n");
    }
    
    return huff_short;
}