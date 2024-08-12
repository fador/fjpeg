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
#include <cctype>
#include <vector>
#include <cassert>
#include <string>
#include <chrono>
#ifdef _WIN32
#include <winsock2.h>
#pragma comment(lib, "ws2_32.lib")
#endif

#include "fjpeg.h"
#include "fjpeg_bitstream.h"
#include "fjpeg_transquant.h"

//#include "fjpeg_tables.h"

//#define FJPEG_DEBUG_BLOCK 1
//#define FJPEG_DEBUG_COEFF 1



// Generate jpeg header
bool fjpeg_generate_header(fjpeg_bitstream* stream, fjpeg_context* context) {
    
    uint8_t tmp[64];
    // SOI
    stream->writeBits(0xFFD8, 16);

    // APP0
    stream->writeBits(0xFFE0, 16);
    stream->writeBits(16, 16);
    stream->writeBits(0x4A4649, 24); // "JFI"
    stream->writeBits(0x4600, 16); // "F\0"
    stream->writeBits(0x0102, 16); // Version
    stream->writeBits(0x00, 8); // Units
    stream->writeBits(0x0001, 16); // X density
    stream->writeBits(0x0001, 16); // Y density
    stream->writeBits(0x00, 8); // X thumbnail
    stream->writeBits(0x00, 8); // Y thumbnail
    // Thumbnail data
    
    // DQT
    
    stream->writeBits(0xFFDB, 16);
    stream->writeBits(67, 16);
    stream->writeBits(0, 4);
    stream->writeBits(0, 4);    
    
    for(int i = 0; i < 64; i++) {
        tmp[fjpeg_zigzag_8x8[i]] = context->fjpeg_luminance_quantization_table[i];
    }
    for (int i = 0; i < 64; i++) {
        stream->writeBits(tmp[i], 8);
    }

    if(context->channels > 1) {
        // chroma quant
        stream->writeBits(0xFFDB, 16);
        stream->writeBits(67, 16);
        stream->writeBits(0, 4);
        stream->writeBits(1, 4);

        for(int i = 0; i < 64; i++) {
            tmp[fjpeg_zigzag_8x8[i]] = context->fjpeg_chrominance_quantization_table[i];
        }
        for (int i = 0; i < 64; i++) {
            stream->writeBits(tmp[i], 8);
        }
    }
    
   
    // SOF0
    stream->writeBits(0xFFC0, 16);
    stream->writeBits(context->channels==1?11:17, 16);
    stream->writeBits(8, 8); // 8 bits per sample
    stream->writeBits(context->height, 16);
    stream->writeBits(context->width, 16);
    stream->writeBits(context->channels, 8);

    for (int i = 0; i < context->channels; i++) {
        stream->writeBits(i + 1, 8);
        stream->writeBits(context->channels == 1?0x11:(i == 0 ? 0x22 : 0x11), 8); // Sampling factors 4:2:0
        stream->writeBits(i==0?0:1, 8); // Quant table
    }


    // DHT
    stream->writeBits(0xFFC4, 16); // Huffman tables
    stream->writeBits(31+179, 16); // Length
    stream->writeBits(0, 4); // DC
    stream->writeBits(0, 4); // Table ID

    for (int i = 0; i < 16; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_dc.bits[i], 8);
    }

    for (int i = 0; i < 12; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_dc.val[i], 8);
    }

    //stream->writeBits(0xFFC4, 16); // Huffman tables
    //stream->writeBits(181, 16); // Length
    stream->writeBits(1, 4); // AC
    stream->writeBits(0, 4); // Table ID

    for (int i = 0; i < 16; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_ac.bits[i], 8);
    }

    for (int i = 0; i < 162; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_ac.val[i], 8);
    }

    if(context->channels > 1) {
        stream->writeBits(0xFFC4, 16); // Huffman tables
        stream->writeBits(210, 16); // Length

        stream->writeBits(0, 4); // DC
        stream->writeBits(1, 4); // Table ID

        for (int i = 0; i < 16; i++) {
            stream->writeBits(context->fjpeg_short_huffman_chroma_dc.bits[i], 8);
        }

        for (int i = 0; i < 12; i++) {
            stream->writeBits(context->fjpeg_short_huffman_chroma_dc.val[i], 8);
        }

        stream->writeBits(1, 4); // AC
        stream->writeBits(1, 4); // Table ID

        for (int i = 0; i < 16; i++) {
            stream->writeBits(context->fjpeg_short_huffman_chroma_ac.bits[i], 8);
        }

        for (int i = 0; i < 162; i++) {
            stream->writeBits(context->fjpeg_short_huffman_chroma_ac.val[i], 8);
        }
    }

    // COM
    stream->writeBits(0xFFFE, 16);
    stream->writeBits((uint32_t)(strlen("FJPEG ")+strlen(FJPEG_VERSION)) + 2, 16);

    for(int i = 0; i < strlen("FJPEG "); i++) {
        stream->writeBits("FJPEG "[i], 8);
    }

    for(int i = 0; i < strlen(FJPEG_VERSION); i++) {
        stream->writeBits(FJPEG_VERSION[i], 8);
    }

    // SOS
    stream->writeBits(0xFFDA, 16);
    stream->writeBits(context->channels==1?8:12, 16); // Length
    stream->writeBits(context->channels, 8);

    for (int i = 0; i < context->channels; i++) {
        stream->writeBits(i + 1, 8); // Component ID
        stream->writeBits(i == 0 ? 0 : 0x11, 8); // Huffman table
    }

    stream->writeBits(0, 8); // DCT coeff start
    stream->writeBits(0x3F, 8); // DCT coeff end
    stream->writeBits(0, 8); // Successive Approximation

    // Entropy coded huffman data
    // Luma from context->fjpeg_ydct
    // Chroma from context->fjpeg_cbdct and context->fjpeg_crdct
    int last_dc_coeff[3] = {0, 0, 0};
    fjpeg_coeff_t dct_block[64];
    stream->avoidFF = true;
    const int inc_xy = context->channels==1?8:16;
    const int max_uv = context->channels==1?1:2;
    
    for(int y = 0; y < context->height; y+=inc_xy) {
        for(int x = 0; x < context->width; x+=inc_xy) {

            for(int v = 0; v < max_uv; v++) {
                for(int u = 0; u < max_uv; u++) {
                    #ifdef FJPEG_DEBUG_BLOCK
                    printf("Encoding block %dx%d + %dx%d\n", x, y, u*8, v*8);
                    #endif
                    fjpeg_extract_coeff_8x8(context, dct_block, x+u*8, y+v*8, 0);
                    // Print out the block
                    #ifdef FJPEG_DEBUG_BLOCK
                    for(int i = 0; i < 64; i++) {
                        printf("%3d ", (int16_t)(dct_block[i]+0.5f));
                        if((i+1)%8 == 0) printf("\r\n");
                    }
                    #endif                    
                    last_dc_coeff[0] = fjpeg_entropy_encode_block(stream, context, dct_block, 0, last_dc_coeff[0]);
                }
            }

            if(context->channels==3) {
                fjpeg_extract_coeff_8x8(context, dct_block, x>>1, y>>1, 1);
                last_dc_coeff[1] = fjpeg_entropy_encode_block(stream, context, dct_block, 1, last_dc_coeff[1]);

                fjpeg_extract_coeff_8x8(context, dct_block, x>>1, y>>1, 2);
                last_dc_coeff[2] = fjpeg_entropy_encode_block(stream, context, dct_block, 2, last_dc_coeff[2]);
            }
        }
    }
    stream->avoidFF = false;


    // EOI
    stream->writeBits(0xFFD9, 16);

    stream->flushToFile();

    return true;
}


void fjpeg_print_usage() {
    printf("Usage: fjpeg [options]\r\n");
    printf("Example: fjpeg -i input.yuv -q 50 -r 1280x720 -o output.jpg\r\n");
    printf("Options:\r\n");
    printf("  -i <input_filename>  input YUV file\r\n");
    printf("  -q <quality>  Set quality factor (1-100)\r\n");
    printf("  -r <width>x<height>  Set resolution\r\n");
    printf("  -o <output_filename>  Output JPEG file\r\n");
    printf("  -h  Show help\r\n");
}
