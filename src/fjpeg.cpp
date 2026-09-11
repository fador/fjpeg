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

bool fjpeg_read_headers(fjpeg_bitstream* stream, fjpeg_context* context) {
    // Read headers
    uint16_t marker = stream->readBits(16);
    if(marker != 0xFFD8) {
        fprintf(stderr, "Error: Not a JPEG file\n");
        return false;
    }

    uint8_t bits[16];
    fjpeg_short_huffman_table_t* huffman_table = nullptr;
    fjpeg_huffman_table_t* huffman_table_long = nullptr;
    uint32_t type;
    uint16_t version;
    uint8_t xthumb;
    uint8_t ythumb;
    uint8_t tmp;
    uint8_t table;
    std::vector<uint8_t> data;
    while(true) {
        marker = stream->readBits(16);
        fprintf(stderr, "Marker: %04X\n", marker);
        if(marker == 0xFFD9) {
            break;
        }

        uint16_t length = stream->readBits(16);
        switch(marker) {
            case 0xFFE0: // APP0
                type = stream->readBits(32); // "JFIF"
                stream->readBits(8); // "\0"
                version = stream->readBits(16); // Version
                stream->readBits(8); // Units
                stream->readBits(16); // X density
                stream->readBits(16); // Y density
                xthumb = stream->readBits(8); // X thumbnail
                ythumb = stream->readBits(8); // Y thumbnail                
                if(xthumb > 0 && ythumb > 0) {
                    stream->readBytes(xthumb * ythumb * 3);
                }
                break;
            case 0xFFDB: // DQT
                while(length > 2) {
                    
                    tmp = stream->readBits(8);
                    int table = tmp & 0x0f;
                    for(int j = 0; j < 64; j++) {
                        context->fjpeg_luminance_quantization_table[j] = stream->readBits(8);
                    }
                }                
                break;
            case 0xFFC0: // SOF0
                context->channels = stream->readBits(8);
                context->height = stream->readBits(16);
                context->width = stream->readBits(16);
                context->bpp = stream->readBits(8); // Bits per sample
                for(int i = 0; i < context->channels; i++) {
                    context->component_id[i] = stream->readBits(8); // Component ID
                    context->sampling_factors[i] = stream->readBits(8); // Sampling factors
                    context->quant_table[i] = stream->readBits(8); // Quant table
                }
                break;
            case 0xFFC4: // DHT
                while(length > 2) {
                    type = stream->readBits(4);
                    table = stream->readBits(4);
                    
                    for(int i = 0; i < 16; i++) {
                        bits[i] = stream->readBits(8);
                    }
                    int total = 0;
                    for(int i = 0; i < 16; i++) {
                        total += bits[i];
                    }
                    
                    if(type == 0) {
                        huffman_table = &context->fjpeg_short_huffman_luma_dc;
                        huffman_table_long = context->fjpeg_huffman_luma_dc;
                    } else if(type == 1) {
                        huffman_table = &context->fjpeg_short_huffman_luma_ac;
                        huffman_table_long = context->fjpeg_huffman_luma_ac;
                    } else if(type == 16) {
                        huffman_table = &context->fjpeg_short_huffman_chroma_dc;
                        huffman_table_long = context->fjpeg_huffman_chroma_dc;
                    } else if(type == 17) {
                        huffman_table = &context->fjpeg_short_huffman_chroma_ac;
                        huffman_table_long = context->fjpeg_huffman_chroma_ac;
                    }
                    if(huffman_table) {
                        for(int i = 0; i < total; i++) {
                            huffman_table->bits[i] = bits[i];
                            huffman_table->val[i] = stream->readBits(8);                        
                        }
                    }
                    fjpeg_generate_tables(huffman_table_long, huffman_table);
                }
                break;
            case 0xFFDA: // SOS
                stream->readBits(8); // Length
                context->channels = stream->readBits(8);
                for(int i = 0; i < context->channels; i++) {
                    stream->readBits(8); // Component ID
                    stream->readBits(8); // Huffman table
                }
                stream->readBits(8); // DCT coeff start
                stream->readBits(8); // DCT coeff end
                stream->readBits(8); // Successive Approximation
                break;
            case 0xFFFE: // COM
                data = stream->readBytes(length-2);
                fprintf(stderr, "Comment: %s\n", std::string(data.begin(), data.end()).c_str());
                break;
            default:
                fprintf(stderr, "Error: Unsupported marker %04X\n", marker);
                return false;
        }
    }


}

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

    // Calculate huffman statistics
    int last_dc_coeff[3] = {0, 0, 0};
    fjpeg_coeff_t dct_block[64];
    const int inc_xy = context->channels==1?8:16;
    const int max_uv = context->channels==1?1:2;
    fjpeg_huffman_statistics_t huff_stats;

    memset(&huff_stats, 0, sizeof(fjpeg_huffman_statistics_t));

    for(int y = 0; y < context->padded_height; y+=inc_xy) {
        for(int x = 0; x < context->padded_width; x+=inc_xy) {
            for(int v = 0; v < max_uv; v++) {
                for(int u = 0; u < max_uv; u++) {
                    fjpeg_extract_coeff_8x8(context, dct_block, x+u*8, y+v*8, 0);
                    last_dc_coeff[0] = fjpeg_entropy_stats(&huff_stats, context, dct_block, 0, last_dc_coeff[0]);                     
                }
            }
            if(context->channels==3) {
                fjpeg_extract_coeff_8x8(context, dct_block, x>>1, y>>1, 1);   
                last_dc_coeff[1] = fjpeg_entropy_stats(&huff_stats, context, dct_block, 1, last_dc_coeff[1]);

                fjpeg_extract_coeff_8x8(context, dct_block, x>>1, y>>1, 2);
                last_dc_coeff[2] = fjpeg_entropy_stats(&huff_stats, context, dct_block, 2, last_dc_coeff[2]);
            }
        }
    }
    // Generate optimal, length-limited Huffman tables for every table from the
    // gathered statistics instead of using the generic default tables.
    memset(context->fjpeg_huffman_luma_dc, 0, sizeof(context->fjpeg_huffman_luma_dc));
    memset(context->fjpeg_huffman_luma_ac, 0, sizeof(context->fjpeg_huffman_luma_ac));
    memset(context->fjpeg_huffman_chroma_dc, 0, sizeof(context->fjpeg_huffman_chroma_dc));
    memset(context->fjpeg_huffman_chroma_ac, 0, sizeof(context->fjpeg_huffman_chroma_ac));

    context->fjpeg_short_huffman_luma_dc = fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_luma_dc, huff_stats.luma_dc, 12);
    context->fjpeg_short_huffman_luma_ac = fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_luma_ac, huff_stats.luma_ac, 256);
    context->fjpeg_short_huffman_chroma_dc = fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_chroma_dc, huff_stats.chroma_dc, 12);
    context->fjpeg_short_huffman_chroma_ac = fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_chroma_ac, huff_stats.chroma_ac, 256);

    int luma_dc_count = 0;
    int luma_ac_count = 0;
    int chroma_dc_count = 0;
    int chroma_ac_count = 0;
    for (int i = 0; i < 16; i++) {
        luma_dc_count += context->fjpeg_short_huffman_luma_dc.bits[i];
        luma_ac_count += context->fjpeg_short_huffman_luma_ac.bits[i];
        chroma_dc_count += context->fjpeg_short_huffman_chroma_dc.bits[i];
        chroma_ac_count += context->fjpeg_short_huffman_chroma_ac.bits[i];
    }

    // DHT: luma DC and luma AC share one segment.
    stream->writeBits(0xFFC4, 16); // Huffman tables
    stream->writeBits(2 + (1 + 16 + luma_dc_count) + (1 + 16 + luma_ac_count), 16); // Length
    stream->writeBits(0, 4); // DC
    stream->writeBits(0, 4); // Table ID

    for (int i = 0; i < 16; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_dc.bits[i], 8);
    }

    for (int i = 0; i < luma_dc_count; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_dc.val[i], 8);
    }

    stream->writeBits(1, 4); // AC
    stream->writeBits(0, 4); // Table ID

    for (int i = 0; i < 16; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_ac.bits[i], 8);
    }

    for (int i = 0; i < luma_ac_count; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_ac.val[i], 8);
    }

    if(context->channels > 1) {
        stream->writeBits(0xFFC4, 16); // Huffman tables
        stream->writeBits(2 + (1 + 16 + chroma_dc_count) + (1 + 16 + chroma_ac_count), 16); // Length

        stream->writeBits(0, 4); // DC
        stream->writeBits(1, 4); // Table ID

        for (int i = 0; i < 16; i++) {
            stream->writeBits(context->fjpeg_short_huffman_chroma_dc.bits[i], 8);
        }

        for (int i = 0; i < chroma_dc_count; i++) {
            stream->writeBits(context->fjpeg_short_huffman_chroma_dc.val[i], 8);
        }

        stream->writeBits(1, 4); // AC
        stream->writeBits(1, 4); // Table ID

        for (int i = 0; i < 16; i++) {
            stream->writeBits(context->fjpeg_short_huffman_chroma_ac.bits[i], 8);
        }

        for (int i = 0; i < chroma_ac_count; i++) {
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

    stream->avoidFF = true;
    memset(&last_dc_coeff, 0, sizeof(last_dc_coeff));
    
    for(int y = 0; y < context->padded_height; y+=inc_xy) {
        for(int x = 0; x < context->padded_width; x+=inc_xy) {

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
    stream->padToByte();
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
    printf("  -d  Decode JPEG file\r\n");
    printf("  -h  Show help\r\n");
}
