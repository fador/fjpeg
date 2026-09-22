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

// Gather Huffman symbol statistics from the quantized (integer-valued) blocks.
static void fjpeg_gather_stats(fjpeg_context* context, fjpeg_huffman_statistics_t* huff_stats) {
    int last_dc_coeff[3] = {0, 0, 0};
    fjpeg_coeff_t dct_block[64];
    const int inc_xy = context->channels==1?8:16;
    const int max_uv = context->channels==1?1:2;

    memset(huff_stats, 0, sizeof(fjpeg_huffman_statistics_t));

    for(int y = 0; y < context->padded_height; y+=inc_xy) {
        for(int x = 0; x < context->padded_width; x+=inc_xy) {
            for(int v = 0; v < max_uv; v++) {
                for(int u = 0; u < max_uv; u++) {
                    fjpeg_extract_coeff_8x8(context, dct_block, x+u*8, y+v*8, 0);
                    last_dc_coeff[0] = fjpeg_entropy_stats(huff_stats, context, dct_block, 0, last_dc_coeff[0]);
                }
            }
            if(context->channels==3) {
                fjpeg_extract_coeff_8x8(context, dct_block, x>>1, y>>1, 1);
                last_dc_coeff[1] = fjpeg_entropy_stats(huff_stats, context, dct_block, 1, last_dc_coeff[1]);

                fjpeg_extract_coeff_8x8(context, dct_block, x>>1, y>>1, 2);
                last_dc_coeff[2] = fjpeg_entropy_stats(huff_stats, context, dct_block, 2, last_dc_coeff[2]);
            }
        }
    }
}

// Build all four optimal, length-limited Huffman tables from the statistics.
static void fjpeg_generate_huffman_tables(fjpeg_context* context, fjpeg_huffman_statistics_t* huff_stats) {
    memset(context->fjpeg_huffman_luma_dc, 0, sizeof(context->fjpeg_huffman_luma_dc));
    memset(context->fjpeg_huffman_luma_ac, 0, sizeof(context->fjpeg_huffman_luma_ac));
    memset(context->fjpeg_huffman_chroma_dc, 0, sizeof(context->fjpeg_huffman_chroma_dc));
    memset(context->fjpeg_huffman_chroma_ac, 0, sizeof(context->fjpeg_huffman_chroma_ac));

    context->fjpeg_short_huffman_luma_dc = fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_luma_dc, huff_stats->luma_dc, 12);
    context->fjpeg_short_huffman_luma_ac = fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_luma_ac, huff_stats->luma_ac, 256);
    context->fjpeg_short_huffman_chroma_dc = fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_chroma_dc, huff_stats->chroma_dc, 12);
    context->fjpeg_short_huffman_chroma_ac = fjpeg_generate_huffman_from_stats(context->fjpeg_huffman_chroma_ac, huff_stats->chroma_ac, 256);
}

// Re-quantize every block with the rate-distortion trellis, writing the chosen
// integer coefficients back into the DCT arrays.
static void fjpeg_run_trellis(fjpeg_context* context) {
    fjpeg_coeff_t in_block[64];
    fjpeg_coeff_t out_block[64];
    const int inc_xy = context->channels==1?8:16;
    const int max_uv = context->channels==1?1:2;

    for(int y = 0; y < context->padded_height; y+=inc_xy) {
        for(int x = 0; x < context->padded_width; x+=inc_xy) {
            for(int v = 0; v < max_uv; v++) {
                for(int u = 0; u < max_uv; u++) {
                    fjpeg_extract_coeff_8x8(context, in_block, x+u*8, y+v*8, 0);
                    fjpeg_trellis_quant_block(context, in_block, context->fjpeg_huffman_luma_ac, 0, out_block);
                    fjpeg_store_coeff_8x8(context, out_block, x+u*8, y+v*8, 0);
                }
            }
            if(context->channels==3) {
                fjpeg_extract_coeff_8x8(context, in_block, x>>1, y>>1, 1);
                fjpeg_trellis_quant_block(context, in_block, context->fjpeg_huffman_chroma_ac, 1, out_block);
                fjpeg_store_coeff_8x8(context, out_block, x>>1, y>>1, 1);

                fjpeg_extract_coeff_8x8(context, in_block, x>>1, y>>1, 2);
                fjpeg_trellis_quant_block(context, in_block, context->fjpeg_huffman_chroma_ac, 2, out_block);
                fjpeg_store_coeff_8x8(context, out_block, x>>1, y>>1, 2);
            }
        }
    }
}

// Apply rate-distortion optimized re-quantization if enabled. The current
// Huffman tables are built from the plain quantization and used as the rate
// model; the caller is expected to re-gather statistics afterwards.
void fjpeg_trellis_optimize(fjpeg_context* context) {
    if (context->trellis_lambda <= 0.0f) {
        return;
    }
    fjpeg_huffman_statistics_t stats;
    fjpeg_gather_stats(context, &stats);
    fjpeg_generate_huffman_tables(context, &stats);
    fjpeg_run_trellis(context);
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

    // Optional rate-distortion optimized re-quantization, then gather the
    // Huffman statistics and generate the optimal tables.
    fjpeg_trellis_optimize(context);

    fjpeg_huffman_statistics_t huff_stats;
    fjpeg_gather_stats(context, &huff_stats);
    fjpeg_generate_huffman_tables(context, &huff_stats);

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

    int last_dc_coeff[3] = {0, 0, 0};
    fjpeg_coeff_t dct_block[64];
    const int inc_xy = context->channels==1?8:16;
    const int max_uv = context->channels==1?1:2;

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
    printf("  -t  Enable rate-distortion optimized (trellis) quantization\r\n");
    printf("  -p  Write a progressive JPEG (spectral selection)\r\n");
    printf("  -l <lambda>  Set trellis Lagrange multiplier (default 0.01)\r\n");
    printf("  -d  Decode JPEG file\r\n");
    printf("  -h  Show help\r\n");
}
