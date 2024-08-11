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
//#include "fjpeg_tables.h"

//#define FJPEG_DEBUG_BLOCK 1
//#define FJPEG_DEBUG_COEFF 1

fjpeg_coeff_t* fjpeg_dct8x8(fjpeg_pixel_t* block, fjpeg_coeff_t* out) {
    for (int v = 0; v < FJPEG_BLOCK_SIZE; v++) {
        for (int u = 0; u < FJPEG_BLOCK_SIZE; u++) {
            float sum = 0.0f;
            float cu = (u == 0) ? 1.0f / sqrtf(2.f) : 1.0f;  // Scaling factors
            float cv = (v == 0) ? 1.0f / sqrtf(2.f) : 1.0f;
            for (int x = 0; x < FJPEG_BLOCK_SIZE; x++) {
                for (int y = 0; y < FJPEG_BLOCK_SIZE; y++) {                    
                    sum += (((float)block[y*FJPEG_BLOCK_SIZE+x])-128.f) * cosf((2.f * x + 1.f) * u * M_PI / 16.f) * cosf((2.f * y + 1.f) * v * M_PI / 16.f);
                }
            }
            out[v*FJPEG_BLOCK_SIZE+u] = 0.25f * sum  * cu * cv;  // Apply constant factor
        }
    }
    return out;
}

fjpeg_pixel_t* fjpeg_idct8x8(fjpeg_coeff_t* block, fjpeg_pixel_t* out) {
    for (int y = 0; y < FJPEG_BLOCK_SIZE; y++) {
        for (int x = 0; x < FJPEG_BLOCK_SIZE; x++) {
            float sum = 0.0f;
            for (int v = 0; v < FJPEG_BLOCK_SIZE; v++) {
                for (int u = 0; u < FJPEG_BLOCK_SIZE; u++) {
                    float cu = (u == 0) ? 1.0f / sqrtf(2.f) : 1.0f;  // Scaling factors
                    float cv = (v == 0) ? 1.0f / sqrtf(2.f) : 1.0f;
                    sum += cu * cv * block[v*FJPEG_BLOCK_SIZE+u] * cosf((2.f * x + 1.f) * u * M_PI / 16.f) * cosf((2.f * y + 1.f) * v * M_PI / 16.f);
                }
            }
            out[y*FJPEG_BLOCK_SIZE+x] = (fjpeg_pixel_t)((0.25f * sum)+128.0f);  // Apply constant factor
        }
    }
    return out;
}

fjpeg_coeff_t* fjpeg_quant8x8(fjpeg_context* context, fjpeg_coeff_t* input, fjpeg_coeff_t *output, int table) {

    uint8_t *quant_table = table == 0 ? context->fjpeg_luminance_quantization_table : context->fjpeg_chrominance_quantization_table;

    for (int i = 0; i < 64; i++) {
        output[i] = input[i] / (float)quant_table[i];
    }

    return output;
}

fjpeg_coeff_t* fjpeg_dequant8x8(fjpeg_context* context, fjpeg_coeff_t* input, fjpeg_coeff_t *output, int table) {

    uint8_t *quant_table = table == 0 ? context->fjpeg_luminance_quantization_table : context->fjpeg_chrominance_quantization_table;

    for (int i = 0; i < 64; i++) {
        output[i] = input[i] * quant_table[i];
    }

    return output;
}

fjpeg_coeff_t* fjpeg_zigzag8x8(fjpeg_coeff_t* block, fjpeg_coeff_t* out) {

    for (int i = 0; i < 64; i++) {
        out[fjpeg_zigzag_8x8[i]] = block[i];
    }

    return out;
}

fjpeg_coeff_t* fjpeg_izigzag8x8(fjpeg_coeff_t* block, fjpeg_coeff_t* out) {

    for (int i = 0; i < 64; i++) {
        out[i] = block[fjpeg_zigzag_8x8[i]];
    }

    return out;
}

fjpeg_pixel_t* fjpeg_extract_8x8(fjpeg_context* context, fjpeg_pixel_t* output, int x, int y, int channel) {    

    fjpeg_pixel_t* image = channel==0?context->fjpeg_y:channel==1?context->fjpeg_cb:context->fjpeg_cr;
    const int input_width = channel==0?context->width:context->width/2;

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            output[j * 8 + i] = image[(y + j) * input_width  + (x + i)];
        }
    }
    return output;
}

fjpeg_coeff_t* fjpeg_extract_coeff_8x8(fjpeg_context* context, fjpeg_coeff_t* output, int x, int y, int channel) {    

    fjpeg_coeff_t* image = channel==0?context->fjpeg_ydct:channel==1?context->fjpeg_cbdct:context->fjpeg_crdct;
    const int input_width = channel==0?context->width:context->width/2;

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            output[j * 8 + i] = image[(y + j) * input_width  + (x + i)];
        }
    }
    return output;
}

bool fjpeg_store_coeff_8x8(fjpeg_context* context, fjpeg_coeff_t* input, int x, int y, int channel) {

    fjpeg_coeff_t* image = channel==0?context->fjpeg_ydct:channel==1?context->fjpeg_cbdct:context->fjpeg_crdct;
    const int input_width = channel==0?context->width:context->width/2;

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            image[(y + j) * input_width + (x + i)] = input[j * 8 + i];
        }
    }
    return true;
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
    if(i < FJPEG_BLOCK_SIZE*FJPEG_BLOCK_SIZE) {
        // EOB
        #ifdef FJPEG_DEBUG_COEFF
        printf("EOB\r\n");
        #endif
        stream->writeBits(huff_ac[0x00].code, huff_ac[0x00].len); // EOB
    }

    return last_dc_coeff;
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
    stream->writeBits(strlen("FJPEG ")+strlen(FJPEG_VERSION) + 2, 16);

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

static const char *fjpeg_version() {
    return FJPEG_VERSION;
}


void print_usage() {
    printf("Usage: fjpeg [options]\r\n");
    printf("Example: fjpeg -i input.yuv -q 50 -r 1280x720 -o output.jpg\r\n");
    printf("Options:\r\n");
    printf("  -i <input_filename>  input YUV file\r\n");
    printf("  -q <quality>  Set quality factor (1-100)\r\n");
    printf("  -r <width>x<height>  Set resolution\r\n");
    printf("  -o <output_filename>  Output JPEG file\r\n");
    printf("  -h  Show help\r\n");
}

int main(int argc, char** argv) {
    printf("FJPEG %s\n", fjpeg_version());
    
    std::string input_filename;
    std::string output_filename;
    int quality = 50;
    int width = 0;
    int height = 0;

    // Parse filename, quality and resolution
    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-q") == 0) {
            if(i+1 < argc) {
                quality = atoi(argv[i+1]);
                if(quality < 1 || quality > 100) {
                    fprintf(stderr, "Error: Invalid quality value\n");
                    return 1;
                }
            } else {
                fprintf(stderr, "Error: Missing quality value\n");
                return 1;
            }
            i++;
        }
        else if(strcmp(argv[i], "-i") == 0) {
            if(i+1 < argc) {
                if(strlen(argv[i+1]) > 255) {
                    fprintf(stderr, "Error: Input filename too long\n");
                    return 1;
                }
                input_filename = argv[i+1];
            } else {
                fprintf(stderr, "Error: Missing input filename\n");
                return 1;
            }
            i++;
        }
        else if(strcmp(argv[i], "-r") == 0) {
            if(i+1 < argc) {
                if(sscanf(argv[i+1], "%dx%d", &width, &height) != 2) {
                    fprintf(stderr, "Error: Invalid resolution\n");
                    return 1;
                }
                if(width < 1 || height < 1) {
                    fprintf(stderr, "Error: Invalid resolution\n");
                    return 1;
                }
            } else {
                fprintf(stderr, "Error: Missing resolution\n");
                return 1;
            }
            i++;
        }
        else if(strcmp(argv[i], "-o") == 0) {
            if(i+1 < argc) {
                if(strlen(argv[i+1]) > 255) {
                    fprintf(stderr, "Error: Output filename too long\n");
                    return 1;
                }
                output_filename = argv[i+1];
            } else {
                fprintf(stderr, "Error: Missing output filename\n");
                return 1;
            }
            i++;
        }
        else if(strcmp(argv[i], "-h") == 0) {
            print_usage();
            return 0;
        }
    }

    if(input_filename.empty()) {
        fprintf(stderr, "Error: Missing input filename\n");
        print_usage();
        return 1;
    }
    if(output_filename.empty()) {
        fprintf(stderr, "Error: Missing output filename\n");
        print_usage();
        return 1;
    }
    if(width == 0 || height == 0) {
        fprintf(stderr, "Error: Missing resolution\n");
        print_usage();
        return 1;
    }

    // Time measurement
    int64_t time_input_read_ms = 0;
    int64_t time_dct_quant_ms = 0;
    int64_t time_header_ms = 0;


    fjpeg_context* context = new fjpeg_context();

    context->setQuality(quality);

    // Calculate time
    auto start = std::chrono::high_resolution_clock::now();
    context->readInput(input_filename.c_str(), width, height);
    auto end = std::chrono::high_resolution_clock::now();
    time_input_read_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    context->channels = 3;

    fjpeg_pixel_t cur_block[64];

    fjpeg_coeff_t dct_block[64];
    fjpeg_coeff_t dct_block2[64];

    fjpeg_pixel_t* image = new fjpeg_pixel_t[1280*720];

    start = std::chrono::high_resolution_clock::now();

    for(int y = 0; y < context->height; y+=8) {
        for(int x = 0; x < context->width; x+=8) {
            fjpeg_extract_8x8(context, cur_block, x, y, 0);
            fjpeg_dct8x8(cur_block, dct_block);
            fjpeg_quant8x8(context, dct_block,dct_block2, 0);
            fjpeg_zigzag8x8(dct_block2, dct_block);
            fjpeg_store_coeff_8x8(context, dct_block, x, y, 0);

            fjpeg_izigzag8x8(dct_block, dct_block2);
            fjpeg_dequant8x8(context, dct_block2, dct_block, 0);
            fjpeg_idct8x8(dct_block, cur_block);
            for (int j = 0; j < 8; j++) {
                for (int i = 0; i < 8; i++) {
                    image[(y + j) * context->width + (x + i)] = cur_block[j * 8 + i];
                }
            }
        }
    }

    if(context->channels == 3) {
        for(int y = 0; y < context->height/2; y+=8) {
            for(int x = 0; x < context->width/2; x+=8) {
                fjpeg_extract_8x8(context, cur_block, x, y, 1);
                fjpeg_dct8x8(cur_block, dct_block);
                fjpeg_quant8x8(context, dct_block,dct_block2, 1);
                fjpeg_zigzag8x8(dct_block2, dct_block);
                fjpeg_store_coeff_8x8(context, dct_block, x, y, 1);
            }
        }

        for(int y = 0; y < context->height/2; y+=8) {
            for(int x = 0; x < context->width/2; x+=8) {
                fjpeg_extract_8x8(context, cur_block, x, y, 2);
                fjpeg_dct8x8(cur_block, dct_block);
                fjpeg_quant8x8(context, dct_block,dct_block2, 2);
                fjpeg_zigzag8x8(dct_block2, dct_block);
                fjpeg_store_coeff_8x8(context, dct_block, x, y, 2);
            }
        }
    }

    end = std::chrono::high_resolution_clock::now();
    time_dct_quant_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

#ifdef FJPEG_DEBUG_DCT_BLOCK
    FILE* dct_out = fopen("dct.yuv", "wb");
    for(int y = 0; y < 720; y++) {
        for(int x = 0; x < 1280; x++) {
            uint8_t val = FJPEG_CLAMP(128+context->fjpeg_ydct[y*1280+x], 0, 255);
            fwrite(&val, 1,1, dct_out);
        }
    }
    fclose(dct_out);


    FILE* idct_out = fopen("idct.yuv", "wb");
    for(int y = 0; y < 720; y++) {
        for(int x = 0; x < 1280; x++) {
            uint8_t val = FJPEG_CLAMP(image[y*1280+x], 0, 255);
            fwrite(&val, 1,1, idct_out);
        }
    }
    fclose(idct_out);
#endif


#ifdef FJPEG_TEST_DCT
    //fjpeg_pixel_t cur_block[64];
    //fjpeg_coeff_t dct_block[64];
    fjpeg_pixel_t cur_block2[64] = {
        52, 55, 61, 66, 70, 61, 64, 73,
        63, 59, 55, 90, 109, 85, 69, 72,
        62, 59, 68, 113, 144, 104, 66, 73,
        63, 58, 71, 122, 154, 106, 70, 69,
        67, 61, 68, 104, 126, 88, 68, 70,
        79, 65, 60, 70, 77, 68, 58, 75,
        85, 71, 64, 59, 55, 61, 65, 83,
        87, 79, 69, 68, 65, 76, 78, 94
    };
    
    fjpeg_dct8x8(cur_block2, dct_block);
    for(int i = 0; i < 64; i++) {
        printf("%.2f ", dct_block[i]);
        if((i+1)%8 == 0) printf("\r\n");
    }
#endif

    FILE *fp = fopen(output_filename.c_str(), "wb");
    if(!fp) {
        fprintf(stderr, "Error: Unable to open output file\n");
        return 1;
    }
    fjpeg_bitstream* stream = new fjpeg_bitstream(fp);

    start = std::chrono::high_resolution_clock::now();
    fjpeg_generate_header(stream, context);
    end = std::chrono::high_resolution_clock::now();
    time_header_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    int file_size = ftell(fp);

    fclose(fp);
    
    printf("Time: Input read %d ms, DCT/Quant %d ms, Header %d ms\r\n", (int)time_input_read_ms, (int)time_dct_quant_ms, (int)time_header_ms);
    printf("Input size: %d bytes\r\n", context->width*context->height*3/2);
    printf("Output size: %d bytes\r\n", file_size);

    delete stream;
    delete [] image;

    return 0;

}