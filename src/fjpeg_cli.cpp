
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

#include "fjpeg.h"
#include "fjpeg_bitstream.h"
#include "fjpeg_transquant.h"

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
            fjpeg_print_usage();
            return 0;
        }
    }

    if(input_filename.empty()) {
        fprintf(stderr, "Error: Missing input filename\n");
        fjpeg_print_usage();
        return 1;
    }
    if(output_filename.empty()) {
        fprintf(stderr, "Error: Missing output filename\n");
        fjpeg_print_usage();
        return 1;
    }
    if(width == 0 || height == 0) {
        fprintf(stderr, "Error: Missing resolution\n");
        fjpeg_print_usage();
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


    #ifdef FJPEG_DEBUG_DCT_BLOCK
    fjpeg_pixel_t* image = new fjpeg_pixel_t[1280*720];
    #endif

    start = std::chrono::high_resolution_clock::now();
    fjpeg_transquant_input(context);

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
    #ifdef FJPEG_DEBUG_DCT_BLOCK
    delete [] image;
    #endif

    return 0;

}