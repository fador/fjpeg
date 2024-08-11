#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cctype>
#include <vector>
#include <cassert>
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

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            output[j * 8 + i] = image[(y + j) * context->width  + (x + i)];
        }
    }
    return output;
}

fjpeg_coeff_t* fjpeg_extract_coeff_8x8(fjpeg_context* context, fjpeg_coeff_t* output, int x, int y, int channel) {    

    fjpeg_coeff_t* image = channel==0?context->fjpeg_ydct:channel==1?context->fjpeg_cbdct:context->fjpeg_crdct;

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            output[j * 8 + i] = image[(y + j) * context->width  + (x + i)];
        }
    }
    return output;
}

bool fjpeg_store_coeff_8x8(fjpeg_context* context, fjpeg_coeff_t* input, int x, int y, int channel) {

    fjpeg_coeff_t* image = channel==0?context->fjpeg_ydct:channel==1?context->fjpeg_cbdct:context->fjpeg_crdct;

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            image[(y + j) * context->width + (x + i)] = input[j * 8 + i];
        }
    }
    return true;
}

// Function to encode a single block of quantized DCT coefficients
int fjpeg_entropy_encode_block(fjpeg_bitstream* stream, fjpeg_context* context, fjpeg_coeff_t* block, int channel, int last_dc) {
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
    printf("Coeff %d\r\n", coeff);
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

    printf("Writing DC coeff %d size %d huff len %d %d\n", orig_diff, size, context->fjpeg_huffman_luma_dc[size].len, context->fjpeg_huffman_luma_dc[size].code);
    stream->writeBits(context->fjpeg_huffman_luma_dc[size].code, context->fjpeg_huffman_luma_dc[size].len); // Write size code
    if(size != 0) {
        diff = orig_diff;
        if(sign) {
            diff = (1 << size) + diff - 1;
        }
        //stream->writeBits(sign, 1); // Write sign bit
        printf("Writing DC coeff %d size %d\n", diff, size);
        stream->writeBits(diff, size); // Write diff value
    }

    if(last_coeff == 0) {
        // All zero block
        stream->writeBits(context->fjpeg_huffman_luma_ac[0x00].code, context->fjpeg_huffman_luma_ac[0x00].len); // EOB
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
                printf("ZRL\r\n");
                stream->writeBits(context->fjpeg_huffman_luma_ac[0xF0].code, context->fjpeg_huffman_luma_ac[0xF0].len); // ZRL
                run_length = 0;
            }
        }
        printf("Run length %d\r\n", run_length);
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
        printf("Writing AC coeff %d size %d huff len %d %d\n", orig_coeff, size, context->fjpeg_huffman_luma_ac[(run_length << 4) + size].len, context->fjpeg_huffman_luma_ac[(run_length << 4) + size].code);
        stream->writeBits(context->fjpeg_huffman_luma_ac[(run_length << 4) + size].code,
                            context->fjpeg_huffman_luma_ac[(run_length << 4) + size].len); // Write run-length/size code
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
        printf("EOB\r\n");
        stream->writeBits(context->fjpeg_huffman_luma_ac[0x00].code, context->fjpeg_huffman_luma_ac[0x00].len); // EOB
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
        stream->writeBits(1, 4);
        stream->writeBits(0, 4);

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
    for(int y = 0; y < context->height; y+=8) {
        for(int x = 0; x < context->width; x+=8) {
            printf("Encoding block %d %d\n", x, y);
            fjpeg_extract_coeff_8x8(context, dct_block, x, y, 0);
            // Print out the block
            for(int i = 0; i < 64; i++) {
                printf("%3d ", (int16_t)(dct_block[i]+0.5f));
                if((i+1)%8 == 0) printf("\r\n");
            }
            last_dc_coeff[0] = fjpeg_entropy_encode_block(stream, context, dct_block, 0, last_dc_coeff[0]);
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



int main(void) {
    printf("FJPEG %s\n", fjpeg_version());
    
    fjpeg_context* context = new fjpeg_context();

    context->setQuality(50);

    context->readInput("G:\\WORK\\sequences\\KristenAndSara_1280x720_60.yuv", 1280, 720);

    context->channels = 1;

    fjpeg_pixel_t cur_block[64];

    fjpeg_coeff_t dct_block[64];
    fjpeg_coeff_t dct_block2[64];

    fjpeg_pixel_t* image = new fjpeg_pixel_t[1280*720];

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

                fjpeg_extract_8x8(context, cur_block, x, y, 2);
                fjpeg_dct8x8(cur_block, dct_block);
                fjpeg_quant8x8(context, dct_block,dct_block2, 1);
                fjpeg_zigzag8x8(dct_block2, dct_block);
                fjpeg_store_coeff_8x8(context, dct_block, x, y, 2);
            }
        }

        for(int y = 0; y < 720; y+=8) {
            for(int x = 0; x < 1280; x+=8) {
                fjpeg_extract_8x8(context, cur_block, x, y, 2);
                fjpeg_dct8x8(cur_block, dct_block);
                fjpeg_quant8x8(context, dct_block,dct_block2, 2);
                fjpeg_zigzag8x8(dct_block2, dct_block);
                fjpeg_store_coeff_8x8(context, dct_block, x, y, 2);

                fjpeg_izigzag8x8(dct_block, dct_block2);
                fjpeg_dequant8x8(context, dct_block2, dct_block, 2);
                fjpeg_idct8x8(dct_block, cur_block);
                for (int j = 0; j < 8; j++) {
                    for (int i = 0; i < 8; i++) {
                        image[(y + j) * context->width + (x + i)] = cur_block[j * 8 + i];
                    }
                }
            }
        }
    }

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

    FILE *fp = fopen("test.jpg", "wb");
    fjpeg_bitstream* stream = new fjpeg_bitstream(fp);

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


    fjpeg_generate_header(stream, context);

    fclose(fp);
    delete stream;
    delete [] image;

    return 0;

}