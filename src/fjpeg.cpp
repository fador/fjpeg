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

// Bitstream handling
typedef struct fjpeg_bitstream {
    uint32_t current;
    int offset;
    FILE *fp;
    std::vector<uint8_t> buffer;

    fjpeg_bitstream(FILE *fp) : current(0), offset(0), fp(fp) {}

    void writeBits(uint32_t input, int bits) {
        assert(bits > 0 && bits <= 24);

        // Is there space in the current buffer
        if (offset + bits >= 32) {
            // Flush the buffer
            while(offset >= 8) {
                buffer.push_back(current & 0xff);
                current >>= 8;
                offset -= 8;
            }
        }
        // Write the remaining bits
        current <<= bits;
        current |= input;
        offset += bits;
    }

    void flushToFile() {
        if (offset > 0) {
            buffer.push_back(current & 0xFF);
        }
        fwrite(buffer.data(), 1, buffer.size(), fp);
        buffer.clear();
    }

    ~fjpeg_bitstream() {
        flushToFile();
    }
} fjpeg_bitstream_t;


fjpeg_coeff_t* fjpeg_dct8x8(fjpeg_pixel_t* block, fjpeg_coeff_t* out) {
    fjpeg_coeff_t tmp[64];

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            float sum = 0.0f;
            for (int v = 0; v < 8; v++) {
                for (int u = 0; u < 8; u++) {
                    float a = (u == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
                    float b = (v == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
                    sum += a * b * block[v * 8 + u] * cosf((2.0f * x + 1.0f) * u * M_PI / 16.0f) * cosf((2.0f * y + 1.0f) * v * M_PI / 16.0f);
                }
            }
            tmp[y * 8 + x] = sum / 4.0f;
        }
    }

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            float sum = 0.0f;
            for (int v = 0; v < 8; v++) {
                for (int u = 0; u < 8; u++) {
                    float a = (u == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
                    float b = (v == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
                    sum += a * b * tmp[v * 8 + u] * cosf((2.0f * x + 1.0f) * u * M_PI / 16.0f) * cosf((2.0f * y + 1.0f) * v * M_PI / 16.0f);
                }
            }
            out[y * 8 + x] = sum / 4.0f;//FJPEG_CLAMP(sum / 4.0f + 128.0f, 0, 255);
        }
    }

    return out;
}

fjpeg_pixel_t* fjpeg_idct8x8(fjpeg_coeff_t* block, fjpeg_pixel_t* out) {
    fjpeg_coeff_t tmp[64];

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            float sum = 0.0f;
            for (int v = 0; v < 8; v++) {
                for (int u = 0; u < 8; u++) {
                    float a = (u == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
                    float b = (v == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
                    sum += a * b * (block[v * 8 + u] - 128.0f) * cosf((2.0f * x + 1.0f) * u * M_PI / 16.0f) * cosf((2.0f * y + 1.0f) * v * M_PI / 16.0f);
                }
            }
            tmp[y * 8 + x] = sum / 4.0f;
        }
    }

    for (int y = 0; y < 8; y++) {
        for (int x = 0; x < 8; x++) {
            float sum = 0.0f;
            for (int v = 0; v < 8; v++) {
                for (int u = 0; u < 8; u++) {
                    float a = (u == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
                    float b = (v == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
                    sum += a * b * tmp[v * 8 + u] * cosf((2.0f * x + 1.0f) * u * M_PI / 16.0f) * cosf((2.0f * y + 1.0f) * v * M_PI / 16.0f);
                }
            }
            out[y * 8 + x] = FJPEG_CLAMP(sum / 4.0f, 0, 255);
        }
    }

    return out;
}

fjpeg_coeff_t* fjpeg_quant8x8(fjpeg_context* context, fjpeg_coeff_t* input, fjpeg_coeff_t *output, int table) {

    uint8_t *quant_table = table == 0 ? context->fjpeg_luminance_quantization_table : context->fjpeg_chrominance_quantization_table;

    for (int i = 0; i < 64; i++) {
        output[i] = input[i] / quant_table[i];
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

bool fjpeg_store_coeff_8x8(fjpeg_context* context, fjpeg_coeff_t* input, int x, int y, int channel) {

    fjpeg_coeff_t* image = channel==0?context->fjpeg_ydct:channel==1?context->fjpeg_cbdct:context->fjpeg_crdct;

    for (int j = 0; j < 8; j++) {
        for (int i = 0; i < 8; i++) {
            image[(y + j) * context->width + (x + i)] = input[j * 8 + i];
        }
    }
    return true;
}

// Generate jpeg header
bool fjpeg_generate_header(fjpeg_bitstream* stream, fjpeg_context* context) {
    
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
    
    for (int i = 0; i < 64; i++) {
        stream->writeBits(context->fjpeg_luminance_quantization_table[fjpeg_zigzag_8x8[i]], 8);
    }

    if(context->channels > 1) {
        // chroma quant
        stream->writeBits(0xFFDB, 16);
        stream->writeBits(67, 16);
        stream->writeBits(1, 4);
        stream->writeBits(0, 4);

        for (int i = 0; i < 64; i++) {
            stream->writeBits(context->fjpeg_chrominance_quantization_table[fjpeg_zigzag_8x8[i]], 8);
        }
    }
    
   
    // SOF0
    stream->writeBits(0xFFC0, 16);
    stream->writeBits(17, 16);
    stream->writeBits(8, 8); // 8 bits per sample
    stream->writeBits(context->height, 16);
    stream->writeBits(context->width, 16);
    stream->writeBits(context->channels, 8);

    for (int i = 0; i < context->channels; i++) {
        stream->writeBits(i + 1, 8);
        stream->writeBits(i == 0 ? 0x22 : 0x11, 8); // Sampling factors 4:2:0
        stream->writeBits(i==0?0:1, 8); // Quant table
    }


    // DHT
    stream->writeBits(0xFFC4, 16); // Huffman tables
    stream->writeBits(0x01A2, 16); // Length
    stream->writeBits(0, 4); // DC
    stream->writeBits(0, 4); // Table ID

    for (int i = 0; i < 16; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_dc.bits[i], 8);
    }

    for (int i = 0; i < 12; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_dc.val[i], 8);
    }

    stream->writeBits(0x10, 4); // AC
    stream->writeBits(1, 4); // Table ID

    for (int i = 0; i < 16; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_ac.bits[i], 8);
    }

    for (int i = 0; i < 162; i++) {
        stream->writeBits(context->fjpeg_short_huffman_luma_ac.val[i], 8);
    }

    if(context->channels > 1) {
        stream->writeBits(0xFFC4, 16); // Huffman tables
        stream->writeBits(0x01A2, 16); // Length

        stream->writeBits(0, 4); // DC
        stream->writeBits(1, 4); // Table ID

        for (int i = 0; i < 16; i++) {
            stream->writeBits(context->fjpeg_short_huffman_chroma_dc.bits[i], 8);
        }

        for (int i = 0; i < 12; i++) {
            stream->writeBits(context->fjpeg_short_huffman_chroma_dc.val[i], 8);
        }

        stream->writeBits(0x11, 4); // AC
        stream->writeBits(1, 4); // Table ID

        for (int i = 0; i < 16; i++) {
            stream->writeBits(context->fjpeg_short_huffman_chroma_ac.bits[i], 8);
        }

        for (int i = 0; i < 162; i++) {
            stream->writeBits(context->fjpeg_short_huffman_chroma_ac.val[i], 8);
        }
    }

    // SOS
    stream->writeBits(0xFFDA, 16);
    stream->writeBits(12, 16);
    stream->writeBits(context->channels, 8);

    for (int i = 0; i < context->channels; i++) {
        stream->writeBits(i + 1, 8);
        stream->writeBits(i == 0 ? 0 : 0x11, 8);
    }

    stream->writeBits(0, 8);
    stream->writeBits(0x3F, 8);
    stream->writeBits(0, 8);

    // COM
    stream->writeBits(0xFFFE, 16);
    stream->writeBits(strlen("FJPEG ")+strlen(FJPEG_VERSION), 16);

    for(int i = 0; i < strlen("FJPEG "); i++) {
        stream->writeBits("FJPEG "[i], 8);
    }

    for(int i = 0; i < strlen(FJPEG_VERSION); i++) {
        stream->writeBits(FJPEG_VERSION[i], 8);
    }

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

    context->readInput("G:\\WORK\\sequences\\KristenAndSara_1280x720_60.yuv", 1280, 720);

    context->channels = 1;

    fjpeg_pixel_t cur_block[64];
    fjpeg_pixel_t block[64];

    fjpeg_coeff_t dct_block[64];
    fjpeg_coeff_t dct_block2[64];


    for(int y = 0; y < 720; y+=8) {
        for(int x = 0; x < 1280; x+=8) {
            fjpeg_extract_8x8(context, cur_block, x, y, 0);
            fjpeg_dct8x8(cur_block, dct_block);
            fjpeg_quant8x8(context, dct_block,dct_block2, 0);
            fjpeg_zigzag8x8(dct_block2, dct_block);
            fjpeg_store_coeff_8x8(context, dct_block, x, y, 0);
        }
    }

    FILE *fp = fopen("test.jpg", "wb");
    fjpeg_bitstream* stream = new fjpeg_bitstream(fp);

    fjpeg_generate_header(stream, context);

    return 0;

}