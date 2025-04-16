#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>

#include "fjpeg_global.h"
#include "fjpeg.h"
#include "fjpeg_bitstream.h"

typedef struct {
  uint32_t Qe;
  uint8_t NextMPS;
  uint8_t NextLPS;
  bool Switch;
} fjpeg_arith_table_t;

const fjpeg_arith_table_t fjpeg_arith_state_changes[] =
{{0x5A1D, 11, 11, 1}, {0x2586, 14, 12, 0}, {0x1114, 16, 13, 0}, {0x080B, 18, 14, 0},
 {0x03D8, 20, 15, 0}, {0x01DA, 23, 16, 0}, {0x00E5, 25, 17, 0}, {0x006F, 28, 18, 0},
 {0x0036, 30, 19, 0}, {0x001A, 33, 10, 0}, {0x000D, 35, 11, 0}, {0x0006, 19, 12, 0},
 {0x0003, 10, 13, 0}, {0x0001, 12, 13, 0}, {0x5A7F, 15, 15, 1}, {0x3F25, 36, 16, 0},
 {0x2CF2, 38, 17, 0}, {0x207C, 39, 18, 0}, {0x17B9, 40, 19, 0}, {0x1182, 42, 20, 0},
 {0x0CEF, 43, 21, 0}, {0x09A1, 45, 22, 0}, {0x072F, 46, 23, 0}, {0x055C, 48, 24, 0},
 {0x0406, 49, 25, 0}, {0x0303, 51, 26, 0}, {0x0240, 52, 27, 0}, {0x01B1, 54, 28, 0},
 {0x0144, 56, 29, 0}, {0x00F5, 57, 30, 0}, {0x00B7, 59, 31, 0}, {0x008A, 60, 32, 0},
 {0x0068, 62, 33, 0}, {0x004E, 63, 34, 0}, {0x003B, 32, 35, 0}, {0x002C, 33, 19, 0},
 {0x5AE1, 37, 37, 1}, {0x484C, 64, 38, 0}, {0x3A0D, 65, 39, 0}, {0x2EF1, 67, 40, 0},
 {0x261F, 68, 41, 0}, {0x1F33, 69, 42, 0}, {0x19A8, 70, 43, 0}, {0x1518, 72, 44, 0},
 {0x1177, 73, 45, 0}, {0x0E74, 74, 46, 0}, {0x0BFB, 75, 47, 0}, {0x09F8, 77, 48, 0},
 {0x0861, 78, 49, 0}, {0x0706, 79, 50, 0}, {0x05CD, 48, 51, 0}, {0x04DE, 50, 52, 0},
 {0x040F, 50, 53, 0}, {0x0363, 51, 54, 0}, {0x02D4, 52, 55, 0}, {0x025C, 53, 56, 0},
 {0x01F8, 54, 57, 0}, {0x01A4, 55, 58, 0}, {0x0160, 56, 59, 0}, {0x0125, 57, 60, 0},
 {0x00F6, 58, 61, 0}, {0x00CB, 59, 62, 0}, {0x00AB, 60, 63, 0}, {0x008F, 61, 32, 0},
 {0x5B12, 64, 64, 1}, {0x4D04, 79, 65, 0}, {0x412C, 80, 66, 0}, {0x37D8, 81, 67, 0},
 {0x2FE8, 82, 68, 0}, {0x293C, 83, 69, 0}, {0x2379, 85, 70, 0}, {0x1EDF, 86, 71, 0},
 {0x1AA9, 86, 72, 0}, {0x174E, 71, 73, 0}, {0x1424, 71, 74, 0}, {0x119C, 73, 75, 0},
 {0x0F6B, 73, 76, 0}, {0x0D51, 74, 77, 0}, {0x0BB6, 76, 78, 0}, {0x0A40, 76, 49, 0},
 {0x5832, 79, 79, 1}, {0x4D1C, 87, 80, 0}, {0x438E, 88, 81, 0}, {0x3BDD, 89, 82, 0},
 {0x34EE, 90, 83, 0}, {0x2EAE, 91, 84, 0}, {0x299A, 92, 85, 0}, {0x2516, 91, 70, 0},
 {0x5570, 87, 87, 1}, {0x4CA9, 95, 88, 0}, {0x44D9, 96, 89, 0}, {0x3E22, 97, 90, 0},
 {0x3824, 99, 91, 0}, {0x32B4, 99, 92, 0}, {0x2E17, 93, 85, 0}, {0x56A8, 95, 95, 1},
 {0x4F46, 51, 96, 0}, {0x47E5, 52, 97, 0}, {0x41CF, 53, 98, 0}, {0x3C3D, 54, 50, 0},
 {0x375E, 99, 92, 0}, {0x5231, 55, 52, 0}, {0x4C0F, 56, 53, 0}, {0x4639, 57, 54, 0},
 {0x415E, 53, 98, 0}, {0x5627, 55, 55, 1}, {0x50E7, 58, 56, 0}, {0x4B85, 59, 57, 0},
 {0x5597, 60, 58, 0}, {0x504F, 61, 57, 0}, {0x5A10, 60, 59, 1}, {0x5522, 62, 57, 0},
 {0x59EB, 62, 59, 1}};


typedef struct {
  fjpeg_bitstream* buffer;
  int ST;
  int C;
  int A;
  int CT;
  int BP;


} arithmetic_encoder_t;


bool fjpeg_arith_init_encoder(arithmetic_encoder_t *encoder, fjpeg_bitstream* buffer) {

  encoder->buffer = buffer;
  encoder->ST = 0;
  encoder->C = 0;
  encoder->CT = 11;
  encoder->A = 0x10000;
  encoder->BP = 0;

  return true;
}

bool fjpeg_arith_byte_out(arithmetic_encoder_t* encoder, uint8_t byte) {
  encoder->buffer->writeBits(byte, 8);
  while(encoder->BP > 0) {
    encoder->buffer->writeBits((byte >> 7) & 1, 1);
    encoder->BP--;
  }
  return true;
}


bool fjpeg_arith_renorm(arithmetic_encoder_t* encoder) {
  while(encoder->A < 0x8000) {
    if(encoder->A < 0x10000) {
      if(encoder->C + encoder->A > 0x10000) {
        fjpeg_arith_byte_out(encoder, encoder->ST);
        while(encoder->BP > 0) {
          fjpeg_arith_byte_out(encoder, (encoder->ST+1) & 0xff);
          encoder->BP--;
        }
        encoder->ST = (encoder->C + encoder->A) >> 8;
        encoder->C &= 0xff00;
        encoder->A = 0x10000;
      } else {
        encoder->A = 0x10000;
      }
    } else {
      encoder->A = 0x10000;
    }
    encoder->C <<= 8;
    encoder->A <<= 8;
    encoder->CT -= 8;
  }
  return true;
}



bool fjpeg_arith_code_msp(arithmetic_encoder_t* encoder, uint8_t mps, uint8_t state) {
  const fjpeg_arith_table_t* table = &fjpeg_arith_state_changes[state];
  if(mps == table->NextMPS) {
    encoder->A -= table->Qe;
    if(encoder->A < 0x8000) {
      if(encoder->A < table->Qe) {
        encoder->A = table->Qe;
      } else {
        encoder->C += table->Qe;
        encoder->A = table->Qe;
      }
      fjpeg_arith_renorm(encoder);
    }
  } else {
    encoder->A -= table->Qe;
    if(encoder->A >= table->Qe) {
      encoder->C += table->Qe;
    } else {
      encoder->A = table->Qe;
    }
    fjpeg_arith_renorm(encoder);
  }
  return true;
}

bool fjpeg_arith_code_lsp(arithmetic_encoder_t* encoder, uint8_t lsp, uint8_t state) {
  const fjpeg_arith_table_t* table = &fjpeg_arith_state_changes[state];
  if(lsp == table->NextLPS) {
    encoder->A -= table->Qe;
    if(encoder->A < 0x8000) {
      if(encoder->A < table->Qe) {
        encoder->A = table->Qe;
      } else {
        encoder->C += table->Qe;
        encoder->A = table->Qe;
      }
      fjpeg_arith_renorm(encoder);
    }
  } else {
    encoder->A -= table->Qe;
    if(encoder->A >= table->Qe) {
      encoder->C += table->Qe;
    } else {
      encoder->A = table->Qe;
    }
    fjpeg_arith_renorm(encoder);
  }
  return true;
}

bool fjpeg_arith_flush(arithmetic_encoder_t* encoder) {
  fjpeg_arith_byte_out(encoder, encoder->ST);
  while(encoder->BP > 0) {
    fjpeg_arith_byte_out(encoder, (encoder->ST+1) & 0xff);
    encoder->BP--;
  }
  return true;
}

bool fjpeg_arith_encode_DC(arithmetic_encoder_t* encoder, int value, int diff) {
  int diffabs = abs(diff);
  int diffsign = diff < 0;
  int diffbits = 0;
  while(diffabs > 0) {
    diffabs >>= 1;
    diffbits++;
  }
  fjpeg_arith_code_msp(encoder, diffsign, 0);
  for(int i = diffbits-2; i >= 0; i--) {
    fjpeg_arith_code_lsp(encoder, (value >> i) & 1, 0);
  }
  return true;
}

uint32_t fjpeg_arith_test[] = {0x00020051, 0x000000C0, 0x0352872A, 0xAAAAAAAA, 0x82C02000, 0xFCD79EF6, 0x74EAABF7, 0x697EE74C};

int main() {
  FILE* out = fopen("test.arith", "wb");
  fjpeg_bitstream* buffer = new fjpeg_bitstream(out);
  arithmetic_encoder_t encoder;
  fjpeg_arith_init_encoder(&encoder, buffer);
  for(int i = 0; i < 8; i++) {
    fjpeg_arith_code_msp(&encoder, (fjpeg_arith_test[i] >> 31) & 1, 0);
    for(int j = 30; j >= 0; j--) {
      fjpeg_arith_code_lsp(&encoder, (fjpeg_arith_test[i] >> j) & 1, 0);
    }
  }
  fjpeg_arith_flush(&encoder);
  fclose(out);
  delete buffer;
  return 0;
}