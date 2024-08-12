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
#pragma once

#include "fjpeg.h"

bool fjpeg_store_coeff_8x8(fjpeg_context* context, fjpeg_coeff_t* input, int x, int y, int channel);
fjpeg_coeff_t* fjpeg_extract_coeff_8x8(fjpeg_context* context, fjpeg_coeff_t* output, int x, int y, int channel);
fjpeg_pixel_t* fjpeg_extract_8x8(fjpeg_context* context, fjpeg_pixel_t* output, int x, int y, int channel);
fjpeg_coeff_t* fjpeg_izigzag8x8(fjpeg_coeff_t* block, fjpeg_coeff_t* out);
fjpeg_coeff_t* fjpeg_zigzag8x8(fjpeg_coeff_t* block, fjpeg_coeff_t* out);

fjpeg_coeff_t* fjpeg_dequant8x8(fjpeg_context* context, fjpeg_coeff_t* input, fjpeg_coeff_t *output, int table);
fjpeg_coeff_t* fjpeg_quant8x8(fjpeg_context* context, fjpeg_coeff_t* input, fjpeg_coeff_t *output, int table);
fjpeg_pixel_t* fjpeg_idct8x8(fjpeg_context* context, fjpeg_coeff_t* block, fjpeg_pixel_t* out);
fjpeg_coeff_t* fjpeg_dct8x8(fjpeg_context* context, fjpeg_pixel_t* block, fjpeg_coeff_t* out);

bool fjpeg_transquant_input(fjpeg_context* context);