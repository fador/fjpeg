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
#include "fjpeg_bitstream.h"

// Standard ITU-T T.81 Annex H spatial predictor function
static inline int fjpeg_lossless_predict(int tl, int t, int l, int predictor) {
    switch (predictor) {
        case 0: return 0;
        case 1: return l;
        case 2: return t;
        case 3: return tl;
        case 4: return l + t - tl;
        case 5: return l + ((t - tl) >> 1);
        case 6: return t + ((l - tl) >> 1);
        case 7:
        default: return (l + t) >> 1;
    }
}

// Generate ITU-T T.81 Annex H Lossless JPEG (SOF3).
// If predictor is 0, evaluates candidate predictors (1..7) and picks the one
// minimizing total compressed bitstream size.
bool fjpeg_generate_lossless(fjpeg_bitstream* stream, fjpeg_context* context, int predictor);
