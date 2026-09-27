**FJPEG**

**Overview**

FJPEG is a compact JPEG encoder written in C++, designed for educational
purposes and experimentation. It reads raw YUV 4:2:0 input and produces a
baseline JPEG file, with a clear, self-contained implementation of every
stage of the pipeline.

**Key Features**

* Baseline (sequential, 8-bit), extended sequential (12-bit, SOF1), progressive (SOF2), and arithmetic-coded sequential (SOF9) JPEG encoding
* Full pipeline: raw YUV input → FDCT → quantization → zigzag → run-length
  coding → Huffman or QM-coder arithmetic entropy coding → bitstream with 0xFF byte stuffing
* Per-image optimal Huffman tables (DC and AC, luma and chroma), generated
  from a first-pass statistics scan
* Arithmetic coding mode: ITU-T T.81 Annex F/D (SOF9 + DAC) with adaptive QM-coder binary probability estimation state machine (8-bit and 12-bit)
* Optional rate-distortion optimized (trellis) quantization of the AC
  coefficients
* Progressive mode: separate DC (with successive approximation) and
  per-component AC spectral-selection scans
* Lossless JPEG mode: ITU-T T.81 Annex H (SOF3) spatial DPCM with predictors
  1-7 and automatic rate-distortion predictor selection (8-bit and 12-bit)
* 12-bit sample precision support for Lossless (SOF3), Extended Sequential DCT (SOF1), and Arithmetic Coding (SOF9)
* Standard ITU-T T.81 quantization tables with libjpeg-compatible quality
  scaling, plus a tunable quantization deadzone
* 4:2:0 chroma subsampling and a separable (2-pass) FDCT
* Baseline, progressive, arithmetic, 12-bit, and lossless JPEG decoding back to raw YUV 4:2:0 (8-bit or 16-bit words)

**Building**

Requires CMake 3.12+ and a C++11 compiler.

```bash
cmake -S . -B build
cmake --build build --config Release
```

The CLI binary is produced as `build/<config>/fjpeg` (`fjpeg.exe` on Windows).
On single-config generators (Make/Ninja) it is `build/fjpeg`.

**Usage**

```bash
# input.yuv must be raw YUV 4:2:0, width*height Y bytes
# followed by (width*height)/4 Cb and (width*height)/4 Cr bytes
./fjpeg -i input.yuv -r 1280x720 -q 70 -o output.jpg

# arithmetic coding (SOF9 / DAC) for 20-30% smaller files at identical quality
./fjpeg -i input.yuv -r 1280x720 -a -q 70 -o output_arith.jpg

# lossless compression with automatic predictor selection
./fjpeg -i input.yuv -r 1280x720 -lossless -o output.jpg

# 12-bit extended sequential compression from 16-bit little-endian raw YUV
./fjpeg -i input_12bit.yuv -r 1280x720 -b 12 -q 70 -o output_12bit.jpg

# decode baseline, progressive, arithmetic, 12-bit, or lossless JPEG back to raw YUV 4:2:0
./fjpeg -d -i input.jpg -o output.yuv
```

The decoder accepts sequential baseline (SOF0), extended sequential (SOF1), progressive (SOF2), lossless (SOF3),
and arithmetic sequential (SOF9) JPEGs. The resolution and bit depth are taken from the JPEG headers; the output
is raw YUV with each component stored at its own sampling resolution (Y plane, then Cb, then Cr for 4:2:0).

Options:

```
-i <input_filename>    input raw YUV 4:2:0 file
-r <width>x<height>    frame resolution
-q <quality>           quality factor (1-100)
-o <output_filename>   output JPEG file
-a, -arith             write an arithmetic-coded JPEG (ITU-T T.81 SOF9 / DAC)
-t                     enable rate-distortion optimized (trellis) quantization
-l <lambda>            trellis Lagrange multiplier (default 0.007)
-p                     write a progressive JPEG
-lossless, -ll         write a lossless JPEG (ITU-T T.81 SOF3 DPCM)
-pred <1-7>            select lossless predictor (1-7, default 0=auto best)
-b <8|12>              sample bit depth (8 or 12, default 8)
-d                     decode an existing JPEG (baseline, progressive, arithmetic, 12-bit, or lossless)
-h                     show help
```

**Source Layout**

| File | Purpose |
| --- | --- |
| `src/fjpeg.cpp` | JPEG header generation, encode loop, CLI usage text |
| `src/fjpeg.h` | `fjpeg_context`, quality scaling, input loading |
| `src/fjpeg_transquant.cpp` | block extraction, FDCT/IDCT, quantization, zigzag, trellis |
| `src/fjpeg_huffman.cpp` | statistics pass, optimal Huffman generation, entropy coding |
| `src/fjpeg_progressive.cpp` | progressive frame/scan layout and scan entropy coding |
| `src/fjpeg_lossless.h` | lossless JPEG declarations and predictor helper |
| `src/fjpeg_lossless.cpp` | ITU-T T.81 SOF3 lossless DPCM encoding and auto-prediction |
| `src/fjpeg_bitstream.h` | bit reader/writer, 0xFF stuffing, file flushing |
| `src/fjpeg_global.h` | types, default quantization tables, zigzag tables |
| `src/fjpeg_huffman.h` | default Huffman tables and statistics struct |
| `src/fjpeg_decode.cpp` | JPEG header parsing, entropy decoding (Huffman & arithmetic), IDCT, YUV output |
| `src/fjpeg_cli.cpp` | command line parsing and program flow |
| `src/fjpeg_arith.h` | ITU-T T.81 Table D.2 probability estimation tables and QM-coder definitions |
| `src/fjpeg_arith.cpp` | ITU-T T.81 SOF9 / DAC arithmetic sequential DCT encoder |

**Compression Performance**

The encoder performs a statistics pass over all DCT blocks before writing the
entropy-coded data, then builds Huffman tables tailored to the actual image.
Verified improvements, each measured with a rate-distortion sweep (BD-rate
reduction at matched PSNR) over several natural and synthetic images:

1. **Optimal Huffman tables for all four tables.** Previously only the luma DC
   table was generated from the statistics; the luma AC, chroma DC and chroma
   AC tables used the generic defaults. All four are now generated with a
   correct, length-limited (≤ 16 bits) Huffman construction based on JPEG
   Annex K.2, and the DHT segments are written dynamically to match the
   generated symbol counts. Together with the bitstream fixes below this gave
   ~10-24% smaller files at equal PSNR.
2. **Bitstream tail fixes.** The final partial byte was written misaligned and
   the last few bits were dropped, and the EOI marker was appended without
   first padding the entropy segment to a byte boundary (so `0xFFD9` never
   appeared literally in the file). Both are fixed: the tail is flushed
   correctly, the entropy segment is padded with 1-bits, and the file now ends
   with a proper EOI marker.
3. **Symmetric coefficient rounding.** `(int)(x + 0.5f)` truncates towards zero
   for negative `x`, so every negative quantized coefficient was shrunk by
   about one level. Using round-half-away-from-zero removed the bias: a
   **17-33% BD-rate** improvement.
4. **Standard luminance quantization table.** Several high-frequency entries of
   the built-in table had an extra leading digit (e.g. `24 → 124`), over
   quantizing detail. The standard ITU-T T.81 table gives an **8-11% BD-rate**
   improvement on natural images (neutral on very smooth content).
5. **Standard libjpeg quality scaling.** The linear `(100 - quality)` ramp was
   replaced by the standard `5000/quality` / `200 - 2*quality` curve. It is
   rate-distortion neutral but quality values below 50 now reach roughly half
   the previous minimum bitrate.
6. **Quantization deadzone.** `FJPEG_QUANT_DEADZONE` (0.65) widens the zero bin,
   dropping near-zero coefficients for a **3-5.5% BD-rate** gain.
7. **Separable FDCT.** The transform is computed as two 1-D passes, reducing
   work from O(N⁴) to O(2·N³) per block (4096 → 1024 multiplies).
8. **Trellis quantization (opt-in, `-t`).** A dynamic program over the 63 AC
   coefficients jointly chooses zero/nonzero, the magnitude, and the end-of-block
   position, minimizing squared error plus `lambda` times the Huffman+VLI rate.
   The distortion is weighted by the squared quantization step and the multiplier
   is scaled by the mean step so the tradeoff is quality independent. Using the
   first-pass Huffman tables as the rate model and rebuilding them afterwards
   gives a mean **−3.8% BD-rate** across natural images (−3.1% to −5.1%),
   including 720p.

On a synthetic gradient/texture 640×480 frame (byte sizes; PSNR in dB against
the source), original vs. the current encoder:

| Quality | Original | Current | Size change | PSNR (orig → new) |
| ---: | ---: | ---: | ---: | --- |
| 25 | 10226 | 7456 | −27.1% | 37.41 → 42.18 |
| 50 | 12839 | 11739 | −8.6% | 39.09 → 46.59 |
| 75 | 17430 | 16574 | −4.9% | 44.91 → 49.69 |
| 90 | 27992 | 26394 | −5.7% | 50.21 → 52.70 |

The per-quality size change is no longer the headline number because the
quantization fixes deliberately spend some of the saved bits on accuracy. The
correct measure is rate-distortion: against the original encoder the current
one achieves about **−40% to −50% BD-rate** on natural images (the original
also suffered from the rounding and tail bugs). At 1280×720 / quality 75 the
separable FDCT reduced the DCT/quantization stage from ~76 ms to ~9 ms on the
test machine.

**Progressive JPEG**

Passing `-p` writes a progressive (SOF2) JPEG instead of the baseline file. The
encoder implements both spectral selection and AC successive approximation (ITU-T T.81 Annex G)
with block-level EOB run-length coding. The default scan script is:

1. DC, all components, first scan (`Ah=0, Al=1`)
2. DC refinement (`Ah=1, Al=0`)
3. Luma AC initial scan (`Ss=1..63, Ah=0, Al=1`)
4. Luma AC refinement (`Ss=1..63, Ah=1, Al=0`)
5. Cb AC initial scan (`Ss=1..63, Ah=0, Al=1`)
6. Cb AC refinement (`Ss=1..63, Ah=1, Al=0`)
7. Cr AC initial scan (`Ss=1..63, Ah=0, Al=1`)
8. Cr AC refinement (`Ss=1..63, Ah=1, Al=0`)

Every scan carries its own optimal Huffman tables, written in a DHT segment
immediately before the scan. Initial AC scans accumulate consecutive zero blocks into
EOB runs (`0x00..0xE0`), and refinement scans encode newly non-zero coefficients
along with point-transform refinement bits for previously non-zero coefficients.

Decoded progressive files are **byte-for-byte identical** to the corresponding
baseline decode, verified across images and qualities with both the internal decoder,
libjpeg (Pillow), and ffmpeg. With AC successive approximation and EOB runs,
progressive mode compresses natural images smaller than baseline (typically −2% to −3%
file size reduction, and up to −9% BD-rate when combined with trellis quantization).

**Arithmetic Coding (ITU-T T.81 SOF9 / DAC)**

Passing `-a` or `-arith` enables standard JPEG arithmetic entropy coding instead of
Huffman coding:

* **Standard QM-Coder State Machine:** Uses the 113-state adaptive binary probability
  estimation state machine (ITU-T T.81 Table D.2) with compact packed state transitions.
* **Conditioning (DAC):** Emits standard DAC (`0xFFCC`) markers defining DC conditioning
  parameters ($L=0, U=1$) and AC conditioning parameter ($K=5$).
* **Full Interoperability:** Generates standard SOF9 (`0xFFC9`) frames compliant with
  libjpeg, libjpeg-turbo, and Pillow. The decoder supports sequential arithmetic scans
  with full 8-bit and 12-bit IDCT reconstruction.
* **Tremendous Rate-Distortion Gains:**
  * `test_natural.yuv`: **−22.90% BD-rate** vs baseline Huffman (−25.16% with trellis)
  * `test_edges.yuv`: **−31.49% BD-rate** vs baseline Huffman (−31.62% with trellis)
  * `test_grad_texture.yuv`: **−20.32% BD-rate** vs baseline Huffman
  * `12-bit DCT`: −14.2% size at $q=50$, −12.3% size at $q=70$, −9.2% size at $q=90$ at bit-for-bit identical PSNR.

**Trellis & Rate-Distortion Optimizations**

Enabled with `-t` (with default Lagrange multiplier `-l 0.007`):

* **Joint DC Trellis Optimization:** A dynamic programming (Viterbi) search over the
  sequence of DC coefficients jointly minimizes spatial reconstruction distortion and the
  differential DC Huffman + VLI code rate across each component.
* **AC Trellis Quantization:** Comprehensive candidate evaluation bracketing `nearest`,
  `floor`, `ceil`, `nearest - 1`, and size-category boundaries `(1 << s) - 1`, with
  consistent Lagrangian rate scaling on EOB and ZRL symbols.
* **Two-Pass Rate Model Adaptation:** Re-estimates Huffman symbol distributions from
  trellis-quantized coefficients and runs a refinement pass with the updated rate model.
* **Unbiased DC Rounding:** Standard half-away-from-zero rounding is used for DC coefficients
  to prevent deadzone distortion on block brightness.

**Further Improvement Opportunities**

* **Integer / AAN fast DCT.** An integer or AAN-scaled transform would remove
  the remaining floating-point cost and improve numerical determinism.
* **Restart markers.** Restarts improve error resilience for both baseline and
  progressive streams.
* **Multithreading.** The per-block DCT and the statistics/entropy passes are
  embarrassingly parallel.
* **Decoder color handling.** The decoder emits raw YUV 4:2:0 only; a color
  conversion and an image container format would make it directly viewable.

**Limitations**

* Sequential baseline (8-bit), extended sequential (12-bit), progressive, arithmetic sequential (8-bit and 12-bit), and lossless (8-bit and 12-bit) JPEG output and decoding are fully supported.
* Raw YUV 4:2:0 input only (8-bit bytes or 16-bit words for 12-bit); there is no color-space conversion or file-format
  handling. Width and height must be even; edge blocks for non-MCU-aligned
  dimensions are handled by replicating the last row/column.
* Error handling is minimal, as befits an educational implementation.

**License**

This project is licensed under the 2-Clause BSD License.

**Acknowledgments**

* The code structure draws inspiration from various JPEG resources and
  tutorials, and the Huffman table generation follows the procedure in
  ITU-T T.81 (JPEG), Annex K.2.
