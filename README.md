**FJPEG**

**Overview**

FJPEG is a compact JPEG encoder written in C++, designed for educational
purposes and experimentation. It reads raw YUV 4:2:0 input and produces a
baseline JPEG file, with a clear, self-contained implementation of every
stage of the pipeline.

**Key Features**

* Baseline (sequential, 8-bit) JPEG encoding
* Full pipeline: raw YUV input → FDCT → quantization → zigzag → run-length
  coding → Huffman entropy coding → bitstream with 0xFF byte stuffing
* Per-image optimal Huffman tables (DC and AC, luma and chroma), generated
  from a first-pass statistics scan
* Standard ITU-T T.81 quantization tables with libjpeg-compatible quality
  scaling, plus a tunable quantization deadzone
* 4:2:0 chroma subsampling and a separable (2-pass) FDCT
* Optional arithmetic-coding experiment (`fjpeg_arith`)
* Header parsing / decoding is under development

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
```

Options:

```
-i <input_filename>    input raw YUV 4:2:0 file
-r <width>x<height>    frame resolution
-q <quality>           quality factor (1-100)
-o <output_filename>   output JPEG file
-d                     decode an existing JPEG (work in progress)
-h                     show help
```

**Source Layout**

| File | Purpose |
| --- | --- |
| `src/fjpeg.cpp` | JPEG header generation/parsing, encode loop, CLI usage text |
| `src/fjpeg.h` | `fjpeg_context`, quality scaling, input loading |
| `src/fjpeg_transquant.cpp` | block extraction, FDCT/IDCT, quantization, zigzag |
| `src/fjpeg_huffman.cpp` | statistics pass, optimal Huffman generation, entropy coding |
| `src/fjpeg_bitstream.h` | bit reader/writer, 0xFF stuffing, file flushing |
| `src/fjpeg_global.h` | types, default quantization tables, zigzag tables |
| `src/fjpeg_huffman.h` | default Huffman tables and statistics struct |
| `src/fjpeg_cli.cpp` | command line parsing and program flow |
| `src/fjpeg_arith.cpp` | experimental arithmetic coder (standalone test) |

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

**Further Improvement Opportunities**

* **Rate-distortion optimized quantization.** A trellis search over the AC
  coefficients that accounts for run-length and code cost typically gains a few
  more percent.
* **Integer / AAN fast DCT.** An integer or AAN-scaled transform would remove
  the remaining floating-point cost and improve numerical determinism.
* **Progressive JPEG and restart markers.** Progressive scans and restarts
  improve error resilience and can improve rate-distortion at low bitrates.
* **Multithreading.** The per-block DCT and the statistics/entropy passes are
  embarrassingly parallel.
* **Decoder.** Complete `fjpeg_read_headers` (it currently returns without a
  value on some paths) and implement entropy decoding.

**Limitations**

* Baseline sequential JPEG only; no progressive, lossless, or arithmetic
  JPEG output.
* Raw YUV 4:2:0 input only; there is no color-space conversion or file-format
  handling. Width and height must be even; edge blocks for non-MCU-aligned
  dimensions are handled by replicating the last row/column.
* Error handling is minimal, as befits an educational implementation.

**License**

This project is licensed under the 2-Clause BSD License.

**Acknowledgments**

* The code structure draws inspiration from various JPEG resources and
  tutorials, and the Huffman table generation follows the procedure in
  ITU-T T.81 (JPEG), Annex K.2.
