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
Recent improvements:

1. **Optimal Huffman tables for all four tables.** Previously only the luma DC
   table was generated from the statistics; the luma AC, chroma DC and chroma
   AC tables used the generic defaults. All four are now generated with a
   correct, length-limited (≤ 16 bits) Huffman construction based on JPEG
   Annex K.2, and the DHT segments are written dynamically to match the
   generated symbol counts.
2. **Bitstream tail fixes.** The final partial byte was written misaligned and
   the last few bits were dropped, and the EOI marker was appended without
   first padding the entropy segment to a byte boundary (so `0xFFD9` never
   appeared literally in the file). Both are fixed: the tail is flushed
   correctly, the entropy segment is padded with 1-bits, and the file now ends
   with a proper EOI marker.
3. **Separable FDCT.** The transform is computed as two 1-D passes, reducing
   work from O(N⁴) to O(2·N³) per block (4096 → 1024 multiplies).

Measured on a synthetic gradient/texture 640×480 frame (byte sizes; PSNR in
dB against the source), original vs. improved:

| Quality | Original | Improved | Size change | PSNR (orig → new) |
| ---: | ---: | ---: | ---: | --- |
| 25 | 10226 | 7786 | −23.9% | 37.41 → 37.41 |
| 50 | 12839 | 10577 | −17.6% | 39.09 → 40.56 |
| 75 | 17430 | 15000 | −13.9% | 44.91 → 44.92 |
| 90 | 27992 | 25001 | −10.7% | 50.21 → 50.22 |
| **Total** | **68487** | **58364** | **−14.8%** | — |

At 1280×720 / quality 75 the same change reduced the file from 51071 to 44499
bytes (−12.9%) with no measurable PSNR loss. The separable FDCT reduced the
720p DCT/quantization stage from ~76 ms to ~9 ms on the test machine.

**Further Improvement Opportunities**

* **Progressively optimized quantization tables.** The tables used are the
  classic JPEG tables plus a simple `(100 - quality)` scaling. Using the
  standard libjpeg scaling and validating the luma table would give a more
  predictable quality/size curve. (Note: several high-frequency entries in
  `fjpeg_default_luma_quant_table` differ from the standard Annex K table,
  e.g. 24 → 124, which makes the default table more aggressive than intended.)
* **Integer / AAN fast DCT.** An integer or AAN-scaled transform would remove
  the remaining floating-point cost and improve numerical determinism.
* **Progressive JPEG and restart markers.** Progressive scans and restarts
  improve error resilience and can improve rate-distortion at low bitrates.
* **Per-component statistics.** The statistics pass is currently scalar; it
  could be vectorized for a further speed-up.
* **Decoder.** Complete `fjpeg_read_headers` (it currently returns without a
  value on some paths) and implement entropy decoding.

**Limitations**

* Baseline sequential JPEG only; no progressive, lossless, or arithmetic
  JPEG output.
* Raw YUV 4:2:0 input only; there is no color-space conversion or file-format
  handling. Width and height should be multiples of 16 (the 4:2:0 MCU size);
  edge blocks are not padded.
* Error handling is minimal, as befits an educational implementation.

**License**

This project is licensed under the 2-Clause BSD License.

**Acknowledgments**

* The code structure draws inspiration from various JPEG resources and
  tutorials, and the Huffman table generation follows the procedure in
  ITU-T T.81 (JPEG), Annex K.2.
