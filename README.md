**FJPEG - A JPEG Encoder in C++**

**Overview**

FJPEG is a compact JPEG encoder implementation in C++, designed for educational purposes and experimentation. 

**Key Features**

* Core JPEG encoding functionalities:
    * DCT (Discrete Cosine Transform)
    * Quantization
    * Huffman encoding
* Supports baseline JPEG compression
* Educational focus: Clear code structure for understanding JPEG principles.

**Usage**

1. **Compile:**
   ```bash
   cmake .
   make
   ```

2. **Run:**
   ```bash
   ./fjpeg
   ```
   This will generate a "test.jpg" file in the current directory based on a default YUV input.

**Understanding the Code**

The code is well-commented, making it suitable for learning about JPEG compression techniques. Key sections include:

* **DCT and IDCT:** Functions for performing the forward and inverse Discrete Cosine Transforms.
* **Quantization and Dequantization:** Functions to apply quantization and dequantization based on quality settings.
* **Huffman Coding:** Functions to encode DCT coefficients using Huffman tables.
* **Bitstream Handling:**  A `fjpeg_bitstream` struct manages writing bits to the output file.
* **Header Generation:**  The `fjpeg_generate_header` function creates the JPEG file header.

**Limitations**

* **Basic Features:**  FJPEG implements baseline JPEG compression, suitable for learning, but lacks advanced features found in full-fledged JPEG libraries.
* **Error Handling:** Robust error handling is minimal in this educational implementation.

**License**

This project is licensed under the 2-Clause BSD License.

**Acknowledgments**

* The code structure draws inspiration from various JPEG resources and tutorials.

