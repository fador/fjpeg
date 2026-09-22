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
#include <cstdint>
#include <cmath>
#include <vector>

#include "fjpeg.h"
#include "fjpeg_bitstream.h"
#include "fjpeg_transquant.h"

// ---------------------------------------------------------------------------
// Byte/bit reader
//
// Header segments are read as plain big-endian bytes. Entropy-coded data is
// read bit by bit with 0xFF/0x00 byte unstuffing. When a marker is encountered
// while reading entropy bits, the reader records it and fills with 1-bits, as
// required by the specification; the marker is then returned by next_marker().
// ---------------------------------------------------------------------------
class fjpeg_byte_reader {
public:
    explicit fjpeg_byte_reader(FILE* f) : fp(f), cur(0), bit_cnt(0), pending_marker(-1) {}

    int read_u8() {
        int c = fgetc(fp);
        return c < 0 ? 0 : c;
    }

    int read_u16() {
        int hi = read_u8();
        int lo = read_u8();
        return (hi << 8) | lo;
    }

    void skip(size_t n) {
        for (size_t i = 0; i < n; i++) {
            if (fgetc(fp) < 0) break;
        }
    }

    int read_bit() {
        if (pending_marker >= 0) {
            return 1;
        }
        if (bit_cnt == 0) {
            int c = fgetc(fp);
            if (c < 0) {
                pending_marker = 0xFFD9;
                return 0;
            }
            if (c == 0xFF) {
                int c2 = fgetc(fp);
                while (c2 == 0xFF) {
                    c2 = fgetc(fp);
                }
                if (c2 < 0) {
                    pending_marker = 0xFFD9;
                    return 0;
                }
                if (c2 == 0x00) {
                    c = 0xFF; // stuffed byte
                } else {
                    pending_marker = (0xFF << 8) | (c2 & 0xFF);
                    return 1;
                }
            }
            cur = c;
            bit_cnt = 8;
        }
        bit_cnt--;
        return (cur >> bit_cnt) & 1;
    }

    int read_bits(int n) {
        int v = 0;
        for (int i = 0; i < n; i++) {
            v = (v << 1) | read_bit();
        }
        return v;
    }

    void align() {
        bit_cnt = 0;
    }

    // Scan forward for the next marker, resynchronising over any stray entropy
    // bytes. Returns the 16-bit marker or -1 at end of file.
    int next_marker() {
        if (pending_marker >= 0) {
            int m = pending_marker;
            pending_marker = -1;
            bit_cnt = 0;
            return m;
        }
        bit_cnt = 0;
        for (;;) {
            int c = fgetc(fp);
            if (c < 0) {
                return -1;
            }
            if (c != 0xFF) {
                continue;
            }
            do {
                c = fgetc(fp);
            } while (c == 0xFF);
            if (c < 0) {
                return -1;
            }
            if (c == 0x00) {
                continue; // stuffed byte
            }
            return (0xFF << 8) | (c & 0xFF);
        }
    }

    // Consume a restart marker at an MCU boundary.
    void restart() {
        bit_cnt = 0;
        if (pending_marker >= 0) {
            pending_marker = -1;
            return;
        }
        next_marker();
    }

private:
    FILE* fp;
    int cur;
    int bit_cnt;
    int pending_marker;
};

// ---------------------------------------------------------------------------
// Canonical Huffman decoder
// ---------------------------------------------------------------------------
struct fjpeg_huff_decoder {
    int mincode[17];
    int maxcode[18];
    int valptr[17];
    uint8_t huffval[256];
    bool valid;

    fjpeg_huff_decoder() : valid(false) {
        memset(mincode, 0, sizeof(mincode));
        memset(maxcode, 0, sizeof(maxcode));
        memset(valptr, 0, sizeof(valptr));
        memset(huffval, 0, sizeof(huffval));
    }

    void build(const uint8_t bits[16], const uint8_t* values, int count) {
        int code = 0;
        int k = 0;
        for (int i = 1; i <= 16; i++) {
            if (bits[i - 1]) {
                valptr[i] = k;
                mincode[i] = code;
                code += bits[i - 1];
                maxcode[i] = code - 1;
                k += bits[i - 1];
            } else {
                valptr[i] = 0;
                mincode[i] = 0;
                maxcode[i] = -1;
            }
            code <<= 1;
        }
        maxcode[17] = 0x7FFFFFFF;
        if (count > 256) count = 256;
        for (int i = 0; i < count; i++) {
            huffval[i] = values[i];
        }
        valid = true;
    }

    int decode(fjpeg_byte_reader& reader) const {
        int code = reader.read_bit();
        for (int i = 1; i <= 16; i++) {
            if (code <= maxcode[i]) {
                int idx = valptr[i] + code - mincode[i];
                if (idx < 0 || idx >= 256) return 0;
                return huffval[idx];
            }
            code = (code << 1) | reader.read_bit();
        }
        return 0;
    }
};

// ---------------------------------------------------------------------------
// Decoded component
// ---------------------------------------------------------------------------
struct fjpeg_dec_component {
    int id;
    int h;
    int v;
    int quant_id;
    int dc_table;
    int ac_table;
    int width;          // component width in samples
    int height;         // component height in samples
    int std_blocks_w;   // conformant (non-interleaved) block counts
    int std_blocks_h;
    int blocks_w;       // MCU-padded block grid (matches the encoder)
    int blocks_h;
    std::vector<int16_t> coeff;   // blocks_w * blocks_h * 64, zigzag order
    std::vector<uint8_t> plane;   // decoded samples, stride = blocks_w * 8
    int dc_pred;
    int eob_run;

    fjpeg_dec_component()
        : id(0), h(1), v(1), quant_id(0), dc_table(0), ac_table(0),
          width(0), height(0), std_blocks_w(0), std_blocks_h(0),
          blocks_w(0), blocks_h(0), dc_pred(0), eob_run(0) {}
};

// JPEG VLI sign extension.
static inline int fjpeg_extend(int v, int n) {
    if (n == 0) return 0;
    if (v < (1 << (n - 1))) {
        v -= (1 << n) - 1;
    }
    return v;
}

static inline int fjpeg_zigzag_scan(int k) {
    // Inverse of fjpeg_zigzag_8x8: maps a zigzag position to a natural index.
    for (int n = 0; n < 64; n++) {
        if (fjpeg_zigzag_8x8[n] == k) return n;
    }
    return k;
}

// ---------------------------------------------------------------------------
// Decoder
// ---------------------------------------------------------------------------
class fjpeg_decoder {
public:
    fjpeg_decoder()
        : reader(nullptr), width(0), height(0), max_h(1), max_v(1), num_comps(0),
          progressive(false), restart_interval(0), mcus_per_row(0), mcus_per_col(0) {
        memset(quant_tables, 0, sizeof(quant_tables));
    }

    bool decode(FILE* fp) {
        reader = new fjpeg_byte_reader(fp);

        if (reader->next_marker() != 0xFFD8) {
            fprintf(stderr, "Error: Not a JPEG file\n");
            return false;
        }

        for (;;) {
            int marker = reader->next_marker();
            if (marker < 0) break;          // end of file
            if (marker == 0xFFD9) break;    // EOI

            // Markers without a length field.
            if (marker == 0xFF01) continue; // TEM
            if (marker >= 0xFFD0 && marker <= 0xFFD7) continue; // RSTn

            int length = reader->read_u16();
            if (length < 2) {
                fprintf(stderr, "Error: Invalid marker length\n");
                return false;
            }

            switch (marker) {
                case 0xFFDB: if (!parse_dqt(length)) return false; break;
                case 0xFFC0: if (!parse_sof(length, false)) return false; break;
                case 0xFFC1: if (!parse_sof(length, false)) return false; break;
                case 0xFFC2: if (!parse_sof(length, true)) return false; break;
                case 0xFFC4: if (!parse_dht(length)) return false; break;
                case 0xFFDD: restart_interval = reader->read_u16(); break;
                case 0xFFDA: if (!parse_sos(length)) return false; break;
                default:
                    reader->skip(length - 2);
                    break;
            }
        }

        return true;
    }

    bool write_yuv(const char* filename) {
        if (num_comps < 1) return false;
        FILE* fp = fopen(filename, "wb");
        if (!fp) {
            fprintf(stderr, "Error: Unable to open output file\n");
            return false;
        }

        // Decode every block into its component plane.
        for (int ci = 0; ci < num_comps; ci++) {
            fjpeg_dec_component& c = comps[ci];
            const uint8_t* q = quant_tables[c.quant_id];
            int stride = c.blocks_w * 8;
            for (int by = 0; by < c.blocks_h; by++) {
                for (int bx = 0; bx < c.blocks_w; bx++) {
                    const int16_t* zz = &c.coeff[(by * c.blocks_w + bx) * 64];
                    fjpeg_coeff_t zz_float[64];
                    fjpeg_coeff_t natural[64];
                    fjpeg_coeff_t dq[64];
                    for (int i = 0; i < 64; i++) {
                        zz_float[i] = (fjpeg_coeff_t)zz[i];
                    }
                    fjpeg_izigzag8x8(zz_float, natural);
                    for (int i = 0; i < 64; i++) {
                        dq[i] = natural[i] * (fjpeg_coeff_t)q[i];
                    }
                    fjpeg_pixel_t pixels[64];
                    fjpeg_idct8x8(&ctx, dq, pixels);
                    uint8_t* dst = &c.plane[(by * 8) * stride + bx * 8];
                    for (int j = 0; j < 8; j++) {
                        memcpy(dst + j * stride, pixels + j * 8, 8);
                    }
                }
            }
        }

        // Write the planes as raw YUV, each cropped to its true dimensions.
        for (int ci = 0; ci < num_comps; ci++) {
            fjpeg_dec_component& c = comps[ci];
            int stride = c.blocks_w * 8;
            for (int y = 0; y < c.height; y++) {
                fwrite(&c.plane[(size_t)y * stride], 1, c.width, fp);
            }
        }

        fclose(fp);
        return true;
    }

    int get_width() const { return width; }
    int get_height() const { return height; }
    int get_num_components() const { return num_comps; }

private:
    int find_component(int id) const {
        for (int i = 0; i < num_comps; i++) {
            if (comps[i].id == id) return i;
        }
        return -1;
    }

    bool parse_dqt(int length) {
        int remaining = length - 2;
        while (remaining > 0) {
            int pq_tq = reader->read_u8();
            remaining--;
            int pq = pq_tq >> 4;
            int tq = pq_tq & 0x0F;
            if (tq > 3) return false;

            if (pq == 0) {
                if (remaining < 64) return false;
                for (int i = 0; i < 64; i++) {
                    int v = reader->read_u8();
                    quant_tables[tq][fjpeg_zigzag_scan(i)] = (uint8_t)v;
                }
                remaining -= 64;
            } else {
                if (remaining < 128) return false;
                for (int i = 0; i < 64; i++) {
                    int hi = reader->read_u8();
                    int lo = reader->read_u8();
                    int v = (hi << 8) | lo;
                    quant_tables[tq][fjpeg_zigzag_scan(i)] = (uint8_t)(v > 255 ? 255 : v);
                }
                remaining -= 128;
            }
        }
        return true;
    }

    bool parse_dht(int length) {
        int remaining = length - 2;
        while (remaining > 0) {
            int tc_th = reader->read_u8();
            remaining--;
            int tc = tc_th >> 4;
            int th = tc_th & 0x0F;
            if (th > 3 || tc > 1) return false;
            if (remaining < 16) return false;

            uint8_t bits[16];
            int total = 0;
            for (int i = 0; i < 16; i++) {
                bits[i] = (uint8_t)reader->read_u8();
                total += bits[i];
            }
            remaining -= 16;
            if (remaining < total) return false;

            uint8_t values[256];
            int n = total > 256 ? 256 : total;
            for (int i = 0; i < n; i++) {
                values[i] = (uint8_t)reader->read_u8();
            }
            // Consume any excess (malformed) entries.
            for (int i = n; i < total; i++) {
                reader->read_u8();
            }
            remaining -= total;

            if (tc == 0) {
                dc_tables[th].build(bits, values, n);
            } else {
                ac_tables[th].build(bits, values, n);
            }
        }
        return true;
    }

    bool parse_sof(int length, bool prog) {
        (void)length;
        progressive = prog;
        int precision = reader->read_u8();
        if (precision != 8) {
            fprintf(stderr, "Error: Only 8-bit samples are supported\n");
            return false;
        }
        height = reader->read_u16();
        width = reader->read_u16();
        num_comps = reader->read_u8();
        if (num_comps < 1 || num_comps > 4) {
            fprintf(stderr, "Error: Unsupported component count %d\n", num_comps);
            return false;
        }

        comps.resize(num_comps);
        max_h = 1;
        max_v = 1;
        for (int i = 0; i < num_comps; i++) {
            comps[i].id = reader->read_u8();
            int hv = reader->read_u8();
            comps[i].h = (hv >> 4) & 0x0F;
            comps[i].v = hv & 0x0F;
            comps[i].quant_id = reader->read_u8();
            if (comps[i].h < 1) comps[i].h = 1;
            if (comps[i].v < 1) comps[i].v = 1;
            if (comps[i].quant_id > 3) comps[i].quant_id = 0;
            max_h = FJPEG_MAX(max_h, comps[i].h);
            max_v = FJPEG_MAX(max_v, comps[i].v);
        }

        mcus_per_row = (width + 8 * max_h - 1) / (8 * max_h);
        mcus_per_col = (height + 8 * max_v - 1) / (8 * max_v);

        for (int i = 0; i < num_comps; i++) {
            fjpeg_dec_component& c = comps[i];
            c.width = (width * c.h + max_h - 1) / max_h;
            c.height = (height * c.v + max_v - 1) / max_v;
            c.std_blocks_w = (c.width + 7) / 8;
            c.std_blocks_h = (c.height + 7) / 8;
            c.blocks_w = mcus_per_row * c.h;
            c.blocks_h = mcus_per_col * c.v;
            c.coeff.assign((size_t)c.blocks_w * c.blocks_h * 64, 0);
            c.plane.assign((size_t)c.blocks_w * 8 * c.blocks_h * 8, 0);
            c.dc_pred = 0;
            c.eob_run = 0;
        }
        return true;
    }

    bool parse_sos(int length) {
        (void)length;
        int ns = reader->read_u8();
        if (ns < 1 || ns > 4) {
            fprintf(stderr, "Error: Invalid scan component count\n");
            return false;
        }

        int scan_comp[4];
        int dc_sel[4];
        int ac_sel[4];
        for (int i = 0; i < ns; i++) {
            int cs = reader->read_u8();
            int t = reader->read_u8();
            int idx = find_component(cs);
            if (idx < 0) {
                fprintf(stderr, "Error: Unknown component id %d\n", cs);
                return false;
            }
            scan_comp[i] = idx;
            dc_sel[i] = (t >> 4) & 0x0F;
            ac_sel[i] = t & 0x0F;
            comps[idx].dc_table = dc_sel[i];
            comps[idx].ac_table = ac_sel[i];
            comps[idx].dc_pred = 0;
            comps[idx].eob_run = 0;
        }

        int ss = reader->read_u8();
        int se = reader->read_u8();
        int a = reader->read_u8();
        int ah = (a >> 4) & 0x0F;
        int al = a & 0x0F;

        if (!progressive) {
            return decode_sequential(scan_comp, dc_sel, ac_sel, ns);
        }
        if (ss == 0) {
            return decode_progressive_dc(scan_comp, dc_sel, ns, ah, al);
        }
        return decode_progressive_ac(scan_comp, ns, ss, se, ah, al);
    }

    // -----------------------------------------------------------------------
    // Sequential (baseline) decoding
    // -----------------------------------------------------------------------
    bool decode_sequential(int* scan_comp, int* dc_sel, int* ac_sel, int ns) {
        if (ns > 1) {
            int mcu = 0;
            for (int my = 0; my < mcus_per_col; my++) {
                for (int mx = 0; mx < mcus_per_row; mx++) {
                    if (restart_interval > 0 && mcu > 0 && (mcu % restart_interval) == 0) {
                        reader->restart();
                        for (int s = 0; s < ns; s++) comps[scan_comp[s]].dc_pred = 0;
                    }
                    for (int s = 0; s < ns; s++) {
                        fjpeg_dec_component& c = comps[scan_comp[s]];
                        for (int v = 0; v < c.v; v++) {
                            for (int h = 0; h < c.h; h++) {
                                decode_block_seq(c, mx * c.h + h, my * c.v + v,
                                                 dc_sel[s], ac_sel[s]);
                            }
                        }
                    }
                    mcu++;
                }
            }
        } else {
            fjpeg_dec_component& c = comps[scan_comp[0]];
            int mcu = 0;
            for (int by = 0; by < c.std_blocks_h; by++) {
                for (int bx = 0; bx < c.std_blocks_w; bx++) {
                    if (restart_interval > 0 && mcu > 0 && (mcu % restart_interval) == 0) {
                        reader->restart();
                        c.dc_pred = 0;
                    }
                    decode_block_seq(c, bx, by, dc_sel[0], ac_sel[0]);
                    mcu++;
                }
            }
        }
        return true;
    }

    void decode_block_seq(fjpeg_dec_component& c, int bx, int by, int dc_sel, int ac_sel) {
        int16_t* blk = &c.coeff[((size_t)by * c.blocks_w + bx) * 64];
        memset(blk, 0, 64 * sizeof(int16_t));

        int t = dc_tables[dc_sel].decode(*reader);
        int diff = 0;
        if (t) {
            diff = fjpeg_extend(reader->read_bits(t), t);
        }
        c.dc_pred += diff;
        blk[0] = (int16_t)c.dc_pred;

        int k = 1;
        while (k < 64) {
            int rs = ac_tables[ac_sel].decode(*reader);
            int r = rs >> 4;
            int s = rs & 0x0F;
            if (s == 0) {
                if (r == 15) {
                    k += 16;
                    continue;
                }
                break; // EOB
            }
            k += r;
            if (k > 63) break;
            blk[k] = (int16_t)fjpeg_extend(reader->read_bits(s), s);
            k++;
        }
    }

    // -----------------------------------------------------------------------
    // Progressive decoding
    // -----------------------------------------------------------------------
    bool decode_progressive_dc(int* scan_comp, int* dc_sel, int ns, int ah, int al) {
        if (ah == 0) {
            for_each_block(scan_comp, ns, [&](fjpeg_dec_component& c, int bx, int by, int table_sel) {
                int16_t* blk = &c.coeff[((size_t)by * c.blocks_w + bx) * 64];
                int t = dc_tables[table_sel].decode(*reader);
                int diff = 0;
                if (t) {
                    diff = fjpeg_extend(reader->read_bits(t), t);
                }
                c.dc_pred += diff;
                blk[0] = (int16_t)(c.dc_pred << al);
            }, dc_sel);
        } else {
            for_each_block(scan_comp, ns, [&](fjpeg_dec_component& c, int bx, int by, int table_sel) {
                (void)table_sel;
                int16_t* blk = &c.coeff[((size_t)by * c.blocks_w + bx) * 64];
                if (reader->read_bit()) {
                    blk[0] |= (int16_t)(1 << al);
                }
            }, dc_sel);
        }
        return true;
    }

    bool decode_progressive_ac(int* scan_comp, int ns, int ss, int se, int ah, int al) {
        if (ns != 1) {
            fprintf(stderr, "Error: Progressive AC scans must be non-interleaved\n");
            return false;
        }
        fjpeg_dec_component& c = comps[scan_comp[0]];
        if (ah == 0) {
            decode_ac_first(c, ss, se, al);
        } else {
            decode_ac_refine(c, ss, se, al);
        }
        return true;
    }

    void decode_ac_first(fjpeg_dec_component& c, int ss, int se, int al) {
        int table_sel = c.ac_table;
        for (int by = 0; by < c.std_blocks_h; by++) {
            for (int bx = 0; bx < c.std_blocks_w; bx++) {
                if (c.eob_run > 0) {
                    c.eob_run--;
                    continue;
                }
                int16_t* blk = &c.coeff[((size_t)by * c.blocks_w + bx) * 64];
                int k = ss;
                while (k <= se) {
                    int rs = ac_tables[table_sel].decode(*reader);
                    int r = rs >> 4;
                    int s = rs & 0x0F;
                    if (s) {
                        k += r;
                        if (k > se) break;
                        blk[k] = (int16_t)(fjpeg_extend(reader->read_bits(s), s) << al);
                        k++;
                    } else {
                        if (r == 15) {
                            k += 16; // ZRL
                        } else {
                            int run = 1 << r;
                            if (r) run += reader->read_bits(r);
                            c.eob_run = run - 1;
                            break;
                        }
                    }
                }
            }
        }
    }

    void decode_ac_refine(fjpeg_dec_component& c, int ss, int se, int al) {
        int table_sel = c.ac_table;
        const int p1 = 1 << al;
        const int m1 = -p1;

        for (int by = 0; by < c.std_blocks_h; by++) {
            for (int bx = 0; bx < c.std_blocks_w; bx++) {
                int16_t* blk = &c.coeff[((size_t)by * c.blocks_w + bx) * 64];
                int k = ss;

                if (c.eob_run == 0) {
                    for (; k <= se; k++) {
                        int rs = ac_tables[table_sel].decode(*reader);
                        int r = rs >> 4;
                        int s = rs & 0x0F;
                        if (s) {
                            // A newly non-zero coefficient, always +/- 1 * p1.
                            int newval = reader->read_bit() ? p1 : m1;
                            do {
                                if (blk[k] != 0) {
                                    if (reader->read_bit()) {
                                        if ((blk[k] & p1) == 0) {
                                            blk[k] += (blk[k] >= 0) ? p1 : m1;
                                        }
                                    }
                                } else {
                                    if (--r < 0) break;
                                }
                                k++;
                            } while (k <= se);
                            if (k > se) break;
                            blk[k] = (int16_t)newval;
                            continue; // for loop advances k
                        }
                        if (r != 15) {
                            int run = 1 << r;
                            if (r) run += reader->read_bits(r);
                            c.eob_run = run;
                            break; // EOB: rest handled below
                        }
                        // ZRL: skip 16 zero coefficients while refining non-zeros.
                        do {
                            if (blk[k] != 0) {
                                if (reader->read_bit()) {
                                    if ((blk[k] & p1) == 0) {
                                        blk[k] += (blk[k] >= 0) ? p1 : m1;
                                    }
                                }
                            } else {
                                if (--r < 0) break;
                            }
                            k++;
                        } while (k <= se);
                    }
                }

                if (c.eob_run > 0) {
                    for (; k <= se; k++) {
                        if (blk[k] != 0) {
                            if (reader->read_bit()) {
                                if ((blk[k] & p1) == 0) {
                                    blk[k] += (blk[k] >= 0) ? p1 : m1;
                                }
                            }
                        }
                    }
                    c.eob_run--;
                }
            }
        }
    }

    // Visit every block of a scan. Interleaved scans follow the MCU grid; a
    // non-interleaved scan visits the component's own conformant block grid.
    template <typename F>
    void for_each_block(int* scan_comp, int ns, F fn, int* table_sel) {
        if (ns > 1) {
            for (int my = 0; my < mcus_per_col; my++) {
                for (int mx = 0; mx < mcus_per_row; mx++) {
                    for (int s = 0; s < ns; s++) {
                        fjpeg_dec_component& c = comps[scan_comp[s]];
                        for (int v = 0; v < c.v; v++) {
                            for (int h = 0; h < c.h; h++) {
                                fn(c, mx * c.h + h, my * c.v + v, table_sel[s]);
                            }
                        }
                    }
                }
            }
        } else {
            fjpeg_dec_component& c = comps[scan_comp[0]];
            for (int by = 0; by < c.std_blocks_h; by++) {
                for (int bx = 0; bx < c.std_blocks_w; bx++) {
                    fn(c, bx, by, table_sel[0]);
                }
            }
        }
    }

    fjpeg_byte_reader* reader;
    fjpeg_context ctx;
    std::vector<fjpeg_dec_component> comps;
    fjpeg_huff_decoder dc_tables[4];
    fjpeg_huff_decoder ac_tables[4];
    uint8_t quant_tables[4][64];
    int width;
    int height;
    int max_h;
    int max_v;
    int num_comps;
    bool progressive;
    int restart_interval;
    int mcus_per_row;
    int mcus_per_col;
};

bool fjpeg_decode_file(const char* input_filename, const char* output_filename) {
    FILE* fp = fopen(input_filename, "rb");
    if (!fp) {
        fprintf(stderr, "Error: Unable to open input file\n");
        return false;
    }

    fjpeg_decoder decoder;
    bool ok = decoder.decode(fp);
    fclose(fp);

    if (!ok) {
        return false;
    }

    printf("Decoded %dx%d, %d component(s)\n",
           decoder.get_width(), decoder.get_height(), decoder.get_num_components());

    return decoder.write_yuv(output_filename);
}
