#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Zigzag 掃描表
const int zigzag_index[64] = {
     0,  1,  8, 16,  9,  2,  3, 10,
    17, 24, 32, 25, 18, 11,  4,  5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13,  6,  7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63
};

// --- 定義 BMP 結構 ---
#pragma pack(push, 1)
typedef struct {
    uint16_t type; uint32_t size; uint16_t res1, res2; uint32_t offset;
} BMPHeader;
typedef struct {
    uint32_t size; int32_t width, height; uint16_t planes, bpp;
    uint32_t compression, image_size; int32_t x_ppm, y_ppm;
    uint32_t colors_used, colors_important;
} BMPInfoHeader;
#pragma pack(pop)

// --- 標準量化表 ---
const int Q_Luminance[8][8] = {
    {16, 11, 10, 16, 24, 40, 51, 61}, {12, 12, 14, 19, 26, 58, 60, 55},
    {14, 13, 16, 24, 40, 57, 69, 56}, {14, 17, 22, 29, 51, 87, 80, 62},
    {18, 22, 37, 56, 68, 109,103, 77}, {24, 35, 55, 64, 81, 104,113, 92},
    {49, 64, 78, 87, 103,121,120,101}, {72, 92, 95, 98, 112,100,103, 99}
};

const int Q_Chrominance[8][8] = {
    {17, 18, 24, 47, 99, 99, 99, 99}, {18, 21, 26, 66, 99, 99, 99, 99},
    {24, 26, 56, 99, 99, 99, 99, 99}, {47, 66, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99}, {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99}, {99, 99, 99, 99, 99, 99, 99, 99}
};

// --- Mode 3 Helper Structures ---
typedef struct {
    FILE *fp;
    uint8_t buffer;
    int bit_count;
} BitWriter;

typedef struct HNode {
    int symbol;
    int freq;
    struct HNode *left, *right;
} HNode;

typedef struct {
    int len;
    uint32_t code;
    // ★★★ 關鍵修正：加大到 256 以支援超大圖片造成的深層 Huffman Tree ★★★
    char str[256]; 
} HCode;

// Global Variables
HCode global_codebook[256];
int global_freq[256];

// --- Helper Functions ---

void init_bitwriter(BitWriter *bw, FILE *fp) {
    bw->fp = fp;
    bw->buffer = 0;
    bw->bit_count = 0;
}

void write_bits(BitWriter *bw, uint32_t bits, int len) {
    for (int i = len - 1; i >= 0; i--) {
        int bit = (bits >> i) & 1;
        bw->buffer = (bw->buffer << 1) | bit;
        bw->bit_count++;
        if (bw->bit_count == 8) {
            fwrite(&bw->buffer, 1, 1, bw->fp);
            bw->buffer = 0;
            bw->bit_count = 0;
        }
    }
}

void flush_bits(BitWriter *bw) {
    if (bw->bit_count > 0) {
        bw->buffer = bw->buffer << (8 - bw->bit_count);
        fwrite(&bw->buffer, 1, 1, bw->fp);
    }
}

void rgb_to_ycbcr(uint8_t r, uint8_t g, uint8_t b, float *y, float *cb, float *cr) {
    *y  =  0.299f * r + 0.587f * g + 0.114f * b;
    *cb = -0.1687f * r - 0.3313f * g + 0.5f * b + 128.0f;
    *cr =  0.5f * r - 0.4187f * g - 0.0813f * b + 128.0f;
}

void perform_dct(float input[8][8], float output[8][8]) {
    for (int u = 0; u < 8; u++) {
        for (int v = 0; v < 8; v++) {
            float sum = 0.0f;
            float cu = (u == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
            float cv = (v == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
            for (int x = 0; x < 8; x++) {
                for (int y = 0; y < 8; y++) {
                    sum += input[x][y] * cosf((2 * x + 1) * u * M_PI / 16.0f) * cosf((2 * y + 1) * v * M_PI / 16.0f);
                }
            }
            output[u][v] = 0.25f * cu * cv * sum;
        }
    }
}

int get_category(int16_t val) {
    if (val == 0) return 0;
    val = abs(val);
    int cat = 0;
    while (val > 0) { val >>= 1; cat++; }
    return cat;
}

uint32_t get_vli(int16_t val, int cat) {
    if (val > 0) return val;
    return val + (1 << cat) - 1;
}

void generate_codes_recursive(HNode *node, uint32_t code, int len) {
    if (!node) return;
    if (!node->left && !node->right) {
        global_codebook[node->symbol].len = len;
        global_codebook[node->symbol].code = code;
        for(int i=0; i<len; i++) 
            global_codebook[node->symbol].str[i] = ((code >> (len-1-i)) & 1) ? '1' : '0';
        global_codebook[node->symbol].str[len] = '\0';
        return;
    }
    generate_codes_recursive(node->left, (code << 1) | 0, len + 1);
    generate_codes_recursive(node->right, (code << 1) | 1, len + 1);
}

void build_huffman_tree() {
    HNode *pool[256];
    int count = 0;
    for (int i = 0; i < 256; i++) {
        if (global_freq[i] > 0) {
            HNode *n = malloc(sizeof(HNode));
            n->symbol = i;
            n->freq = global_freq[i];
            n->left = n->right = NULL;
            pool[count++] = n;
        }
    }
    while (count > 1) {
        int m1 = -1, m2 = -1;
        for (int i = 0; i < count; i++) {
            if (m1 == -1 || pool[i]->freq < pool[m1]->freq) {
                m2 = m1; m1 = i;
            } else if (m2 == -1 || pool[i]->freq < pool[m2]->freq) {
                m2 = i;
            }
        }
        HNode *parent = malloc(sizeof(HNode));
        parent->symbol = -1;
        parent->freq = pool[m1]->freq + pool[m2]->freq;
        parent->left = pool[m1];
        parent->right = pool[m2];
        pool[m1] = parent;
        pool[m2] = pool[count - 1];
        count--;
    }
    if (count > 0) generate_codes_recursive(pool[0], 0, 0);
}

void write_qt_txt(const char *filename, const int qt[8][8]) {
    FILE *fp = fopen(filename, "w");
    if(!fp) return;
    for(int i=0; i<8; i++) {
        for(int j=0; j<8; j++) fprintf(fp, "%d%c", qt[i][j], (j==7?'\n':' '));
    }
    fclose(fp);
}

// --- Mode 0 Encoder ---
void mode_0_encoder(char *input_bmp, char *r_txt, char *g_txt, char *b_txt, char *dim_txt) {
    FILE *fp = fopen(input_bmp, "rb");
    if (!fp) { printf("Error opening input BMP\n"); exit(1); }
    BMPHeader header; BMPInfoHeader info;
    fread(&header, sizeof(BMPHeader), 1, fp); fread(&info, sizeof(BMPInfoHeader), 1, fp);
    int width = info.width, height = abs(info.height); 
    int padding = (4 - (width * 3) % 4) % 4; 
    uint8_t **R = (uint8_t **)malloc(height * sizeof(uint8_t *));
    uint8_t **G = (uint8_t **)malloc(height * sizeof(uint8_t *));
    uint8_t **B = (uint8_t **)malloc(height * sizeof(uint8_t *));
    for (int i = 0; i < height; i++) {
        R[i] = malloc(width); G[i] = malloc(width); B[i] = malloc(width);
    }
    fseek(fp, header.offset, SEEK_SET);
    for (int i = 0; i < height; i++) {
        int target_row = (info.height > 0) ? (height - 1 - i) : i;
        for (int j = 0; j < width; j++) {
            uint8_t bgr[3]; fread(bgr, 3, 1, fp);
            B[target_row][j] = bgr[0]; G[target_row][j] = bgr[1]; R[target_row][j] = bgr[2];
        }
        fseek(fp, padding, SEEK_CUR); 
    }
    fclose(fp);
    FILE *fr = fopen(r_txt, "w"), *fg = fopen(g_txt, "w"), *fb = fopen(b_txt, "w");
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            fprintf(fr, "%d%c", R[i][j], (j == width - 1 ? '\n' : ' '));
            fprintf(fg, "%d%c", G[i][j], (j == width - 1 ? '\n' : ' '));
            fprintf(fb, "%d%c", B[i][j], (j == width - 1 ? '\n' : ' '));
        }
    }
    fclose(fr); fclose(fg); fclose(fb);
    FILE *fd = fopen(dim_txt, "w"); fprintf(fd, "%d %d\n", width, height); fclose(fd);
    for (int i = 0; i < height; i++) { free(R[i]); free(G[i]); free(B[i]); }
    free(R); free(G); free(B);
    printf("Mode 0: BMP to TXT conversion completed.\n");
}

// --- Mode 1 Encoder ---
void mode_1_encoder(int argc, char *argv[]) {
    FILE *fp_in = fopen(argv[2], "rb");
    if (!fp_in) exit(1);
    BMPHeader header; BMPInfoHeader info;
    fread(&header, sizeof(BMPHeader), 1, fp_in); fread(&info, sizeof(BMPInfoHeader), 1, fp_in);
    int width = info.width, height = abs(info.height);
    int padding = (4 - (width * 3) % 4) % 4;
    // 使用 malloc 處理大檔案
    uint8_t *rgb_data = (uint8_t*)malloc(width * height * 3);
    if (!rgb_data) { printf("Memory allocation failed!\n"); exit(1); }

    fseek(fp_in, header.offset, SEEK_SET);
    for (int i = 0; i < height; i++) {
        int row = (info.height > 0) ? (height - 1 - i) : i;
        fread(rgb_data + row * width * 3, 3, width, fp_in);
        fseek(fp_in, padding, SEEK_CUR);
    }
    fclose(fp_in);
    FILE *f_dim = fopen(argv[6], "w"); fprintf(f_dim, "%d %d\n", width, height); fclose(f_dim);
    write_qt_txt(argv[3], Q_Luminance); write_qt_txt(argv[4], Q_Chrominance); write_qt_txt(argv[5], Q_Chrominance);
    FILE *f_qY = fopen(argv[7], "wb"), *f_qCb = fopen(argv[8], "wb"), *f_qCr = fopen(argv[9], "wb");
    FILE *f_eY = fopen(argv[10], "wb"), *f_eCb = fopen(argv[11], "wb"), *f_eCr = fopen(argv[12], "wb");
    
    // 使用 calloc 避免 stack overflow，改用 heap
    double *signal_pow = calloc(3 * 64, sizeof(double));
    double *noise_pow = calloc(3 * 64, sizeof(double));

    for (int r = 0; r < height; r += 8) {
        for (int c = 0; c < width; c += 8) {
            float blk[3][8][8];
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 8; j++) {
                    int idx = ((r + i < height ? r+i : height-1) * width + (c + j < width ? c+j : width-1)) * 3;
                    float y, cb, cr;
                    rgb_to_ycbcr(rgb_data[idx+2], rgb_data[idx+1], rgb_data[idx], &y, &cb, &cr);
                    blk[0][i][j] = y - 128.0f; blk[1][i][j] = cb - 128.0f; blk[2][i][j] = cr - 128.0f;
                }
            }
            FILE *f_q[3] = {f_qY, f_qCb, f_qCr};
            FILE *f_e[3] = {f_eY, f_eCb, f_eCr};
            const int (*qts[3])[8] = {Q_Luminance, Q_Chrominance, Q_Chrominance};
            for (int ch = 0; ch < 3; ch++) {
                float dct_out[8][8];
                perform_dct(blk[ch], dct_out);
                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 8; j++) {
                        int k = i * 8 + j;
                        int16_t qF = (int16_t)roundf(dct_out[i][j] / qts[ch][i][j]);
                        float eF = dct_out[i][j] - (qF * qts[ch][i][j]);
                        fwrite(&qF, sizeof(int16_t), 1, f_q[ch]);
                        fwrite(&eF, sizeof(float), 1, f_e[ch]);
                        signal_pow[ch * 64 + k] += dct_out[i][j] * dct_out[i][j];
                        noise_pow[ch * 64 + k] += eF * eF;
                    }
                }
            }
        }
    }
    fclose(f_qY); fclose(f_qCb); fclose(f_qCr); fclose(f_eY); fclose(f_eCb); fclose(f_eCr); free(rgb_data);
    
    char *ch_names[3] = {"Y", "Cb", "Cr"};
    for (int ch = 0; ch < 3; ch++) {
        printf("SQNR (dB) for Channel %s:\n", ch_names[ch]);
        for (int i = 0; i < 64; i++) {
            double s = signal_pow[ch * 64 + i], n = noise_pow[ch * 64 + i];
            printf("%6.2f ", (n == 0) ? 99.99 : (s == 0 ? 0.0 : 10.0 * log10(s / n)));
            if ((i + 1) % 8 == 0) printf("\n");
        }
        printf("\n");
    }
    free(signal_pow); free(noise_pow);
}

// --- Mode 2 Encoder ---
void mode_2_encoder(int argc, char *argv[]) {
    int is_binary = (strcmp(argv[3], "binary") == 0);
    FILE *fp_in = fopen(argv[2], "rb"); if(!fp_in) exit(1);
    BMPHeader header; BMPInfoHeader info;
    fread(&header, sizeof(header), 1, fp_in); fread(&info, sizeof(info), 1, fp_in);
    int width = info.width, height = abs(info.height);
    int padding = (4 - (width * 3) % 4) % 4;
    uint8_t *rgb_data = malloc(width * height * 3);
    fseek(fp_in, header.offset, SEEK_SET);
    for(int i=0; i<height; i++) {
        int row = (info.height>0) ? (height-1-i) : i;
        fread(rgb_data + row*width*3, 3, width, fp_in);
        fseek(fp_in, padding, SEEK_CUR);
    }
    fclose(fp_in);

    FILE *fout = fopen(argv[4], is_binary ? "wb" : "w");
    if (is_binary) { fwrite(&width, sizeof(int), 1, fout); fwrite(&height, sizeof(int), 1, fout); }
    else fprintf(fout, "%d %d\n", width, height);
    
    int prev_dc[3] = {0};
    for (int r = 0; r < height; r += 8) {
        for (int c = 0; c < width; c += 8) {
            float blk[3][8][8];
            for (int i=0; i<8; i++) {
                for (int j=0; j<8; j++) {
                    int idx = ((r+i < height ? r+i : height-1)*width + (c+j < width ? c+j : width-1))*3;
                    float y, cb, cr;
                    rgb_to_ycbcr(rgb_data[idx+2], rgb_data[idx+1], rgb_data[idx], &y, &cb, &cr);
                    blk[0][i][j] = y - 128.0f; blk[1][i][j] = cb - 128.0f; blk[2][i][j] = cr - 128.0f;
                }
            }
            const int (*qts[3])[8] = {Q_Luminance, Q_Chrominance, Q_Chrominance};
            for (int ch = 0; ch < 3; ch++) {
                float dct[8][8]; perform_dct(blk[ch], dct);
                int16_t zz[64];
                for (int k = 0; k < 64; k++) {
                    int u = zigzag_index[k] / 8; int v = zigzag_index[k] % 8;
                    zz[k] = (int16_t)roundf(dct[u][v] / qts[ch][u][v]);
                }
                int16_t diff = zz[0] - prev_dc[ch]; prev_dc[ch] = zz[0];
                if (is_binary) fwrite(&diff, sizeof(int16_t), 1, fout);
                else fprintf(fout, "(%d,%d,%d) 0 %d ", r/8, c/8, ch, diff);
                
                int zero_run = 0;
                for (int k = 1; k < 64; k++) {
                    if (zz[k] == 0) zero_run++;
                    else {
                        if (is_binary) { 
                            uint8_t run = zero_run; fwrite(&run, 1, 1, fout); fwrite(&zz[k], 2, 1, fout); 
                        } else fprintf(fout, "%d %d ", zero_run, zz[k]);
                        zero_run = 0;
                    }
                }
                if (is_binary) { uint8_t z=0; int16_t v=0; fwrite(&z,1,1,fout); fwrite(&v,2,1,fout); }
                else fprintf(fout, "\n");
            }
        }
    }
    long sz = ftell(fout); fclose(fout); free(rgb_data);
    printf("Mode 2 Encoded. Size: %ld\n", sz);
}

// --- Mode 3 Encoder ---
void mode_3_encoder(int argc, char *argv[]) {
    int is_binary = (strcmp(argv[3], "binary") == 0);
    char *cb_file = argv[4];
    char *hf_file = argv[5];

    FILE *fp_in = fopen(argv[2], "rb");
    if(!fp_in) exit(1);
    BMPHeader header; BMPInfoHeader info;
    fread(&header, sizeof(header), 1, fp_in); fread(&info, sizeof(info), 1, fp_in);
    int width = info.width, height = abs(info.height);
    int padding = (4 - (width * 3) % 4) % 4;
    uint8_t *rgb_data = malloc(width * height * 3);
    fseek(fp_in, header.offset, SEEK_SET);
    for(int i=0; i<height; i++) {
        int row = (info.height>0) ? (height-1-i) : i;
        fread(rgb_data + row*width*3, 3, width, fp_in);
        fseek(fp_in, padding, SEEK_CUR);
    }
    fclose(fp_in);

    memset(global_freq, 0, sizeof(global_freq));
    memset(global_codebook, 0, sizeof(global_codebook));

    FILE *fout = NULL;
    BitWriter bw;

    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1) {
            build_huffman_tree();
            FILE *fcb = fopen(cb_file, "w");
            fprintf(fcb, "Symbol\tFreq\tCode\n");
            for(int i=0; i<256; i++) {
                if(global_freq[i] > 0) fprintf(fcb, "0x%02X\t%d\t%s\n", i, global_freq[i], global_codebook[i].str);
            }
            fclose(fcb);
            fout = fopen(hf_file, is_binary ? "wb" : "w");
            if (is_binary) {
                fwrite(&width, sizeof(int), 1, fout); fwrite(&height, sizeof(int), 1, fout);
                init_bitwriter(&bw, fout);
            } else fprintf(fout, "%d %d\n", width, height);
        }

        int prev_dc[3] = {0};
        const int (*qts[3])[8] = {Q_Luminance, Q_Chrominance, Q_Chrominance};

        for (int r = 0; r < height; r += 8) {
            for (int c = 0; c < width; c += 8) {
                float blk[3][8][8];
                for (int i=0; i<8; i++) {
                    for (int j=0; j<8; j++) {
                        int idx = ((r+i < height ? r+i : height-1)*width + (c+j < width ? c+j : width-1))*3;
                        float y, cb, cr;
                        rgb_to_ycbcr(rgb_data[idx+2], rgb_data[idx+1], rgb_data[idx], &y, &cb, &cr);
                        blk[0][i][j] = y - 128.0f; blk[1][i][j] = cb - 128.0f; blk[2][i][j] = cr - 128.0f;
                    }
                }

                for (int ch = 0; ch < 3; ch++) {
                    float dct[8][8];
                    perform_dct(blk[ch], dct);
                    int16_t zz[64];
                    for(int k=0; k<64; k++) {
                        int u=zigzag_index[k]/8, v=zigzag_index[k]%8;
                        zz[k] = (int16_t)roundf(dct[u][v]/qts[ch][u][v]);
                    }

                    int16_t diff = zz[0] - prev_dc[ch];
                    prev_dc[ch] = zz[0];
                    int dc_cat = get_category(diff);
                    int dc_symbol = dc_cat; 

                    if (pass == 0) global_freq[dc_symbol]++;
                    else if (is_binary) {
                        write_bits(&bw, global_codebook[dc_symbol].code, global_codebook[dc_symbol].len);
                        write_bits(&bw, get_vli(diff, dc_cat), dc_cat);
                    } else fprintf(fout, "%s %d ", global_codebook[dc_symbol].str, diff);

                    int zero_run = 0;
                    for (int k = 1; k < 64; k++) {
                        if (zz[k] == 0) zero_run++;
                        else {
                            while (zero_run >= 16) {
                                int zrl = 0xF0;
                                if (pass == 0) global_freq[zrl]++;
                                else if (is_binary) write_bits(&bw, global_codebook[zrl].code, global_codebook[zrl].len);
                                else fprintf(fout, "%s ", global_codebook[zrl].str);
                                zero_run -= 16;
                            }
                            int ac_cat = get_category(zz[k]);
                            int ac_sym = (zero_run << 4) | ac_cat;
                            if (pass == 0) global_freq[ac_sym]++;
                            else if (is_binary) {
                                write_bits(&bw, global_codebook[ac_sym].code, global_codebook[ac_sym].len);
                                write_bits(&bw, get_vli(zz[k], ac_cat), ac_cat);
                            } else fprintf(fout, "%s %d ", global_codebook[ac_sym].str, zz[k]);
                            zero_run = 0;
                        }
                    }
                    if (zero_run > 0) {
                        int eob = 0x00;
                        if (pass == 0) global_freq[eob]++;
                        else if (is_binary) write_bits(&bw, global_codebook[eob].code, global_codebook[eob].len);
                        else fprintf(fout, "%s ", global_codebook[eob].str);
                    }
                    if (pass == 1 && !is_binary) fprintf(fout, "\n");
                }
            }
        }
    }

    if (is_binary && fout) flush_bits(&bw);
    if (fout) fclose(fout);
    free(rgb_data);
    printf("Mode 3 Encoding Complete.\n");
}

// --- Main ---
int main(int argc, char *argv[]) {
    if (argc < 2) { printf("Usage: %s <mode> ...\n", argv[0]); return 1; }
    int mode = atoi(argv[1]);
    if (mode == 0) mode_0_encoder(argv[2], argv[3], argv[4], argv[5], argv[6]);
    else if (mode == 1) mode_1_encoder(argc, argv);
    else if (mode == 2) mode_2_encoder(argc, argv);
    else if (mode == 3) mode_3_encoder(argc, argv);
    else printf("Unknown mode %d\n", mode);
    return 0;
}