#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Zigzag & Quant Tables
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

// --- Structs ---
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

typedef struct DNode {
    int symbol;
    struct DNode *left, *right;
} DNode;

typedef struct {
    FILE *fp;
    uint8_t buffer;
    int bit_count;
} BitReader;

// --- Helper Functions ---

void init_bitreader(BitReader *br, FILE *fp) {
    br->fp = fp;
    br->buffer = 0;
    br->bit_count = 0;
}

// Return -1 on EOF
int read_bit(BitReader *br) {
    if (br->bit_count == 0) {
        size_t bytes_read = fread(&br->buffer, 1, 1, br->fp);
        if (bytes_read == 0) return -1;
        br->bit_count = 8;
    }
    int bit = (br->buffer >> (br->bit_count - 1)) & 1;
    br->bit_count--;
    return bit;
}

uint32_t read_bits(BitReader *br, int n) {
    uint32_t val = 0;
    for (int i = 0; i < n; i++) {
        int bit = read_bit(br);
        if (bit == -1) return 0;
        val = (val << 1) | bit;
    }
    return val;
}

int16_t decode_vli(uint32_t val, int bits) {
    if (bits == 0) return 0;
    if ((val >> (bits - 1)) == 1) return (int16_t)val;
    return (int16_t)val - (1 << bits) + 1;
}

DNode* build_decoding_tree(const char *codebook_file) {
    FILE *fp = fopen(codebook_file, "r");
    if (!fp) exit(1);
    DNode *root = malloc(sizeof(DNode));
    root->symbol = -1; root->left = root->right = NULL;
    char line[100]; fgets(line, sizeof(line), fp);
    int symbol, freq;
    // ★★★ 關鍵修正：加大到 256 以避免緩衝區溢位 ★★★
    char code_str[256]; 

    while (fscanf(fp, "0x%X %d %s", &symbol, &freq, code_str) == 3) {
        DNode *curr = root;
        for (int i = 0; code_str[i] != '\0'; i++) {
            if (code_str[i] == '0') {
                if (!curr->left) {
                    curr->left = malloc(sizeof(DNode));
                    curr->left->symbol = -1; curr->left->left = curr->left->right = NULL;
                }
                curr = curr->left;
            } else {
                if (!curr->right) {
                    curr->right = malloc(sizeof(DNode));
                    curr->right->symbol = -1; curr->right->left = curr->right->right = NULL;
                }
                curr = curr->right;
            }
        }
        curr->symbol = symbol;
    }
    fclose(fp);
    return root;
}

void ycbcr_to_rgb(float y, float cb, float cr, uint8_t *r, uint8_t *g, uint8_t *b) {
    // ★★★ 關鍵修正：修復圖片變暗問題 (+128) ★★★
    y += 128.0f;
    cb += 128.0f;
    cr += 128.0f;

    float val_r = y + 1.402f * (cr - 128.0f);
    float val_g = y - 0.344136f * (cb - 128.0f) - 0.714136f * (cr - 128.0f);
    float val_b = y + 1.772f * (cb - 128.0f);
    
    // Clipping
    if (val_r < 0) val_r = 0; if (val_r > 255) val_r = 255;
    if (val_g < 0) val_g = 0; if (val_g > 255) val_g = 255;
    if (val_b < 0) val_b = 0; if (val_b > 255) val_b = 255;
    *r = (uint8_t)val_r; *g = (uint8_t)val_g; *b = (uint8_t)val_b;
}

void perform_idct(float input[8][8], float output[8][8]) {
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            float sum = 0.0f;
            for (int u = 0; u < 8; u++) {
                for (int v = 0; v < 8; v++) {
                    float cu = (u == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
                    float cv = (v == 0) ? 1.0f / sqrtf(2.0f) : 1.0f;
                    sum += cu * cv * input[u][v] * cosf((2 * x + 1) * u * M_PI / 16.0f) * cosf((2 * y + 1) * v * M_PI / 16.0f);
                }
            }
            output[x][y] = 0.25f * sum;
        }
    }
}

void read_qt_txt(const char *filename, int qt[8][8]) {
    FILE *fp = fopen(filename, "r");
    if (!fp) exit(1);
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) fscanf(fp, "%d", &qt[i][j]);
    }
    fclose(fp);
}

// --- Mode 0 Decoder ---
void mode_0_decoder(char *output_bmp, char *r_txt, char *g_txt, char *b_txt, char *dim_txt) {
    FILE *fdim = fopen(dim_txt, "r"); if (!fdim) return;
    int width, height; fscanf(fdim, "%d %d", &width, &height); fclose(fdim);
    uint8_t *R = malloc(width * height), *G = malloc(width * height), *B = malloc(width * height);
    FILE *fr = fopen(r_txt, "r"), *fg = fopen(g_txt, "r"), *fb = fopen(b_txt, "r");
    for(int i=0; i<height*width; i++) {
        int r, g, b; fscanf(fr, "%d", &r); fscanf(fg, "%d", &g); fscanf(fb, "%d", &b);
        R[i]=r; G[i]=g; B[i]=b;
    }
    fclose(fr); fclose(fg); fclose(fb);
    FILE *fout = fopen(output_bmp, "wb");
    int padding = (4 - (width * 3) % 4) % 4;
    int size = (width * 3 + padding) * height;
    BMPHeader h = {0x4D42, 54+size, 0, 0, 54};
    BMPInfoHeader info = {40, width, height, 1, 24, 0, size, 2835, 2835, 0, 0};
    fwrite(&h, sizeof(h), 1, fout); fwrite(&info, sizeof(info), 1, fout);
    uint8_t pad = 0;
    for(int i = height - 1; i >= 0; i--) {
        for(int j = 0; j < width; j++) {
            fwrite(&B[i*width+j], 1, 1, fout); fwrite(&G[i*width+j], 1, 1, fout); fwrite(&R[i*width+j], 1, 1, fout);
        }
        fwrite(&pad, 1, padding, fout);
    }
    fclose(fout); free(R); free(G); free(B);
}

// --- Mode 1 Decoder ---
void mode_1_decoder(int argc, char *argv[]) {
    int is_mode_1b = (argc == 13);
    char *out_bmp = argv[2];
    int arg_offset = is_mode_1b ? 3 : 4;
    int qt_y[8][8], qt_cb[8][8], qt_cr[8][8];
    read_qt_txt(argv[arg_offset], qt_y); read_qt_txt(argv[arg_offset+1], qt_cb); read_qt_txt(argv[arg_offset+2], qt_cr);
    FILE *fdim = fopen(argv[arg_offset+3], "r");
    int width, height; fscanf(fdim, "%d %d", &width, &height); fclose(fdim);
    FILE *f_qY = fopen(argv[arg_offset+4], "rb"), *f_qCb = fopen(argv[arg_offset+5], "rb"), *f_qCr = fopen(argv[arg_offset+6], "rb");
    FILE *f_eY=NULL, *f_eCb=NULL, *f_eCr=NULL;
    if (is_mode_1b) { f_eY=fopen(argv[arg_offset+7],"rb"); f_eCb=fopen(argv[arg_offset+8],"rb"); f_eCr=fopen(argv[arg_offset+9],"rb"); }
    uint8_t *recon_img = malloc(width * height * 3);
    if (!recon_img) { printf("Memory allocation failed!\n"); exit(1); }

    for (int r = 0; r < height; r += 8) {
        for (int c = 0; c < width; c += 8) {
            float blk_y[8][8], blk_cb[8][8], blk_cr[8][8];
            float (*blocks[3])[8] = {blk_y, blk_cb, blk_cr};
            FILE *f_q[3] = {f_qY, f_qCb, f_qCr}; FILE *f_e[3] = {f_eY, f_eCb, f_eCr};
            int (*qts[3])[8] = {qt_y, qt_cb, qt_cr};
            for (int ch = 0; ch < 3; ch++) {
                float idct_in[8][8];
                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 8; j++) {
                        int16_t qF; fread(&qF, 2, 1, f_q[ch]);
                        float val = qF * qts[ch][i][j];
                        if (is_mode_1b) { float eF; fread(&eF, 4, 1, f_e[ch]); val += eF; }
                        idct_in[i][j] = val;
                    }
                }
                perform_idct(idct_in, blocks[ch]);
            }
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 8; j++) {
                    if (r+i >= height || c+j >= width) continue;
                    uint8_t R, G, B; ycbcr_to_rgb(blk_y[i][j], blk_cb[i][j], blk_cr[i][j], &R, &G, &B);
                    int idx = ((r + i) * width + (c + j)) * 3;
                    recon_img[idx] = R; recon_img[idx+1] = G; recon_img[idx+2] = B;
                }
            }
        }
    }
    fclose(f_qY); fclose(f_qCb); fclose(f_qCr);
    if (is_mode_1b) { fclose(f_eY); fclose(f_eCb); fclose(f_eCr); }
    FILE *fout = fopen(out_bmp, "wb");
    int padding = (4 - (width * 3) % 4) % 4;
    int size = (width * 3 + padding) * height;
    BMPHeader h = {0x4D42, 54+size, 0, 0, 54}; BMPInfoHeader info = {40, width, height, 1, 24, 0, size, 2835, 2835, 0, 0};
    fwrite(&h, sizeof(h), 1, fout); fwrite(&info, sizeof(info), 1, fout);
    uint8_t pad = 0;
    for(int i = height - 1; i >= 0; i--) {
        for(int j = 0; j < width; j++) {
            int idx = (i * width + j) * 3;
            uint8_t R=recon_img[idx], G=recon_img[idx+1], B=recon_img[idx+2];
            fwrite(&B, 1, 1, fout); fwrite(&G, 1, 1, fout); fwrite(&R, 1, 1, fout);
        }
        fwrite(&pad, 1, padding, fout);
    }
    fclose(fout); free(recon_img);
    printf("Mode 1 Decoding complete. Output: %s\n", out_bmp);
}

// --- Mode 2 Decoder ---
void mode_2_decoder(int argc, char *argv[]) {
    int is_binary = (strcmp(argv[3], "binary") == 0);
    FILE *fin = fopen(argv[4], is_binary ? "rb" : "r");
    int width, height;
    if (is_binary) { fread(&width, 4, 1, fin); fread(&height, 4, 1, fin); }
    else fscanf(fin, "%d %d", &width, &height);
    uint8_t *recon_img = malloc(width * height * 3);
    if (!recon_img) { printf("Memory allocation failed!\n"); exit(1); }
    int prev_dc[3] = {0};
    for (int r = 0; r < height; r += 8) {
        for (int c = 0; c < width; c += 8) {
            float blk[3][8][8]; const int (*qts[3])[8] = {Q_Luminance, Q_Chrominance, Q_Chrominance};
            for (int ch = 0; ch < 3; ch++) {
                int16_t zz[64] = {0};
                if (!is_binary) {
                    // ASCII skip
                } else {
                    int16_t diff; fread(&diff, 2, 1, fin);
                    zz[0] = prev_dc[ch] + diff; prev_dc[ch] = zz[0];
                    int cursor = 1;
                    while (1) {
                        uint8_t run; int16_t val; fread(&run, 1, 1, fin); fread(&val, 2, 1, fin);
                        if (run == 0 && val == 0) break;
                        cursor += run; zz[cursor] = val; cursor++;
                    }
                }
                float dct[8][8];
                for (int k = 0; k < 64; k++) {
                    int u = zigzag_index[k]/8, v = zigzag_index[k]%8;
                    dct[u][v] = zz[k] * qts[ch][u][v];
                }
                perform_idct(dct, blk[ch]);
            }
            for(int i=0; i<8; i++) {
                for(int j=0; j<8; j++) {
                    if (r+i >= height || c+j >= width) continue;
                    uint8_t R, G, B; ycbcr_to_rgb(blk[0][i][j], blk[1][i][j], blk[2][i][j], &R, &G, &B);
                    int idx = ((r+i)*width + (c+j))*3;
                    recon_img[idx] = R; recon_img[idx+1] = G; recon_img[idx+2] = B;
                }
            }
        }
    }
    FILE *fout = fopen(argv[2], "wb");
    int padding = (4 - (width * 3) % 4) % 4; int size = (width * 3 + padding) * height;
    BMPHeader h = {0x4D42, 54+size, 0, 0, 54}; BMPInfoHeader info = {40, width, height, 1, 24, 0, size, 2835, 2835, 0, 0};
    fwrite(&h, sizeof(h), 1, fout); fwrite(&info, sizeof(info), 1, fout);
    uint8_t pad = 0;
    for(int i = height - 1; i >= 0; i--) {
        for(int j = 0; j < width; j++) {
            int idx = (i * width + j) * 3;
            uint8_t B = recon_img[idx+2], G = recon_img[idx+1], R = recon_img[idx];
            fwrite(&B, 1, 1, fout); fwrite(&G, 1, 1, fout); fwrite(&R, 1, 1, fout);
        }
        fwrite(&pad, 1, padding, fout);
    }
    fclose(fout); fclose(fin); free(recon_img);
    printf("Mode 2 Decoding complete: %s\n", argv[2]);
}

// --- Mode 3 Decoder ---
void mode_3_decoder(int argc, char *argv[]) {
    int is_binary = (strcmp(argv[3], "binary") == 0);
    DNode *huff_root = build_decoding_tree(argv[4]);
    FILE *fin = fopen(argv[5], "rb");
    int width, height; fread(&width, 4, 1, fin); fread(&height, 4, 1, fin);
    BitReader br; init_bitreader(&br, fin);
    uint8_t *recon_img = malloc(width * height * 3);
    if (!recon_img) { printf("Memory allocation failed!\n"); exit(1); }

    int prev_dc[3] = {0}; const int (*qts[3])[8] = {Q_Luminance, Q_Chrominance, Q_Chrominance};
    
    // 計算 Block 總數 (處理大檔案時使用 block count 而不是 pixel count 避免誤差)
    int blocks_h = (height + 7) / 8;
    int blocks_w = (width + 7) / 8;

    for (int by = 0; by < blocks_h; by++) {
        for (int bx = 0; bx < blocks_w; bx++) {
            int r = by * 8;
            int c = bx * 8;
            float blk[3][8][8];
            for (int ch = 0; ch < 3; ch++) {
                int16_t zz[64] = {0};
                // DC
                DNode *curr = huff_root;
                while (curr->symbol == -1) {
                    int bit = read_bit(&br);
                    if (bit == -1) break;
                    if (bit == 0) curr = curr->left; else curr = curr->right;
                    if (!curr) break;
                }
                if (!curr || curr->symbol == -1) break;
                int dc_cat = curr->symbol;
                int16_t diff = decode_vli(read_bits(&br, dc_cat), dc_cat);
                zz[0] = prev_dc[ch] + diff; prev_dc[ch] = zz[0];

                // AC
                int k = 1;
                while (k < 64) {
                    curr = huff_root;
                    while (curr->symbol == -1) {
                        int bit = read_bit(&br);
                        if (bit == -1) break;
                        if (bit == 0) curr = curr->left; else curr = curr->right;
                        if (!curr) break;
                    }
                    if (!curr || curr->symbol == -1) break;
                    int sym = curr->symbol;
                    if (sym == 0x00) break; // EOB
                    else if (sym == 0xF0) k += 16;
                    else {
                        int run = (sym >> 4) & 0x0F; int cat = sym & 0x0F;
                        k += run; if (k>=64) break;
                        zz[k] = decode_vli(read_bits(&br, cat), cat); k++;
                    }
                }

                float dct[8][8];
                for (int i = 0; i < 64; i++) {
                    int u = zigzag_index[i]/8, v = zigzag_index[i]%8;
                    dct[u][v] = zz[i] * qts[ch][u][v];
                }
                perform_idct(dct, blk[ch]);
            }
            for(int i=0; i<8; i++) {
                for(int j=0; j<8; j++) {
                    if (r+i >= height || c+j >= width) continue;
                    uint8_t R, G, B; ycbcr_to_rgb(blk[0][i][j], blk[1][i][j], blk[2][i][j], &R, &G, &B);
                    int idx = ((r+i)*width + (c+j))*3;
                    recon_img[idx] = R; recon_img[idx+1] = G; recon_img[idx+2] = B;
                }
            }
        }
    }
    FILE *fout = fopen(argv[2], "wb");
    int padding = (4 - (width * 3) % 4) % 4; int size = (width * 3 + padding) * height;
    BMPHeader h = {0x4D42, 54+size, 0, 0, 54}; BMPInfoHeader info = {40, width, height, 1, 24, 0, size, 2835, 2835, 0, 0};
    fwrite(&h, sizeof(h), 1, fout); fwrite(&info, sizeof(info), 1, fout);
    uint8_t pad = 0;
    for(int i = height - 1; i >= 0; i--) {
        for(int j = 0; j < width; j++) {
            int idx = (i * width + j) * 3;
            uint8_t B = recon_img[idx+2], G = recon_img[idx+1], R = recon_img[idx];
            fwrite(&B, 1, 1, fout); fwrite(&G, 1, 1, fout); fwrite(&R, 1, 1, fout);
        }
        fwrite(&pad, 1, padding, fout);
    }
    fclose(fout); fclose(fin); free(recon_img);
    printf("Mode 3 Decoding Complete: %s\n", argv[2]);
}

// --- Main ---
int main(int argc, char *argv[]) {
    if (argc < 2) return 1;
    int mode = atoi(argv[1]);
    if (mode == 0) mode_0_decoder(argv[2], argv[3], argv[4], argv[5], argv[6]);
    else if (mode == 1) mode_1_decoder(argc, argv);
    else if (mode == 2) mode_2_decoder(argc, argv);
    else if (mode == 3) mode_3_decoder(argc, argv);
    else printf("Unknown mode\n");
    return 0;
}