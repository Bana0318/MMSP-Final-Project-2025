#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 將此陣列放在所有函數 (main, mode_1, mode_2) 之前
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

// --- 1. 定義 BMP 結構 (同前) ---
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

// --- 2. 標準量化表 (Standard JPEG Quantization Tables) ---
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
// --- Mode 3 專用結構 ---

// 用來寫入任意長度 bit 的緩衝區
typedef struct {
    FILE *fp;
    uint8_t buffer;
    int bit_count;
} BitWriter;

// 霍夫曼樹節點
typedef struct HNode {
    int symbol; // 0x00 ~ 0xFF (Run << 4 | Size)
    int freq;
    struct HNode *left, *right;
} HNode;

// 霍夫曼編碼表 (存字串方便輸出 codebook.txt，也方便轉換)
typedef struct {
    int len;
    uint32_t code; // 整數形式的 code ( binary )
    char str[32];  // 字串形式的 code "0101"
} HCode;

HCode global_codebook[256]; // 簡化版：所有通道共用一個 Codebook
int global_freq[256];       // 頻率統計


// --- 3. 核心數學運算 ---

// RGB 轉 YCbCr (整數輸入，浮點輸出)
void rgb_to_ycbcr(uint8_t r, uint8_t g, uint8_t b, float *y, float *cb, float *cr) {
    *y  =  0.299f * r + 0.587f * g + 0.114f * b;
    *cb = -0.1687f * r - 0.3313f * g + 0.5f * b + 128.0f;
    *cr =  0.5f * r - 0.4187f * g - 0.0813f * b + 128.0f;
}

// 2D DCT
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

// 寫入量化表 TXT
void write_qt_txt(const char *filename, const int qt[8][8]) {
    FILE *fp = fopen(filename, "w");
    if(!fp) return;
    for(int i=0; i<8; i++) {
        for(int j=0; j<8; j++) fprintf(fp, "%d%c", qt[i][j], (j==7?'\n':' '));
    }
    fclose(fp);
}

// --- 4. Mode 0 (保持不變) ---
void mode_0_encoder(char *input_bmp, char *r_txt, char *g_txt, char *b_txt, char *dim_txt) {
FILE *fp = fopen(input_bmp, "rb");
    if (!fp) {
        perror("Error opening input BMP");
        exit(1);
    }

    BMPHeader header;
    BMPInfoHeader info;
    fread(&header, sizeof(BMPHeader), 1, fp);
    fread(&info, sizeof(BMPInfoHeader), 1, fp);

    int width = info.width;
    int height = abs(info.height); // 高度可能為負值
    int padding = (4 - (width * 3) % 4) % 4; // BMP 每橫列需為 4 bytes 倍數

    // 動態分配記憶體儲存 R, G, B
    uint8_t **R = (uint8_t **)malloc(height * sizeof(uint8_t *));
    uint8_t **G = (uint8_t **)malloc(height * sizeof(uint8_t *));
    uint8_t **B = (uint8_t **)malloc(height * sizeof(uint8_t *));
    for (int i = 0; i < height; i++) {
        R[i] = (uint8_t *)malloc(width);
        G[i] = (uint8_t *)malloc(width);
        B[i] = (uint8_t *)malloc(width);
    }

    // 尋找到像素起始點
    fseek(fp, header.offset, SEEK_SET);

    // 讀取像素資料 (BMP 存儲順序為 BGR，且通常由下往上)
    for (int i = 0; i < height; i++) {
        int target_row = (info.height > 0) ? (height - 1 - i) : i;
        for (int j = 0; j < width; j++) {
            uint8_t bgr[3];
            fread(bgr, 3, 1, fp);
            B[target_row][j] = bgr[0];
            G[target_row][j] = bgr[1];
            R[target_row][j] = bgr[2];
        }
        fseek(fp, padding, SEEK_CUR); // 跳過 Padding bits
    }
    fclose(fp);

    // 輸出 R, G, B txt 檔案
    FILE *fr = fopen(r_txt, "w");
    FILE *fg = fopen(g_txt, "w");
    FILE *fb = fopen(b_txt, "w");
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            fprintf(fr, "%d%c", R[i][j], (j == width - 1 ? '\n' : ' '));
            fprintf(fg, "%d%c", G[i][j], (j == width - 1 ? '\n' : ' '));
            fprintf(fb, "%d%c", B[i][j], (j == width - 1 ? '\n' : ' '));
        }
    }
    fclose(fr); fclose(fg); fclose(fb);

    // 輸出維度資訊
    FILE *fd = fopen(dim_txt, "w");
    fprintf(fd, "%d %d\n", width, height);
    fclose(fd);

    // 釋放記憶體
    for (int i = 0; i < height; i++) {
        free(R[i]); free(G[i]); free(B[i]);
    }
    free(R); free(G); free(B);
    printf("Mode 0 executed (Placeholder).\n"); 
}

// --- 5. Mode 1 實作 ---
void mode_1_encoder(int argc, char *argv[]) {
    // 參數對應:
    // argv[2]: bmp, argv[3-5]: Qt_txt, argv[6]: dim
    // argv[7-9]: qF_raw (short), argv[10-12]: eF_raw (float)

    FILE *fp_in = fopen(argv[2], "rb");
    if (!fp_in) { perror("Open BMP failed"); exit(1); }

    BMPHeader header; BMPInfoHeader info;
    fread(&header, sizeof(BMPHeader), 1, fp_in);
    fread(&info, sizeof(BMPInfoHeader), 1, fp_in);
    
    int width = info.width;
    int height = abs(info.height);
    int padding = (4 - (width * 3) % 4) % 4;
    
    // 讀取 BMP 像素
    uint8_t *rgb_data = (uint8_t*)malloc(width * height * 3);
    fseek(fp_in, header.offset, SEEK_SET);
    
    // 處理 Bottom-up 與 Padding
    for (int i = 0; i < height; i++) {
        int row = (info.height > 0) ? (height - 1 - i) : i;
        fread(rgb_data + row * width * 3, 3, width, fp_in);
        fseek(fp_in, padding, SEEK_CUR);
    }
    fclose(fp_in);

    // 輸出 dim.txt
    FILE *f_dim = fopen(argv[6], "w");
    fprintf(f_dim, "%d %d\n", width, height);
    fclose(f_dim);

    // 輸出量化表 TXT
    write_qt_txt(argv[3], Q_Luminance); // Qt_Y
    write_qt_txt(argv[4], Q_Chrominance); // Qt_Cb
    write_qt_txt(argv[5], Q_Chrominance); // Qt_Cr

    // 開啟輸出檔案 (qF: binary short, eF: binary float)
    FILE *f_qY = fopen(argv[7], "wb"), *f_qCb = fopen(argv[8], "wb"), *f_qCr = fopen(argv[9], "wb");
    FILE *f_eY = fopen(argv[10], "wb"), *f_eCb = fopen(argv[11], "wb"), *f_eCr = fopen(argv[12], "wb");

    // 用來計算 SQNR 的累積能量 (Signal Power, Noise Power)
    // index 0:Y, 1:Cb, 2:Cr. 每個通道有 64 個頻率成分
    double signal_pow[3][64] = {0};
    double noise_pow[3][64] = {0};

    // 處理每個 8x8 Block
    // 注意：如果有邊緣不足 8 的情況，通常 JPEG 會填充，但這裡假設圖檔是 8 的倍數或直接處理
    for (int r = 0; r < height; r += 8) {
        for (int c = 0; c < width; c += 8) {
            float blk_y[8][8], blk_cb[8][8], blk_cr[8][8];
            
            // 1. 擷取 8x8 並做 RGB -> YCbCr & Level Shift
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 8; j++) {
                    int idx = ((r + i) * width + (c + j)) * 3;
                    // 邊界檢查：如果超出範圍，重複邊緣像素 (簡易處理)
                    if (r+i >= height || c+j >= width) idx = ((height-1)*width + (width-1))*3; 
                    
                    uint8_t B = rgb_data[idx];
                    uint8_t G = rgb_data[idx+1];
                    uint8_t R = rgb_data[idx+2];
                    
                    float y, cb, cr;
                    rgb_to_ycbcr(R, G, B, &y, &cb, &cr);
                    blk_y[i][j] = y - 128.0f;
                    blk_cb[i][j] = cb - 128.0f;
                    blk_cr[i][j] = cr - 128.0f;
                }
            }

            // 定義通道處理陣列以方便迴圈
            float (*blocks[3])[8] = {blk_y, blk_cb, blk_cr};
            FILE *f_q[3] = {f_qY, f_qCb, f_qCr};
            FILE *f_e[3] = {f_eY, f_eCb, f_eCr};
            const int (*qts[3])[8] = {Q_Luminance, Q_Chrominance, Q_Chrominance};

            // 針對 Y, Cb, Cr 三個通道分別處理
            for (int ch = 0; ch < 3; ch++) {
                float dct_out[8][8];
                perform_dct(blocks[ch], dct_out);

                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 8; j++) {
                        int k = i * 8 + j; // 1D index (0-63)
                        
                        // Quantization: qF = round(F / Q)
                        int q_val = Q_Chrominance[i][j]; // Default
                        if (ch == 0) q_val = Q_Luminance[i][j];
                        
                        int16_t qF = (int16_t)roundf(dct_out[i][j] / q_val);
                        
                        // Quantization Error: eF = F - qF * Q
                        float rec_F = qF * q_val;
                        float eF = dct_out[i][j] - rec_F;

                        // 寫入檔案
                        fwrite(&qF, sizeof(int16_t), 1, f_q[ch]);
                        fwrite(&eF, sizeof(float), 1, f_e[ch]);

                        // 累積 SQNR 能量
                        signal_pow[ch][k] += dct_out[i][j] * dct_out[i][j];
                        noise_pow[ch][k] += eF * eF;
                    }
                }
            }
        }
    }

    // 關閉所有檔案
    fclose(f_qY); fclose(f_qCb); fclose(f_qCr);
    fclose(f_eY); fclose(f_eCb); fclose(f_eCr);
    free(rgb_data);

    // 列印 SQNR (3 channels x 64 freqs)
    char *ch_names[3] = {"Y", "Cb", "Cr"};
    for (int ch = 0; ch < 3; ch++) {
        printf("SQNR (dB) for Channel %s:\n", ch_names[ch]);
        for (int i = 0; i < 64; i++) {
            double s = signal_pow[ch][i];
            double n = noise_pow[ch][i];
            double sqnr = 0.0;
            if (n == 0) sqnr = 99.99; // 避免除以零，表示無誤差
            else if (s == 0) sqnr = 0.0; // 無訊號
            else sqnr = 10.0 * log10(s / n);
            
            printf("%6.2f ", sqnr);
            if ((i + 1) % 8 == 0) printf("\n");
        }
        printf("\n");
    }
}
// --- 6. Mode 2: RLE Encoder ---

// 輔助結構：儲存 RLE Pair
typedef struct {
    uint8_t run;
    int16_t value;
} RLEPair;

void mode_2_encoder(int argc, char *argv[]) {
    // 參數: 2 input.bmp [ascii/binary] output_file
    if (argc != 5) {
        printf("Usage mode 2: encoder 2 input.bmp [ascii/binary] output_file\n");
        exit(1);
    }
    
    int is_binary = (strcmp(argv[3], "binary") == 0);
    char *output_file = argv[4];

    // 讀取 BMP (複製自 Mode 1 邏輯，為求簡潔這裡直接寫重點)
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

    // 開啟輸出檔案
    FILE *fout = fopen(output_file, is_binary ? "wb" : "w");
    
    // 寫入 Header (寬高)
    if (is_binary) {
        fwrite(&width, sizeof(int), 1, fout);
        fwrite(&height, sizeof(int), 1, fout);
    } else {
        fprintf(fout, "%d %d\n", width, height);
    }

    // DPCM 狀態變數 (記錄前一個 Block 的 DC 值)
    int prev_dc[3] = {0, 0, 0}; 

    // 處理 Blocks
    for (int r = 0; r < height; r += 8) {
        for (int c = 0; c < width; c += 8) {
            float blk[3][8][8];
            // 擷取 RGB -> YCbCr -> Level Shift (同 Mode 1)
            for (int i=0; i<8; i++) {
                for (int j=0; j<8; j++) {
                    int idx = ((r+i < height ? r+i : height-1)*width + (c+j < width ? c+j : width-1))*3;
                    float y, cb, cr;
                    rgb_to_ycbcr(rgb_data[idx+2], rgb_data[idx+1], rgb_data[idx], &y, &cb, &cr);
                    blk[0][i][j] = y - 128.0f; 
                    blk[1][i][j] = cb - 128.0f; 
                    blk[2][i][j] = cr - 128.0f;
                }
            }

            // 對三個通道分別做 DCT -> Quant -> ZigZag -> RLE
            const int (*qts[3])[8] = {Q_Luminance, Q_Chrominance, Q_Chrominance};
            char *ch_names[3] = {"Y", "Cb", "Cr"};

            for (int ch = 0; ch < 3; ch++) {
                float dct[8][8];
                perform_dct(blk[ch], dct);

                // Quantization & ZigZag
                int16_t zz_block[64];
                for (int k = 0; k < 64; k++) {
                    int u = zigzag_index[k] / 8;
                    int v = zigzag_index[k] % 8;
                    zz_block[k] = (int16_t)roundf(dct[u][v] / qts[ch][u][v]);
                }

                // --- DPCM for DC (Index 0) ---
                int16_t dc_val = zz_block[0];
                int16_t dc_diff = dc_val - prev_dc[ch];
                prev_dc[ch] = dc_val; // 更新

                // --- RLE for AC (Index 1~63) ---
                RLEPair ac_rle[64];
                int rle_count = 0;
                int zero_run = 0;

                for (int k = 1; k < 64; k++) {
                    if (zz_block[k] == 0) {
                        zero_run++;
                    } else {
                        ac_rle[rle_count].run = zero_run;
                        ac_rle[rle_count].value = zz_block[k];
                        rle_count++;
                        zero_run = 0;
                    }
                }
                // 注意：最後剩下的 zero_run 不需要存 (EOB 概念)

                // --- 寫入檔案 ---
                if (!is_binary) {
                    // ASCII: ($m,$n, Channel) DC_Diff run val run val ...
                    // 這裡我們把 DC 當作 (0, diff) 輸出
                    fprintf(fout, "(%d,%d, %s) 0 %d ", r/8, c/8, ch_names[ch], dc_diff);
                    for (int k = 0; k < rle_count; k++) {
                        fprintf(fout, "%d %d ", ac_rle[k].run, ac_rle[k].value);
                    }
                    fprintf(fout, "\n");
                } else {
                    // Binary: [DC_Diff(2B)] [AC_Pairs...] [EOB(0,0)]
                    fwrite(&dc_diff, sizeof(int16_t), 1, fout);
                    for (int k = 0; k < rle_count; k++) {
                        fwrite(&ac_rle[k].run, sizeof(uint8_t), 1, fout);
                        fwrite(&ac_rle[k].value, sizeof(int16_t), 1, fout);
                    }
                    // 寫入 EOB (0, 0) 作為區塊結束標記
                    uint8_t eob_run = 0; int16_t eob_val = 0;
                    fwrite(&eob_run, sizeof(uint8_t), 1, fout);
                    fwrite(&eob_val, sizeof(int16_t), 1, fout);
                }
            }
        }
    }
    
    // 計算壓縮率
    long original_size = 54 + width * height * 3;
    long compressed_size = ftell(fout);
    fclose(fout);
    free(rgb_data);

    if (is_binary) {
        printf("Compression Ratio: %.2f%% (Size: %ld / %ld)\n", 
               (double)compressed_size / original_size * 100.0, compressed_size, original_size);
    } else {
        printf("ASCII RLE file generated: %s\n", output_file);
    }
}

// --- BitWriter 實作 ---
void init_bitwriter(BitWriter *bw, FILE *fp) {
    bw->fp = fp;
    bw->buffer = 0;
    bw->bit_count = 0;
}

void write_bits(BitWriter *bw, uint32_t bits, int len) {
    // 從高位寫到低位
    for (int i = len - 1; i >= 0; i--) {
        int bit = (bits >> i) & 1;
        bw->buffer = (bw->buffer << 1) | bit;
        bw->bit_count++;
        if (bw->bit_count == 8) {
            // 如果是 JPEG 標準，0xFF 後面要補 0x00 (Byte Stuffing)
            // 這裡為了作業簡化，我們直接寫入 (或者依照教授要求)
            // 建議：為了避免跟 Marker 衝突，通常需要 stuffing，但作業可先忽略
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

// --- JPEG 數值編碼 (Category & VLI) ---
// 計算數值需要幾個 bits (Category/Size)
int get_category(int16_t val) {
    if (val == 0) return 0;
    val = abs(val);
    int cat = 0;
    while (val > 0) { val >>= 1; cat++; }
    return cat;
}

// 取得數值的實際位元 (Variable Length Integer)
// 正數就是原碼，負數是反碼 (One's Complement)
uint32_t get_vli(int16_t val, int cat) {
    if (val > 0) return val;
    return val + (1 << cat) - 1; // 負數轉換
}

// --- Huffman Tree 實作 (簡化版) ---
// 遞迴產生編碼
void generate_codes_recursive(HNode *node, uint32_t code, int len) {
    if (!node) return;
    if (!node->left && !node->right) { // 葉節點
        global_codebook[node->symbol].len = len;
        global_codebook[node->symbol].code = code;
        // 生成字串方便 debug/output
        for(int i=0; i<len; i++) 
            global_codebook[node->symbol].str[i] = ((code >> (len-1-i)) & 1) ? '1' : '0';
        global_codebook[node->symbol].str[len] = '\0';
        return;
    }
    generate_codes_recursive(node->left, (code << 1) | 0, len + 1);
    generate_codes_recursive(node->right, (code << 1) | 1, len + 1);
}

// 建立 Huffman Tree 並生成 Codebook
void build_huffman_tree() {
    HNode *pool[256];
    int count = 0;
    
    // 1. 初始化葉節點
    for (int i = 0; i < 256; i++) {
        if (global_freq[i] > 0) {
            HNode *n = malloc(sizeof(HNode));
            n->symbol = i;
            n->freq = global_freq[i];
            n->left = n->right = NULL;
            pool[count++] = n;
        }
    }
    
    // 2. 建樹 (每次找兩個最小的合併)
    while (count > 1) {
        // 尋找最小的兩個 (min1 < min2)
        int m1 = -1, m2 = -1;
        for (int i = 0; i < count; i++) {
            if (m1 == -1 || pool[i]->freq < pool[m1]->freq) {
                m2 = m1; m1 = i;
            } else if (m2 == -1 || pool[i]->freq < pool[m2]->freq) {
                m2 = i;
            }
        }
        
        // 合併
        HNode *parent = malloc(sizeof(HNode));
        parent->symbol = -1; // Internal node
        parent->freq = pool[m1]->freq + pool[m2]->freq;
        parent->left = pool[m1];
        parent->right = pool[m2];
        
        // 移除這兩個，加入 parent
        // 簡單作法：把 m1 換成 parent，把最後一個搬到 m2
        pool[m1] = parent;
        pool[m2] = pool[count - 1];
        count--;
    }
    
    // 3. 產生編碼
    if (count > 0) {
        generate_codes_recursive(pool[0], 0, 0);
    }
}

// 釋放 Tree 記憶體 (略，作業可由 OS 回收)

void mode_3_encoder(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Usage mode 3: encoder 3 input.bmp [ascii/binary] codebook.txt huffman_code.bin\n");
        exit(1);
    }
    
    int is_binary = (strcmp(argv[3], "binary") == 0);
    char *cb_file = argv[4];
    char *hf_file = argv[5];

    // --- 讀取 BMP (標準起手式) ---
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

    // 初始化頻率表
    memset(global_freq, 0, sizeof(global_freq));
    memset(global_codebook, 0, sizeof(global_codebook));

    // --- 準備 Two-Pass 迴圈 ---
    // pass=0: 統計頻率, pass=1: 真正寫入
    FILE *fout = NULL;
    BitWriter bw;

    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1) {
            // Pass 1 結束，建立 Codebook 並輸出
            build_huffman_tree();
            
            // 輸出 Codebook.txt
            FILE *fcb = fopen(cb_file, "w");
            fprintf(fcb, "Symbol\tFreq\tCode\n");
            for(int i=0; i<256; i++) {
                if(global_freq[i] > 0) {
                    fprintf(fcb, "0x%02X\t%d\t%s\n", i, global_freq[i], global_codebook[i].str);
                }
            }
            fclose(fcb);

            // 準備 Pass 2 寫入
            fout = fopen(hf_file, is_binary ? "wb" : "w");
            if (is_binary) {
                // 寫入 Header: Width, Height
                fwrite(&width, sizeof(int), 1, fout);
                fwrite(&height, sizeof(int), 1, fout);
                init_bitwriter(&bw, fout);
            } else {
                fprintf(fout, "%d %d\n", width, height);
            }
        }

        // DPCM 重置
        int prev_dc[3] = {0};
        const int (*qts[3])[8] = {Q_Luminance, Q_Chrominance, Q_Chrominance};

        for (int r = 0; r < height; r += 8) {
            for (int c = 0; c < width; c += 8) {
                float blk[3][8][8];
                // 擷取 Pixel -> YCbCr -> Level Shift
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

                    // --- Huffman Encoding Logic ---
                    
                    // 1. DC Component
                    int16_t diff = zz[0] - prev_dc[ch];
                    prev_dc[ch] = zz[0];
                    int dc_cat = get_category(diff);
                    int dc_symbol = (0x00) | dc_cat; // DC Run is always 0

                    if (pass == 0) {
                        global_freq[dc_symbol]++;
                    } else {
                        // 寫入 DC Symbol Code
                        if (is_binary) {
                            write_bits(&bw, global_codebook[dc_symbol].code, global_codebook[dc_symbol].len);
                            write_bits(&bw, get_vli(diff, dc_cat), dc_cat);
                        } else {
                            // ASCII Output (簡化)
                            fprintf(fout, "%s %d ", global_codebook[dc_symbol].str, diff);
                        }
                    }

                    // 2. AC Components
                    int zero_run = 0;
                    for (int k = 1; k < 64; k++) {
                        if (zz[k] == 0) {
                            zero_run++;
                        } else {
                            while (zero_run >= 16) {
                                // ZRL (Zero Run Length) for run >= 16
                                int zrl_sym = 0xF0; 
                                if (pass == 0) global_freq[zrl_sym]++;
                                else if (is_binary) write_bits(&bw, global_codebook[zrl_sym].code, global_codebook[zrl_sym].len);
                                else fprintf(fout, "%s ", global_codebook[zrl_sym].str);
                                zero_run -= 16;
                            }
                            
                            int ac_cat = get_category(zz[k]);
                            int ac_symbol = (zero_run << 4) | ac_cat;
                            
                            if (pass == 0) {
                                global_freq[ac_symbol]++;
                            } else {
                                if (is_binary) {
                                    write_bits(&bw, global_codebook[ac_symbol].code, global_codebook[ac_symbol].len);
                                    write_bits(&bw, get_vli(zz[k], ac_cat), ac_cat);
                                } else {
                                    fprintf(fout, "%s %d ", global_codebook[ac_symbol].str, zz[k]);
                                }
                            }
                            zero_run = 0;
                        }
                    }

                    // EOB (End of Block) if zeros remain
                    if (zero_run > 0) {
                        int eob_sym = 0x00;
                        if (pass == 0) global_freq[eob_sym]++;
                        else if (is_binary) write_bits(&bw, global_codebook[eob_sym].code, global_codebook[eob_sym].len);
                        else fprintf(fout, "%s ", global_codebook[eob_sym].str);
                    }
                    
                    if (pass == 1 && !is_binary) fprintf(fout, "\n"); // ASCII 換行
                }
            }
        }
    } // End of Two-Pass Loop

    if (is_binary && fout) flush_bits(&bw);
    if (fout) fclose(fout);
    free(rgb_data);
    printf("Mode 3 Encoding Complete.\n");
}

// --- 6. Main 函數 ---
int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s <mode> [params...]\n", argv[0]);
        return 1;
    }

    
    int mode = atoi(argv[1]);

    if (mode == 0) {
        if (argc != 7) {
            printf("Usage mode 0: encoder 0 input.bmp R.txt G.txt B.txt dim.txt\n");
            return 1;
        }
        mode_0_encoder(argv[2], argv[3], argv[4], argv[5], argv[6]);
    } 
    else if (mode == 1) {
        // encoder 1 Kimberly.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw ... (total 13 args)
        if (argc != 13) {
            printf("Usage mode 1: encoder 1 input.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw\n");
            return 1;
        }
        mode_1_encoder(argc, argv);
    }
    else if (mode == 2) {
        mode_2_encoder(argc, argv);
    }
    else if (mode == 3) {
        mode_3_encoder(argc, argv);
    }
    else {
        printf("Mode %d not implemented yet.\n", mode);
    }

    return 0;
}