#!/bin/bash

# 設定遇到錯誤立即停止
set -e

# 定義顏色 (讓輸出比較好看)
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo "========================================"
echo "      MMSP 2025 JPEG Project Workflow   "
echo "========================================"

# --- 1. 準備階段 ---
echo -e "${GREEN}[Step 1] Preparing environment...${NC}"

# 下載圖片 (如果不存在)
if [ ! -f "Kimberly.bmp" ]; then
    echo "Downloading Kimberly.bmp..."
    curl -o Kimberly.bmp https://raw.githubusercontent.com/cychiang-ntpu/mmsp-2025-final-jpeg/main/Kimberly.bmp
fi

# 清除舊的執行檔與輸出檔
rm -f encoder decoder *.txt *.bin *.raw *.bmp
# 注意：保留剛剛下載的 Kimberly.bmp (但其他 .bmp 清除)
# 為了避免誤刪剛下載的圖，我們把 Kimberly.bmp 備份一下，或者只刪除特定 pattern
# 這裡簡單處理：只刪除 Output 開頭的 bmp
rm -f Res*.bmp QRes*.bmp

# 編譯程式
echo "Compiling C code..."
gcc encoder.c -o encoder -lm
gcc decoder.c -o decoder -lm

if [ -f "encoder" ] && [ -f "decoder" ]; then
    echo "Compilation Successful."
else
    echo -e "${RED}Compilation Failed!${NC}"
    exit 1
fi

# --- 2. 測試 Mode 0 (BMP 讀寫驗證) ---
echo ""
echo "========================================"
echo -e "${GREEN}[Step 2] Testing Mode 0 (BMP I/O Check)${NC}"
echo "========================================"

# 執行 Encoder Mode 0
./encoder 0 Kimberly.bmp R.txt G.txt B.txt dim.txt

# 執行 Decoder Mode 0
./decoder 0 ResKimberly_0.bmp R.txt G.txt B.txt dim.txt

# 驗證
if cmp -s Kimberly.bmp ResKimberly_0.bmp; then
    echo -e "Mode 0 Check: ${GREEN}PASS${NC} (Files are identical)"
else
    echo -e "Mode 0 Check: ${RED}FAIL${NC} (Files differ)"
    exit 1
fi


# --- 3. 測試 Mode 1 (DCT & Quantization) ---
echo ""
echo "========================================"
echo -e "${GREEN}[Step 3] Testing Mode 1 (DCT & Quantization)${NC}"
echo "========================================"

# 定義變數名稱以保持整潔
QT_Y="Qt_Y.txt"
QT_CB="Qt_Cb.txt"
QT_CR="Qt_Cr.txt"
DIM="dim.txt"
QF_Y="qF_Y.raw"
QF_CB="qF_Cb.raw"
QF_CR="qF_Cr.raw"
EF_Y="eF_Y.raw"
EF_CB="eF_Cb.raw"
EF_CR="eF_Cr.raw"

# 執行 Encoder Mode 1
# 這會產生 qF (量化係數) 和 eF (量化誤差)
echo "Running Encoder Mode 1..."
./encoder 1 Kimberly.bmp $QT_Y $QT_CB $QT_CR $DIM \
  $QF_Y $QF_CB $QF_CR \
  $EF_Y $EF_CB $EF_CR

# --- Mode 1(a): 失真解碼 (Lossy) ---
echo "Running Decoder Mode 1(a) [Lossy]..."
./decoder 1 QResKimberly_1a.bmp Kimberly.bmp \
  $QT_Y $QT_CB $QT_CR $DIM \
  $QF_Y $QF_CB $QF_CR

# --- Mode 1(b): 無損解碼 (Lossless, using Error files) ---
echo "Running Decoder Mode 1(b) [Lossless]..."
./decoder 1 ResKimberly_1b.bmp \
  $QT_Y $QT_CB $QT_CR $DIM \
  $QF_Y $QF_CB $QF_CR \
  $EF_Y $EF_CB $EF_CR

# 驗證 Mode 1(b) 是否完美還原
if cmp -s Kimberly.bmp ResKimberly_1b.bmp; then
    echo -e "Mode 1(b) Lossless Check: ${GREEN}PASS${NC}"
else
    echo -e "Mode 1(b) Lossless Check: ${RED}FAIL${NC}"
    exit 1
fi


# --- 4. 測試 Mode 2 (RLE Compression) ---
echo ""
echo "========================================"
echo -e "${GREEN}[Step 4] Testing Mode 2 (RLE Compression)${NC}"
echo "========================================"

# Encoder Mode 2 (Binary Output)
echo "Running Encoder Mode 2..."
./encoder 2 Kimberly.bmp binary rle_code.bin

# Decoder Mode 2
echo "Running Decoder Mode 2..."
./decoder 2 QResKimberly_Mode2.bmp binary rle_code.bin

# 驗證 Mode 2
# 理論上：Mode 2 (DCT->Quant->RLE->DeRLE->DeQuant->IDCT) 的結果
# 應該要跟 Mode 1a (DCT->Quant->DeQuant->IDCT) 的結果 完全一樣
if cmp -s QResKimberly_1a.bmp QResKimberly_Mode2.bmp; then
    echo -e "Mode 2 Consistency Check (vs Mode 1a): ${GREEN}PASS${NC}"
else
    echo -e "Mode 2 Consistency Check: ${RED}FAIL${NC}"
    echo "Note: Mode 2 output should match Mode 1a output."
    exit 1
fi


# --- 5. 測試 Mode 3 (Huffman Compression) ---
echo ""
echo "========================================"
echo -e "${GREEN}[Step 5] Testing Mode 3 (Huffman Compression)${NC}"
echo "========================================"

# Encoder Mode 3 (Binary Output)
echo "Running Encoder Mode 3..."
./encoder 3 Kimberly.bmp binary codebook.txt huffman_code.bin

# Decoder Mode 3
echo "Running Decoder Mode 3..."
./decoder 3 QResKimberly_Mode3.bmp binary codebook.txt huffman_code.bin

# 驗證 Mode 3
# Huffman 是無失真編碼 (針對 Quantized data 而言)，
# 所以 Mode 3 的還原圖應該要跟 Mode 2 一模一樣
if cmp -s QResKimberly_Mode2.bmp QResKimberly_Mode3.bmp; then
    echo -e "Mode 3 Consistency Check (vs Mode 2): ${GREEN}PASS${NC}"
else
    echo -e "Mode 3 Consistency Check: ${RED}FAIL${NC}"
    exit 1
fi

echo ""
echo "========================================"
echo -e "${GREEN}ALL TESTS PASSED SUCCESSFULLY!${NC}"
echo "Generated files:"
ls -lh huffman_code.bin rle_code.bin Kimberly.bmp QResKimberly_Mode3.bmp
echo "========================================"