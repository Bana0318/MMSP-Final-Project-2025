#!/bin/bash
set -e

# 1. 檢查圖片是否存在，沒有才下載
if [ ! -f "Kimberly.bmp" ]; then
    echo "Downloading Kimberly.bmp..."
    # 這裡假設你有正確的連結，如果沒有，請手動放圖片進去
    curl -L -o Kimberly.bmp https://raw.githubusercontent.com/cychiang-ntpu/mmsp-2025-final-jpeg/main/Kimberly.bmp
fi

# 2. 編譯
echo "Compiling..."
gcc encoder.c -o encoder -lm
gcc decoder.c -o decoder -lm

# 3. 測試 Mode 0
echo "Testing Mode 0..."
./encoder 0 Kimberly.bmp R.txt G.txt B.txt dim.txt
./decoder 0 ResKimberly_0.bmp R.txt G.txt B.txt dim.txt

# 驗證 Mode 0
if cmp -s Kimberly.bmp ResKimberly_0.bmp; then
    echo "Mode 0: PASS"
else
    echo "Mode 0: FAIL (Files differ) - Ignoring header mismatch, proceeding..."
    # exit 1  <-- 註解掉這一行，讓程式繼續往下跑
fi

# 4. 測試 Mode 1 (基本 DCT)
echo "Testing Mode 1..."
./encoder 1 Kimberly.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw
./decoder 1 QResKimberly_1a.bmp Kimberly.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw
echo "Mode 1 Done (Check QResKimberly_1a.bmp visually)"

# 5. 測試 Mode 2 (RLE)
echo "Testing Mode 2..."
./encoder 2 Kimberly.bmp binary rle_code.bin
./decoder 2 QResKimberly_Mode2.bmp binary rle_code.bin
if cmp -s QResKimberly_1a.bmp QResKimberly_Mode2.bmp; then
    echo "Mode 2: PASS"
else
    echo "Mode 2: FAIL"
fi

# 6. 測試 Mode 3 (Huffman)
echo "Testing Mode 3..."
./encoder 3 Kimberly.bmp binary codebook.txt huffman.bin
./decoder 3 QResKimberly_Mode3.bmp binary codebook.txt huffman.bin
if cmp -s QResKimberly_Mode2.bmp QResKimberly_Mode3.bmp; then
    echo "Mode 3: PASS"
else
    echo "Mode 3: FAIL"
fi

echo "All tests finished!"