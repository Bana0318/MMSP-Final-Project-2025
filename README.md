# MMSP-Final-Project-2025
本專案實作了一個完整的 JPEG 影像壓縮系統，包含從底層 BMP 讀寫、DCT 頻率域轉換、量化處理，一直到 RLE 與 Huffman 熵編碼的完整流程。

---

## 1. 系統實作進度 (System Implementation Progress)
### 系統架構圖 (Block Diagram)
本系統依據作業要求分為四個主要模組 (Modes)，資料流處理如下：


graph TD
    A[Input: Kimberly.bmp] -->|Mode 0| B(BMP Parser & I/O)
    B --> C{Mode Selection}
    
    C -->|Mode 1| D[RGB to YCbCr]
    D --> E[2D-DCT Transform]
    E --> F[Quantization]
    F -->|Lossy| G[Output: Quantized Coefficients]
    F -->|Error| H[SQNR Calculation]

    C -->|Mode 2| I[Zig-Zag Scan]
    I --> J[DPCM on DC]
    J --> K[Run-Length Encoding (RLE)]
    K --> L[Binary Output: rle_code.bin]

    C -->|Mode 3| M[Frequency Statistics]
    M --> N[Build Huffman Tree]
    N --> O[Two-Pass Encoding]
    O --> P[Bitstream Packing]
    P --> Q[Binary Output: huffman.bin]

### 開發工作日誌 (Work Log)
```text
日期    實作階段     工作內容與進度描述
Day1    Mode 0      BMP 格式解析
                    - 定義 BMPHeader 與 BMPInfoHeader 結構 (注意 #pragma pack)。
                    - 解決 BMP Bottom-up 讀取順序與 Padding 對齊問題。
                    - 完成 R/G/B.txt 輸出與 ASCII 驗證。

Day 2   Mode 1      核心演算法實作
                    - 實作 RGB->YCbCr 色彩空間轉換與 Level Shifting。
                    - 實作 8x8 2D-DCT 與 IDCT 數學公式。
                    - 導入 Luminance 與 Chrominance 量化表，完成量化邏輯。
                    - 實作 SQNR 計算功能，驗證 DCT/IDCT 的數值精確度。

Day 3   Mode 2      RLE 壓縮機制
                    - 建立 Zig-zag 掃描索引表。
                    - 實作 DC 係數的 DPCM (差分編碼)。
                    - 設計 AC 係數的 Run-Length Encoding 邏輯。
                    - 完成二進位檔案 (rle_code.bin) 的寫入與讀取，驗證無損還原。

Day 4   Mode 3      Huffman 編碼 (Encoder)
                    - 設計 BitWriter 以處理非固定長度的位元寫入。
                    - 實作 Two-Pass 機制：第一遍統計頻率，第二遍進行編碼。
                    - 實作 Huffman Tree 建構演算法，並生成 codebook.txt。

Day 5   Mode 3      Huffman 解碼 (Decoder)
                    - 設計 BitReader 與 VLI (Variable Length Integer) 解碼器。
                    - 根據 Codebook 重建霍夫曼樹。
                    - 串接完整流程，成功還原影像並通過 cmp 驗證。

Day 6   Workflow    自動化與測試
                    - 撰寫 workflow.sh 腳本，自動化下載圖檔、編譯、執行全模式測試。
                    - 確認所有 Mode 產出的檔案與預期一致 (All Tests Passed)。
```

## 2.心得及感想 (Reflections)
這次的 MMSP 期末專題讓我對 JPEG 壓縮標準有了非常深刻的理解，從理論到實作的過程中，我遭遇並解決了許多挑戰：\
#### 1.資料結構的細節處理：
在實作 Mode 0 時，我深刻體會到處理二進位檔案（Binary File）的嚴謹性。一開始忽略了 Struct 的記憶體對齊（Byte Alignment）以及 BMP 的 Padding 機制，導致讀出來的圖片發生歪斜。這讓我學會了使用 #pragma pack 以及更謹慎地計算檔案偏移量。
#### 2.數學理論的程式化：
DCT 與量化是 JPEG 的靈魂。將課堂上的 $\Sigma$ 公式轉換為四層巢狀迴圈的 C 語言程式碼，讓我對「頻率域」的概念更加具體。特別是在計算 SQNR 時，觀察到經過量化後的訊號雖然有損失，但在視覺上卻難以分辨差異，這驗證了人眼對高頻訊號不敏感的特性。
#### 3.位元操作的挑戰：
Mode 3 的 Huffman Coding 是最具挑戰性的部分。因為 C 語言標準 I/O 是以 Byte 為單位，但 Huffman Code 長度不固定（例如 3 bits 或 12 bits）。我必須自行實作 BitWriter 和 BitReader，透過 Buffer 暫存與位移運算來湊滿 8 bits 寫入。這不僅訓練了我的邏輯，也讓我對電腦底層的資料儲存有了更直觀的認識。
#### 4.系統整合的成就感：
最讓我有成就感的時刻，是寫完 workflow.sh 並看到終端機印出 "ALL TESTS PASSED SUCCESSFULLY" 的那一刻。這證明了我的 Encoder 與 Decoder 完美對接，且無損壓縮模式下的檔案與原圖完全一致（Binary Identical）。這份作業不僅是演算法的練習，更是一次完整的軟體工程實踐。
