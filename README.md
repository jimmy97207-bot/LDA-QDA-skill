# LDA／QDA 區別分析 Skill

Claude Skill：針對已知類別的表格資料（CSV），以 R 完成「變數介紹 → 資料檢查 → 分組 EDA 與假設診斷 → 模型可估計性 → 重複分層交叉驗證比較 LDA／QDA → 繁體中文報告」的區別分析流程，並輸出 Markdown 報告、CSV 表格、模型檔與圖檔。

---

## 功能特色

- **半自動決策**：Y 由使用者指定；AI 自動提出 X 清單，使用者確認後才執行分析。資料乾淨時，確認完 Y、X 就一路跑完；其餘只在遇到關卡（疑似同一對象多筆、有空格、某類少於 5 筆、兩種方法都算不出來）時詢問或停止，並用白話說明原因。
- **先確認 Y 與 X**：`--plan` 列出 Y 候選與預設 X 清單，使用者只需確認或修改；名稱是編號、非數值、常數的欄位自動排除，數值像編號的欄位保留並請使用者確認，不自行刪除；讀取欄位說明檔產生變數介紹表。
- **資料檢查只標記不修改**：重複列、類別比例、$n_k\le p$ 的類別、組內零變異、X 是否滿秩、組內 1.5×IQR 離群列。
- **假設診斷**：各組 Mardia 多元偏態與峰態、Box's M 共變異數同質性；以 base R 直接計算，結果與 MVN 6.3／biotools 4.3 一致。報告附兩者的合併解讀表；Mardia 顯著時自動加註。
- **公平比較**：LDA、QDA 使用相同的 10 組分層切分；先驗機率與模型都只用訓練折估計；逐折檢查可估計性，任何一折失敗就標為「驗證不完整」，不參加比較。
- **用考試成績選方法**：Mardia、Box's M 只作診斷參考；選模依分層 5-fold 重複 10 次的平均 Macro recall，較高者勝出，完全相同時選 LDA；另報告逐次勝出次數說明穩定度；選定後用全部資料建立最終模型。
- **固定、不改研究問題**：固定分層 5-fold × 10 次，不開放調整；某類少於 5 筆就停止；不合併或刪除類別、不自行改模型。
- **缺失值不洩漏**：X 有空格時，使用者同意補值後，只用每份練習資料計算中位數再補值（`--impute median`），50 輪的中位數都有紀錄。
- **兩種報告**：`report.html`（R 直接寫出，一定會有）與 `report.qmd`（含 R 程式碼區塊，Render 時重新計算）。

## 評估指標

| | 內容 |
|---|---|
| 主要指標（選模） | Macro recall（各類別召回率等權平均，即平衡正確率） |
| 其他報告指標 | 正確率、Macro-F1、平均答對筆數與範圍、逐次勝出次數、各類別 recall／precision／F1、混淆矩陣（列百分比） |
| 切分方式 | 分層 5-fold 重複 10 次（seed 2026），兩模型共用；固定，不開放調整 |
| 預測來源 | 折外預測（每次重複中每筆觀測恰好被預測一次）；表格為 10 次平均，各類別表現與混淆矩陣為 10 次合併 |

## 檔案結構

```
├── README.md
├── .gitignore                               # 不上傳的暫存檔（驗收輸出、Quarto 暫存等）
├── discriminant-analysis.zip                # Skill 本體（上傳到 Claude Skill）
├── discriminant-analysis/                   # Skill 本體（安裝時使用此資料夾）
│   ├── SKILL.md                             # Skill 主說明（Claude 讀取的入口）
│   ├── requirements.txt                     # R 與套件版本
│   ├── load_project.R                       # 在 R／RStudio 載入分析函式
│   ├── scripts/
│   │   ├── lda_qda_pipeline.R               # 主程式：命令列一鍵執行完整流程
│   │   └── run_analysis.R                   # 核心函式：模型、診斷、圖、報告、輸出
│   ├── references/
│   │   ├── decision-tree.md                 # 資料問題與模型可行性的判斷規則
│   │   ├── lda-qda-guide.md                 # 方法比較、公式、輸出解讀、參考資料
│   │   ├── evaluation.md                    # 驗證設計、指標與模型比較規則
│   │   ├── pipeline_output_iris.txt         # 主程式在 iris 上的實際輸出
│   │   └── pipeline_output_iris_missing.txt # 主程式在 iris_missing 上的實際輸出
│   ├── assets/
│   │   ├── iris.csv                         # demo 資料
│   │   ├── iris_dictionary.csv              # 欄位說明
│   │   ├── iris_missing.csv                 # 第二個 demo：挖掉 5% 的 X
│   │   └── iris_missing_dictionary.csv
│   ├── demo_iris_missing/                   # 有缺失值 demo 的分析輸出（含 imputation.csv）
│   ├── demo_iris/                           # iris demo 的分析輸出
│   │   ├── report.html                      #   繁中報告（直接產生）
│   │   ├── report.qmd                       #   含 R 程式碼的報告原始檔
│   │   ├── report-quarto.html               #   report.qmd Render 的結果
│   │   ├── report.md                        #   Markdown 版報告
│   │   ├── report.txt                       #   主控台輸出
│   │   ├── summary.json                     #   關鍵數字
│   │   ├── figures/                         #   6 張圖
│   │   └── *.csv、*.rds                     #   表格、模型與設定
│   └── tests/
│       └── validate.R                       # 程式驗收（136 項檢查）
└── materials/                               # Skill 說明報告
    ├── make-skill-tutorial.qmd              # 上台報告原始檔（Quarto）
    ├── make-skill-tutorial.html             # 上台報告 Render 輸出
    └── VALIDATION.md                        # 驗收紀錄
```

- **`discriminant-analysis.zip`**：Skill 本體，安裝到 Claude 時使用。
- **`materials/`**：說明報告，介紹流程設計、決策邏輯與 iris 示範結果；`.qmd` 可用 Quarto 重新產生 `.html`。

## 程式說明與執行流程

主程式為 **`lda_qda_pipeline.R`**，完整流程都從這裡執行；計算都在 `run_analysis.R` 的函式內。

| 程式 | 用途 | 輸入 | 輸出 |
|---|---|---|---|
| **`lda_qda_pipeline.R`** | 命令列入口：`--plan` 提出 Y、X 清單；變數介紹、完整分析與輸出 | 任何 CSV + `--target` | `lda_qda_<檔名>/` 資料夾 |
| `run_analysis.R` | 核心函式，可在 R 中個別呼叫 | 資料框或 CSV 路徑 | R 物件；`save_analysis()`／`save_figures()` 寫檔 |
| `tests/validate.R` | 驗收：數值正確性、讀寫、重現性、失敗處理、命令列 | 內建 iris | `validation/<時間>/` |

```
assets/iris.csv (+ iris_dictionary.csv)
                │
                ▼
   lda_qda_pipeline.R --plan             （先列出 Y 候選與 X 清單，由使用者確認）
                │
                ▼
   lda_qda_pipeline.R --target Species
   1. 讀檔與欄位說明、套用確認過的 X（名稱是編號／非數值／常數的欄位排除）、變數介紹
   2. 資料檢查（只標記）：重複、類別比例、n_k ≤ p、零變異、滿秩、IQR 離群
   3. 分組描述統計；各組 Mardia、Box's M
   4. 檢查 LDA／QDA 在全資料與每個訓練折是否可估計
   5. 相同切分的分層 5-fold × 10 次：折外預測、正確率、Macro recall、Macro-F1、每次結果
   6. 選方法：平均 Macro recall 較高者勝出，相同選 LDA；用全部資料建立最終模型
   7. 輸出 report.html、report.qmd、report.md、report.txt、summary.json、figures/、CSV、RDS
                │
                ▼
            demo_iris/
```

## 安裝與使用

### 作為 Claude Skill

將 `discriminant-analysis/` 資料夾放進 Claude 的 skills 目錄，或在 Claude 的 Skill 設定中上傳 `discriminant-analysis.zip`。

### 直接執行

```bash
# 只需要 R（>= 4.3）與內建的 MASS；不需要另外安裝套件
cd discriminant-analysis

# 執行 demo（demo_iris/ 已附上跑好的結果；輸出資料夾已存在時會停止，所以用新名稱）
Rscript scripts/lda_qda_pipeline.R assets/iris.csv --target Species --out my_iris_result

# 分析自己的資料
Rscript scripts/lda_qda_pipeline.R my_data.csv --target group --dict my_data_dictionary.csv --out lda_qda_my_data

# 驗收
Rscript tests/validate.R
```

> **可重現性**：10 組切分都由 `--seed` 決定（預設 2026）。相同資料、參數與 R／MASS 版本可得到相同結果。

| 參數 | 說明 |
|---|---|
| `--target` | 分類欄位 Y（必填；`--plan` 時可省略） |
| `--plan` | 只列出 Y 候選、預設 X 與要確認的欄位，不分析 |
| `--x` | 逗號分隔的預測變數；省略時自動選數值欄位，排除名稱是編號的欄位；數值像編號的欄位保留並提示確認 |
| `--drop` | 自動選 X 時額外排除的欄位 |
| `--dict` | 欄位說明檔；未指定時自動尋找 `<檔名>_dictionary.csv` |
| `--seed` | 亂數種子（預設 2026） |
| `--impute` | `none`（預設，X 有空格就停止）或 `median`（每份練習資料內補中位數） |
| `--out` | 輸出資料夾（預設 `lda_qda_<檔名>`）；已存在時停止，避免覆蓋 |

欄位說明檔格式（CSV，UTF-8）：`欄位,說明,單位,角色`；`欄位 = __dataset__` 的那一列寫整份資料的描述。

資料要求：每列一筆觀測、第一列為欄名；Y 至少兩類、不能有空格；X 為數值、沒有無限值（有空格時須加 `--impute median`）；每一類至少 5 筆。程式不會自行刪除或轉換資料，不符合時會停止並說明原因。

## 示範結果

| 資料 | 診斷 | LDA 平均 Macro recall（答對） | QDA | 逐次勝出（LDA／QDA／相同） | 選擇 |
|---|---|---|---|---|---|
| iris（150 筆、3 類、4 個 X） | Mardia 皆未拒絕；Box's M p = 3.35 × 10⁻²⁰ | 0.9793（146.9，146–147） | 0.9733（146.0，145–147） | 6／0／4 | LDA |
| iris_missing（同上，挖掉 5% 的 X，`--impute median`） | 以 120 筆完整觀測計算 | 0.9400（141.0，139–143） | 0.9553（143.3，140–145） | 1／8／1 | QDA |

第一個 demo 示範 Box's M 拒絕共變異數同質，QDA 卻沒有分得比較準（10 次中一次都沒贏）；第二個 demo 示範缺失值關卡與各折內的中位數補值，QDA 多數時候較高，但也有 1 次反過來。

## 主要參考資料

- Fisher, R. A. (1936). The use of multiple measurements in taxonomic problems. *Annals of Eugenics*, 7(2), 179–188.
- Box, G. E. P. (1949). A general distribution theory for a class of likelihood criteria. *Biometrika*, 36(3/4), 317–346.
- Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with applications. *Biometrika*, 57(3), 519–530.
- Johnson, R. A., & Wichern, D. W. (2007). *Applied Multivariate Statistical Analysis* (6th ed.). Pearson.
- Hastie, T., Tibshirani, R., & Friedman, J. (2009). *The Elements of Statistical Learning* (2nd ed.). Springer.
- Venables, W. N., & Ripley, B. D. (2002). *Modern Applied Statistics with S* (4th ed.). Springer.

**資料來源**：iris 取自 R `datasets::iris`；資料採集歸於 Anderson（1935），並用於 Fisher（1936）的分析。
