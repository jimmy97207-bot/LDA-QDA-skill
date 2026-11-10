---
name: discriminant-analysis
description: 用 R 做 LDA／QDA 區別分析：以多個數值變數判斷已知類別，含資料檢查、Mardia 與 Box's M 診斷、重複交叉驗證比較兩法，輸出繁體中文報告。使用者提到區別分析、判別分析、LDA、QDA、Fisher 判別時使用；不適用於主題模型的 LDA。
---

# LDA／QDA 區別分析 Skill

協助使用者以 LDA／QDA 完成已知類別的**分類**，固定流程：

**① 使用者指定 Y → ② AI 提出 X 清單請使用者確認，並檢查資料 → ③ 執行分析程式 → ④ 資料能直接用嗎 → ⑤ 檢查資料假設（Mardia、Box's M，只作診斷） → ⑥ LDA、QDA 算得出來嗎 → ⑦ 分層 5-fold 重複 10 次比較兩種方法 → ⑧ 用全部資料建立最終模型 → ⑨ 交付報告（report.html ＋ report.qmd）**

採半自動：Y 由使用者指定；AI 自動提出 X 清單，使用者確認後才執行分析（使用者已指定 X 時不再問）。資料乾淨時，確認完 Y 和 X 就一路跑完；其餘只在遇到下方關卡 C–F 時詢問或停止。設計原則是**固定且可解釋的流程**：AI 不合併或刪除類別、不自行改模型、不改變研究問題；做不到就停下來說明原因。已提供的資訊不重複詢問。

## 環境準備（每次開始前）

這個 skill 用 R 執行。開始分析前先檢查：

```bash
Rscript -e 'cat(R.version.string, as.character(packageVersion("MASS")), "\n")'
```

- 有輸出 R 與 MASS 版本：直接進行。
- 找不到 R：先告訴使用者「需要先安裝 R，約一兩分鐘」，再嘗試安裝（Debian／Ubuntu 環境：`apt-get install -y r-base-core r-cran-mass`，必要時加 `sudo`）；中文欄位或檔名出錯時，執行指令前加 `LC_ALL=C.UTF-8`。
- 安裝失敗（例如沒有網路或沒有權限）：不要改用其他語言自行重寫分析，也不要假裝有結果；明確告訴使用者目前環境無法執行 R，建議改在已安裝 R（>= 4.3）的電腦上用 Claude Code 執行。
- 要產生 `report-quarto.html` 才需要 Quarto 與 R 套件 knitr、rmarkdown；沒有時只交付 `report.html` 與 `report.qmd`，並說明 `.qmd` 可在 RStudio 按 Render。

## 快速索引

| 需要什麼 | 去哪裡 |
|---|---|
| 一支跑完 ③–⑨ 並輸出報告、表格、圖（`--plan` 供 ①② 確認 Y、X） | `scripts/lda_qda_pipeline.R`（命令列） |
| 在 R 裡逐步呼叫的函式 | `scripts/run_analysis.R`（`run_pipeline()`、`save_analysis()`、`save_figures()`） |
| 資料問題、何時詢問、模型不可行時怎麼辦 | `references/decision-tree.md` |
| LDA／QDA 原理、公式、R 輸出解讀、常見問題 | `references/lda-qda-guide.md` |
| 切分、指標、模型比較與報告規則 | `references/evaluation.md` |
| demo 的實際主控台輸出 | `references/pipeline_output_iris.txt`、`references/pipeline_output_iris_missing.txt` |
| demo 資料與欄位說明 | `assets/` |
| 程式驗收 | `tests/validate.R` |

## 何時使用 / 不使用

- 使用：Y 是至少兩類的已知類別、X 是數值變數、觀測彼此獨立；想知道「哪些變數能區分類別」或「新樣本該分到哪一類」。
- 先說明範圍、不直接套用：同一受試者多筆、群組、時間序列（本版只有隨機分層 5-fold 重複 10 次，不做分組驗證）；X 以類別變數為主（本版只接受數值 X）。
- 不使用：無監督分群（沒有已知類別）、文字主題模型（Latent Dirichlet Allocation）。

## 停下來問使用者或停止的關卡（A–F）

怎麼問：

- 用白話說明「為什麼要問」，不要只丟專有名詞；開工前的 A、B、C 合併成一次詢問，其餘關卡一次只問一個。
- 需要選擇時給選項，**第一個放建議選項並標「（建議）」**；有可點選的選項工具（如 AskUserQuestion）時就用它。
- 開工前的 A、B、C 一次問完；D 在程式停止後才問；E、F 不給修改資料的選項，直接停止並說明原因。
- 選完後回到主流程重新執行；停止時說明原因，並告訴使用者需要什麼條件才能繼續。
- AI **不**合併或刪除類別、**不**自行刪變數或改模型、**不**刪除有空格的列——這些都會改變原本的研究問題或樣本。

| 關卡 | 什麼時候 | 為什麼（對使用者說的話） | 怎麼處理 |
|---|---|---|---|
| A | 沒指定要分類的欄位（Y） | 要分什麼類是研究問題本身，只能由你決定 | **詢問**，不猜測；列出 `--plan` 找到的候選欄位讓使用者點選，也可自己輸入 |
| B | 確認用來判斷的欄位（X） | 有些欄位不該拿來用，例如編號、結果發生後才知道的資訊 | 先列出預設清單請使用者**確認或修改**（見下方問法）；名稱是編號、文字、常數的欄位自動排除；**數值像編號的欄位保留並標出來請使用者確認**，不自行刪除 |
| C | 疑似同一對象有多筆：腳本偵測到編號欄位有重複值、有日期／時間欄位，或使用者描述資料是多次紀錄 | 同一個人可能同時出現在練習題和考試題，分數會灌水；正確做法是分組驗證，本版不支援 | **暫停**，先確認研究單位是「對象」還是「每一筆紀錄」；確認每一列是獨立個體後才繼續 |
| D | 資料有空格（缺失值） | 模型沒辦法處理空白 | **X 有空格**：先報告空格數量與位置，再問：① 補中位數：`--impute median`，**只用每份練習資料計算中位數**，再補該輪的練習與考試資料（建議） ② 停止分析。不在切分前用全部資料補值，那會讓考試資料洩漏到練習裡。**Y 有空格**：**停止**並列出是哪幾列；Y 是類別答案，不補、不猜，也不自行刪除，由使用者確認後再執行 |
| E | 某一類少於 5 筆 | 分層 5-fold 需要每一類至少 5 筆，才能每一份考試都有這一類 | **停止**並說明哪一類不足；不減少折數、不合併或刪除類別 |
| F | LDA、QDA 都算不出來 | 通常是某個欄位可以由其他欄位算出來（例如同時有公分和公尺，或總分＝各科相加），或總筆數相對欄位數太少，導致共變異數矩陣無法計算 | **停止**並說明原因（例如哪些欄位線性相依）；由使用者決定要不要換一組 X 重新執行 |

A、B 的問法：先執行 `Rscript scripts/lda_qda_pipeline.R data.csv --plan`（不分析，只列出 Y 候選、預設 X、排除原因與要確認的欄位；使用者選了別的 Y 時加 `--target <Y>` 再跑一次），再一次問兩題，不要把所有欄位逐一列成複選：

> 1. 要分類的是哪一欄？（列出 2–4 個候選，第一個標「建議」；可自己輸入）
> 2. 我打算用這 N 個欄位：……；已排除：……（原因）；⚠ 請確認：○○ 每筆都不同而且都是整數，看起來像編號，它是一般的測量值嗎？
>    - 照這樣用（建議）
>    - 我要修改（輸入要加入或排除的欄位）

使用者修改時用 `--x` 或 `--drop` 帶入。

C 的問法：「資料裡同一個對象（例如同一位病人、同一台機器）是不是有好幾筆？」答「不是」就繼續。腳本的偵測只是提示，沒有編號欄位時偵測不到，所以使用者描述資料來源時也要留意。

另外也要在開工前確認：分析目的（解釋類別差異或預測新樣本，影響報告重點）、有沒有欄位說明（有就用 `--dict` 帶入；沒有就依欄位名稱推測並標明是推測）。

## 不問、自動處理並在報告註明

| 情況 | 自動怎麼做 |
|---|---|
| 只有 QDA 算不出來（某類 $n_k\le p$ 或某折的組內共變異數奇異） | 只評估 LDA，報告寫明 QDA 不可行的原因 |
| 常數欄位、名稱是編號的欄位 | 自動排除在 X 之外，報告列出排除原因 |
| 某個方法在任何一次重複的任何一折失敗 | 標為「驗證不完整」，不參加比較、不提供分數；只評估另一個方法 |
| 重複資料、特別大或特別小的值（組內 1.5×IQR） | 只標記，全部保留 |
| 各類別數量差很多 | 以「每一類各自答對多少」的平均（Macro recall）為主，不只看總正確率 |
| Mardia、Box's M 結果 | 寫進報告當參考，不用來選方法；報告附兩者的合併解讀表，並標出本資料屬於哪一格 |
| Mardia 顯著（某類 p < 0.05） | 報告自動加註：該類不太像多元常態，模型前提可能不成立、表現可能受影響；後驗機率與 Box's M 要保守解讀；可另外比較不需要常態假設的方法 |
| 不是分類問題（Y 是連續數值）、Y 只有一類 | 說明 LDA／QDA 不適用，結束 |

## 一鍵執行

優先用腳本產生結果，再依結果解說；**不要憑記憶寫數字**。

```bash
# 在 skill 資料夾內執行；只需要 R 與 MASS
Rscript scripts/lda_qda_pipeline.R data.csv --plan                 # 關卡 A、B：只列出 Y 候選與預設 X
Rscript scripts/lda_qda_pipeline.R data.csv --target <Y> [--x a,b,c] [--drop id] [--dict data_dictionary.csv] [--seed 2026] [--impute none|median] [--out 資料夾]
```

| 參數 | 說明 |
|---|---|
| `--target` | 分類欄位 Y（必填；`--plan` 時可省略） |
| `--plan` | 只列出 Y 候選、預設 X、排除原因與要確認的欄位，不分析（關卡 A、B 使用） |
| `--x` | 逗號分隔的預測變數；省略時自動選數值欄位，排除名稱是編號的欄位；數值像編號的欄位保留並提示確認 |
| `--drop` | 自動選 X 時額外排除的欄位 |
| `--dict` | 欄位說明 CSV（`欄位,說明,單位,角色`；`__dataset__` 列寫資料描述）；省略時自動找 `<檔名>_dictionary.csv` |
| `--seed` | 亂數種子（預設 2026）。驗證方式固定為分層 5-fold 重複 10 次，沒有 `--folds`、`--repeats` 參數 |
| `--impute` | `none`（預設；X 有空格就停止）或 `median`（在每份練習資料內算中位數補值）；只有在使用者於關卡 D 選擇補值後才用 `median` |
| `--out` | 輸出資料夾（預設 `lda_qda_<檔名>`）；**已存在就停止**，避免覆蓋舊結果 |

輸出：

| 檔案 | 用途 |
|---|---|
| `report.html` | 繁中報告（R 直接寫出、圖片內嵌，不需要 Quarto，**一定會有**）——交給使用者的主要檔案 |
| `report.qmd` | 同內容但含 R 程式碼區塊；在 RStudio 按 Render 會重新讀資料、重新計算、重新畫圖，輸出 `report-quarto.html`（需要 Quarto、knitr、rmarkdown） |
| `input_data.csv`、`run_analysis.R` | 給 `report.qmd` 用的資料與程式副本，讓輸出資料夾可以獨立重現 |
| `report.md`、`report.txt`、`summary.json` | Markdown 版報告、主控台紀錄、關鍵數字 |
| `figures/*.png`（6 張）、各表 CSV、`analysis.rds`／`final_model.rds` | 圖、表格、模型 |

環境有 Quarto 時，可在輸出資料夾執行 `quarto render report.qmd` 一併產生 `report-quarto.html`；沒有 Quarto 時交付 `report.html` 與 `report.qmd`，並告訴使用者在 RStudio 按 Render 即可。

在 R 裡使用（RStudio 互動或需要客製時）：

```r
source("load_project.R", encoding = "UTF-8")   # 只載入函式，不會自動分析
result <- run_pipeline("assets/iris.csv", y = "Species",
                       x = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width"),
                       seed = 2026)   # 固定 5 折 × 10 次
result$metrics; result$diagnostics; result$data_checks
cat(result$report_lines, sep = "\n")
save_analysis(result, "output/iris-001")                       # 表格、模型、report.md
save_figures(result, read.csv("assets/iris.csv"), "output/iris-001/figures")
write_report_html(result, "output/iris-001")                   # report.html（圖片取自 figures/）
```

## 報告章節與要向使用者說明的重點

1. **資料與設定**：資料說明、Y 與各類別筆數、X 清單、驗證方式、主要指標、先驗機率；變數介紹表（說明、單位、角色、型態、缺失數）。
2. **資料檢查**（只標記、不修改）：重複列、類別比例、$n_k\le p$ 的類別、組內零變異、X 是否滿秩、組內 1.5×IQR 離群列數。觸發問題時依 `decision-tree.md` 區分「可繼續但要揭露」「需使用者決定」「無法建模」。
3. **模型驗證表現**：LDA／QDA 的狀態、10 次重複的平均正確率、平均 Macro recall（含標準差）、平均 Macro-F1、平均答對筆數與範圍；每一次重複的結果與勝出次數；各類別 recall／precision／F1 與混淆矩陣（10 次合併，以百分比呈現）。
4. **假設診斷**：各組 Mardia 偏態與峰態、Box's M（df = $p(p+1)(K-1)/2$）；兩者的合併解讀表（見下）與自動註記。
5. **分組描述統計**、**圖檔**、**警告與錯誤**、**範圍與限制**。

腳本遇到 Y 有空格、X 有空格（未加 `--impute median`）、非有限值、非數值 X、類別數 < 2、某類少於 5 筆時會**停止並說明原因**，不會自行刪除或插補；依上方關卡 D、E 的方式處理，取得使用者決定後再跑。

## 向使用者解說時務必講到

- **資料與目標**：Y 是什麼、有幾類各幾筆、用了哪些 X、排除了哪些欄位與原因。
- **資料檢查**：重複與離群只標記、全部保留；相同測量值不能單憑資料表認定為同一個體。
- **診斷不是選模開關**：Mardia 不顯著只是「未拒絕多元常態」，不是已證明；Box's M 對偏離常態很敏感，顯著時不一定代表共變異數真的差很多，也**不等於 QDA 會預測得比較好**——QDA 多估 $(K-1)p(p+1)/2$ 個參數，可能抵消好處。以相同切分的折外表現比較。
- **驗證設計**：分層 5-fold 重複 10 次（每次重新隨機分組，共 50 輪練習／考試）；兩模型每一次都用同一組切分；先驗與模型都只用訓練折估計；這組 CV 同時用來比較，分數不是獨立測試表現。次數固定，不開放調整。
- **選擇規則**：比較兩者 10 次重複的**平均 Macro recall，較高者勝出；完全相同時選較簡單的 LDA**。方法必須 10 次全部完成才參加比較。報告寫出兩者平均答對幾筆（含範圍）與逐次勝出次數。
- **最終模型**：選定後，用全部分析資料重新配適該模型（`final_model.rds`；有補值時也用全部資料的中位數）。交叉驗證分數只用來比較，不宣稱是獨立測試表現。
- **不確定性**：用逐次勝出次數說明穩定度——落後的方法也贏了好幾次時，要明說兩者實際表現接近，**不代表較高者真的比較好**，不能宣稱統計優越性；勝出次數只描述穩定度，不是檢定。
- **解釋界線**：LDA 有判別軸（最多 $\min(p,K-1)$ 個），QDA 沒有共同線性判別軸；判別係數大小不等於變數重要性；組間差異不是因果。
- **未完成就說未完成**：模型不可估計、折失敗、診斷無法計算時，保留原因，不產生替代數字，也不只用成功折排名。

## Mardia 與 Box's M 的合併解讀

報告會附這張表並標出本資料屬於哪一格（`diagnostic_notes()`）。兩個檢定都只是診斷，**不用來選方法**。

| Mardia | Box's M | 解讀 |
|---|---|---|
| 不顯著 | 不顯著 | 未發現明顯違反假設；仍依交叉驗證結果選方法 |
| 不顯著 | 顯著 | 各類散開方式可能不同，QDA 有機會較好 |
| 顯著 | 不顯著 | 資料不像鐘形，算出的機率可能不準 |
| 顯著 | 顯著 | Box's M 可能只是受非鐘形影響，不能直接下結論 |

任何一類的 Mardia 偏態或峰態 p < 0.05 時，報告自動加註受影響的類別，說明：交叉驗證仍可評估 LDA／QDA 的表現，但模型的前提可能不成立、表現可能因此受影響；算出的機率和 Box's M 要保守解讀；可另外比較不需要常態假設的方法（例如羅吉斯迴歸、決策樹），超出本 skill 範圍。AI 最後總結時要把這點告訴使用者，但不自行執行其他方法，也不自行轉換變數。

## 分析原則

- 所有數值、圖表與模型結果必須由程式實際產生。
- 保留原始資料，記錄分析使用的欄位、設定與處理方式。
- 不刪除觀測值、不合併或刪除類別、不轉換變數；X 的缺失值只在使用者同意後，以各訓練折中位數補值；Y 的缺失值不補，停止並詢問。
- 不依單一假設檢定的 p 值直接決定 LDA 或 QDA。
- 不把訓練資料（再代入）的正確率當成泛化能力。
- 模型或驗證失敗時，保留原因，不虛構結果或默默省略失敗折。
- 不擅自改用其他分類方法；QDA 不可行時可提及正則化區別分析，但需使用者同意擴充範圍。

## 模型可行性決策

```
X 線性相依／合併共變異數奇異   ──► LDA、QDA 都不可估計 → 停止，說明原因（關卡 F）
某類 n_k ≤ p 或組內共變異數奇異 ──► QDA 不可估計 → 只評估 LDA，說明 QDA 不可行原因
兩者皆可估計                    ──► 相同切分做 10 次分層 5-fold
某方法任一次任一折失敗          ──► 標為「驗證不完整」，不參加比較，只評估另一個
平均 Macro recall 較高者勝出；完全相同 ──► 選 LDA
選定後 ──► 用全部資料重新配適，存成 final_model.rds
落後者也贏了好幾次 ──► 照規則選，但明說兩者表現接近
```

## LDA／QDA 一頁重點

| | LDA | QDA |
|---|---|---|
| 假設 | 各類 $N_p(\boldsymbol\mu_k,\boldsymbol\Sigma)$，共用 $\boldsymbol\Sigma$ | 各類 $N_p(\boldsymbol\mu_k,\boldsymbol\Sigma_k)$ |
| 判別函數 | $\delta_k=\mathbf x^\top\boldsymbol\Sigma^{-1}\boldsymbol\mu_k-\tfrac12\boldsymbol\mu_k^\top\boldsymbol\Sigma^{-1}\boldsymbol\mu_k+\log\pi_k$ | $\delta_k=-\tfrac12\log\lvert\boldsymbol\Sigma_k\rvert-\tfrac12(\mathbf x-\boldsymbol\mu_k)^\top\boldsymbol\Sigma_k^{-1}(\mathbf x-\boldsymbol\mu_k)+\log\pi_k$ |
| 邊界 | 超平面 | 二次曲面 |
| 參數個數 | $Kp+p(p+1)/2$ | $Kp+Kp(p+1)/2$ |
| 最常見的坑 | 拿判別係數絕對值當重要性 | 小類別 $n_k\le p$ 無法估計；訓練正確率偏樂觀 |

分類規則都是 $\hat y=\arg\max_k\delta_k(\mathbf x)$（Bayes 規則：後驗機率最大）。LDA／QDA 對仿射轉換不變，不需要先標準化。完整推導與 Mardia、Box's M 公式見 `references/lda-qda-guide.md`。

## Demo

```bash
Rscript scripts/lda_qda_pipeline.R assets/iris.csv --target Species --out my_iris_result
Rscript scripts/lda_qda_pipeline.R assets/iris_missing.csv --target Species --impute median --out my_iris_missing_result
```

這兩個 demo 跑好的完整結果放在 GitHub repo 的 `demo_iris/`、`demo_iris_missing/`（不包含在 skill 的 zip 裡）；skill 內可參考 `references/pipeline_output_iris.txt`、`references/pipeline_output_iris_missing.txt` 的主控台輸出。輸出資料夾已存在時程式會停止，所以重跑要用新的名稱。

| 資料 | 要分類什麼 | 筆數 | X | 示範什麼 |
|---|---|---|---|---|
| `assets/iris.csv` | 鳶尾花品種 `Species`（setosa／versicolor／virginica 各 50） | 150 | 花萼長、寬，花瓣長、寬（cm） | 乾淨資料：確認 Y、X 後一路跑完；R `datasets::iris`（Anderson 1935；Fisher 1936） |
| `assets/iris_missing.csv` | 同上 | 150 | 同上，隨機挖掉 5%（30 格，每列最多 1 格） | 關卡 D：不加 `--impute` 會停止；使用者選補值後，在每份練習資料內補中位數 |

iris 實際結果（`references/pipeline_output_iris.txt`）：

| 項目 | 結果 |
|---|---|
| 資料檢查 | 1 筆完全重複列、13 列組內 IQR 離群標記，全部保留；類別平衡；X 滿秩 |
| Mardia | 三組偏態 p = 0.157–0.194、峰態 p = 0.195–0.879，皆未拒絕多元常態 |
| Box's M | $\chi^2$ = 140.94，df = 20，p = 3.35 × 10⁻²⁰，拒絕共變異數同質 |
| 合併解讀 | Mardia 不顯著、Box's M 顯著：各類散開方式可能不同，QDA 有機會較好，但要看驗證分數 |
| LDA（5-fold × 10 次，seed 2026） | 平均 Macro recall 0.9793；平均答對 146.9／150（範圍 146–147） |
| QDA（相同切分） | 平均 Macro recall 0.9733；平均答對 146.0／150（範圍 145–147） |
| 逐次比較 | LDA 較高 6 次、QDA 較高 0 次、相同 4 次 |
| 判別軸 | LD1 解釋 99.1% 組間變異；setosa 完全分開，誤判都在 versicolor／virginica 之間 |
| 選擇 | LDA（`lda_better`）；最終模型以全部 150 筆重新配適 |

解說時：Box's M 拒絕共變異數同質，但 QDA 沒有因此分得比較準；10 次中 QDA 一次都沒有贏過 LDA。依規則選 LDA，但平均只差 0.9 筆，不能說 LDA 統計上較佳。

iris_missing 實際結果（`references/pipeline_output_iris_missing.txt`）：30 格空格分散在 4 個欄位；LDA 平均答對 141.0／150（139–143，平均 Macro recall 0.9400）、QDA 143.3／150（140–145，0.9553）；QDA 較高 8 次、LDA 1 次、相同 1 次。依規則選 QDA（`qda_better`），報告說明「多數重複由 QDA 較高，但也有 1 次反過來」。50 輪各自使用的中位數記錄在 `imputation.csv`。

## 檔案索引

- `scripts/lda_qda_pipeline.R`：命令列入口；讀 CSV 與欄位說明、選 X、呼叫 `run_pipeline()`，輸出 `report.html`、`report.qmd`、表格、圖、`report.txt`、`summary.json`。
- `scripts/run_analysis.R`：核心函式——`run_analysis()`（LDA／QDA 與重複分層 CV）、`diagnose_data()`（Mardia、Box's M、資料檢查）、`save_figures()`、`make_report()`、`selection_text()`、`diagnostic_notes()`、`write_report_html()`、`write_report_qmd()`、`run_pipeline()`、`save_analysis()`。
- `load_project.R`：在 R／RStudio 中載入上述函式（不會自動執行分析）。
- `references/decision-tree.md`：資料問題與模型可行性的判斷規則。
- `references/lda-qda-guide.md`：方法比較、公式、R 輸出解讀、常見問題、參考資料。
- `references/evaluation.md`：驗證設計、評估指標、模型比較與報告規則。
- `references/pipeline_output_iris.txt`、`references/pipeline_output_iris_missing.txt`：兩個 demo 的主控台輸出。
- `assets/iris.csv`、`assets/iris_missing.csv` 與各自的 `*_dictionary.csv`：demo 資料與欄位說明。
- `tests/validate.R`：程式驗收；從 skill 資料夾執行 `Rscript tests/validate.R`。
- `requirements.txt`：R 與套件版本。
