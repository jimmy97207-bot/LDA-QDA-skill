# LDA／QDA 分析報告

整體狀態：完成。
模型分析：完成；假設診斷：完成。

## 資料與設定

- 資料說明：R datasets::iris 原樣匯出；150 朵鳶尾花，三個品種各 50 朵（Anderson 1935 採集，Fisher 1936 分析）
- 分類變數：Species
- 預測變數：Sepal.Length、Sepal.Width、Petal.Length、Petal.Width
- 樣本數：150
- 驗證方式：重複分層 5-fold 交叉驗證 10 次（共 50 輪練習／考試）；隨機種子：2026
- 主要指標：Macro recall（各類別召回率的等權平均，即平衡正確率）。
- 先驗機率：由各訓練折的類別比例估計。
- 資料處理：未刪除觀測值、插補缺失值或轉換變數。

| class | n |
| --- | --- |
| setosa |    50 |
| versicolor |    50 |
| virginica |    50 |

### 變數介紹

| column | description | unit | role | type | missing |
| --- | --- | --- | --- | --- | --- |
| Sepal.Length | 花萼長度 | cm | X（預測變數） | numeric |     0 |
| Sepal.Width | 花萼寬度 | cm | X（預測變數） | numeric |     0 |
| Petal.Length | 花瓣長度 | cm | X（預測變數） | numeric |     0 |
| Petal.Width | 花瓣寬度 | cm | X（預測變數） | numeric |     0 |
| Species | 鳶尾花品種（setosa、versicolor、virginica） |  | Y（分類變數） | character |     0 |

## 資料檢查

以下只做標記與說明，不刪除或轉換任何資料（缺失值處理另見下一節）。

| check | value | note |
| --- | --- | --- |
| observations |   150 | 3 classes, 4 predictors |
| duplicate_rows |     1 | Identical rows on Y and X; kept, not assumed to be the same unit. |
| class_balance_ratio |     1 |  |
| classes_with_n_le_p |     0 |  |
| zero_within_class_variance |     0 |  |
| full_rank_X |     1 |  |
| iqr_outlier_rows |    13 | Within-class 1.5 IQR rule per variable; flagged only, all kept. |

## 模型驗證表現

| 模型 | 模型狀態 | 預測涵蓋率 | 平均正確率 | 平均 Macro recall | Macro recall 標準差 | 平均 Macro-F1 | 平均答對筆數 | 最少答對 | 最多答對 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LDA | 驗證完整 | 100.00% | 97.93% | 97.93% | 0.21% | 97.93% | 146.9 |   146 |   147 |
| QDA | 驗證完整 | 100.00% | 97.33% | 97.33% | 0.54% | 97.33% | 146.0 |   145 |   147 |

10 次重複的結果：LDA 每次平均答對 146.9／150 筆（範圍 146–147，平均 Macro recall 97.93%）；QDA 每次平均答對 146.0／150 筆（範圍 145–147，平均 Macro recall 97.33%）。逐次比較 Macro recall：LDA 較高 6 次、QDA 較高 0 次、相同 4 次。LDA 平均較高，選 LDA。10 次中沒有一次由 QDA 較高；平均差 0.9 筆。
選定後，用全部 150 筆資料重新配適 LDA，作為最終模型（final_model.rds）。
選擇規則：兩種方法用相同的切分，做 10 次分層 5-fold 交叉驗證（每次重新隨機分組），比較 10 次的平均 Macro recall，較高者勝出；完全相同時選較簡單的 LDA。逐次勝出次數只用來說明結果穩不穩定。Mardia、Box's M 只作診斷參考，不用來選方法。

表中分數是 10 次重複的平均；每次重複中，每筆觀測恰好被預測一次。必須完成全部重複的模型才參加比較，驗證不完整的模型不提供分數。
這組交叉驗證同時用於比較模型，分數不是獨立測試表現。

## 每一次重複的結果

每一次重複都把資料重新隨機分成 5 份；下表看得出分數隨切分變動多少，以及哪個方法較常勝出。

| 第幾次 | LDA 答對 | LDA Macro recall | QDA 答對 | QDA Macro recall | 較高 |
| --- | --- | --- | --- | --- | --- |
|     1 |   147 |  0.98 |   146 | 0.9733 | LDA |
|     2 |   147 |  0.98 |   147 |  0.98 | 相同 |
|     3 |   147 |  0.98 |   145 | 0.9667 | LDA |
|     4 |   147 |  0.98 |   147 |  0.98 | 相同 |
|     5 |   147 |  0.98 |   147 |  0.98 | 相同 |
|     6 |   147 |  0.98 |   146 | 0.9733 | LDA |
|     7 |   146 | 0.9733 |   146 | 0.9733 | 相同 |
|     8 |   147 |  0.98 |   145 | 0.9667 | LDA |
|     9 |   147 |  0.98 |   145 | 0.9667 | LDA |
|    10 |   147 |  0.98 |   146 | 0.9733 | LDA |

## 各類別表現

把 10 次重複的預測合併計算（每筆觀測被預測 10 次）。

| model | class | recall | precision | f1 |
| --- | --- | --- | --- | --- |
| LDA | setosa | 100.00% | 100.00% | 100.00% |
| LDA | versicolor | 96.00% | 97.76% | 96.87% |
| LDA | virginica | 97.80% | 96.07% | 96.93% |
| QDA | setosa | 100.00% | 100.00% | 100.00% |
| QDA | versicolor | 93.80% | 98.12% | 95.91% |
| QDA | virginica | 98.20% | 94.06% | 96.09% |

混淆矩陣以列百分比呈現：每一列是一個真實類別，數字是被預測成各類別的比例。

### 混淆矩陣：LDA

| 真實類別 \ 預測類別 | setosa | versicolor | virginica |
| --- | --- | --- | --- |
| setosa | 100.0% | 0.0% | 0.0% |
| versicolor | 0.0% | 96.0% | 4.0% |
| virginica | 0.0% | 2.2% | 97.8% |

### 混淆矩陣：QDA

| 真實類別 \ 預測類別 | setosa | versicolor | virginica |
| --- | --- | --- | --- |
| setosa | 100.0% | 0.0% | 0.0% |
| versicolor | 0.0% | 93.8% | 6.2% |
| virginica | 0.0% | 1.8% | 98.2% |

## 假設診斷

| 診斷 | 組別 | 樣本數 | 統計量 | p 值 | 狀態 | 備註 |
| --- | --- | --- | --- | --- | --- | --- |
| Mardia Skewness | setosa |    50 | 25.66 | 0.1772 | 已計算 |  |
| Mardia Kurtosis | setosa |    50 | 1.295 | 0.1953 | 已計算 |  |
| Mardia Skewness | versicolor |    50 | 25.19 | 0.1944 | 已計算 |  |
| Mardia Kurtosis | versicolor |    50 | -0.5719 | 0.5674 | 已計算 |  |
| Mardia Skewness | virginica |    50 | 26.27 | 0.1571 | 已計算 |  |
| Mardia Kurtosis | virginica |    50 | 0.1526 | 0.8787 | 已計算 |  |
| Box M | ALL |   150 | 140.9 | 3.352e-20 | 已計算 | df = 20 |

檢定不顯著不代表假設已被證明；無法估計時，不產生替代 p 值。
Box's M 對非常態很敏感，且樣本大時幾乎必然顯著。
診斷結果用於說明適用性，不直接依單一 p 值決定 LDA／QDA。

### Mardia 與 Box's M 的合併解讀

| Mardia | Box's M | 解讀 | 本資料 |
| --- | --- | --- | --- |
| 不顯著 | 不顯著 | 未發現明顯違反假設；仍依交叉驗證結果選方法。 |  |
| 不顯著 | 顯著 | 各類散開方式可能不同，QDA 有機會較好。 | ← 本資料 |
| 顯著 | 不顯著 | 資料不像鐘形，算出的機率可能不準。 |  |
| 顯著 | 顯著 | Box's M 可能只是受非鐘形影響，不能直接下結論。 |  |

合併解讀：Mardia 不顯著、Box's M 顯著——各類散開方式可能不同，QDA 有機會較好。

## 分組描述統計

| class | variable | n | mean | sd | min | median | max |
| --- | --- | --- | --- | --- | --- | --- | --- |
| setosa | Sepal.Length |    50 | 5.006 | 0.3525 |   4.3 |     5 |   5.8 |
| setosa | Sepal.Width |    50 | 3.428 | 0.3791 |   2.3 |   3.4 |   4.4 |
| setosa | Petal.Length |    50 | 1.462 | 0.1737 |     1 |   1.5 |   1.9 |
| setosa | Petal.Width |    50 | 0.246 | 0.1054 |   0.1 |   0.2 |   0.6 |
| versicolor | Sepal.Length |    50 | 5.936 | 0.5162 |   4.9 |   5.9 |     7 |
| versicolor | Sepal.Width |    50 |  2.77 | 0.3138 |     2 |   2.8 |   3.4 |
| versicolor | Petal.Length |    50 |  4.26 | 0.4699 |     3 |  4.35 |   5.1 |
| versicolor | Petal.Width |    50 | 1.326 | 0.1978 |     1 |   1.3 |   1.8 |
| virginica | Sepal.Length |    50 | 6.588 | 0.6359 |   4.9 |   6.5 |   7.9 |
| virginica | Sepal.Width |    50 | 2.974 | 0.3225 |   2.2 |     3 |   3.8 |
| virginica | Petal.Length |    50 | 5.552 | 0.5519 |   4.5 |  5.55 |   6.9 |
| virginica | Petal.Width |    50 | 2.026 | 0.2747 |   1.4 |     2 |   2.5 |

## 圖檔

- `figures/01_class_counts.png`
- `figures/02_boxplots.png`
- `figures/03_qq_by_class.png`
- `figures/04_pairs.png`
- `figures/05_lda_scores.png`
- `figures/06_confusion_matrix.png`

LDA 判別分數圖以全部資料配適，僅供描述類別分離，不代表驗證表現。

## 模型警告與錯誤

模型執行過程沒有記錄到警告或錯誤；診斷未完成原因另見診斷表。

## 執行範圍與限制

本版處理獨立觀測與數值型預測變數，使用重複分層交叉驗證（固定 5 折、10 次）。
各類別診斷有多個檢定，p 值以診斷用途呈現，未進行多重比較校正。

R：R version 4.3.3 (2024-02-29)；MASS：7.3.60.0.1；診斷：base R (Mardia 1970; Box's M chi-square approximation)
