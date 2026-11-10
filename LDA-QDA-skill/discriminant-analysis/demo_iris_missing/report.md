# LDA／QDA 分析報告

整體狀態：完成。
模型分析：完成；假設診斷：完成。

## 資料與設定

- 資料說明：iris_missing：R datasets::iris 隨機挖掉 5% 的 X（30 格，每列最多 1 格，完全隨機），用來示範缺失值關卡與各訓練折內的中位數補值
- 分類變數：Species
- 預測變數：Sepal.Length、Sepal.Width、Petal.Length、Petal.Width
- 樣本數：150
- 驗證方式：重複分層 5-fold 交叉驗證 10 次（共 50 輪練習／考試）；隨機種子：2026
- 主要指標：Macro recall（各類別召回率的等權平均，即平衡正確率）。
- 先驗機率：由各訓練折的類別比例估計。
- 資料處理：未刪除觀測值或轉換變數；缺失值依使用者同意在各練習資料內補中位數。

| class | n |
| --- | --- |
| setosa |    50 |
| versicolor |    50 |
| virginica |    50 |

### 變數介紹

| column | description | unit | role | type | missing |
| --- | --- | --- | --- | --- | --- |
| Sepal.Length | 花萼長度 | cm | X（預測變數） | numeric |    13 |
| Sepal.Width | 花萼寬度 | cm | X（預測變數） | numeric |     5 |
| Petal.Length | 花瓣長度 | cm | X（預測變數） | numeric |     5 |
| Petal.Width | 花瓣寬度 | cm | X（預測變數） | numeric |     7 |
| Species | 鳶尾花品種（setosa、versicolor、virginica） |  | Y（分類變數） | character |     0 |

## 資料檢查

以下只做標記與說明，不刪除或轉換任何資料（缺失值處理另見下一節）。

| check | value | note |
| --- | --- | --- |
| missing_cells |    30 | Sepal.Length=13, Sepal.Width=5, Petal.Length=5, Petal.Width=7 |
| rows_with_missing |    30 | Diagnostics, summaries and figures use the 120 complete rows; models use fold-wise median imputation. |
| complete_observations |   120 | 3 classes, 4 predictors |
| duplicate_rows |     1 | Identical rows on Y and X; kept, not assumed to be the same unit. |
| class_balance_ratio | 1.051 |  |
| classes_with_n_le_p |     0 |  |
| zero_within_class_variance |     0 |  |
| full_rank_X |     1 |  |
| iqr_outlier_rows |    13 | Within-class 1.5 IQR rule per variable; flagged only, all kept. |

## 缺失值處理

X 共有 30 個空格。依使用者決定，在每一份練習資料內計算各欄中位數，再用這個中位數補該輪的練習與考試資料；考試資料不參與計算，避免資料洩漏。最終模型以全部資料的中位數補值。
假設診斷、描述統計與探索圖（01–05）只使用沒有空格的完整觀測；混淆矩陣（06）來自交叉驗證，包含全部觀測。
下表摘要 50 輪練習資料的中位數；每一輪的完整紀錄見 `imputation.csv`。

| 變數 | 中位數最小 | 中位數最大 | 每輪補入練習格數 | 每輪補入考試格數 |
| --- | --- | --- | --- | --- |
| Sepal.Length |   5.8 |  5.95 |  10.4 |   2.6 |
| Sepal.Width |     3 |     3 |     4 |     1 |
| Petal.Length |   4.2 |   4.5 |     4 |     1 |
| Petal.Width |   1.3 |   1.4 |   5.6 |   1.4 |

## 模型驗證表現

| 模型 | 模型狀態 | 預測涵蓋率 | 平均正確率 | 平均 Macro recall | Macro recall 標準差 | 平均 Macro-F1 | 平均答對筆數 | 最少答對 | 最多答對 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LDA | 驗證完整 | 100.00% | 94.00% | 94.00% | 0.77% | 94.04% | 141.0 |   139 |   143 |
| QDA | 驗證完整 | 100.00% | 95.53% | 95.53% | 1.09% | 95.55% | 143.3 |   140 |   145 |

10 次重複的結果：LDA 每次平均答對 141.0／150 筆（範圍 139–143，平均 Macro recall 94.00%）；QDA 每次平均答對 143.3／150 筆（範圍 140–145，平均 Macro recall 95.53%）。逐次比較 Macro recall：LDA 較高 1 次、QDA 較高 8 次、相同 1 次。QDA 平均較高，選 QDA。多數重複由 QDA 較高，但也有 1 次反過來；平均差 2.3 筆。
選定後，用全部 150 筆資料重新配適 QDA，作為最終模型（final_model.rds）。
選擇規則：兩種方法用相同的切分，做 10 次分層 5-fold 交叉驗證（每次重新隨機分組），比較 10 次的平均 Macro recall，較高者勝出；完全相同時選較簡單的 LDA。逐次勝出次數只用來說明結果穩不穩定。Mardia、Box's M 只作診斷參考，不用來選方法。

表中分數是 10 次重複的平均；每次重複中，每筆觀測恰好被預測一次。必須完成全部重複的模型才參加比較，驗證不完整的模型不提供分數。
這組交叉驗證同時用於比較模型，分數不是獨立測試表現。

## 每一次重複的結果

每一次重複都把資料重新隨機分成 5 份；下表看得出分數隨切分變動多少，以及哪個方法較常勝出。

| 第幾次 | LDA 答對 | LDA Macro recall | QDA 答對 | QDA Macro recall | 較高 |
| --- | --- | --- | --- | --- | --- |
|     1 |   143 | 0.9533 |   144 |  0.96 | QDA |
|     2 |   142 | 0.9467 |   144 |  0.96 | QDA |
|     3 |   141 |  0.94 |   141 |  0.94 | 相同 |
|     4 |   139 | 0.9267 |   143 | 0.9533 | QDA |
|     5 |   141 |  0.94 |   145 | 0.9667 | QDA |
|     6 |   141 |  0.94 |   144 |  0.96 | QDA |
|     7 |   140 | 0.9333 |   143 | 0.9533 | QDA |
|     8 |   141 |  0.94 |   140 | 0.9333 | LDA |
|     9 |   142 | 0.9467 |   144 |  0.96 | QDA |
|    10 |   140 | 0.9333 |   145 | 0.9667 | QDA |

## 各類別表現

把 10 次重複的預測合併計算（每筆觀測被預測 10 次）。

| model | class | recall | precision | f1 |
| --- | --- | --- | --- | --- |
| LDA | setosa | 96.80% | 100.00% | 98.37% |
| LDA | versicolor | 94.20% | 88.53% | 91.28% |
| LDA | virginica | 91.00% | 94.01% | 92.48% |
| QDA | setosa | 98.80% | 100.00% | 99.40% |
| QDA | versicolor | 91.20% | 96.41% | 93.73% |
| QDA | virginica | 96.60% | 90.62% | 93.51% |

混淆矩陣以列百分比呈現：每一列是一個真實類別，數字是被預測成各類別的比例。

### 混淆矩陣：LDA

| 真實類別 \ 預測類別 | setosa | versicolor | virginica |
| --- | --- | --- | --- |
| setosa | 96.8% | 3.2% | 0.0% |
| versicolor | 0.0% | 94.2% | 5.8% |
| virginica | 0.0% | 9.0% | 91.0% |

### 混淆矩陣：QDA

| 真實類別 \ 預測類別 | setosa | versicolor | virginica |
| --- | --- | --- | --- |
| setosa | 98.8% | 0.0% | 1.2% |
| versicolor | 0.0% | 91.2% | 8.8% |
| virginica | 0.0% | 3.4% | 96.6% |

## 假設診斷

| 診斷 | 組別 | 樣本數 | 統計量 | p 值 | 狀態 | 備註 |
| --- | --- | --- | --- | --- | --- | --- |
| Mardia Skewness | setosa |    39 | 16.45 | 0.6881 | 已計算 |  |
| Mardia Kurtosis | setosa |    39 | 0.2634 | 0.7923 | 已計算 |  |
| Mardia Skewness | versicolor |    40 | 18.59 | 0.5483 | 已計算 |  |
| Mardia Kurtosis | versicolor |    40 | -1.034 | 0.3013 | 已計算 |  |
| Mardia Skewness | virginica |    41 | 27.41 | 0.1241 | 已計算 |  |
| Mardia Kurtosis | virginica |    41 | 0.124 | 0.9013 | 已計算 |  |
| Box M | ALL |   120 | 111.8 | 9.254e-15 | 已計算 | df = 20 |

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
| setosa | Sepal.Length |    39 | 5.038 | 0.3689 |   4.3 |   5.1 |   5.8 |
| setosa | Sepal.Width |    39 | 3.451 | 0.3501 |   2.9 |   3.4 |   4.4 |
| setosa | Petal.Length |    39 | 1.451 | 0.1745 |     1 |   1.5 |   1.9 |
| setosa | Petal.Width |    39 | 0.2436 | 0.09678 |   0.1 |   0.2 |   0.5 |
| versicolor | Sepal.Length |    40 | 6.015 | 0.4969 |     5 |     6 |     7 |
| versicolor | Sepal.Width |    40 | 2.817 | 0.3046 |   2.2 |   2.9 |   3.4 |
| versicolor | Petal.Length |    40 | 4.325 | 0.4331 |   3.3 |   4.4 |   5.1 |
| versicolor | Petal.Width |    40 | 1.353 |   0.2 |     1 |   1.4 |   1.8 |
| virginica | Sepal.Length |    41 | 6.522 | 0.6605 |   4.9 |   6.4 |   7.9 |
| virginica | Sepal.Width |    41 | 2.944 | 0.3294 |   2.2 |     3 |   3.8 |
| virginica | Petal.Length |    41 | 5.493 | 0.5698 |   4.5 |   5.4 |   6.9 |
| virginica | Petal.Width |    41 | 2.005 | 0.2738 |   1.4 |     2 |   2.5 |

## 圖檔

- `figures/01_class_counts.png`
- `figures/02_boxplots.png`
- `figures/03_qq_by_class.png`
- `figures/04_pairs.png`
- `figures/05_lda_scores.png`
- `figures/06_confusion_matrix.png`

LDA 判別分數圖以完整觀測配適，僅供描述類別分離，不代表驗證表現。

## 模型警告與錯誤

模型執行過程沒有記錄到警告或錯誤；診斷未完成原因另見診斷表。

## 執行範圍與限制

本版處理獨立觀測與數值型預測變數，使用重複分層交叉驗證（固定 5 折、10 次）。
各類別診斷有多個檢定，p 值以診斷用途呈現，未進行多重比較校正。

R：R version 4.3.3 (2024-02-29)；MASS：7.3.60.0.1；診斷：base R (Mardia 1970; Box's M chi-square approximation)
