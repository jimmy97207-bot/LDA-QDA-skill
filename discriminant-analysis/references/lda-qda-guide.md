# LDA／QDA 方法指南

需要確認 LDA／QDA 的假設、公式、R 函式及輸出解讀時閱讀。

## 1. 方法比較

以下以傳統的多變量常態模型與一般共變異數估計為範圍。

| 項目           | LDA 線性區別分析               | QDA 二次區別分析               |
| -------------- | ------------------------------ | ------------------------------ |
| 分析目的       | 根據預測變數判斷已知類別       | 相同                           |
| 組內分布假設   | 各類別內為多變量常態           | 各類別內為多變量常態           |
| 共變異數假設   | 各類別共用一個共變異數矩陣     | 各類別各自估計共變異數矩陣     |
| 類別間決策邊界 | 線性（超平面）                 | 一般為二次曲面                 |
| 參數個數       | $Kp + p(p+1)/2$                | $Kp + K\,p(p+1)/2$             |
| 資料需求       | 合併組內共變異數可逆：$n-K \ge p$ | 每一組都要 $n_k > p$ 且 $\mathbf S_k$ 可逆 |
| R 函式         | `MASS::lda()`                  | `MASS::qda()`                  |

常態假設針對「各類別內的 X 聯合分布」，不是要求所有資料混在一起呈常態；各變數分別接近常態，也不能單獨證明多變量常態。

## 2. 核心公式

### Bayes 分類規則

先驗機率 $\pi_k$、類別條件密度 $f_k(\mathbf x)$。在各種誤判成本相同時，期望誤判成本最小的規則是把 $\mathbf x$ 分到後驗機率最大的類別：

$$
P(Y=k\mid\mathbf x)=\frac{\pi_k f_k(\mathbf x)}{\sum_{l=1}^{K}\pi_l f_l(\mathbf x)},\qquad
\hat y(\mathbf x)=\arg\max_k\ \pi_k f_k(\mathbf x).
$$

### QDA

$\mathbf X\mid Y=k\sim N_p(\boldsymbol\mu_k,\boldsymbol\Sigma_k)$，取對數並去掉與 $k$ 無關的常數：

$$
\delta_k^{Q}(\mathbf x)=-\tfrac12\log|\boldsymbol\Sigma_k|-\tfrac12(\mathbf x-\boldsymbol\mu_k)^\top\boldsymbol\Sigma_k^{-1}(\mathbf x-\boldsymbol\mu_k)+\log\pi_k .
$$

### LDA

再假設 $\boldsymbol\Sigma_k=\boldsymbol\Sigma$。$-\tfrac12\log|\boldsymbol\Sigma|$ 與 $-\tfrac12\mathbf x^\top\boldsymbol\Sigma^{-1}\mathbf x$ 對所有類別相同而消去：

$$
\delta_k^{L}(\mathbf x)=\mathbf x^\top\boldsymbol\Sigma^{-1}\boldsymbol\mu_k-\tfrac12\boldsymbol\mu_k^\top\boldsymbol\Sigma^{-1}\boldsymbol\mu_k+\log\pi_k .
$$

### 估計量（`method = "moment"`）

$$
\hat\pi_k=\frac{n_k}{n},\quad
\hat{\boldsymbol\mu}_k=\bar{\mathbf x}_k,\quad
\mathbf S_k=\frac{1}{n_k-1}\sum_{i:y_i=k}(\mathbf x_i-\bar{\mathbf x}_k)(\mathbf x_i-\bar{\mathbf x}_k)^\top,\quad
\mathbf S_{\text{pooled}}=\frac{\sum_k (n_k-1)\mathbf S_k}{n-K}.
$$

本 skill 的先驗機率在每一個訓練折內重新估計。

### Fisher 判別方向

令 $\mathbf B$ 為組間散佈矩陣、$\mathbf W$ 為組內散佈矩陣：

$$
\max_{\mathbf a}\frac{\mathbf a^\top\mathbf B\mathbf a}{\mathbf a^\top\mathbf W\mathbf a}
\ \Longrightarrow\ \mathbf W^{-1}\mathbf B\,\mathbf a=\lambda\mathbf a .
$$

$\operatorname{rank}(\mathbf B)\le K-1$，所以最多 $\min(p,K-1)$ 個判別軸。兩類時方向為 $\mathbf a\propto\mathbf W^{-1}(\bar{\mathbf x}_1-\bar{\mathbf x}_2)$。`lda()` 輸出的 `scaling` 就是這些軸，`svd^2/sum(svd^2)` 是各軸的相對判別資訊比例（各軸的組間對組內變異比占總和的比例；iris：LD1 99.1%），不是分類正確率，也不是原始資料的變異解釋率。

### 仿射不變性

LDA 與 QDA 的分類結果對 $\mathbf x\mapsto\mathbf A\mathbf x+\mathbf b$（$\mathbf A$ 可逆）不變，所以**不需要先標準化**；標準化只會改變判別係數的尺度，便於比較係數大小。

### 診斷檢定（`scripts/run_analysis.R` 內建，與 MVN 6.3／biotools 4.3 結果一致）

Mardia（1970），各組分別計算，$\mathbf S$ 用除以 $n$ 的最大概似估計，$d_{ij}=(\mathbf x_i-\bar{\mathbf x})^\top\mathbf S^{-1}(\mathbf x_j-\bar{\mathbf x})$：

$$
b_{1,p}=\frac{1}{n^2}\sum_{i,j}d_{ij}^3,\quad \frac{n}{6}b_{1,p}\ \dot\sim\ \chi^2_{p(p+1)(p+2)/6};\qquad
b_{2,p}=\frac1n\sum_i d_{ii}^2,\quad \frac{b_{2,p}-p(p+2)}{\sqrt{8p(p+2)/n}}\ \dot\sim\ N(0,1).
$$

$n<20$ 時偏態統計量乘上小樣本校正 $k=\frac{(p+1)(n+1)(n+3)}{n[(n+1)(p+1)-6]}$。

Box's M（1949），$H_0:\boldsymbol\Sigma_1=\cdots=\boldsymbol\Sigma_K$：

$$
M=(n-K)\ln|\mathbf S_{\text{pooled}}|-\sum_k(n_k-1)\ln|\mathbf S_k|,\quad
c=\frac{2p^2+3p-1}{6(p+1)(K-1)}\Big(\sum_k\frac{1}{n_k-1}-\frac{1}{n-K}\Big),
$$

$$
M(1-c)\ \dot\sim\ \chi^2_{p(p+1)(K-1)/2}.
$$

自由度 $p(p+1)(K-1)/2$ 正好是 QDA 比 LDA 多估的參數個數。Box's M 對非常態非常敏感，樣本量增加時可能偵測到較小的差異，應結合驗證表現解讀；它說明「共變異數不同」，不說明「多估這些參數後預測會比較好」。

### Mardia 與 Box's M 的合併解讀

兩個檢定要一起看，因為 Box's M 的結論依賴常態假設。報告由 `diagnostic_notes()` 產生這張表並標出本資料所在的格子；兩者都不用來選方法。

| Mardia | Box's M | 解讀 |
|---|---|---|
| 不顯著 | 不顯著 | 未發現明顯違反假設；仍依交叉驗證結果選方法 |
| 不顯著 | 顯著 | 各類散開方式可能不同，QDA 有機會較好 |
| 顯著 | 不顯著 | 資料不像鐘形，算出的機率可能不準 |
| 顯著 | 顯著 | Box's M 可能只是受非鐘形影響，不能直接下結論 |

「Mardia 顯著」指任何一類的偏態或峰態 $p<0.05$。此時報告自動加註受影響的類別：交叉驗證仍可評估 LDA／QDA 的表現，但模型前提可能不成立、表現可能因此受影響；後驗機率與 Box's M 要保守解讀；可另外比較不需要常態假設的方法（羅吉斯迴歸、決策樹）。iris 落在第二格（Mardia 皆不顯著、Box's M 顯著），但 QDA 的驗證分數並沒有比 LDA 好。

## 3. R 模型設定

- 將 Y 明確設為分類因子，保存類別名稱與順序。
- 第一版以數值型 X 為主，不直接把類別編碼當成連續量。
- 使用 `method = "moment"` 作為一般估計設定。
- 明確記錄缺失值處理，不讓程式默默刪除觀測值；補值只能在各訓練折內估計（`impute = "median"`），見 `evaluation.md` §8。
- 未另行指定時，先驗機率使用訓練資料的類別比例。
- 自訂先驗機率必須依 Y 的因子層級排列，並說明依據（本版程式尚未開放此參數）。

先驗機率會影響分類；在 MASS 的 LDA 實作中，也可能影響判別軸的方向。

## 4. 主要輸出與解讀

| 輸出                 | 意義                   | 解讀界線                         |
| -------------------- | ---------------------- | -------------------------------- |
| `means`              | 各類別的預測變數平均數 | 描述組間差異，不代表因果效果     |
| `prior`              | 模型採用的類別先驗機率 | 需與類別順序一起記錄             |
| 預測的 `class`       | 模型分派的類別         | 必須與真實類別比較，才能評估表現 |
| 預測的 `posterior`   | 模型估計的類別後驗機率 | 依賴模型與先驗設定，不保證已校準 |
| LDA 的 `scaling`     | 線性判別軸係數         | 不能只依係數絕對值排列變數重要性 |
| LDA 預測的 `x`       | 樣本在判別軸上的分數   | 可用於呈現類別分離情形           |

兩類別時只有 LD1，不強行產生 LD1–LD2 圖（`save_figures()` 會改畫 LD1 密度圖）。
QDA 沒有與 LDA 相同的一組共同線性判別軸；其 `scaling` 是各組的轉換矩陣，不能當成 LDA 判別係數解讀。

## 5. 常見問題與處理

| 問題                       | 處理原則                             |
| -------------------------- | ------------------------------------ |
| 變數在組內沒有變異         | 指出欄位與影響，提出處理建議         |
| 共變異數奇異或變數線性相依 | 記錄問題，不直接宣稱模型有效         |
| QDA 某組無法估計（$n_k\le p$） | 說明受影響類別，保留 LDA 等其他可行分析 |
| 類別樣本數過少             | 檢查每個驗證訓練折，不只檢查完整資料 |
| 圖上類別分得很開           | 仍需驗證，不能用圖形取代分類評估     |
| QDA 訓練正確率較高         | 不直接判定較好，需比較驗證表現       |
| Box's M 顯著               | 不等於應改用 QDA；比較兩者在重複交叉驗證的表現 |
| Mardia 顯著                | 分類仍可用，後驗機率可能不準；Box's M 要保守解讀（報告自動加註） |

QDA 不可行時，可向使用者提及正則化區別分析（Friedman, 1989）：$\hat{\boldsymbol\Sigma}_k(\alpha)=\alpha\mathbf S_k+(1-\alpha)\mathbf S_{\text{pooled}}$，把 LDA（$\alpha=0$）和 QDA（$\alpha=1$）連成一族。本版不實作，需使用者同意擴充範圍。

## 6. 驗證介面的注意事項

`lda(..., CV = TRUE)` 與 `qda(..., CV = TRUE)` 提供留一交叉驗證，不是一般的 k-fold 交叉驗證。內建交叉驗證在估計先驗機率時使用整份輸入資料的類別比例；若流程需要每折重新估計先驗、插補、標準化或選變數，應建立完整的重抽樣流程（本 skill 的 `run_analysis()` 即如此）。

詳細的切分方式、評估指標與模型比較規則，另列於 `evaluation.md`。

## 7. 參考資料

- Fisher, R. A. (1936). The use of multiple measurements in taxonomic problems. *Annals of Eugenics*, 7(2), 179–188.
- Box, G. E. P. (1949). A general distribution theory for a class of likelihood criteria. *Biometrika*, 36(3/4), 317–346.
- Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with applications. *Biometrika*, 57(3), 519–530.
- Friedman, J. H. (1989). Regularized discriminant analysis. *JASA*, 84(405), 165–175.
- Johnson, R. A., & Wichern, D. W. (2007). *Applied Multivariate Statistical Analysis* (6th ed.), Ch. 11. Pearson.
- Hastie, T., Tibshirani, R., & Friedman, J. (2009). *The Elements of Statistical Learning* (2nd ed.), §4.3. Springer.
- Venables, W. N., & Ripley, B. D. (2002). *Modern Applied Statistics with S* (4th ed.). Springer.（MASS 套件）
- 函式文件：[MASS::lda](https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/lda.html)、[MASS::qda](https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/qda.html)、[predict.lda](https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/predict.lda.html)。
