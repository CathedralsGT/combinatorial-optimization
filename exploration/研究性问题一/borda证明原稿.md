# 研究型问题一

首先注意到题目中给出的 Borda Count 与最小化 Spearman rho distance 之间的等价性（证明详见 The MedRank algorithm and Spearman's footrule versus the Borda count and Spearman's rank correlation，Bednay et al.），知 $h_{borda}$ 相当于直接对各分量求均方最小值得到一个包含于 $\mathbb{R}^n$ 的向量；（然后再采用任意的平局打破规则即可得到 Borda Count 对应的重排，由先前的结论可知这一重排满足了最小 Spearman rho distance）

我们记 
$$h_{borda}(\Sigma) = m = \frac{1}{k} \sum_{l=1}^k \sigma^l \in \mathbb{R}^n$$
并令 $\beta = ind(m)$ 为由 $m$ 诱导出的排列；
同时记 $\pi^*$ 为最小化 Spearman Footrule 的最优解，并记该最小值为 
$$OPT_F = \sum_{l=1}^k || \pi^* - \sigma^l ||_1$$
则有
$$ D_F(\beta, \Sigma) = \sum_{l=1}^k ||\beta - \sigma^l||_1$$
$$ \leq k||\beta-m||_1 + \sum_{l=1}^k||m-\sigma^l||_1$$
由 An Algorithm View of Voting 中的 lemma 2.4（诱导排列的最优性） 得
$$k||\beta-m||_1 \leq k||\pi^*-m||_1$$
$$ = ||k\pi^*-km||_1$$
$$ = ||k\pi^* - \sum_{l=1}^k \sigma^l||_1$$
$$ \leq \sum_{l=1}^k ||\pi^* - \sigma^l||_1 = OPT_F$$

再由 An Algorithm View of Voting 中的 lemma 4.5 得
$$\sum_{k=1}^l ||m-\sigma^l|| \leq 2OPT_F$$
得
$$ D_F(\beta, \Sigma) \leq 3OPT_F , \,\, \forall \Sigma \in \mathbb{S}_n$$
因此当 A 代表 Borda Count 且采用 F 距离时，$r_A = 3$