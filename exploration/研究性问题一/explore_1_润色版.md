# 研究性问题一：Borda 法在 F 距离下的近似比

题目给出了 Borda 法在 Spearman footrule 距离（以下简称 F 距离）下的已知估计 $3\le r_A\le 4$。下面利用均值向量的结构与诱导排列的最优性，将上界改进为 $3$，从而与已知下界一致。这里的 $r_A$ 指对所有问题规模与输入取最坏情形的近似比。

设 $\Sigma=(\sigma^1,\ldots,\sigma^k)\in(S_n)^k$ 为给定的 $k$ 个排列，其中 $S_n$ 表示 $\{1,\ldots,n\}$ 的所有排列构成的集合，$\sigma^\ell(i)$ 表示对象 $i$ 在第 $\ell$ 个排列中的名次。Borda 法先计算各对象的平均名次，得到均值向量

$$
m=h_{\mathrm{Borda}}(\Sigma)
=\frac{1}{k}\sum_{\ell=1}^k\sigma^\ell\in\mathbb{R}^n,
$$

再按 $m$ 的各分量从小到大赋予名次，得到 Borda 排列 $\beta=\operatorname{ind}(m)\in S_n$。具体而言，$m_i<m_j$ 蕴含 $\beta(i)<\beta(j)$；若平均名次相同，则可采用任意平局打破规则。应当区分这两个对象：$m$ 是实向量，而 $\beta$ 才是算法最终输出的排列。

Borda 法与 Spearman $\rho$ 距离的最小化具有一致性。这里的 $\rho$ 距离采用题目中的定义，即平方欧氏距离 $d_\rho(u,v)=\lVert u-v\rVert_2^2$。均值向量 $m$ 是实向量空间中总平方距离的唯一最优解，因为对任意 $x\in\mathbb{R}^n$，有

$$
\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2
=k\lVert x-m\rVert_2^2
+\sum_{\ell=1}^k\lVert m-\sigma^\ell\rVert_2^2.
$$

由文献 [1] 的引理 2.4，诱导排列 $\beta$ 在所有排列中到 $m$ 的平方欧氏距离最小，因此它也最小化排列空间中的总 $\rho$ 距离。这一结论对任意平局打破规则均成立，亦可参见文献 [2] 的命题 7.1。对于 F 距离，Borda 排列一般不再是最优解，但仍可利用同一均值向量控制其近似比。

记 $D_F(\pi,\Sigma)=\sum_{\ell=1}^k\lVert\pi-\sigma^\ell\rVert_1$，并取最优排列 $\pi^*\in\arg\min_{\pi\in S_n}D_F(\pi,\Sigma)$，其目标值记为 $\mathrm{OPT}_F=D_F(\pi^*,\Sigma)$。由三角不等式，Borda 排列的目标值满足

$$
D_F(\beta,\Sigma)
\le k\lVert\beta-m\rVert_1
+\sum_{\ell=1}^k\lVert m-\sigma^\ell\rVert_1.
$$

下面分别估计右侧两项。首先，由文献 [1] 的引理 2.4 在 $p=1$ 时的结论，$\beta$ 也是所有排列中到 $m$ 的 $L_1$ 距离最小的排列。再利用 $m$ 的均值表示，可得

$$
\begin{aligned}
k\lVert\beta-m\rVert_1
&\le k\lVert\pi^*-m\rVert_1\\
&=\left\lVert\sum_{\ell=1}^k(\pi^*-\sigma^\ell)\right\rVert_1\\
&\le\sum_{\ell=1}^k\lVert\pi^*-\sigma^\ell\rVert_1
=\mathrm{OPT}_F.
\end{aligned}
$$

其次，对每个输入排列应用三角不等式，并利用上式对 $k\lVert\pi^*-m\rVert_1$ 的估计，得到

$$
\begin{aligned}
\sum_{\ell=1}^k\lVert m-\sigma^\ell\rVert_1
&\le\sum_{\ell=1}^k
\bigl(\lVert m-\pi^*\rVert_1
+\lVert\pi^*-\sigma^\ell\rVert_1\bigr)\\
&=k\lVert m-\pi^*\rVert_1+\mathrm{OPT}_F\\
&\le 2\mathrm{OPT}_F.
\end{aligned}
$$

这一估计亦可由文献 [1] 的引理 4.5 直接得到。合并上述两项估计，有 $D_F(\beta,\Sigma)\le 3\mathrm{OPT}_F$。该不等式对任意 $n,k$、任意输入 $\Sigma\in(S_n)^k$ 及任意平局打破规则均成立，因此 $r_A\le 3$。当 $\mathrm{OPT}_F=0$ 时，上式同时保证 Borda 排列的目标值为 $0$；计算近似比时，只需对 $\mathrm{OPT}_F>0$ 的输入取上确界。

结合题目给出的已知下界 $r_A\ge 3$，最终得到

$$
\boxed{r_A=3.}
$$

因此，Borda 法在 F 距离下的上界可以由 $4$ 改进为 $3$，并由此确定其最坏情形近似比。上述改进的关键在于：直接利用 $m$ 是输入排列的均值这一事实，将诱导排列与均值向量之间的误差控制在 $\mathrm{OPT}_F$ 以内。

---

参考文献：

[1] Ronald Fagin, Ravi Kumar, Mohammad Mahdian, D. Sivakumar, Erik Vee. *An Algorithmic View of Voting*. SIAM Journal on Discrete Mathematics, 30(4): 1978–1996, 2016. DOI: 10.1137/15M1046915. 本文使用引理 2.4 与引理 4.5。

[2] Dezső Bednay, Balázs Fleiner, Attila Tasnádi. *The MedRank algorithm and Spearman’s footrule versus the Borda count and Spearman’s rank correlation*. Optimization, 2026. DOI: 10.1080/02331934.2026.2649825. 关于 Borda 法与总 $\rho$ 距离最小化的关系，参见命题 7.1。
