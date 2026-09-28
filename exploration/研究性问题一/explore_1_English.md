# Research Problem I: The Approximation Ratio of the Borda Method under the F Distance

The problem statement gives the known bounds $`3\le r_A\le 4`$ for the Borda method under the Spearman footrule distance (hereafter referred to as the F distance). Below, we use the structure of the mean vector and the optimality of induced permutations to improve the upper bound to $`3`$, matching the known lower bound. Here, $`r_A`$ denotes the worst-case approximation ratio over all problem sizes and inputs.

Let $`\Sigma=(\sigma^1,\ldots,\sigma^k)\in(S_n)^k`$ be the given collection of $`k`$ permutations, where $`S_n`$ is the set of all permutations of $`\{1,\ldots,n\}`$ and $`\sigma^\ell(i)`$ denotes the rank of item $`i`$ in the $`\ell`$ th permutation. The Borda method first computes the average rank of each item, yielding the mean vector

```math
m=h_{\mathrm{Borda}}(\Sigma)
=\frac{1}{k}\sum_{\ell=1}^k\sigma^\ell\in\mathbb{R}^n,
```

It then assigns ranks in increasing order of the components of $`m`$, producing the Borda permutation $`\beta=\mathrm{ind}(m)\in S_n`$. More precisely, $`m_i\lt m_j`$ implies $`\beta(i)\lt\beta(j)`$; ties in average rank may be resolved by any tie-breaking rule. These two objects should be distinguished: $`m`$ is a real vector, whereas $`\beta`$ is the permutation ultimately returned by the algorithm.

The Borda method is equivalent to minimizing the total Spearman $`\rho`$ distance. Here, the $`\rho`$ distance follows the definition in the problem statement, namely the squared Euclidean distance $`d_\rho(u,v)=\lVert u-v\rVert_2^2`$. The mean vector $`m`$ is the unique minimizer of the total squared distance over the space of real vectors, since, for every $`x\in\mathbb{R}^n`$,

```math
\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2
=k\lVert x-m\rVert_2^2
+\sum_{\ell=1}^k\lVert m-\sigma^\ell\rVert_2^2.
```

By Lemma 2.4 in [1], the induced permutation $`\beta`$ minimizes the squared Euclidean distance to $`m`$ among all permutations. It therefore also minimizes the total $`\rho`$ distance over the space of permutations. This conclusion holds for any tie-breaking rule; see also Proposition 7.1 in [2]. Under the F distance, the Borda permutation is generally no longer optimal, but the same mean vector can still be used to bound its approximation ratio.

Let $`D_F(\pi,\Sigma)=\sum_{\ell=1}^k\lVert\pi-\sigma^\ell\rVert_1`$, choose an optimal permutation $`\pi^{*}\in\arg\min_{\pi\in S_n}D_F(\pi,\Sigma)`$, and denote its objective value by $`\mathrm{OPT}_F=D_F(\pi^{*},\Sigma)`$. By the triangle inequality, the objective value of the Borda permutation satisfies

```math
D_F(\beta,\Sigma)
\le k\lVert\beta-m\rVert_1
+\sum_{\ell=1}^k\lVert m-\sigma^\ell\rVert_1.
```

We bound the two terms on the right-hand side separately. First, by the case $`p=1`$ of Lemma 2.4 in [1], $`\beta`$ also minimizes the $`L_1`$ distance to $`m`$ among all permutations. Using the expression for $`m`$ as the mean of the input permutations, we obtain

```math
\begin{aligned}
k\lVert\beta-m\rVert_1
&\le k\lVert\pi^{*}-m\rVert_1\\
&=\left\lVert\sum_{\ell=1}^k(\pi^{*}-\sigma^\ell)\right\rVert_1\\
&\le\sum_{\ell=1}^k\lVert\pi^{*}-\sigma^\ell\rVert_1
=\mathrm{OPT}_F.
\end{aligned}
```

Next, applying the triangle inequality to each input permutation and using the bound on $`k\lVert\pi^{*}-m\rVert_1`$ established above gives

```math
\begin{aligned}
\sum_{\ell=1}^k\lVert m-\sigma^\ell\rVert_1
&\le\sum_{\ell=1}^k
\bigl(\lVert m-\pi^{*}\rVert_1
+\lVert\pi^{*}-\sigma^\ell\rVert_1\bigr)\\
&=k\lVert m-\pi^{*}\rVert_1+\mathrm{OPT}_F\\
&\le 2\mathrm{OPT}_F.
\end{aligned}
```

This estimate also follows directly from Lemma 4.5 in [1]. Combining the two bounds yields $`D_F(\beta,\Sigma)\le 3\mathrm{OPT}_F`$. This inequality holds for all $`n,k`$, every input $`\Sigma\in(S_n)^k`$, and any tie-breaking rule. Hence, $`r_A\le 3`$. When $`\mathrm{OPT}_F=0`$, the same inequality ensures that the objective value of the Borda permutation is also $`0`$; in defining the approximation ratio, we need only take the supremum over inputs with $`\mathrm{OPT}_F\gt 0`$.

Combining this upper bound with the known lower bound $`r_A\ge 3`$ given in the problem statement, we conclude that

```math
\boxed{r_A=3.}
```

Thus, the upper bound for the Borda method under the F distance can be improved from $`4`$ to $`3`$, determining its worst-case approximation ratio. The key to this improvement is to use directly the fact that $`m`$ is the mean of the input permutations, thereby bounding the error between the induced permutation and the mean vector by $`\mathrm{OPT}_F`$.

---

References:

[1] Ronald Fagin, Ravi Kumar, Mohammad Mahdian, D. Sivakumar, Erik Vee. *An Algorithmic View of Voting*. SIAM Journal on Discrete Mathematics, 30(4): 1978–1996, 2016. DOI: 10.1137/15M1046915. This proof uses Lemmas 2.4 and 4.5.

[2] Dezső Bednay, Balázs Fleiner, Attila Tasnádi. *The MedRank algorithm and Spearman’s footrule versus the Borda count and Spearman’s rank correlation*. Optimization, 2026. DOI: 10.1080/02331934.2026.2649825. For the relationship between the Borda method and minimization of the total $`\rho`$ distance, see Proposition 7.1.
