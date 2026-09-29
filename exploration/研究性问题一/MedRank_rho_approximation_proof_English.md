# Research Problem 1: An Upper Bound on the Approximation Ratio of MedRank under the ρ Distance

Using the mean decomposition of squared distances, a pairing estimate for medians, and the optimality of induced permutations, we prove that MedRank has approximation ratio at most $`3+2\sqrt{2}`$ under the Spearman $`\rho`$ distance. Here the distance is defined as the squared Euclidean distance, $`d_\rho(u,v)=\lVert u-v\rVert_2^2`$. This is an upper bound; the argument does not establish that the bound is attained.

Let $`\Sigma=(\sigma^1,\ldots,\sigma^k)\in(S_n)^k`$ be the input permutations, where $`\sigma^\ell(i)`$ denotes the rank of object $`i`$ in the $`\ell`$th permutation. Write $`x=\frac1k\sum_{\ell=1}^k\sigma^\ell\in\mathbb R^n`$ for the mean vector and $`h\in\mathbb R^n`$ for the vector of coordinatewise medians. The output of MedRank is $`\mu=\operatorname{ind}(h)\in S_n`$: ranks are assigned in increasing order of the components of $`h`$, with ties broken arbitrarily. Consistently with the experiments in this directory, we use the upper median when the number of inputs is even. In fact, the proof below allows any choice within the median interval.

Define $`D_\rho(\pi,\Sigma)=\sum_{\ell=1}^k\lVert\pi-\sigma^\ell\rVert_2^2`$, choose $`\pi^*\in\arg\min_{\pi\in S_n}D_\rho(\pi,\Sigma)`$, and write $`\mathrm{OPT}_\rho=D_\rho(\pi^*,\Sigma)`$. We will prove that

```math
D_\rho(\mu,\Sigma)\le(3+2\sqrt2)\,\mathrm{OPT}_\rho.
```

First, for any $`y\in\mathbb R^n`$, expanding the squares and using $`\sum_{\ell=1}^k(x-\sigma^\ell)=0`$ gives

```math
\begin{aligned}
\sum_{\ell=1}^k\lVert y-\sigma^\ell\rVert_2^2
&=\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2
  +k\lVert y-x\rVert_2^2
  +2\left\langle y-x,\sum_{\ell=1}^k(x-\sigma^\ell)\right\rangle\\
&=\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2
  +k\lVert y-x\rVert_2^2.
\end{aligned}
```

This is the mean decomposition of squared distances. Substituting $`y=\pi^*`$ and $`y=\mu`$, respectively, yields

```math
\begin{aligned}
\mathrm{OPT}_\rho
&=\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2
  +k\lVert\pi^*-x\rVert_2^2,\\
D_\rho(\mu,\Sigma)
&=\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2
  +k\lVert\mu-x\rVert_2^2.
\end{aligned}
```

We now invoke Lemma 2.4 in Section 2.3 of Fagin et al. [1], which establishes the optimality of induced permutations: for any real vector $`z\in\mathbb R^n`$ and any $`p\ge1`$, its induced permutation minimizes $`\sum_{i=1}^n|\pi(i)-z_i|^p`$ over all permutations, with arbitrary tie breaking allowed. We use the case $`p=2`$; taking square roots of the nonnegative squared distances also gives optimality for the Euclidean distance. In particular, the permutation induced by the mean vector is the Borda permutation. The mean decomposition above shows that it minimizes the $`\rho`$ objective, so it can be chosen as $`\pi^*`$.

Applying the optimality of induced permutations to $`z=h`$ gives $`\lVert\mu-h\rVert_2\le\lVert\pi^*-h\rVert_2`$. Applying the triangle inequality twice, we obtain

```math
\begin{aligned}
\lVert\mu-x\rVert_2
&\le\lVert\mu-h\rVert_2+\lVert h-x\rVert_2\\
&\le\lVert\pi^*-h\rVert_2+\lVert h-x\rVert_2\\
&\le\lVert\pi^*-x\rVert_2+2\lVert h-x\rVert_2.
\end{aligned}
```

It remains to control the distance between the median vector and the mean vector. Pairing the input values on opposite sides of each coordinatewise median gives

```math
\sum_{\ell=1}^k\lVert h-\sigma^\ell\rVert_2^2
\le2\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2.
```

The details are provided in Appendix A. This pairing argument also appears in the proof of Lemma 4.3 and equation (4.1) of [1]. Here we take the comparison point to be the real vector $`x`$; the estimate does not require the comparison point to be a permutation.

Taking $`y=h`$ in the mean decomposition, we therefore obtain

```math
\begin{aligned}
\lVert h-x\rVert_2^2
&=\frac1k\left(
\sum_{\ell=1}^k\lVert h-\sigma^\ell\rVert_2^2
-\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2\right)\\
&\le\frac1k\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2.
\end{aligned}
```

To simplify the final estimate, set $`a=\lVert\pi^*-x\rVert_2^2`$ and $`b=\frac1k\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2`$. Then $`a,b\ge0`$ and $`\mathrm{OPT}_\rho=k(a+b)`$. The preceding distance estimates give $`\lVert\mu-x\rVert_2\le\sqrt a+2\sqrt b`$, and hence

```math
D_\rho(\mu,\Sigma)
=kb+k\lVert\mu-x\rVert_2^2
\le k\bigl[b+(\sqrt a+2\sqrt b)^2\bigr]
=k(a+4\sqrt{ab}+5b).
```

If $`\mathrm{OPT}_\rho>0`$, then $`a+b>0`$, so division is valid. Setting $`v=(\sqrt a,\sqrt b)^{\mathsf T}\ne0`$, we recognize the resulting expression as the Rayleigh quotient of a symmetric matrix:

```math
\frac{D_\rho(\mu,\Sigma)}{\mathrm{OPT}_\rho}
\le\frac{a+4\sqrt{ab}+5b}{a+b}
=\frac{v^{\mathsf T}
\begin{pmatrix}1&2\\2&5\end{pmatrix}v}
{v^{\mathsf T}v}
\le\lambda_{\max}\!\begin{pmatrix}1&2\\2&5\end{pmatrix}.
```

The characteristic polynomial of this matrix is $`(1-\lambda)(5-\lambda)-4=\lambda^2-6\lambda+1`$. Its eigenvalues are $`3\pm2\sqrt2`$, so its largest eigenvalue is $`3+2\sqrt2`$. If $`\mathrm{OPT}_\rho=0`$, then $`a=b=0`$, and the preceding bound on the objective directly gives $`D_\rho(\mu,\Sigma)=0`$, without division. Thus, for every input and every choice of tie breaking,

```math
\boxed{D_\rho(\mu,\Sigma)\le(3+2\sqrt2)\,\mathrm{OPT}_\rho,
\qquad r_{\mathrm{MedRank},\rho}\le3+2\sqrt2.}
```

Here $`r_{\mathrm{MedRank},\rho}`$ denotes the worst-case approximation ratio, defined as the supremum over inputs with positive optimal objective value. The key is to decompose the total squared distance into a common dispersion term and a squared distance to the mean, and then control the latter using the median pairing estimate. Retaining the cross term $`4\sqrt{ab}`$ allows the Rayleigh quotient to optimize the contributions of the two nonnegative quantities together.

## Appendix A: The Median Pairing Estimate

Fix a coordinate $`i`$ and arrange its $`k`$ input ranks in nondecreasing order as $`t_1\le\cdots\le t_k`$. Set $`m=h(i)`$ and $`c=x(i)`$, and pair $`t_j`$ with $`t_{k+1-j}`$ for $`1\le j\le\lfloor k/2\rfloor`$. By the definition of a median, the endpoints $`u=t_j`$ and $`v=t_{k+1-j}`$ of every pair satisfy $`u\le m\le v`$. Consequently,

```math
\begin{aligned}
(u-m)^2+(v-m)^2
&\le\bigl((m-u)+(v-m)\bigr)^2\\
&=(v-u)^2\\
&=\bigl((v-c)-(u-c)\bigr)^2\\
&\le2\bigl((u-c)^2+(v-c)^2\bigr).
\end{aligned}
```

If $`k`$ is even, every input rank belongs to a pair. If $`k`$ is odd, the remaining middle value is $`m=t_{(k+1)/2}`$: its squared distance to the median is zero, whereas its squared distance to the mean is nonnegative. Summing over all pairs therefore gives $`\sum_{\ell=1}^k(h(i)-\sigma^\ell(i))^2\le2\sum_{\ell=1}^k(x(i)-\sigma^\ell(i))^2`$ in both cases. Summing over all coordinates $`i`$ yields

```math
\sum_{\ell=1}^k\lVert h-\sigma^\ell\rVert_2^2
\le2\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2.
```

## References

[1] Ronald Fagin, Ravi Kumar, Mohammad Mahdian, D. Sivakumar, Erik Vee. *An Algorithmic View of Voting*. SIAM Journal on Discrete Mathematics, **30**(4): 1978–1996, 2016. DOI: [10.1137/15M1046915](https://doi.org/10.1137/15M1046915). [Full text](https://s3.us.cloud-object-storage.appdomain.cloud/res-files/500-sidma16.pdf). We use Lemma 2.4 in Section 2.3 (statement on p. 1982 and proof on p. 1983), as well as the proof of Lemma 4.3 and equation (4.1) in Section 4.1 (p. 1985).
