# 研究性问题一：MedRank 法在 ρ 距离下的近似比上界

下面利用平方距离的均值分解、中位数的配对估计以及诱导排列的最优性，证明 MedRank 法在 Spearman $`\rho`$ 距离下满足近似比上界 $`3+2\sqrt{2}`$。这里采用平方欧氏距离的定义，即 $`d_\rho(u,v)=\lVert u-v\rVert_2^2`$；所得结论是一个上界，并不意味着该上界一定能够达到。

设 $`\Sigma=(\sigma^1,\ldots,\sigma^k)\in(S_n)^k`$ 为给定的输入排列，$`\sigma^\ell(i)`$ 表示对象 $`i`$ 在第 $`\ell`$ 个排列中的名次。记均值向量为 $`x=\frac1k\sum_{\ell=1}^k\sigma^\ell\in\mathbb R^n`$，逐坐标中位数向量为 $`h\in\mathbb R^n`$，MedRank 的输出为 $`\mu=\mathrm{ind}(h)\in S_n`$，即按 $`h`$ 的分量从小到大赋予名次，并在分量相同时任意打破平局。与本目录中的实验约定一致，偶数个输入时取上中位数；下面的证明事实上适用于中位数区间内的任意选择。

定义 $`D_\rho(\pi,\Sigma)=\sum_{\ell=1}^k\lVert\pi-\sigma^\ell\rVert_2^2`$，并取 $`\pi^*\in\arg\min_{\pi\in S_n}D_\rho(\pi,\Sigma)`$，记 $`\mathrm{OPT}_\rho=D_\rho(\pi^*,\Sigma)`$。我们要证明

```math
D_\rho(\mu,\Sigma)\le(3+2\sqrt2)\,\mathrm{OPT}_\rho.
```

首先，对任意 $`y\in\mathbb R^n`$，由 $`\sum_{\ell=1}^k(x-\sigma^\ell)=0`$，展开平方可得

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

这就是平方距离的均值分解。分别取 $`y=\pi^*`$ 和 $`y=\mu`$，得到

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

接下来引用 Fagin 等人在文献 [1] 第 2.3 节中证明的引理 2.4（induced permutations are optimal，诱导排列的最优性）：对任意实向量 $`z\in\mathbb R^n`$ 及任意 $`p\ge1`$，其诱导排列在所有排列中最小化 $`\sum_{i=1}^n|\pi(i)-z_i|^p`$，且允许任意平局打破规则。这里使用 $`p=2`$ 的情形；对非负平方距离取平方根后，同样得到欧氏距离的最优性。特别地，均值向量的诱导排列即 Borda 排列，由上述均值分解可知它是 $`\rho`$ 目标的最优解，故可以将其取作 $`\pi^*`$。

将诱导排列的最优性用于 $`z=h`$，有 $`\lVert\mu-h\rVert_2\le\lVert\pi^*-h\rVert_2`$。再两次使用三角不等式，得到

```math
\begin{aligned}
\lVert\mu-x\rVert_2
&\le\lVert\mu-h\rVert_2+\lVert h-x\rVert_2\\
&\le\lVert\pi^*-h\rVert_2+\lVert h-x\rVert_2\\
&\le\lVert\pi^*-x\rVert_2+2\lVert h-x\rVert_2.
\end{aligned}
```

因此，剩下的关键是控制中位数向量与均值向量之间的距离。利用中位数两侧的输入值逐对配对，可以证明

```math
\sum_{\ell=1}^k\lVert h-\sigma^\ell\rVert_2^2
\le2\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2.
```

详细推导见附录 A。这一配对方法也见于文献 [1] 引理 4.3 的证明及其中的式 (4.1)；这里将比较点取为实向量 $`x`$，所用估计并不要求比较点是排列。

在均值分解中取 $`y=h`$，于是

```math
\begin{aligned}
\lVert h-x\rVert_2^2
&=\frac1k\left(
\sum_{\ell=1}^k\lVert h-\sigma^\ell\rVert_2^2
-\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2\right)\\
&\le\frac1k\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2.
\end{aligned}
```

为了简化最后的估计，记 $`a=\lVert\pi^*-x\rVert_2^2`$、$`b=\frac1k\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2`$，则 $`a,b\ge0`$，并且 $`\mathrm{OPT}_\rho=k(a+b)`$。由前面的距离估计，$`\lVert\mu-x\rVert_2\le\sqrt a+2\sqrt b`$，因此

```math
D_\rho(\mu,\Sigma)
=kb+k\lVert\mu-x\rVert_2^2
\le k\bigl[b+(\sqrt a+2\sqrt b)^2\bigr]
=k(a+4\sqrt{ab}+5b).
```

若 $`\mathrm{OPT}_\rho>0`$，则 $`a+b>0`$，故可以相除。令 $`v=(\sqrt a,\sqrt b)^{\mathsf T}\ne0`$，右侧恰好是一个对称矩阵的 Rayleigh 商：

```math
\frac{D_\rho(\mu,\Sigma)}{\mathrm{OPT}_\rho}
\le\frac{a+4\sqrt{ab}+5b}{a+b}
=\frac{v^{\mathsf T}
\begin{pmatrix}1&2\\2&5\end{pmatrix}v}
{v^{\mathsf T}v}
\le\lambda_{\max}\!\begin{pmatrix}1&2\\2&5\end{pmatrix}.
```

矩阵的特征多项式为 $`(1-\lambda)(5-\lambda)-4=\lambda^2-6\lambda+1`$，两个特征值为 $`3\pm2\sqrt2`$，所以最大特征值为 $`3+2\sqrt2`$。若 $`\mathrm{OPT}_\rho=0`$，则 $`a=b=0`$，上面的目标值估计直接给出 $`D_\rho(\mu,\Sigma)=0`$，无需进行除法。因此，对所有输入以及任意平局打破规则均有

```math
\boxed{D_\rho(\mu,\Sigma)\le(3+2\sqrt2)\,\mathrm{OPT}_\rho,
\qquad r_{\mathrm{MedRank},\rho}\le3+2\sqrt2.}
```

这里 $`r_{\mathrm{MedRank},\rho}`$ 是在最优值为正的输入上取上确界得到的最坏情形近似比。证明的核心是将总平方距离分成公共的离散程度项与到均值的距离项，再通过中位数配对估计控制后者；最后保留交叉项 $`4\sqrt{ab}`$，利用 Rayleigh 商统一优化两个非负量的贡献。

## 附录 A：中位数配对估计

固定一个坐标 $`i`$，将该坐标的 $`k`$ 个输入名次按非降顺序记为 $`t_1\le\cdots\le t_k`$。令 $`m=h(i)`$、$`c=x(i)`$，并将 $`t_j`$ 与 $`t_{k+1-j}`$ 配成一对，其中 $`1\le j\le\lfloor k/2\rfloor`$。由中位数的定义，每一对的两个端点 $`u=t_j`$、$`v=t_{k+1-j}`$ 均满足 $`u\le m\le v`$，从而

```math
\begin{aligned}
(u-m)^2+(v-m)^2
&\le\bigl((m-u)+(v-m)\bigr)^2\\
&=(v-u)^2\\
&=\bigl((v-c)-(u-c)\bigr)^2\\
&\le2\bigl((u-c)^2+(v-c)^2\bigr).
\end{aligned}
```

当 $`k`$ 为偶数时，所有输入名次恰好配完；当 $`k`$ 为奇数时，剩下的中间项就是 $`m=t_{(k+1)/2}`$，其到中位数的平方距离为零，而到均值的平方距离非负。因此，对各对求和后，两种情形都给出 $`\sum_{\ell=1}^k(h(i)-\sigma^\ell(i))^2\le2\sum_{\ell=1}^k(x(i)-\sigma^\ell(i))^2`$。再对坐标 $`i`$ 求和，得到

```math
\sum_{\ell=1}^k\lVert h-\sigma^\ell\rVert_2^2
\le2\sum_{\ell=1}^k\lVert x-\sigma^\ell\rVert_2^2.
```

## 参考文献

[1] Ronald Fagin, Ravi Kumar, Mohammad Mahdian, D. Sivakumar, Erik Vee. *An Algorithmic View of Voting*. SIAM Journal on Discrete Mathematics, **30**(4): 1978–1996, 2016. DOI: [10.1137/15M1046915](https://doi.org/10.1137/15M1046915). [论文全文](https://s3.us.cloud-object-storage.appdomain.cloud/res-files/500-sidma16.pdf)。本文引用第 2.3 节引理 2.4（陈述位于第 1982 页，证明位于第 1983 页），以及第 4.1 节引理 4.3 的证明与式 (4.1)（第 1985 页）。
