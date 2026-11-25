# 第 5 章 紧算子和 Fredholm 算子：零基础期末急救与作业详解

> **适用对象**：完全没学过泛函分析、需要快速打通核心逻辑并稳拿期末与作业分的本科生  
> **对照教材**：《泛函分析》（科学出版社“十四五”普通高等教育本科规划教材，第 5 章，教材第 182–205 页）  
> **覆盖作业**：习题 5 第 **5、6、14、15、20、22** 题（共 6 题，附教材笔误勘误与严密排雷）

---

## 零、零基础前置急救包（第 1–4 章核心概念 5 分钟速通）

**什么是泛函分析？** 一句话总结：**把线性代数推广到无穷维，再加上数学分析里的“极限与收敛”。**  
在线性代数里，向量是有限维数组 $(x_1, \dots, x_n)$，线性变换是 $n \times n$ 矩阵；在泛函分析里，“向量”可以是一个无穷数列或一个连续函数，“矩阵”就升级成了作用在函数/数列空间上的**线性算子**。

下面按“第 5 章真的会用到什么”提炼 4 组必背前置知识，每条都标明了它在第 5 章的直接用处。

### 0.1 度量空间、完备性、列紧性与 Arzelà-Ascoli 定理

1. **度量空间 $(X, d)$**：就是定义了“两点之间距离 $d(x, y)$”的集合，只要距离满足非负性、对称性和三角不等式，就能讨论序列收敛 $x_n \to x$（即 $d(x_n, x) \to 0$）。
2. **完备性（Cauchy 列必收敛）**：如果空间里任何“越挤越近”的 Cauchy 序列（$\lim_{m,n\to\infty} d(x_m, x_n) = 0$）都能在空间内部找到极限点（空间没有“针孔漏洞”），就称这个空间是**完备的**。
   - 反例：有理数 $\mathbb{Q}$ 不完备（$1, 1.4, 1.41, \dots$ 趋向于 $\sqrt{2} \notin \mathbb{Q}$）；实数 $\mathbb{R}$ 完备。
3. **列紧集（相对紧集）**：子集 $E \subset X$ 称为**列紧的**，如果 $E$ 中**任意一个序列** $\{x_n\} \subset E$ 都能抽出一个在 $X$ 中收敛的子列 $\{x_{n_k}\}$。（注意：子列的极限允许落在 $E$ 的闭包 $\overline{E}$ 中，不要求一定在 $E$ 内。）
   - **高数类比**：在有限维空间 $\mathbb{R}^n$ 里，Bolzano-Weierstrass 定理告诉我们“有界集 $\iff$ 列紧集”。**但在无穷维空间里，有界集通常不再列紧！** 这正是引入“紧算子”的全部动机。
4. **Arzelà-Ascoli 定理（判定连续函数集列紧的头号神器 → 习题 5 第 14 题直接要用！）**：
   设 $C([a, b])$ 是闭区间 $[a, b]$ 上全体连续函数组成的空间，配最大值范数 $\|f\|_\infty = \max_{x \in [a, b]} |f(x)|$。子集 $\mathcal{F} \subset C([a, b])$ 是**列紧集**，当且仅当它同时满足：
   - **一致有界**：存在统一常数 $M > 0$，使对一切 $f \in \mathcal{F}$ 和一切 $x \in [a, b]$，都有 $|f(x)| \le M$；
   - **等度连续**：对任意 $\varepsilon > 0$，存在统一的 $\delta > 0$（只依赖于 $\varepsilon$，**与 $f \in \mathcal{F}$ 及点 $x_1, x_2$ 均无关**），使得只要 $|x_1 - x_2| < \delta$，就对一切 $f \in \mathcal{F}$ 同时有 $|f(x_1) - f(x_2)| < \varepsilon$。

### 0.2 赋范线性空间、Banach 空间、四大基础定理与 Riesz 引理

1. **赋范线性空间与 Banach 空间**：
   - **范数 $\|x\|$** 就是“向量的长度”，由它自然诱导距离 $d(x, y) = \|x - y\|$。
   - **Banach 空间** = **完备的赋范线性空间**。期末最常见的三类 Banach 空间：

| 空间记号 | 元素是什么 | 范数定义 | 备注 |
| :--- | :--- | :--- | :--- |
| $C([0,1])$ | $[0,1]$ 上的连续函数 $f$ | $\lVert f\rVert_\infty = \max_{x \in [0,1]} \lvert f(x)\rvert$ | 收敛即“一致收敛”，非自反 |
| $\ell^p$ ($1 \le p < \infty$) | 满足 $\sum_{n=1}^\infty \lvert x_n\rvert^p < \infty$ 的数列 | $\lVert x\rVert_p = \left(\sum_{n=1}^\infty \lvert x_n\rvert^p\right)^{1/p}$ | $p=2$ 是 Hilbert 空间；$p=1$ 非自反 |
| $L^p(\Omega)$ ($1 \le p < \infty$) | $\int_\Omega \lvert f\rvert^p < \infty$ 的可测函数（几乎处处相等视为同一元素） | $\lVert f\rVert_p = \left(\int_\Omega \lvert f(x)\rvert^p \mathrm{d}x\right)^{1/p}$ | $p=2$ 是 Hilbert 空间 |

2. **有界线性算子 $\mathcal{L}(X, Y)$ 与算子范数**：
   - 设 $X, Y$ 是赋范线性空间，线性映射 $A: X \to Y$ 若满足存在常数 $M \ge 0$ 使 $\|Ax\| \le M\|x\|$（$\forall x \in X$），则称 $A$ 为**有界线性算子**。对线性算子而言，**有界 $\iff$ 连续**。
   - 全体有界线性算子记为 $\mathcal{L}(X, Y)$（$X = Y$ 时简记为 $\mathcal{L}(X)$），其**算子范数**定义为：
     $$\|A\| = \sup_{\|x\| \le 1} \|Ax\| = \sup_{x \neq 0} \frac{\|Ax\|}{\|x\|}.$$
     （直觉：$A$ 对单位球最大的“拉伸倍数”。）若 $Y$ 是 Banach 空间，则 $\mathcal{L}(X, Y)$ 也是 Banach 空间。
   - 常用记号：核空间 $\operatorname{Ker} A = \{x \in X : Ax = 0\}$（总为闭子空间），值域 $\mathcal{R}(A) = \{Ax : x \in X\}$（**不一定闭！**）。
3. **Banach 空间三大核心定理（第 5 章多处用到！）**：
   - **开映射定理**：若 $X, Y$ 为 Banach 空间，$A \in \mathcal{L}(X, Y)$ 是**满射**，则 $A$ 把开集映成开集（特别地，$A$ 把 $X$ 的单位开球映成的集合必然包含 $Y$ 中某个以原点为心的开球 $B_Y(0, r)$）。→ **第 6 题 4 行秒杀证明的核心！**
   - **Banach 逆算子定理**：若 $X, Y$ 为 Banach 空间，$A \in \mathcal{L}(X, Y)$ 是**双射**（既单又满），则逆算子 $A^{-1}$ 自动有界，即 $A^{-1} \in \mathcal{L}(Y, X)$。→ **定理 5.5、定理 5.13、第 5/6 题核心！**
   - **闭图像定理**：若 $X, Y$ 为 Banach 空间，线性算子 $A: X \to Y$ 的图像 $\{(x, Ax) : x \in X\}$ 在 $X \times Y$ 中是闭集，则 $A$ 自动有界。
4. **Riesz 引理与无穷维的根本特征（第 5、6 题与引理 5.1 的核心！）**：
   - **Riesz 引理**：设 $M$ 是赋范线性空间 $X$ 的**真闭子空间**（$M \subsetneq X$）。对任意 $0 < \theta < 1$（例如 $\theta = \frac{1}{2}$），存在单位向量 $x_\theta \in X$（$\|x_\theta\| = 1$）使得它到 $M$ 的距离至少为 $\theta$：
     $$\operatorname{dist}(x_\theta, M) = \inf_{y \in M} \|x_\theta - y\| \ge \theta.$$
   - **最关键推论（闭单位球列紧 $\iff$ 有限维）**：赋范线性空间 $X$ 是**有限维的**，当且仅当其闭单位球 $B_X = \{x \in X : \|x\| \le 1\}$ 是**列紧集**。
     - 证明直觉：若 $\dim X = \infty$，任取单位向量 $x_1$；令 $M_1 = \operatorname{span}\{x_1\}$（有限维必真闭），由 Riesz 引理取单位向量 $x_2$ 使 $\|x_2 - x_1\| \ge \frac{1}{2}$；再令 $M_2 = \operatorname{span}\{x_1, x_2\}$ 取 $x_3$ 使到 $M_2$ 距离 $\ge \frac{1}{2}$……如此得单位向量列 $\{x_n\}$，任意两项距离 $\|x_m - x_n\| \ge \frac{1}{2}$ ($m \neq n$)，绝不可能有收敛子列！

### 0.3 Hilbert 空间、标准正交基与共轭算子

1. **内积与 Hilbert 空间**：
   - **Hilbert 空间** $H$ = 带有内积 $(\cdot, \cdot)$、且由诱导范数 $\|x\| = \sqrt{(x, x)}$ 构成**完备空间**的线性空间（即“带角度与垂直概念的无穷维欧氏空间”）。
   - **本教材内积约定**（务必看准）：**第一变元线性、第二变元共轭线性**，即对任意 $\alpha, \beta \in \mathbb{C}$：
     $$(\alpha x + \beta y, z) = \alpha(x, z) + \beta(y, z), \qquad (x, \alpha y + \beta z) = \bar{\alpha}(x, y) + \bar{\beta}(x, z), \qquad (y, x) = \overline{(x, y)}.$$
   - **Cauchy-Schwarz 不等式**：$|(x, y)| \le \|x\|\|y\|$。
2. **正交投影定理（第 15 题核心！）**：
   - 设 $M$ 是 Hilbert 空间 $H$ 的**闭**线性子空间，记其正交补 $M^\perp = \{y \in H : (y, x) = 0, \, \forall x \in M\}$。则 $M^\perp$ 也是闭子空间，且 $H$ 具有唯一的正交直和分解：
     $$H = M \oplus M^\perp,$$
     即任意 $x \in H$ 都可唯一写成 $x = x_1 + x_2$（$x_1 \in M, x_2 \in M^\perp$），且满足勾股定理 $\|x\|^2 = \|x_1\|^2 + \|x_2\|^2$。
3. **标准正交基与 Parseval 等式（5.2 节与第 15 题核心！）**：
   - 若 $\{e_n\}_{n=1}^\infty$ 是可分 Hilbert 空间 $H$ 的一组**标准正交基**（两两正交、单位长度、且极大/张成稠密子空间），则对任意 $x \in H$ 都有：
     - **Fourier 展开**：$x = \sum_{n=1}^\infty (x, e_n) e_n$；
     - **Parseval 等式（无穷维勾股定理）**：
       $$\|x\|^2 = \sum_{n=1}^\infty |(x, e_n)|^2.$$
     - 若仅满足标准正交（未必极大），则把等号换成 $\le$，称为 **Bessel 不等式**。
4. **Hilbert 空间上的共轭算子 $A^*$**：
   - 对 $A \in \mathcal{L}(H, G)$，唯一存在有界线性算子 $A^* \in \mathcal{L}(G, H)$ 使得对一切 $x \in H, y \in G$ 都有：
     $$(Ax, y)_G = (x, A^* y)_H.$$
   - 性质：$\|A^*\| = \|A\|$，$(A^*)^* = A$，$(AB)^* = B^* A^*$。（矩阵类比：复矩阵的共轭转置 $A^H = \overline{A}^T$。）

### 0.4 对偶空间、弱收敛、商空间与余维数

1. **对偶空间 $X'$ 与 Banach 空间共轭算子**：
   - 设 $X$ 为赋范线性空间，所有连续线性泛函（即值域在数域 $\mathbb{C}$ 或 $\mathbb{R}$ 中的有界线性算子）构成的空间 $X' = \mathcal{L}(X, \mathbb{C})$ 称为 $X$ 的**对偶空间**（或共轭空间）。
   - 对 $A \in \mathcal{L}(X, Y)$，其**Banach 空间共轭算子** $A^* \in \mathcal{L}(Y', X')$ 定义为：对任意 $y' \in Y'$，$(A^* y')(x) = y'(Ax)$（$\forall x \in X$）。同样有 $\|A^*\| = \|A\|$。
2. **强收敛 vs. 弱收敛（定理 5.1 与第 6 题反例用到！）**：
   - **强收敛**（按范数收敛）：$\lim_{n\to\infty} \|x_n - x_0\| = 0$，记作 $x_n \to x_0$。
   - **弱收敛**：对**每一个**有界线性泛函 $f \in X'$，都有数列的极限 $\lim_{n\to\infty} f(x_n) = f(x_0)$（在 Hilbert 空间中由 Riesz 表示定理等价于：对一切 $y \in H$，$(x_n, y) \to (x_0, y)$），常记作 $x_n \rightharpoonup x_0$。
   - **直觉对比**：在 $\ell^2$ 中，标准正交基 $\{e_n\}$ 满足对任意 $y = (y_1, y_2, \dots) \in \ell^2$ 都有 $(e_n, y) = \bar{y}_n \to 0$，所以 $e_n \rightharpoonup 0$（弱收敛到 $0$）；但 $\|e_n - 0\| = 1 \not\to 0$，所以**不强收敛**！
   - **重要特例**：强收敛必推弱收敛；反之一般不成立。但 $\ell^1$ 空间具有特殊的 **Schur 性质**：在 $\ell^1$ 中，**弱收敛与强收敛等价**（教材第 183 页注记 5.1）。
   - **自反空间**：若 $X$ 到其二次对偶 $X''$ 的自然嵌入是满射，则称 $X$ 为**自反 Banach 空间**。所有 Hilbert 空间和 $\ell^p, L^p$（$1 < p < \infty$）都自反；$C([0,1])$ 和 $\ell^1$ 不自反。自反空间的有界闭球弱列紧（任意有界列必有弱收敛子列）。
3. **商空间 $X / M$ 与余维数 $\operatorname{codim} M$（5.3 节 Fredholm 算子与第 6、22 题必考！）**：
   - 设 $M$ 是线性空间 $X$ 的线性子空间。把相差一个 $M$ 中元素的向量视为“同一类”，即定义等价类（陪集）$\tilde{x} = x + M = \{x + m : m \in M\}$。所有等价类组成的线性空间称为**商空间**，记作 $X / M$。
   - **几何直觉**：把 $\mathbb{R}^3$ 模掉 $xy$-平面 $M$，就是把所有水平高度相同的点捏成一个点，剩下的空间 $\mathbb{R}^3 / M$ 本质上就是 1 维的 $z$-轴！
   - **商范数**：当 $M$ 是 Banach 空间 $X$ 的**闭**子空间时，在 $X / M$ 上定义范数（几何上即点 $x$ 到闭子空间 $M$ 的距离）：
     $$\|\tilde{x}\| = \inf_{m \in M} \|x + m\| = \inf_{z \in \tilde{x}} \|z\|,$$
     此时 $X / M$ 成为一个 **Banach 空间**，且自然商映射 $q: X \to X/M$（$q(x) = \tilde{x}$）是连续线性满射（$\|q\| \le 1$）。
   - **余维数（codimension）**：定义子空间 $M$ 在 $X$ 中的余维数为商空间的维数：
     $$\operatorname{codim} M \triangleq \dim(X / M).$$
     **一句话理解**：$\operatorname{codim} M$ 就是“为了把 $M$ 补全成整个空间 $X$，还需要额外添加多少个线性无关的向量”。例如对算子 $A: X \to Y$，$\operatorname{codim} \mathcal{R}(A)$ 就是值域 $\mathcal{R}(A)$ 距离铺满 $Y$ 还“差了几个维度”。

---

## 一、第 5 章核心知识零基础精讲

本章恒设 $X, Y$ 为 Banach 空间（教材第 182 页约定）。整章围绕两条主线展开：
- **紧算子**：无穷维空间里“脾气最像有限维矩阵”的算子（把庞大松散的无穷维单位球“压扁、拧紧”成列紧集）；
- **Fredholm 算子**：无穷维空间里“几乎可逆”的算子（只在有限维上丢失信息、只在有限维上没覆盖到）。

### 5.1 紧算子 $\mathcal{C}(X, Y)$

#### 1. 定义与直觉（定义 5.1，教材第 182 页）

> **定义 5.1**：设 $A \in \mathcal{L}(X, Y)$。如果 $A$ 把 $X$ 中每个有界集都映射为 $Y$ 中的列紧集，则称 $A$ 为**紧算子**（或**全连续算子**）。$X$ 到 $Y$ 的全体紧算子记为 $\mathcal{C}(X, Y)$；当 $X = Y$ 时简记为 $\mathcal{C}(X)$。

- **等价序列刻画（做题最常用！）**：$A \in \mathcal{L}(X, Y)$ 是紧算子 $\iff$ 对 $X$ 中**任意有界序列** $\{x_n\}$，其像序列 $\{Ax_n\}$ 在 $Y$ 中**必有收敛子列**。
- **反面典型（例 5.2，教材第 183 页）**：在无穷维空间（如 $\ell^2$）上，**恒等算子 $I$ 不是紧算子**！因为标准正交基 $\{e_n\}$ 有界（$\|e_n\|=1$），但 $\|Ie_j - Ie_k\| = \|e_j - e_k\| = \sqrt{2}$（$j \neq k$），根本抽不出收敛子列。
- **正面典型（例 5.1，教材第 182–183 页）**：若 $K(s, t)$ 在 $[0,1] \times [0,1]$ 上连续，则 $C([0,1])$ 上的 Fredholm 积分算子
  $$(Ax)(s) = \int_0^1 K(s, t) x(t) \mathrm{d}t$$
  是紧算子（证明靠 Arzelà-Ascoli 定理：先证一致有界，再利用 $K$ 的一致连续性证等度连续）。**直觉**：积分有“平滑/抚平毛刺”的作用，所以积分算子常常是紧算子；微分算子会放大震荡，通常是无界且不紧的。

#### 2. 紧算子的核心性质（定理 5.1–5.8，教材第 183–187 页）

| 教材编号 | 精准结论 | 零基础直觉与考试用途 |
| :--------------- | :----------------------------------- | :----------------------------- |
| **定理 5.1**                                       | $A \in \mathcal{C}(X, Y)$ 且 $x_n \rightharpoonup x_0$（弱收敛） $\implies Ax_n \to Ax_0$（强收敛）。<br>（注记 5.1：$X$ 自反时逆命题成立；非自反时如 $\ell^1$ 不成立。）                                                                                                                                                                                                                                                                                                  | 紧算子有“变弱为强”的魔力；常用于处理弱收敛极限（如习题 13、16）。                                                                                                                                                                            |
| **定理 5.2**                                       | $A \in \mathcal{C}(X, Y) \implies \mathcal{R}(A)$ 是可分的（有可数稠密子集）。                                                                                                                                                                                                                                                                                                                                                                        | 说明紧算子的像集“体量不会太大”。                                                                                                                                                                                               |
| **定理 5.3**                                       | **(i)** $\alpha_1 A_1 + \alpha_2 A_2 \in \mathcal{C}(X, Y)$；<br>**(ii)** $B_1 \in \mathcal{L}(Y, Z), B_2 \in \mathcal{L}(Z, X) \implies A_1 B_2 \in \mathcal{C}(Z, Y), \; B_1 A_1 \in \mathcal{C}(X, Z)$。                                                                                                                                                                                                                               | **双边理想性质（超高频！）**：紧算子左乘或右乘**任何**有界算子，结果还是紧算子！第 5、20、22 题全靠它。                                                                                                                                                     |
| **定理 5.4**                                       | 若 $A_n \in \mathcal{C}(X, Y)$ 且 $\lVert A_n - A\rVert \to 0$，则 $A \in \mathcal{C}(X, Y)$。                                                                                                                                                                                                                                                                                                                                               | **按算子范数极限封闭**！证明算子紧的核武器：找一列有限秩算子 $A_n$ 逼近 $A$。                                                                                                                                                                  |
| **推论 5.1**                                       | $\mathcal{C}(X, Y)$ 是 $\mathcal{L}(X, Y)$ 的闭线性子空间，本身是 Banach 空间。                                                                                                                                                                                                                                                                                                                                                                        | 由定理 5.3 + 定理 5.4 立得。                                                                                                                                                                                            |
| **定理 5.5**                                       | 若 $X$ 无穷维、$A \in \mathcal{C}(X, Y)$ 且 $A$ 单射，则 $\mathcal{R}(A) \neq Y$。                                                                                                                                                                                                                                                                                                                                                                 | **无穷维紧算子不可能有界可逆！**（若可逆则 $I_X = A^{-1}A$ 紧，矛盾于 Riesz 引理。）→ 对应习题 5 第 5 题。                                                                                                                                         |
| **定理 5.6**                                       | $A \in \mathcal{C}(X), \, \lambda \neq 0 \implies \mathcal{R}(\lambda I - A)$ 是**闭子空间**。                                                                                                                                                                                                                                                                                                                                                | 单个紧算子 $A$ 的值域通常不闭，但加上 $\lambda I$ 后值域一定闭！                                                                                                                                                                       |
| **定理 5.7**<br>*(Schauder)*                       | $A \in \mathcal{C}(X, Y) \implies A^* \in \mathcal{C}(Y', X')$。（事实上互为充要。）                                                                                                                                                                                                                                                                                                                                                               | 紧性在取共轭时保持不变。                                                                                                                                                                                                    |
| **引理 5.1 &**<br>**定理 5.8**<br>*(Riesz-Fredholm)* | 设 $A \in \mathcal{C}(X), \, \lambda \neq 0$，则：<br>(1) $\mathcal{R}(\lambda I - A) = \mathcal{N}'(\operatorname{Ker}(\lambda I - A^*))$；<br>(2) $\mathcal{R}(\lambda I - A^*) = \mathcal{N}(\operatorname{Ker}(\lambda I - A))$；<br>(3) $\dim \operatorname{Ker}(\lambda I - A) = \dim \operatorname{Ker}(\lambda I - A^*)$；<br>(4) $\operatorname{codim} \mathcal{R}(\lambda I - A) = \dim \operatorname{Ker}(\lambda I - A) < \infty$。 | **超级核心结论**：<br>① **Fredholm 择一性**：对 $\lambda \neq 0$，$\lambda I - A$ 单射 $\iff$ 满射 $\iff$ 有界可逆（像极了有限维方阵！）；<br>② 令 $\lambda = 1$ 得 $I - A \in \mathcal{F}(X)$ 且 **$\operatorname{ind}(I - A) = 0$**（第 20 题关键基石！）。 |

> **教材勘误提示**：教材第 187 页定理 5.8(2) 末尾多印了一个“$< \infty$”，正确表述为 $\mathcal{R}(\lambda I - A^*) = \mathcal{N}(\operatorname{Ker}(\lambda I - A))$。

> **符号解码：$\mathcal{N}(M)$ 与 $\mathcal{N}'(N)$ 是什么意思？（教材 4.4 节零化子空间）**
>
> 1. **严格定义**：设 $X$ 是 Banach 空间、$X'$ 是其对偶空间，
>    - 对子集 $M \subset X$，其**零化子空间（Annihilator）**定义为
>      $$\mathcal{N}(M) \triangleq \{f \in X' : f(x) = 0,\; \forall x \in M\} \subset X',$$
>      大白话就是：$X'$ 中所有把 $M$ 里的向量全都“吃成 0”的有界线性泛函 $f$。
>    - 对子集 $N \subset X'$，其**前零化子空间（Pre-annihilator）**定义为
>      $$\mathcal{N}'(N) \triangleq \{x \in X : f(x) = 0,\; \forall f \in N\} \subset X,$$
>      大白话就是：$X$ 中所有被 $N$ 里的泛函全都“吃成 0”的向量 $x$。
> 2. **零基础直观理解（Banach 空间版的正交补 $M^\perp$）**：
>    在 Hilbert 空间或线性代数里，有内积可以写 $(x, y) = 0$ 表示“垂直”；一般 Banach 空间没有内积，就用对偶配对 $f(x) = 0$ 来代替“垂直”——因此 $\mathcal{N}(M)$ 在许多教材里直接记作 $M^\perp$！
> 3. **翻成大白话看懂定理 5.8(1)(2)**：
>    - **(1) $\mathcal{R}(\lambda I - A) = \mathcal{N}'(\operatorname{Ker}(\lambda I - A^*))$**：原方程 $(\lambda I - A)x = y$ 有解 $\iff$ 右端项 $y$ 被共轭齐次方程 $(\lambda I - A^*)f = 0$ 的所有解 $f$ 都作用成 $0$（即 $f(y) = 0$）。
>    - **(2) $\mathcal{R}(\lambda I - A^*) = \mathcal{N}(\operatorname{Ker}(\lambda I - A))$**：共轭方程 $(\lambda I - A^*)f = g$ 有解 $\iff$ 右端项 $g$ 把原齐次方程 $(\lambda I - A)x = 0$ 的所有解 $x$ 都作用成 $0$（即 $g(x) = 0$）。这正是大一线代里 $\operatorname{Im}(B^T) = (\operatorname{Ker} B)^\perp$ 的无穷维推广！

#### 3. 有限秩算子与逼近（定义 5.2、定理 5.9，教材第 189 页）

- **定义 5.2**：设 $A \in \mathcal{L}(X, Y)$，如果 $\dim \mathcal{R}(A) < \infty$，则称 $A$ 为**有限秩算子**。
- **结构与紧性（例 5.3）**：$A$ 是有限秩算子 $\iff$ 存在 $f_1, \dots, f_k \in X'$ 与 $y_1, \dots, y_k \in Y$ 使
  $$Ax = \sum_{j=1}^k f_j(x) y_j.$$
  由于有限维空间中的有界集必列紧，**任何有限秩算子都是紧算子**。
- **定理 5.9（可分 Hilbert 空间上的有限秩逼近）**：设 $A$ 是可分 Hilbert 空间 $H$ 上的紧算子，则存在一列有限秩算子 $\{A_n\}_{n=1}^\infty$ 使得 $\lim_{n\to\infty} \|A_n - A\| = 0$。（具体构造就是截断正交展开 $A_n x = \sum_{j=1}^n (x, e_j) A e_j$。）
  - **重要推论**：结合定理 5.4，在 Hilbert 空间上，**紧算子 $\iff$ 有限秩算子按算子范数的极限**！
  - *延伸说明*：对一般 Banach 空间，Enflo (1973) 证明存在不具逼近性质的 Banach 空间，但本科常见的 $C([0,1]), \ell^p, L^p$ 等空间都具备有限秩逼近性质。

#### 4. 期末实操：证明/否定“$A$ 是紧算子”的四大高频套路

| 题目类型 | 优先武器 | 怎么操作（3 步套路） |
| :--- | :--- | :--- |
| **套路 1**：$C([a,b])$ 上的积分算子 | **Arzelà-Ascoli 定理** | ① 先证 $A$ 线性有界；② 取有界集 $\lVert f\rVert_\infty \le M$，证像集**一致有界**；③ 估计 $\lvert Af(x_1) - Af(x_2)\rvert$ 证**等度连续**（如第 14 题）。 |
| **套路 2**：$\ell^p$ 上的对角/移位/矩阵算子 | **有限秩逼近（定理 5.4）** | ① 构造前 $N$ 项截断算子 $A_N$（显为有限秩 $\implies$ 紧）；② 估计余项算子范数 $\lVert A - A_N\rVert \to 0$；③ 由定理 5.4 立即得 $A$ 紧。 |
| **套路 3**：Hilbert 空间（尤其是 $L^2, \ell^2$）上的积分/矩阵算子 | **Hilbert-Schmidt 判定（定理 5.12(5)）** | 验证 $\sum_{n=1}^\infty \lVert Ae_n\rVert^2 < \infty$（对 $L^2$ 积分算子即验证 $\iint \lvert K(x,y)\rvert^2 \mathrm{d}x\mathrm{d}y < \infty$），立即推出紧性！ |
| **套路 4**：证明算子**不是**紧算子 / 紧算子不可能具备某性质 | **反证法 + Riesz 引理（或直接找反例有界列）** | ① 反设紧性，利用理想性质凑出 $I_X \in \mathcal{C}(X)$（第 5、6 题）；② 或直接构造有界列 $\{x_n\}$ 使 $\lVert Ax_m - Ax_n\rVert \ge c > 0$（例 5.2）。 |

---

### 5.2 Hilbert-Schmidt 算子 $\mathcal{L}^2(H, G)$

**零基础直觉**：在线性代数中，衡量矩阵 $A = (a_{ij})$ 大小的“平方和开根号”叫 **Frobenius 范数** $\|A\|_F = \left(\sum_{i,j} |a_{ij}|^2\right)^{1/2}$。把矩阵的列向量 $Ae_n$ 的长度平方加起来再开方，推广到无穷维 Hilbert 空间，就是 **Hilbert-Schmidt 范数**！

本节默认 $H, G$ 为可分 Hilbert 空间（定义对一般 Hilbert 空间亦可自然推广）。

#### 1. 定义与基无关性（定义 5.4、式 (5.8)，教材第 196 页）

> **定义 5.4**：设 $A$ 是 $H$ 到 $G$ 的线性算子。如果存在 $H$ 上的一组标准正交基 $\{e_n\}_{n=1}^\infty$ 使得
> $$\|A\|_2 \triangleq \left(\sum_{n=1}^\infty \|Ae_n\|^2\right)^{\frac{1}{2}} < \infty,$$
> 则称 $A$ 是 $H$ 到 $G$ 的 **Hilbert-Schmidt 算子**，全体记为 $\mathcal{L}^2(H, G)$；当 $H = G$ 时简记为 $\mathcal{L}^2(H)$。

- **为什么与标准正交基的选取无关？（式 (5.8)）**：设 $\{e_m\}, \{\tilde{e}_m\}$ 为 $H$ 的两组标准正交基，$\{f_n\}$ 为 $G$ 的标准正交基，由 Parseval 等式与 $A^*$ 的定义：
  $$\sum_{m=1}^\infty \|Ae_m\|^2 = \sum_{m,n=1}^\infty |(Ae_m, f_n)|^2 = \sum_{m,n=1}^\infty |(e_m, A^* f_n)|^2 = \sum_{n=1}^\infty \|A^* f_n\|^2 = \sum_{m,n=1}^\infty |(\tilde{e}_m, A^* f_n)|^2 = \sum_{m=1}^\infty \|A\tilde{e}_m\|^2.$$
  这个等式一石二鸟：同时证明了 **$\|A\|_2$ 与基的选取无关** 和 **$\|A\|_2 = \|A^*\|_2$**！
  *(注：教材第 196 页式 (5.8) 中将复内积的模平方 $|(Ae_m, f_n)|^2$ 漏写了绝对值符号，做题时请务必写绝对值。)*

#### 2. Hilbert-Schmidt 算子的六大性质（定理 5.12，教材第 197–198 页）

设 $A \in \mathcal{L}^2(H, G)$，则：
1. **共轭等距**：$\|A\|_2 = \|A^*\|_2$，故 $A^* \in \mathcal{L}^2(G, H)$；
2. **控制算子范数**：$\|A\| \le \|A\|_2$（特别地，HS 算子自动是有界算子）；
3. **双边理想与范数不等式**：若 $B_1 \in \mathcal{L}(H), B_2 \in \mathcal{L}(G)$，则 $AB_1, B_2 A \in \mathcal{L}^2(H, G)$，且
   $$\|B_2 A\|_2 \le \|B_2\|\|A\|_2, \qquad \|AB_1\|_2 \le \|B_1\|\|A\|_2;$$
4. **自身构成 Hilbert 空间**：在定义于 $\mathcal{L}^2(H, G)$ 的内积
   $$(A_1, A_2)_2 \triangleq \sum_{n=1}^\infty (A_1 e_n, A_2 e_n)$$
   及诱导范数 $\|\cdot\|_2$ 下，$\mathcal{L}^2(H, G)$ 是 Hilbert 空间；
5. **必为紧算子**：每个 Hilbert-Schmidt 算子都是**紧算子**，且有限秩算子按 $\|\cdot\|_2$ 范数在 $\mathcal{L}^2(H, G)$ 中稠密；
6. **可分性**：$\mathcal{L}^2(H, G)$ 是可分 Hilbert 空间。

- **必背算子包含链（在无穷维可分 Hilbert 空间上均为严格真包含）**：
  $$\{\text{有限秩算子}\} \;\subsetneq\; \mathcal{L}^2(H, G) \;\subsetneq\; \mathcal{C}(H, G) \;\subsetneq\; \mathcal{L}(H, G).$$
  *(对角算子记忆法：$A e_n = a_n e_n$。有限秩 $\iff$ 仅有限个 $a_n \neq 0$；HS $\iff \sum |a_n|^2 < \infty$；紧 $\iff a_n \to 0$；有界 $\iff \{a_n\}$ 有界。例如 $a_n = 1/\sqrt{n}$ 紧但非 HS，$a_n = 1/n$ 是 HS 但非有限秩！)*

---

### 5.3 Fredholm 算子 $\mathcal{F}(X, Y)$

**零基础直觉**：一个有界线性算子 $A: X \to Y$ 可逆，意味着 $\dim \operatorname{Ker} A = 0$（不压缩任何维度）且 $\operatorname{codim} \mathcal{R}(A) = 0$（完全铺满 $Y$）。Fredholm 算子把“等于 0”放宽到“**有限**”——也就是**“几乎可逆”的算子**！

#### 1. 定义与指标（定义 5.5、定义 5.6，教材第 198–199 页）

> **定义 5.5**：设 $X, Y$ 均为 Banach 空间，$A \in \mathcal{L}(X, Y)$。如果同时满足：
> 1. $\mathcal{R}(A)$ 是 $Y$ 中的**闭**子空间；
> 2. $\dim \operatorname{Ker}(A) < +\infty$（核空间有限维）；
> 3. $\operatorname{codim} \mathcal{R}(A) = \dim(Y / \mathcal{R}(A)) < +\infty$（值域余维有限），
>
> 则称 $A$ 为 **Fredholm 算子**，全体记为 $\mathcal{F}(X, Y)$；当 $X = Y$ 时简记为 $\mathcal{F}(X)$。

> **重要注记（深度理解：为什么条件 (3) $\operatorname{codim}\mathcal{R}(A) < +\infty$ 能自动推出条件 (1) $\mathcal{R}(A)$ 闭？）**：
>
> - **先避开一个直觉误区**：一个普通的有限余维线性子空间**并不一定是闭的**！（反例：在无穷维 Banach 空间中，不连续线性泛函 $f$ 的核 $\operatorname{Ker} f$ 余维数为 $1$，却在空间中处处稠密而不闭。）
> - **这里为何能自动变闭？** 全靠两大前提联手：**$X, Y$ 完备（Banach 空间）** $+$ **$\mathcal{R}(A)$ 是有界线性算子 $A$ 的值域**！既然 $\mathcal{R}(A)$ 距铺满 $Y$ 只差一个有限维小缺口 $W$，我们把 $W$ 补上去就能造出一个 Banach 空间之间的连续双射，再用开映射定理（Banach 逆算子定理）直接逼出 $\mathcal{R}(A)$ 必闭！
> - **5 步完整证明（Kato 引理）**：设 $\operatorname{codim}\mathcal{R}(A) = \dim(Y/\mathcal{R}(A)) = k < +\infty$。
>   1. **消核变单射**：因 $A$ 连续，$\operatorname{Ker}A$ 闭，故商空间 $\tilde{X} \triangleq X/\operatorname{Ker}A$ 是 Banach 空间，诱导算子 $A_0([x]) \triangleq Ax$ 为有界单射，且 $\mathcal{R}(A_0) = \mathcal{R}(A)$。
>   2. **找 $k$ 维补丁 $W$**：取 $Y$ 的 $k$ 维子空间 $W$ 使代数直和 $Y = \mathcal{R}(A) \oplus W$ 成立。
>   3. **拼成直和 Banach 空间**：令 $Z \triangleq \tilde{X} \oplus W$，配范数 $\lVert([x], w)\rVert_Z \triangleq \lVert[x]\rVert_{\tilde{X}} + \lVert w\rVert_Y$。因 $\tilde{X}$ 完备、$W$ 有限维必完备，故 $Z$ 是 Banach 空间。
>   4. **构造连续双射 $T: Z \to Y$**：定义 $T([x], w) \triangleq A_0[x] + w = Ax + w$。显然 $T$ 有界（$\lVert T([x],w)\rVert_Y \le \max(\lVert A_0\rVert, 1)\lVert([x],w)\rVert_Z$）；由 $Y = \mathcal{R}(A) \oplus W$ 与 $A_0$ 单射立得 $T$ 既是单射又是满射。
>   5. **开映射定理收尾**：由 Banach 逆算子定理（开映射定理推论），$T^{-1}: Y \to Z$ 连续，即 $T$ 是拓扑同构（把闭集映成闭集）。因 $\tilde{X} \oplus \{0\}$ 显然是 $Z$ 中的闭子空间，故其像 $\mathcal{R}(A) = T(\tilde{X} \oplus \{0\})$ 必为 $Y$ 中的**闭子空间**！$\blacksquare$
> - **考场提醒**：虽然条件 (3) 能推出条件 (1)，但教材定义 5.5 明确列出了全部 3 条，**闭卷写定义与做常规证明题时仍请按定义 5.5 把这 3 条写全（尤其别漏证“$\mathcal{R}(A)$ 闭”）！**

> **定义 5.6（Fredholm 指标）**：设 $A \in \mathcal{F}(X, Y)$，定义 $A$ 的**指标（index）**为整数：
> $$\operatorname{ind} A = \dim \operatorname{Ker}(A) - \operatorname{codim} \mathcal{R}(A).$$

- **三个必记标杆例子（例 5.8–5.9，教材第 199 页）**：
  1. **可逆算子**（注记 5.3）：$\dim \operatorname{Ker} A = 0, \, \operatorname{codim} \mathcal{R}(A) = 0 \implies A \in \mathcal{F}(X, Y)$ 且 $\operatorname{ind} A = 0$。
  2. **紧算子的恒等扰动 $I - K$**（例 5.8）：对任意 $K \in \mathcal{C}(X)$，由定理 5.6 与定理 5.8(4) 知 $I - K \in \mathcal{F}(X)$ 且 **$\operatorname{ind}(I - K) = 0$**！（第 20 题必用！）
  3. **$\ell^2$ 上的左/右平移算子**（例 5.9，最能体现“无穷维与有限维不同”的例子！）：
     - 左平移 $L(x_1, x_2, x_3, \dots) = (x_2, x_3, \dots)$：丢掉了第 1 个分量，故 $\operatorname{Ker} L = \operatorname{span}\{e_1\}$（1 维），但它是满射 $\mathcal{R}(L) = \ell^2$（闭，余维 0），所以 $L \in \mathcal{F}(\ell^2)$ 且 **$\operatorname{ind} L = 1 - 0 = 1$**！
     - 右平移 $R(x_1, x_2, \dots) = (0, x_1, x_2, \dots) = L^*$：单射（核 0 维），但值域是第 1 个分量为 0 的闭超平面（余维 1），所以 $R \in \mathcal{F}(\ell^2)$ 且 **$\operatorname{ind} R = 0 - 1 = -1$**！
     - 进一步：$\operatorname{ind}(L^n) = n$，$\operatorname{ind}(R^n) = -n$。*(注：有限维方阵的秩-零度定理强制要求 $\dim \operatorname{Ker} A - \operatorname{codim} \mathcal{R}(A) = 0$，只有无穷维才会出现非零指标！)*

#### 2. Fredholm 算子的三大主定理（定理 5.13–5.15，教材第 199–203 页）

1. **定理 5.13（Atkinson 定理 — 判定 Fredholm 算子的超级神器！）**：
   - **(1) 必要性**：若 $A \in \mathcal{F}(X, Y)$，则存在 $S \in \mathcal{L}(Y, X)$（称为 $A$ 的正则化子/拟逆）及紧算子 $A_1 \in \mathcal{C}(X), \, A_2 \in \mathcal{C}(Y)$，使得
     $$SA = I_X - A_1, \qquad AS = I_Y - A_2.$$
     *(注记 5.4：此处的 $S$ 自动也满足 $S \in \mathcal{F}(Y, X)$！)*
   - **(2) 充分性**：若 $A \in \mathcal{L}(X, Y)$，且存在 $S_1, S_2 \in \mathcal{L}(Y, X)$ 及紧算子 $A_1 \in \mathcal{C}(X), \, A_2 \in \mathcal{C}(Y)$，使得
     $$S_1 A = I_X - A_1, \qquad A S_2 = I_Y - A_2,$$
     则 $A \in \mathcal{F}(X, Y)$。
   - **一句话核心思想（注记 5.5）**：**Fredholm 算子 $\iff$ “模掉紧算子之后可逆”！**（即在商代数 $\mathcal{L}(X)/\mathcal{C}(X)$ 中可逆。）只要能凑出“乘完等于 $I$ 减紧算子”，就可以完全绕开繁琐的“核维数 + 余维数 + 值域闭”逐项检验！*(教材勘误：第 200 页定理 5.13(2) 末尾把“则 $A \in \mathcal{F}(X, Y)$”误印成了 $T$。)*
2. **定理 5.14（乘积指标公式 — 对数法则）**：
   设 $X, Y, Z$ 为 Banach 空间，$A_1 \in \mathcal{F}(X, Y), \, A_2 \in \mathcal{F}(Y, Z)$，则 $A_2 A_1 \in \mathcal{F}(X, Z)$ 且
   $$\operatorname{ind}(A_2 A_1) = \operatorname{ind} A_2 + \operatorname{ind} A_1.$$
3. **定理 5.15（小扰动稳定性）**：
   设 $A \in \mathcal{F}(X, Y)$，则存在 $\varepsilon > 0$，使得对一切满足 $\|S\| < \varepsilon$ 的 $S \in \mathcal{L}(X, Y)$，均有 $A + S \in \mathcal{F}(X, Y)$ 且
   $$\operatorname{ind}(A + S) = \operatorname{ind} A.$$
   *(直觉：Fredholm 算子集在 $\mathcal{L}(X, Y)$ 中是开集；指标是整数值连续函数，微小扰动掰不动它！)*
4. **紧扰动不变性（即习题 5 第 20 题）**：
   设 $B \in \mathcal{F}(X, Y)$，$A \in \mathcal{C}(X, Y)$，则 $B + A \in \mathcal{F}(X, Y)$ 且 $\operatorname{ind}(B + A) = \operatorname{ind} B$。*(不管紧扰动 $A$ 的范数有多大，都不改变 Fredholm 性与指标！)*

---

## 二、第 5 章课后作业逐题详解（第 5、6、14、15、20、22 题）

---

### 习题 5 第 5 题（教材第 204 页）

#### 【原题】
> 5. 设 $X$ 为无穷维 Banach 空间，$A \in \mathcal{C}(X)$，证明：$A$ 没有有界逆。

#### 【零基础破题思路】
- **第一反应**：题目出现两个黄金关键词——**“无穷维”** + **“紧算子 $A$”** + **“没有有界逆”**。只要在无穷维空间里证明一个紧算子做不到某件事，**99% 是用反证法导出恒等算子 $I_X$ 是紧算子，再用 Riesz 引理矛盾！**
- **逻辑链条**：如果 $A$ 有定义在全空间上的有界逆 $A^{-1} \in \mathcal{L}(X)$，那么 $I_X = A^{-1} A$。由紧算子的双边理想性质（定理 5.3(ii)），有界算子乘以紧算子还是紧算子，所以 $I_X \in \mathcal{C}(X)$。但 $I_X$ 紧意味着 $X$ 的闭单位球列紧，这与“$X$ 是无穷维空间”（Riesz 引理）直接冲突！
- **更严谨的覆盖**：教材和试卷中“$A$ 没有有界逆”通常指不存在全空间上的 $A^{-1} \in \mathcal{L}(X)$；但有时也指更强的情形——即便 $A$ 不是满射，仅要求 $A$ 单射且在其值域 $\mathcal{R}(A)$ 上有界可逆（即 $A^{-1}: \mathcal{R}(A) \to X$ 有界，等价于存在下界 $\|Ax\| \ge c\|x\|$）。下面我们先给出标准简洁证法，再把“值域上有界逆也不存在”的加强版证法一并补齐，做到滴水不漏！

#### 【规范解答】

**证明（标准版：不存在全空间上的有界逆 $A^{-1} \in \mathcal{L}(X)$）：**

用反证法。假设 $A$ 存在有界逆 $A^{-1} \in \mathcal{L}(X)$，则有
$$I_X = A^{-1} A,$$
其中 $I_X$ 为 $X$ 上的恒等算子。

因为 $A \in \mathcal{C}(X)$ 且 $A^{-1} \in \mathcal{L}(X)$，由紧算子的理想性质（定理 5.3(ii)）可知
$$I_X = A^{-1} A \in \mathcal{C}(X).$$

令 $B_X = \{x \in X : \|x\| \le 1\}$ 为 $X$ 中的闭单位球。由于 $B_X$ 是有界集且 $I_X$ 是紧算子，根据紧算子定义（定义 5.1），像集 $I_X(B_X) = B_X$ 是 $X$ 中的列紧集。

然而，$X$ 是无穷维赋范线性空间，由 Riesz 引理可知 $X$ 的闭单位球 $B_X$ 不是列紧集（否则可利用 Riesz 引理递推构造出单位向量列 $\{x_n\}_{n=1}^\infty \subset B_X$ 满足 $\|x_m - x_n\| \ge \frac{1}{2}$ 对一切 $m \neq n$，它没有任何收敛子列），矛盾！

因此，$A$ 没有有界逆。 $\blacksquare$

---

**补充证明（加强版：即便不要求 $A$ 满射，也不存在定义在 $\mathcal{R}(A)$ 上的有界逆 $A^{-1}: \mathcal{R}(A) \to X$）：**

假设 $A$ 单射且 $A^{-1}: \mathcal{R}(A) \to X$ 是有界线性算子，记 $M = \|A^{-1}\| > 0$，则对一切 $x \in X$ 均有
$$\|x\| = \|A^{-1}(Ax)\| \le M \|Ax\|.$$

任取 $X$ 中的有界序列 $\{x_n\}_{n=1}^\infty$。因为 $A \in \mathcal{C}(X)$，所以像序列 $\{Ax_n\}_{n=1}^\infty$ 在 $X$ 中必有收敛子列 $\{Ax_{n_k}\}_{k=1}^\infty$。特别地，$\{Ax_{n_k}\}$ 是 Cauchy 列。于是对任意 $j, k \in \mathbb{N}$ 有
$$\|x_{n_k} - x_{n_j}\| \le M \|A(x_{n_k} - x_{n_j})\| = M \|Ax_{n_k} - Ax_{n_j}\| \to 0 \qquad (j, k \to \infty).$$

由此可见 $\{x_{n_k}\}_{k=1}^\infty$ 是 $X$ 中的 Cauchy 列。由于 $X$ 是 Banach 空间（完备），$\{x_{n_k}\}$ 在 $X$ 中收敛。这表明 $X$ 中任意有界序列都有收敛子列，即闭单位球 $B_X$ 是列紧集，与 $\dim X = \infty$（Riesz 引理）矛盾！故即使仅限定在 $\mathcal{R}(A)$ 上，$A$ 也没有有界逆。 $\blacksquare$

#### 【考点复盘与避坑提示】
1. **与定理 5.5 的关系**：定理 5.5 证明了“单射紧算子不可能是满射”（实际上同样是通过反设满射 $\Rightarrow$ Banach 逆算子定理得 $A^{-1} \in \mathcal{L}(Y, X)$ $\Rightarrow$ $I_X = A^{-1}A$ 紧，得出矛盾）。本题直接从定义 + 定理 5.3 + Riesz 引理出发，自包含且最稳。
2. **加强版的隐藏考点**：在证明加强版时，**必须用到 $X$ 的完备性**（从 Cauchy 列推出收敛子列）；如果只证标准版，只需 $X$ 是无穷维赋范线性空间即可。
3. **常见误区**：注意题干给的是“无穷维”！有限维空间中每个线性算子都是紧算子，且可逆矩阵显然有有界逆，所以千万别漏写“与 $\dim X = \infty$（Riesz 引理）矛盾”这一步！

---

### 习题 5 第 6 题（教材第 204 页）【重点排雷题！】

#### 【原题】
> 6. 设 $X$ 为无穷维 Banach 空间. 证明：$A \in \mathcal{C}(X)$ 当且仅当 $\mathcal{R}(A)$ 不包含无穷维的子空间.

> **⚠️ 极其重要的严谨性排雷（必看！）**：
> 原题表述存在 **两处实质性的数学瑕疵**，在理解和作答时必须心中有数：
> 1. **漏写了关键限定词“闭”**：如果不要求子空间是**闭**的，那么连必要性 $(\Rightarrow)$ 都是**假命题**！
>    - **反例**：取 $X = \ell^2$，定义紧算子 $A(x_1, x_2, \dots) = \left(x_1, \frac{x_2}{2}, \dots, \frac{x_n}{n}, \dots\right)$（显然是 Hilbert-Schmidt 算子，故 $A \in \mathcal{C}(\ell^2)$）。对每个标准基向量 $e_n$，都有 $e_n = A(ne_n) \in \mathcal{R}(A)$，于是所有有限非零项序列构成的子空间 $c_{00} = \operatorname{span}\{e_1, e_2, \dots\} \subset \mathcal{R}(A)$，它显然是**无穷维**线性子空间（只是在 $\ell^2$ 中不闭）。
>    - **教材自身佐证**：同书下一页（第 205 页）习题 17(ii) 就正确地表述为“如果**闭**线性子空间 $M \subset \mathcal{R}(A)$，则 $\dim M < +\infty$”。
> 2. **补上“闭”字之后，充分性 $(\Leftarrow)$ 仅在 Hilbert 空间上成立，在一般 Banach 空间上仍不成立！**
>    - **必要性 $(\Rightarrow)$**（若 $A \in \mathcal{C}(X)$，则 $\mathcal{R}(A)$ 不包含无穷维**闭**子空间）：对所有 Banach 空间**完全成立**！
>    - **充分性 $(\Leftarrow)$**（若 $A \in \mathcal{L}(X)$ 且 $\mathcal{R}(A)$ 不包含无穷维闭子空间，则 $A \in \mathcal{C}(X)$）：
>      - 在 **Hilbert 空间** 上**成立**；
>      - 在 **一般 Banach 空间** 上**不成立**！（满足“值域不含无穷维闭子空间”的算子在泛函分析中称为“严格奇异算子”，范围严格大于紧算子；下文给出使用教材第 183 页注记 5.1 即可完全看懂的标准反例！）

#### 【零基础破题思路】
1. **怎么证必要性 $(\Rightarrow)$？（核心必考，对应习题 17(ii)！）**
   - 假设存在一个**闭**子空间 $M \subset \mathcal{R}(A)$，我们的目标是证 $\dim M < \infty$，也就是证 **$M$ 的某个闭球是列紧集**（再靠 Riesz 引理秒杀！）。
   - 难点在于：$A$ 未必是单射，所以不能直接写 $A^{-1}$！怎么把 $A$ 变成从某个 Banach 空间到 $M$ 的满射/双射？
   - **核心桥梁**：把 $A$ 的定义域收缩到原像 $X_0 = A^{-1}(M) = \{x \in X : Ax \in M\}$ 上！因为 $A$ 连续、$M$ 闭，所以 $X_0$ 也是闭子空间（即 Banach 空间），且 $A|_{X_0}: X_0 \to M$ 是**连续满射**！
   - 接下来有两条康庄大道：
     - **路子 A（开映射定理，4 行秒杀）**：$A|_{X_0}: X_0 \to M$ 是 Banach 空间间的满射 $\implies$ 它把 $X_0$ 的单位开球映成的集合包含 $M$ 里的某个开球 $B_M(0, r)$ $\implies$ $B_M(0, r)$ 被包在一个列紧集 $A(B_{X_0}(0, 1))$ 里 $\implies$ $M$ 的闭球列紧 $\implies$ $\dim M < \infty$！
     - **路子 B（商空间 + Banach 逆算子定理，与教材风格一致）**：模掉 $\operatorname{Ker} A$，得到双射 $\tilde{A}: X_0 / \operatorname{Ker} A \to M$，由逆算子定理得 $\tilde{A}^{-1}$ 有界，再对 $M$ 的单位球里的任意序列利用商范数拉回 $X_0$ 中找有界代表元，借 $A$ 的紧性抽出收敛子列 $\implies$ $\dim M < \infty$！
2. **怎么理解充分性 $(\Leftarrow)$ 的一般 Banach 空间反例？**
   - 我们需要一个“不是紧算子、但值域里装不下任何无穷维闭子空间”的算子。利用两个“几何性质格格不入”的空间之间的恒等嵌入 $J: \ell^1 \hookrightarrow \ell^2$：
     - $J$ 不紧（$Je_n = e_n$ 间距恒为 $\sqrt{2}$）；
     - 如果 $J(\ell^1)$ 里藏着一个无穷维闭子空间 $M_2 \subset \ell^2$，那 $M_2$ 就会通过 $J^{-1}$ 拓扑同构于 $\ell^1$ 的一个闭子空间 $N$——可 $M_2$ 里有“弱收敛于 $0$ 但不强收敛”的正交列，而 $\ell^1$ 里弱收敛等价于强收敛（Schur 性质，教材第 183 页注记 5.1），直接爆出矛盾！
   - 为凑成 $X \to X$ 的算子，只需把两个空间拼成直和 $X = \ell^1 \oplus \ell^2$，令 $A(x, y) = (0, Jx)$ 即可！

#### 【规范解答】

**一、修正原题为“不含无穷维闭子空间”后的必要性 $(\Rightarrow)$ 严格证明**

> **待证命题**：设 $X$ 为 Banach 空间，$A \in \mathcal{C}(X)$。若 $M$ 是 $X$ 的**闭**线性子空间且 $M \subset \mathcal{R}(A)$，则 $\dim M < +\infty$。

**证法一（商空间 + Banach 逆算子定理，最贴合教材 5.1/5.3 节风格）：**

1. **构造原像子空间**：令
   $$X_0 = A^{-1}(M) = \{x \in X : Ax \in M\}.$$
   因为 $A \in \mathcal{L}(X)$ 连续且 $M$ 是闭子空间，所以 $X_0$ 是 $X$ 的闭线性子空间，从而 $X_0$ 与 $M$ 均为 Banach 空间。
2. **构造商空间上的双射**：显然 $\operatorname{Ker} A \subset X_0$ 且 $\operatorname{Ker} A$ 是 $X_0$ 的闭线性子空间。作商空间 $\tilde{X} = X_0 / \operatorname{Ker} A$，赋予商范数
   $$\|\tilde{x}\| = \inf_{z \in \operatorname{Ker} A} \|x + z\|, \qquad \tilde{x} = x + \operatorname{Ker} A \in \tilde{X},$$
   则 $\tilde{X}$ 是 Banach 空间。定义映射 $\tilde{A}: \tilde{X} \to M$ 为
   $$\tilde{A}(\tilde{x}) = Ax \qquad (x \in \tilde{x}).$$
   - **良定性**：若 $x_1 - x_2 \in \operatorname{Ker} A$，则 $Ax_1 = Ax_2$，故 $\tilde{A}$ 与代表元选取无关。
   - **线性与有界性**：$\tilde{A}$ 显然线性；对任意 $x \in \tilde{x}$ 与 $z \in \operatorname{Ker} A$，有 $\|\tilde{A}\tilde{x}\| = \|A(x + z)\| \le \|A\|\|x + z\|$，对所有 $z \in \operatorname{Ker} A$ 取下确界得 $\|\tilde{A}\tilde{x}\| \le \|A\|\|\tilde{x}\|$，故 $\tilde{A} \in \mathcal{L}(\tilde{X}, M)$。
   - **单射与满射**：若 $\tilde{A}\tilde{x} = 0$，则 $Ax = 0 \implies x \in \operatorname{Ker} A \implies \tilde{x} = \tilde{0}$，故单射；对任意 $y \in M \subset \mathcal{R}(A)$，存在 $x \in X$ 使 $Ax = y \in M$，于是 $x \in X_0$ 且 $\tilde{A}\tilde{x} = y$，故满射。
3. **利用逆算子定理与紧性**：$\tilde{A}$ 是两个 Banach 空间之间的有界线性双射，由 **Banach 逆算子定理**，其逆算子 $\tilde{A}^{-1} \in \mathcal{L}(M, \tilde{X})$ 有界。
4. **证明 $M$ 的闭单位球列紧**：任取 $M$ 的闭单位球中的序列 $\{y_n\}_{n=1}^\infty \subset M$（$\|y_n\| \le 1$）。令 $\tilde{x}_n = \tilde{A}^{-1} y_n \in \tilde{X}$，则
   $$\|\tilde{x}_n\| \le \|\tilde{A}^{-1}\|\|y_n\| \le \|\tilde{A}^{-1}\|.$$
   根据商范数的下确界定义，对每个 $n \in \mathbb{N}$，可选取代表元 $x_n \in \tilde{x}_n \subset X_0$ 使得
   $$\|x_n\| \le \|\tilde{x}_n\| + 1 \le \|\tilde{A}^{-1}\| + 1, \qquad \text{且 } Ax_n = \tilde{A}\tilde{x}_n = y_n.$$
   由此可见 $\{x_n\}_{n=1}^\infty$ 是 $X$ 中的有界序列。因为 $A \in \mathcal{C}(X)$ 是紧算子，序列 $\{Ax_n\} = \{y_n\}$ 在 $X$ 中必有收敛子列 $\{y_{n_k}\}$；又因 $M$ 是闭子空间，该子列的极限仍落在 $M$ 中。
5. **结论**：$M$ 的闭单位球是列紧集。由 **Riesz 引理**，必有 $\dim M < +\infty$。 $\blacksquare$

---

**证法二（开映射定理，4 行极速版）：**

令 $X_0 = A^{-1}(M)$。因为 $A$ 连续且 $M$ 闭，故 $X_0$ 与 $M$ 均为 Banach 空间，且限制算子 $A_0 = A|_{X_0}: X_0 \to M$ 是有界线性满射。记 $U = \{x \in X_0 : \|x\| < 1\}$ 为 $X_0$ 的单位开球，由**开映射定理**，$A_0(U)$ 是 $M$ 中包含 $0$ 的开集，故存在 $r > 0$ 使 $M$ 中的闭球 $B_M(0, r) = \{y \in M : \|y\| \le r\} \subset A_0(U) = A(U)$。因为 $U$ 有界且 $A$ 为紧算子，$A(U)$ 列紧，其子集 $B_M(0, r)$ 亦列紧（且因 $M$ 闭，子列极限必在 $M$ 中），从而 $M$ 的闭单位球列紧。由 **Riesz 引理**，$\dim M < +\infty$。 $\blacksquare$

---

**二、充分性 $(\Leftarrow)$ 的深度分析：一般 Banach 空间反例与 Hilbert 空间说明**

1. **一般 Banach 空间中充分性 $(\Leftarrow)$ 不成立（标准反例）**：
   取 Banach 空间 $X = \ell^1 \oplus \ell^2$，赋予范数 $\|(x, y)\| = \|x\|_1 + \|y\|_2$。记 $J: \ell^1 \to \ell^2$ 为自然恒等嵌入 $Jx = x$（因为 $\|x\|_2 \le \|x\|_1$，故 $J$ 线性有界）。定义 $A \in \mathcal{L}(X)$ 为
   $$A(x, y) = (0, Jx), \qquad (x, y) \in \ell^1 \oplus \ell^2.$$
   - **$A$ 不是紧算子**：取标准单位向量 $u_n = (e_n, 0) \in X$，有 $\|u_n\| = \|e_n\|_1 = 1$，但 $Au_n = (0, e_n)$ 满足 $\|Au_m - Au_n\| = \|e_m - e_n\|_2 = \sqrt{2}$（$m \neq n$），无收敛子列，故 $A \notin \mathcal{C}(X)$。
   - **$\mathcal{R}(A)$ 不包含无穷维闭子空间**：反设 $M \subset \mathcal{R}(A) = \{0\} \oplus J(\ell^1)$ 是无穷维闭子空间，将其投影到第二分量得到 $\ell^2$ 的无穷维闭子空间 $M_2 \subset J(\ell^1)$。令 $N = J^{-1}(M_2) \subset \ell^1$，则 $N$ 是 $\ell^1$ 的闭子空间，且 $J|_N: N \to M_2$ 为 Banach 空间之间的连续双射，由 Banach 逆算子定理，$N$ 与 $M_2$ 拓扑同构。
     作为无穷维 Hilbert 空间 $\ell^2$ 的无穷维闭子空间，$M_2$ 中存在标准正交列 $\{g_n\}_{n=1}^\infty$：由 Bessel 不等式，对任意 $h \in M_2$ 有 $(g_n, h) \to 0$，即 $g_n \rightharpoonup 0$（弱收敛于 $0$），但 $\|g_n\|_2 = 1$ 不强收敛于 $0$。令 $x_n = (J|_N)^{-1} g_n \in N \subset \ell^1$，因拓扑同构保持弱收敛与强收敛，故 $\{x_n\}$ 在 $\ell^1$ 中弱收敛于 $0$ 但不强收敛于 $0$，这与 **$\ell^1$ 的 Schur 性质**（教材第 183 页注记 5.1：$\ell^1$ 中弱收敛等价于强收敛）矛盾！因此 $\mathcal{R}(A)$ 不包含无穷维闭子空间。 $\blacksquare$

2. **当 $X = H$ 为 Hilbert 空间时，充分性 $(\Leftarrow)$ 成立（思路简述）**：
   设 $A \in \mathcal{L}(H)$ 且 $A$ 不是紧算子，则正自伴算子 $|A| = (A^*A)^{1/2}$ 也不是紧算子（因 $\|Ax\| = \||A|x\|$）。利用自伴算子的谱分解（第 6 章），存在 $\delta > 0$ 使谱投影对应的闭子空间 $N = E([\delta, \infty))H$ 为无穷维闭子空间，且对一切 $x \in N$ 满足 $\|Ax\| = \||A|x\| \ge \delta \|x\|$。这种有正下界的有界算子把闭子空间 $N$ 映成闭子空间 $M = A(N) \subset \mathcal{R}(A)$，且与 $N$ 同构，故 $M$ 是无穷维闭子空间。逆否即得充分性！

#### 【考点复盘与避坑提示】
1. **期末考场实战策略（怎么写最稳拿满分？）**：
   - 如果考卷直接照搬原题，先简明写一句：“注：当 $A \in \mathcal{C}(X)$ 时，若不要求子空间闭，取 $\ell^2$ 上 $Ae_n = \frac{1}{n}e_n$ 其值域含无穷维非闭子空间 $c_{00}$，故题意应理解为不含无穷维**闭**子空间。”
   - 然后把上面的 **必要性 $(\Rightarrow)$ 完整证明（证法一或证法二）** 漂漂亮亮写出来——这正是出题人核心想考的知识点（也是习题 17(ii)）。
   - 最后用 3–4 行注明：对于充分性 $(\Leftarrow)$，在 Hilbert 空间情形下成立，而在一般 Banach 空间下不成立（简述 $X = \ell^1 \oplus \ell^2, \, A(x, y) = (0, Jx)$ 反例）。这样无论阅卷老师按哪套标准批改，都能稳拿满分！
2. **高频副产品**：证法一（商空间 + 逆算子定理）和证法二（开映射定理）直接秒杀同页的 **习题 17(ii)**（只需把 $A \in \mathcal{C}(X)$ 换成 $A \in \mathcal{C}(X, Y)$，字都不用改！）。

---

### 习题 5 第 14 题（教材第 205 页）

#### 【原题】
> 14. 定义 $C([0,1])$ 上的算子 $A$ 为
> $$Af(x) = \int_0^{1-x} f(t)\mathrm{d}t, \qquad x \in [0,1], \quad f \in C([0,1]).$$
> 证明：$A$ 是紧算子。

#### 【零基础破题思路】
- **第一反应**：空间是 $C([0,1])$，算子是积分算子 $\implies$ 立刻想到 **Arzelà-Ascoli 定理**！
- **千万别踩坑**：能不能直接说“由教材例 5.1 知积分算子是紧算子”？**绝对不行！** 因为如果把积分上限 $1-x$ 强行改写成 $\int_0^1 K(x, t) f(t)\mathrm{d}t$，核函数 $K(x, t) = \begin{cases} 1, & 0 \le t \le 1-x \\ 0, & 1-x < t \le 1 \end{cases}$ 在对角线 $x + t = 1$ 上是**间断的**，不满足例 5.1 要求的“$K$ 连续”！所以必须亲自老老实实用 Arzelà-Ascoli 定理检验。
- **三步走战略**：
  1. 确认 $A \in \mathcal{L}(C([0,1]))$（把连续函数映成连续函数、线性、有界）；
  2. 对任意有界集 $B = \{f \in C([0,1]) : \|f\|_\infty \le M\}$，证 $A(B)$ **一致有界**；
  3. 估计 $|Af(x_1) - Af(x_2)|$，得到 Lipschitz 估计 $\le M|x_1 - x_2|$，立刻推出 **等度连续**！

#### 【规范解答】

记 $C([0,1])$ 上的范数为一致范数 $\|f\|_\infty = \max_{x \in [0,1]} |f(x)|$。

**第一步：验证 $A \in \mathcal{L}(C([0,1]))$。**

- **良定性（值域落在 $C([0,1])$ 中）**：对任意 $f \in C([0,1])$，变上限积分 $g(u) = \int_0^u f(t)\mathrm{d}t$ 在 $[0,1]$ 上连续（事实上连续可微），复合连续函数 $u = 1 - x$ 后，$Af(x) = g(1 - x) \in C([0,1])$。
- **线性性**：对任意 $f, g \in C([0,1])$ 及常数 $\alpha, \beta$，由定积分的线性性有
  $$A(\alpha f + \beta g)(x) = \int_0^{1-x} (\alpha f(t) + \beta g(t))\mathrm{d}t = \alpha Af(x) + \beta Ag(x).$$
- **有界性**：对任意 $f \in C([0,1])$ 和任意 $x \in [0,1]$，有
  $$|Af(x)| = \left|\int_0^{1-x} f(t)\mathrm{d}t\right| \le \int_0^{1-x} |f(t)|\mathrm{d}t \le (1 - x)\|f\|_\infty \le \|f\|_\infty.$$
  对 $x \in [0,1]$ 取最大值，得 $\|Af\|_\infty \le \|f\|_\infty$。故 $A$ 是 $C([0,1])$ 上的有界线性算子（且 $\|A\| \le 1$）。

**第二步：利用 Arzelà-Ascoli 定理证明 $A$ 将有界集映为列紧集。**

任取 $C([0,1])$ 中的有界子集 $B$，则存在常数 $M > 0$，使得对一切 $f \in B$ 均有 $\|f\|_\infty \le M$。下面验证像集 $A(B) = \{Af : f \in B\}$ 在 $[0,1]$ 上一致有界且等度连续：

1. **一致有界性**：对任意 $f \in B$ 及任意 $x \in [0,1]$，由第一步的不等式立得
   $$|Af(x)| \le \|f\|_\infty \le M.$$
   故 $A(B)$ 在 $[0,1]$ 上一致有界。
2. **等度连续性**：任取 $f \in B$ 及 $x_1, x_2 \in [0,1]$。不妨设 $x_1 \le x_2$（若 $x_1 > x_2$，交换二者角色即可），则 $0 \le 1 - x_2 \le 1 - x_1 \le 1$，从而
   $$|Af(x_1) - Af(x_2)| = \left|\int_0^{1-x_1} f(t)\mathrm{d}t - \int_0^{1-x_2} f(t)\mathrm{d}t\right| = \left|\int_{1-x_2}^{1-x_1} f(t)\mathrm{d}t\right| \le \int_{1-x_2}^{1-x_1} |f(t)|\mathrm{d}t \le M \bigl((1 - x_1) - (1 - x_2)\bigr) = M|x_1 - x_2|.$$
   于是对任意 $\varepsilon > 0$，取 $\delta = \frac{\varepsilon}{M} > 0$（注意 $\delta$ 只依赖于 $\varepsilon$ 和 $M$，与 $f \in B$ 以及点 $x_1, x_2$ 均无关）。当 $|x_1 - x_2| < \delta$ 时，对一切 $f \in B$ 统一有
   $$|Af(x_1) - Af(x_2)| \le M|x_1 - x_2| < M \cdot \frac{\varepsilon}{M} = \varepsilon.$$
   故 $A(B)$ 在 $[0,1]$ 上等度连续。

由 **Arzelà-Ascoli 定理**，像集 $A(B)$ 是 $C([0,1])$ 中的列紧集。根据定义 5.1，$A$ 是 $C([0,1])$ 上的紧算子。 $\blacksquare$

#### 【考点复盘与避坑提示】
1. **等度连续的“一致”二字怎么写才不扣分**：一定要明确点出“$\delta = \varepsilon / M$ **与 $f \in B$ 及 $x_1, x_2 \in [0,1]$ 无关**”！如果 $\delta$ 依赖于具体的函数 $f$，那只证明了单个函数一致连续，没证明一族函数等度连续。
2. **同类题通杀**：无论是 Volterra 算子 $(Vf)(x) = \int_0^x f(t)\mathrm{d}t$，还是带连续核的变上限算子 $\int_0^{\varphi(x)} K(x,t)f(t)\mathrm{d}t$（$\varphi$ 连续），统统照搬这个“一致有界 + 等度连续 + Arzelà-Ascoli”模板！
3. **细节规范**：证明 $A$ 是紧算子时，别忘了先花两行交代 $A$ 是 $C([0,1])$ 上的有界线性算子（良定性 + 线性 + 有界）。

---

### 习题 5 第 15 题（教材第 205 页）

#### 【原题】
> 15. 设 $H, G$ 为 Hilbert 空间，$A \in \mathcal{L}(H, G)$. 证明：如果 $A$ 是有限秩算子，则 $A$ 必是 Hilbert-Schmidt 算子.

#### 【零基础破题思路】
- **回忆定义 5.4**：要证 $A$ 是 Hilbert-Schmidt 算子，只需找一组 $H$ 的标准正交基 $\{e_\alpha\}$（可分情形即 $\{e_n\}_{n=1}^\infty$），验证 $\sum_\alpha \|Ae_\alpha\|^2 < +\infty$ 即可！（且由教材式 (5.8)，只要对**某一组**标准正交基成立，对所有标准正交基都成立且和相等！）
- **灵光一现（证法一：挑一组最聪明的基！）**：
  - $A$ 是有限秩算子，意思是它的像 $\mathcal{R}(A)$ 只有有限维（设为 $m$ 维）——这意味着空间 $H$ 里绝大部分方向都被 $A$ 拍成了 $0$（落进了 $\operatorname{Ker} A$）！
  - 用 Hilbert 空间的正交分解：$H = \operatorname{Ker} A \oplus (\operatorname{Ker} A)^\perp$。
  - 在 $(\operatorname{Ker} A)^\perp$ 上 $A$ 是单射且映满 $\mathcal{R}(A)$，所以 $(\operatorname{Ker} A)^\perp$ 的维数也是有限的 $m$ 维！
  - 于是我们把 $(\operatorname{Ker} A)^\perp$ 的 $m$ 个标准正交基向量，与 $\operatorname{Ker} A$ 的标准正交基拼成 $H$ 的标准正交基——算 $\sum \|Ae_\alpha\|^2$ 时，落在 $\operatorname{Ker} A$ 里的基向量全被 $A$ 映成 $0$，无穷级数瞬间坍缩成只有 $m$ 项的有限和，当然有限！
- **另一条路（证法二：对任意标准正交基直接算，用到 $A^*$ 和 Parseval！）**：
  - 在有限维空间 $\mathcal{R}(A)$ 里取标准正交基 $\{f_1, \dots, f_m\}$，把 $Ax$ 展开成 $\sum_{j=1}^m (Ax, f_j) f_j = \sum_{j=1}^m (x, A^* f_j) f_j$，再用 Parseval 等式把对基向量 $e_\alpha$ 的无穷和一键收成 $\sum_{j=1}^m \|A^* f_j\|^2 < \infty$！

#### 【规范解答】

**证法一（正交分解 + 构造标准正交基，最直观简洁）：**

因为 $A \in \mathcal{L}(H, G)$ 是有限秩算子，故 $m \triangleq \dim \mathcal{R}(A) < +\infty$。若 $m = 0$，则 $A = 0$，显然 $\|A\|_2 = 0 < +\infty$，结论成立。下设 $m \ge 1$。

1. **正交分解与维数确定**：由于 $A$ 连续，其核空间 $\operatorname{Ker} A$ 是 $H$ 的闭线性子空间。由 Hilbert 空间的正交分解定理，有
   $$H = \operatorname{Ker} A \oplus (\operatorname{Ker} A)^\perp.$$
   考察限制映射 $A_1 = A|_{(\operatorname{Ker} A)^\perp}: (\operatorname{Ker} A)^\perp \to \mathcal{R}(A)$：
   - **单射**：若 $x \in (\operatorname{Ker} A)^\perp$ 且 $Ax = 0$，则 $x \in \operatorname{Ker} A \cap (\operatorname{Ker} A)^\perp = \{0\}$；
   - **满射**：任取 $y = Ax \in \mathcal{R}(A)$，将 $x$ 正交分解为 $x = x_1 + x_2$（$x_1 \in \operatorname{Ker} A, \, x_2 \in (\operatorname{Ker} A)^\perp$），则 $A_1 x_2 = A x_2 = A(x_1 + x_2) = Ax = y$。

   因此 $A_1$ 是线性同构，从而
   $$\dim (\operatorname{Ker} A)^\perp = \dim \mathcal{R}(A) = m < +\infty.$$

2. **构造 $H$ 的标准正交基**：在有限维子空间 $(\operatorname{Ker} A)^\perp$ 中取一组标准正交基 $\{e_1, e_2, \dots, e_m\}$；在闭子空间 $\operatorname{Ker} A$ 中任取一组标准正交基 $\{u_\beta\}_{\beta \in \Lambda}$（按教材 5.2 节可分假设，可编号为 $\{e_{m+1}, e_{m+2}, \dots\}$，有限或可数）。
   因为 $(\operatorname{Ker} A)^\perp \perp \operatorname{Ker} A$ 且 $H = \operatorname{Ker} A \oplus (\operatorname{Ker} A)^\perp$，二者合并 $\{e_1, \dots, e_m\} \cup \{u_\beta\}_{\beta \in \Lambda}$ 构成 $H$ 的一组标准正交基。

3. **计算 Hilbert-Schmidt 范数**：由于每个 $u_\beta \in \operatorname{Ker} A$ 均满足 $A u_\beta = 0$，按定义 5.4 有
   $$\|A\|_2^2 = \sum_{j=1}^m \|Ae_j\|^2 + \sum_{\beta \in \Lambda} \|Au_\beta\|^2 = \sum_{j=1}^m \|Ae_j\|^2 \le \sum_{j=1}^m \|A\|^2 \|e_j\|^2 = m \|A\|^2 < +\infty.$$
   因此 $\|A\|_2 \le \sqrt{m}\,\|A\| < +\infty$，即 $A$ 是 Hilbert-Schmidt 算子。 $\blacksquare$

---

**证法二（任意标准正交基 + 共轭算子与 Parseval 等式）：**

设 $m = \dim \mathcal{R}(A) < +\infty$（$m = 0$ 平凡，下设 $m \ge 1$），在有限维子空间 $\mathcal{R}(A) \subset G$ 中取标准正交基 $\{f_1, f_2, \dots, f_m\}$。对任意 $x \in H$，因为 $Ax \in \mathcal{R}(A)$，由有限维正交展开与共轭算子 $A^* \in \mathcal{L}(G, H)$ 的定义得
$$Ax = \sum_{j=1}^m (Ax, f_j)_G f_j = \sum_{j=1}^m (x, A^* f_j)_H f_j.$$
由勾股定理，对一切 $x \in H$ 有
$$\|Ax\|^2 = \sum_{j=1}^m |(x, A^* f_j)_H|^2.$$
任取 $H$ 的一组标准正交基 $\{e_\alpha\}_{\alpha \in \Gamma}$（可分情形即为 $\{e_n\}_{n=1}^\infty$），对所有基向量求和并交换有限和与级数顺序，再利用 **Parseval 等式**：
$$\sum_{\alpha \in \Gamma} \|Ae_\alpha\|^2 = \sum_{\alpha \in \Gamma} \sum_{j=1}^m |(e_\alpha, A^* f_j)_H|^2 = \sum_{j=1}^m \sum_{\alpha \in \Gamma} |(e_\alpha, A^* f_j)_H|^2 = \sum_{j=1}^m \|A^* f_j\|^2 \le m\|A^*\|^2 = m\|A\|^2 < +\infty.$$
由定义 5.4，故 $A$ 是 Hilbert-Schmidt 算子。 $\blacksquare$

#### 【考点复盘与避坑提示】
1. **副产品不等式**：两种证法都顺带给出了非常实用的估计式 $\|A\|_2 \le \sqrt{\operatorname{rank}(A)}\,\|A\|$（和有限维矩阵的 $\|A\|_F \le \sqrt{\operatorname{rank}(A)}\,\|A\|_2$ 完全一致！）。
2. **为何不会受“可分/不可分”限制**：题目没明写“可分”，上面两种证法使用指标集 $\Lambda / \Gamma$（顺带括号注明可分情形就是 $\{e_n\}_{n=1}^\infty$），无论是按教材可分设定还是对一般 Hilbert 空间，都完全无懈可击！
3. **常见误区**：别只说“把 $\mathcal{R}(A)$ 的基扩充一下”，切记定义 5.4 求和的是**定义域 $H$** 的标准正交基，而不是值域 $G$ 的标准正交基！

---

### 习题 5 第 20 题（教材第 205 页）

#### 【原题】
> 20. 设 $X, Y$ 为 Banach 空间，$B \in \mathcal{F}(X, Y)$，$A \in \mathcal{C}(X, Y)$. 证明：
> (i) $B + A \in \mathcal{F}(X, Y)$；
> (ii) $\operatorname{ind}(T + A) = \operatorname{ind} T$.

> **⚠️ 教材印刷笔误勘误**：原题第 (ii) 问中的 $T$ 系原书印刷笔误（把题干里的 $B$ 误印成了 $T$），应为证明 **$\operatorname{ind}(B + A) = \operatorname{ind} B$**。

#### 【零基础破题思路】
- **第 (i) 问（紧扰动不破坏 Fredholm 性）**：
  - 千万不要傻傻地去按定义 5.5 算 $\operatorname{Ker}(B + A)$ 和 $\mathcal{R}(B + A)$——加法之后核和值域都会彻底乱掉！
  - 记住金句：**“看到证明某算子属于 $\mathcal{F}(X, Y)$，首选 Atkinson 定理（定理 5.13）！”**
  - $B$ 是 Fredholm 算子 $\implies$ 存在拟逆 $S$ 使 $SB = I_X - K_1, \, BS = I_Y - K_2$。
  - 把同一个 $S$ 乘到 $B + A$ 上：$S(B + A) = I_X - (K_1 - SA)$，$(B + A)S = I_Y - (K_2 - AS)$。
  - 因为紧算子是双边理想（定理 5.3），$SA$ 和 $AS$ 还是紧的，紧减紧还是紧！于是一步秒杀 $B + A \in \mathcal{F}(X, Y)$！
- **第 (ii) 问（紧扰动不改变指标）**：
  - **神仙 3 行代数法（最推荐！）**：观察刚刚算出的两条式子——**同一个 $S$ 同时充当了 $B$ 和 $B + A$ 的左拟逆！**
    - 由注记 5.4，$S$ 自己也是 Fredholm 算子；
    - 由乘积指标公式（定理 5.14）与恒等减紧算子指标为 $0$（例 5.8 / 定理 5.8(4)）：
      $$\operatorname{ind} S + \operatorname{ind} B = \operatorname{ind}(SB) = \operatorname{ind}(I_X - K_1) = 0,$$
      $$\operatorname{ind} S + \operatorname{ind}(B + A) = \operatorname{ind}(S(B + A)) = \operatorname{ind}(I_X - (K_1 - SA)) = 0.$$
    - 两式相减，$\operatorname{ind} S$ 直接抵消，立得 $\operatorname{ind}(B + A) = \operatorname{ind} B$！（和教材第 203 页证明定理 5.15 的手法一模一样！）
  - **经典同伦法（附注备用）**：也可以构造连续路径 $F(t) = B + tA$（$t \in [0, 1]$），由 (i) 知整条路径都在 $\mathcal{F}(X, Y)$ 里，由定理 5.15 知整数值函数 $t \mapsto \operatorname{ind}(B + tA)$ 在连通区间 $[0, 1]$ 上局部常值（连续），故必为常数，从而 $t=1$ 与 $t=0$ 处相等！

#### 【规范解答】

*(注：第 (ii) 问中 $T$ 应为题干中的 $B$，以下证 $\operatorname{ind}(B + A) = \operatorname{ind} B$。)*

**(i) 证明 $B + A \in \mathcal{F}(X, Y)$：**

因为 $B \in \mathcal{F}(X, Y)$，由 **Atkinson 定理（定理 5.13(1)）**，存在算子 $S \in \mathcal{L}(Y, X)$ 以及紧算子 $K_1 \in \mathcal{C}(X), \, K_2 \in \mathcal{C}(Y)$，使得
$$SB = I_X - K_1, \qquad BS = I_Y - K_2.$$

将 $S$ 分别左乘、右乘 $B + A$，得
$$S(B + A) = SB + SA = I_X - (K_1 - SA),$$
$$(B + A)S = BS + AS = I_Y - (K_2 - AS).$$

由于 $A \in \mathcal{C}(X, Y)$，$S \in \mathcal{L}(Y, X)$，由紧算子的双边理想性质及线性性（定理 5.3），有 $SA \in \mathcal{C}(X)$、$AS \in \mathcal{C}(Y)$，从而
$$\tilde{K}_1 \triangleq K_1 - SA \in \mathcal{C}(X), \qquad \tilde{K}_2 \triangleq K_2 - AS \in \mathcal{C}(Y).$$

再次应用 **Atkinson 定理（定理 5.13(2)）**，立即得到 $B + A \in \mathcal{F}(X, Y)$。

**(ii) 证明 $\operatorname{ind}(B + A) = \operatorname{ind} B$：**

**证法一（Atkinson 定理 + 乘积指标公式，最简洁推荐）：**

承接 (i) 中的算子 $S \in \mathcal{L}(Y, X)$。由于 $BS = I_Y - K_2$ 且 $SB = I_X - K_1$，将 $B$ 视为 $S$ 的正则化子，由定理 5.13(2)（即注记 5.4）可知 $S \in \mathcal{F}(Y, X)$。

另一方面，由定理 5.6 与定理 5.8(4)（即教材例 5.8），对任意紧算子 $K \in \mathcal{C}(X)$，均有 $I_X - K \in \mathcal{F}(X)$ 且 $\operatorname{ind}(I_X - K) = 0$。

分别对乘积 $SB$ 与 $S(B + A)$ 应用 **Fredholm 算子乘积指标公式（定理 5.14）**，可得：、
$$\operatorname{ind} S + \operatorname{ind} B = \operatorname{ind}(SB) = \operatorname{ind}(I_X - K_1) = 0,$$
$$\operatorname{ind} S + \operatorname{ind}(B + A) = \operatorname{ind}\bigl(S(B + A)\bigr) = \operatorname{ind}(I_X - \tilde{K}_1) = 0.$$

两式相减（或直接由两式均给出 $-\operatorname{ind} S$），即得
$$\operatorname{ind}(B + A) = \operatorname{ind} B. \qquad \blacksquare$$

---

**证法二（同伦路径 + 小扰动稳定性，几何直觉法）：**

对每个 $t \in [0, 1]$，令 $F(t) = B + tA$。因为 $tA \in \mathcal{C}(X, Y)$，由 (i) 知 $F(t) \in \mathcal{F}(X, Y)$。定义整数值函数
$$\varphi(t) = \operatorname{ind}(B + tA), \qquad t \in [0, 1].$$

任取 $t_0 \in [0, 1]$。由 **定理 5.15（小扰动指标不变性）**，存在 $\varepsilon > 0$，使得对任意满足 $\|S\| < \varepsilon$ 的 $S \in \mathcal{L}(X, Y)$，都有 $\operatorname{ind}(F(t_0) + S) = \operatorname{ind} F(t_0)$。取 $\delta > 0$ 满足 $\delta \|A\| < \varepsilon$（若 $A = 0$ 则结论平凡，下设 $\|A\| > 0$，取 $\delta = \frac{\varepsilon}{\|A\|}$）。于是当 $t \in [0, 1]$ 且 $|t - t_0| < \delta$ 时，有
$$\|F(t) - F(t_0)\| = |t - t_0| \cdot \|A\| < \varepsilon,$$
从而 $\varphi(t) = \operatorname{ind}\bigl(F(t_0) + (F(t) - F(t_0))\bigr) = \operatorname{ind} F(t_0) = \varphi(t_0)$。

这表明 $\varphi: [0, 1] \to \mathbb{Z}$ 是连通闭区间 $[0, 1]$ 上的连续（局部常值）整数值函数，因此 $\varphi$ 在 $[0, 1]$ 上恒为常数。特别地，
$$\operatorname{ind}(B + A) = \varphi(1) = \varphi(0) = \operatorname{ind} B. \qquad \blacksquare$$

#### 【考点复盘与避坑提示】
1. **两个超级高频定理的黄金组合**：证明指标相等时，**“找公共正则化子 $S$ + 定理 5.14 + $\operatorname{ind}(I - K) = 0$”** 几乎是最强套路！不仅秒杀本题，也是教材证明定理 5.15 的核心。
2. **别忘交代 $S \in \mathcal{F}(Y, X)$**：在证法一中，使用定理 5.14（$\operatorname{ind}(SB) = \operatorname{ind} S + \operatorname{ind} B$）的前提是参与相乘的两个算子**都**必须是 Fredholm 算子！所以务必明确写上一句“由注记 5.4 / 定理 5.13(2) 知 $S \in \mathcal{F}(Y, X)$”，避免逻辑跳步被扣细节分。
3. **直觉升华**：定理 5.15 说“范数很小的任意扰动不改指标”，本题说“任意大的紧扰动不改指标”——二者合起来表明，Fredholm 指标是一个极其稳健的**拓扑不变量**！

---

### 习题 5 第 22 题（教材第 205 页）

#### 【原题】
> 22. 设 $X, Y, Z$ 均为 Banach 空间，$A \in \mathcal{L}(X, Y)$，$B \in \mathcal{L}(Y, Z)$，$BA \in \mathcal{F}(X, Z)$ 且 $\dim \operatorname{Ker} B < +\infty$. 证明：$A \in \mathcal{F}(X, Y)$ 且 $B \in \mathcal{F}(Y, Z)$.

#### 【零基础破题思路】
- **先看条件给谁铺好了路**：
  - 关于 $B$：题目直接给了 $\dim \operatorname{Ker} B < \infty$（条件 2 到手！）；而 $\mathcal{R}(BA) \subset \mathcal{R}(B) \subset Z$，既然小的 $\mathcal{R}(BA)$ 在 $Z$ 中余维有限，那么包着它的更大的 $\mathcal{R}(B)$ 余维只会更小，必有限（条件 3 到手！）。
  - **关键卡点（绝对不能漏！）**：还要证明 **$\mathcal{R}(B)$ 是闭子空间**（条件 1）！怎么证？因为 $\mathcal{R}(BA)$ 是闭子空间且余维有限，所以 $\mathcal{R}(B)$ 只比 $\mathcal{R}(BA)$ 多了一个**有限维子空间 $W$**，即 $\mathcal{R}(B) = \mathcal{R}(BA) + W$！利用经典小引理“**闭子空间 + 有限维子空间必闭**”（用商映射 2 行就能证完！），立刻补齐最后一块拼图，证得 **$B \in \mathcal{F}(Y, Z)$**！
- **再攻克 $A \in \mathcal{F}(X, Y)$**：
  - 手里现在握着两个 Fredholm 算子：$BA \in \mathcal{F}(X, Z)$ 和刚证出的 $B \in \mathcal{F}(Y, Z)$！
  - 想证 $A \in \mathcal{F}(X, Y)$，果断请出 **Atkinson 定理（定理 5.13）**：
    - $BA$ 有双边拟逆 $S_1$：$S_1 BA = I_X - K_1$，$BA S_1 = I_Z - K_2$；
    - $B$ 有左拟逆 $S_2$：$S_2 B = I_Y - K_3$。
  - 怎么凑出 $A$ 的拟逆？看第一个式子 $(S_1 B) A = I_X - K_1$，很自然猜 $R = S_1 B$ 就是 $A$ 的拟逆！左乘已经搞定，那右乘 $AR = A S_1 B$ 怎么办？
  - **神奇恒等式（乘法结合律）**：给 $AR$ 左边乘上 $S_2 B = I_Y - K_3$：
    - 从左边结合：$(S_2 B) AR = (I_Y - K_3) AR = AR - K_3 AR$；
    - 从中间结合：$S_2 (B A S_1) B = S_2 (I_Z - K_2) B = S_2 B - S_2 K_2 B = I_Y - K_3 - S_2 K_2 B$！
    - 两边一相等，移项立马出现 $AR = I_Y - (\text{一堆紧算子})$！由 Atkinson 定理秒杀 $A \in \mathcal{F}(X, Y)$！

#### 【规范解答】

**预备小引理（闭子空间与有限维子空间之和必闭）：**
> 设 $U$ 是 Banach 空间 $Z$ 的闭线性子空间，$W \subset Z$ 是有限维线性子空间，则 $U + W$ 是 $Z$ 的闭线性子空间。
>
> *简证*：考虑连续线性商映射 $q: Z \to Z/U$（$q(z) = z + U$）。因 $U$ 闭，商空间 $Z/U$ 赋予商范数后是赋范线性空间。像集 $q(W)$ 是 $Z/U$ 的有限维线性子空间，而赋范线性空间中的有限维子空间必闭，故 $q(W)$ 在 $Z/U$ 中闭。注意到 $U + W = q^{-1}(q(W))$，由 $q$ 的连续性可知 $U + W$ 是 $Z$ 的闭线性子空间。

---

**第一步：证明 $B \in \mathcal{F}(Y, Z)$。**

按照 Fredholm 算子的定义（定义 5.5），逐一验证以下三条：

1. **$\dim \operatorname{Ker} B < +\infty$**：题设直接给出。
2. **$\operatorname{codim} \mathcal{R}(B) < +\infty$**：
   对任意 $x \in X$，$BAx = B(Ax) \in \mathcal{R}(B)$，故 $\mathcal{R}(BA) \subset \mathcal{R}(B) \subset Z$。定义商空间之间的自然线性映射
   $$\pi: Z / \mathcal{R}(BA) \to Z / \mathcal{R}(B), \qquad \pi\bigl(z + \mathcal{R}(BA)\bigr) = z + \mathcal{R}(B).$$
   若 $z_1 - z_2 \in \mathcal{R}(BA)$，则 $z_1 - z_2 \in \mathcal{R}(B)$，故 $\pi$ 良定，且显然是线性满射。因为 $BA \in \mathcal{F}(X, Z)$，有 $\dim\bigl(Z / \mathcal{R}(BA)\bigr) < +\infty$，从而
   $$\operatorname{codim} \mathcal{R}(B) = \dim\bigl(Z / \mathcal{R}(B)\bigr) \le \dim\bigl(Z / \mathcal{R}(BA)\bigr) = \operatorname{codim} \mathcal{R}(BA) < +\infty.$$
3. **$\mathcal{R}(B)$ 是 $Z$ 的闭子空间**：
   因为 $\mathcal{R}(BA) \subset \mathcal{R}(B) \subset Z$，商空间 $\mathcal{R}(B) / \mathcal{R}(BA)$ 是有限维空间 $Z / \mathcal{R}(BA)$ 的线性子空间，故其维数 $k \triangleq \dim\bigl(\mathcal{R}(B) / \mathcal{R}(BA)\bigr) < +\infty$。
   选取 $w_1, \dots, w_k \in \mathcal{R}(B)$ 使它们的陪集构成 $\mathcal{R}(B) / \mathcal{R}(BA)$ 的一组基，并记有限维子空间 $W = \operatorname{span}\{w_1, \dots, w_k\} \subset \mathcal{R}(B)$（若 $k = 0$ 则取 $W = \{0\}$），于是有
   $$\mathcal{R}(B) = \mathcal{R}(BA) + W.$$
   由 $BA \in \mathcal{F}(X, Z)$ 知 $\mathcal{R}(BA)$ 是 $Z$ 的闭子空间。应用上述预备小引理，立即得出 $\mathcal{R}(B)$ 是 $Z$ 的闭子空间。

综合 (1)(2)(3)，按定义 5.5 已证得 **$B \in \mathcal{F}(Y, Z)$**。

---

**第二步：证明 $A \in \mathcal{F}(X, Y)$。**

因为 $BA \in \mathcal{F}(X, Z)$ 且第一步已证 $B \in \mathcal{F}(Y, Z)$，分别应用 **Atkinson 定理（定理 5.13(1)）**：
- 存在 $S_1 \in \mathcal{L}(Z, X)$ 及紧算子 $K_1 \in \mathcal{C}(X), \, K_2 \in \mathcal{C}(Z)$，使得
  $$S_1 (BA) = I_X - K_1, \qquad (BA) S_1 = I_Z - K_2;$$
- 存在 $S_2 \in \mathcal{L}(Z, Y)$ 及紧算子 $K_3 \in \mathcal{C}(Y)$，使得
  $$S_2 B = I_Y - K_3.$$

令 $R \triangleq S_1 B \in \mathcal{L}(Y, X)$。下面验证 $R$ 满足定理 5.13(2) 的条件：

1. **左乘检验**：
   $$RA = (S_1 B) A = S_1 (BA) = I_X - K_1, \qquad K_1 \in \mathcal{C}(X).$$
2. **右乘检验**：用两种结合方式计算乘积 $S_2 B A S_1 B \in \mathcal{L}(Y)$：
   - 一方面，结合中间的 $BA S_1 = I_Z - K_2$，得
     $$S_2 (BA S_1) B = S_2 (I_Z - K_2) B = S_2 B - S_2 K_2 B = I_Y - K_3 - S_2 K_2 B;$$
   - 另一方面，结合左端的 $S_2 B = I_Y - K_3$ 以及 $S_1 B = R$，得
     $$(S_2 B) A (S_1 B) = (I_Y - K_3) AR = AR - K_3 AR.$$

   令两式相等并移项整理，得
   $$AR = I_Y - \bigl(K_3 + S_2 K_2 B - K_3 AR\bigr).$$
   由紧算子的双边理想性质与线性性（定理 5.3）：$K_3 \in \mathcal{C}(Y)$，$S_2 K_2 B \in \mathcal{C}(Y)$，$K_3 AR \in \mathcal{C}(Y)$，故
   $$\tilde{K} \triangleq K_3 + S_2 K_2 B - K_3 AR \in \mathcal{C}(Y),$$
   即有 $AR = I_Y - \tilde{K}$。

综上，$RA = I_X - K_1$ 且 $AR = I_Y - \tilde{K}$（其中 $K_1 \in \mathcal{C}(X), \, \tilde{K} \in \mathcal{C}(Y)$）。由 **Atkinson 定理（定理 5.13(2)）**，立即得出 **$A \in \mathcal{F}(X, Y)$**。 $\blacksquare$

#### 【考点复盘与避坑提示】
1. **最大的扣分黑洞**：第一步证 $B \in \mathcal{F}(Y, Z)$ 时，很多同学只写了“$\dim \operatorname{Ker} B < \infty$ 且 $\operatorname{codim} \mathcal{R}(B) \le \operatorname{codim} \mathcal{R}(BA) < \infty$”就宣布完毕，**漏证 $\mathcal{R}(B)$ 是闭子空间**！请一定把“$\mathcal{R}(B) = \mathcal{R}(BA) + W$ + 预备小引理”完整写上（教材第 201 页定理 5.13(2) 证明中正使用了完全相同的步骤）。
2. **为什么条件 $\dim \operatorname{Ker} B < +\infty$ 缺一不可？（反例直觉）**：
   如果去掉这个条件，结论不成立！反例：取 $X = Z = \ell^2, \, Y = \ell^2 \oplus \ell^2$，令 $Ax = (x, 0)$，$B(x, y) = x$。此时 $BA = I_{\ell^2} \in \mathcal{F}(\ell^2)$，但 $\operatorname{Ker} B = \{0\} \oplus \ell^2$ 是无穷维的（$B \notin \mathcal{F}$），同时 $\mathcal{R}(A) = \ell^2 \oplus \{0\}$ 在 $Y$ 中余维无穷（$A \notin \mathcal{F}$）。
3. **延伸考点**：本题条件可对称替换为“$BA \in \mathcal{F}(X, Z)$ 且 $\operatorname{codim}\mathcal{R}(A) < +\infty$”，同样能推出 $A \in \mathcal{F}(X, Y)$ 且 $B \in \mathcal{F}(Y, Z)$（证明只需先用 $\operatorname{Ker} A \subset \operatorname{Ker}(BA)$ 与 $\operatorname{codim}\mathcal{R}(A)<\infty$ 证 $A \in \mathcal{F}$，再如法炮制用 Atkinson 定理证 $B \in \mathcal{F}$）。

---

## 三、考前 3 分钟速背清单（全章定理与题型总览）

| 考场遇到什么问题 | 立刻条件反射的工具 | 对应教材与作业 |
| :--- | :--- | :--- |
| 证无穷维紧算子**不可逆** / 某算子**不紧** | 反证法凑出 $I_X \in \mathcal{C}(X)$ + **Riesz 引理** | 定理 5.5、例 5.2、**第 5 题** |
| 证紧算子值域**不含无穷维闭子空间** | 收缩到 $X_0 = A^{-1}(M)$ + **开映射定理**（或商空间 + 逆算子定理）+ Riesz 引理 | **第 6 题**、习题 17(ii) |
| 证 $C([0,1])$ 上积分算子是**紧算子** | **Arzelà-Ascoli 定理**（一致有界 + 等度连续，强调 $\delta$ 与 $f$ 无关） | 例 5.1、**第 14 题** |
| 证 $\ell^p / L^2$ 算子是**紧算子** | 有限秩截断逼近（定理 5.4）或 **Hilbert-Schmidt 判定**（定理 5.12(5)） | 定理 5.4、定理 5.9、定理 5.12 |
| 证有限秩算子是 **Hilbert-Schmidt 算子** | $H = \operatorname{Ker} A \oplus (\operatorname{Ker} A)^\perp$ 凑基，或 Parseval 等式 + $A^*$ | 定义 5.4、**第 15 题** |
| 证算子是 **Fredholm 算子** | 优先 **Atkinson 定理**（定理 5.13：凑 $SA = I - K_1, \, AS = I - K_2$）；按定义证时**切勿漏证值域闭**！ | 定义 5.5、定理 5.13、**第 20/22 题** |
| 算/证 **Fredholm 指标**相等 | **乘积指标公式**（定理 5.14）+ **$\operatorname{ind}(I - K) = 0$**（定理 5.8(4)）+ 公共拟逆 $S$ | 定理 5.14–5.15、**第 20 题** |
