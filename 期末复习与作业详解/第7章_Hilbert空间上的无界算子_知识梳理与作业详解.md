# 第 7 章 Hilbert 空间上的无界算子：零基础知识梳理与作业详解

> **使用说明**：本资料面向**零泛函分析基础**的本科生期末急救与刷题，内容、记号与定理编号严格对齐教材第 7 章（p. 272–294）。所有作业解答均自包含、步骤严谨、零跳步，可直接作为期末考试规范作答模板。

---

## 零基础导读：3 分钟看懂第 7 章到底在学什么

前六章里，我们一直生活在“温室”里：算子 $A \in \mathcal{L}(H)$ 都是**有界**的，而且**定义在整个空间 $H$ 上**。但只要走出纯数学温室、踏进**微分方程**和**量子力学**的真实世界，你遇到的核心算子几乎全都是**无界算子**——它们根本躲不掉！

先看两个让你一秒恍然大悟的“主角”：

1. **求导算子（微分算子 / 量子力学中的动量算子）**：$u(x) \mapsto u'(x)$（或 $-\mathrm{i} u'(x)$）。
   - 取一串振幅始终为 $1$ 的高频正弦波 $u_n(x) = \sin(nx)$（$x \in [0, 2\pi]$），它们的“大小”一致固定（$\|u_n\|_{L^2} = \sqrt{\pi}$），但求导后 $u_n'(x) = n\cos(nx)$，大小被放大了 **$n$ 倍**！随着频率 $n \to \infty$，放大倍数冲向 $+\infty$！
   - 而且，不是每个 $L^2$ 函数（甚至不是每个连续函数）都可导——所以求导算子**不可能定义在全空间上**！
2. **位置算子（乘法算子，教材例 7.1）**：在 $L^2(\mathbb{R})$ 上定义 $(Af)(t) = t f(t)$。
   - 把一个宽度为 $1$、高度为 $1$ 的方波脉冲 $f_n(t) = \chi_{[n, n+1]}(t)$ 往远处平移到 $[n, n+1]$，脉冲本身能量不变（$\|f_n\|_2 = 1$），但乘上 $t$ 之后，幅度暴涨到至少 $n$ 倍（$\|Af_n\|_2 \ge n \to +\infty$）！
   - 而且，像 $f(t) = \frac{1}{1+|t|}$ 这样的函数明明属于 $L^2(\mathbb{R})$，乘上 $t$ 变成 $\frac{t}{1+|t|}$ 却**掉出了 $L^2(\mathbb{R})$**！所以 $A$ 的定义域最大也只能取 $\mathcal{D}(A) = \{f \in L^2(\mathbb{R}) \mid \int_{\mathbb{R}} t^2|f(t)|^2\,\mathrm{d}t < \infty\}$，绝不可能是整个 $L^2(\mathbb{R})$！

进入无界算子的世界后，有三件颠覆直觉、也是期末最爱考的大事：

1. **定义域 $\mathcal{D}(A)$ 是算子的灵魂（不再只是背景板）**：同样的运算公式（比如 $f \mapsto \mathrm{i} f'$），给它配不同的边界条件（即不同的定义域），它可能**连对称都不是**，也可能**仅仅对称但不自伴**，也可能是完美的**自伴算子**！在无界算子里，**必须把一对 $(A, \mathcal{D}(A))$ 当成一个整体**。
2. **“闭算子”接棒“连续算子”**：无界算子都不连续，直接取极限 $\lim Ax_n = A(\lim x_n)$ 会失效；于是数学家退一步，寻找“如果 $x_n \to x$ **并且** $Ax_n$ 自己也收敛到 $y$，那 $x$ 必在定义域内且 $Ax=y$”的**闭算子**——幸运的是，所有共轭算子 $A^*$（特别是自伴算子）自动都是闭算子！
3. **自伴算子依然能像矩阵一样做谱分解**：虽然无界自伴算子的谱不再有界（比如位置算子的谱是整条实轴 $\mathbb{R}$），但利用神奇的 **Cayley 变换** $U = (A-\mathrm{i}I)(A+\mathrm{i}I)^{-1}$，我们可以把整条无界实轴“压扁卷成”单位圆周，把无界自伴算子 $A$ 变成第 6 章已经解决的**有界酉算子 $U$**，借它的谱分解反推回来，就得到了无界自伴算子的谱分解定理 $A = \int_{\mathbb{R}} \lambda\,\mathrm{d}E_\lambda$！

### 有界算子（第 6 章）vs 无界算子（第 7 章）全景对照表

| 对比维度 | 有界线性算子（第 6 章） | 无界线性算子（第 7 章） |
| :--- | :--- | :--- |
| **定义域 $\mathcal{D}(A)$** | 默认是全空间 $H$ | 通常是真子空间，基本要求在 $H$ 中**稠密**（稠定算子） |
| **算子相等 $A = B$** | 只需对所有 $x \in H$ 有 $Ax = Bx$ | 必须同时满足 **$\mathcal{D}(A) = \mathcal{D}(B)$** 且作用公式相同 |
| **共轭算子 $A^*$** | 永远存在、唯一，且 $\mathcal{D}(A^*)=H$ | 必须 $A$ **稠定**才能唯一定义 $A^*$；$\mathcal{D}(A^*)$ 甚至可能只有 $\{0\}$ |
| **对称 vs 自伴** | 全空间对称自动有界且自伴（二者无区别！） | **天壤之别**：对称仅要求 $A \subset A^*$；自伴要求 $A = A^*$（含定义域相等！） |
| **连续性 / 闭性** | 有界 $\iff$ 连续 $\implies$ 闭算子 | 不连续，但核心研究对象（$A^*$、自伴算子）恒为**闭算子** |
| **二阶共轭 $A^{**}$** | 恒有 $A^{**} = A$ | 需 $\mathcal{D}(A^*)$ 稠密才存在；$A$ 稠定闭 $\iff A = A^{**}$ |
| **自伴算子的谱 $\sigma(A)$** | 实轴上的非空**有界闭集**（紧集） | 实轴上的闭集 $\sigma(A) \subset \mathbb{R}$，但**可以无界**（如 $\mathbb{R}$），且无剩余谱 |
| **谱分解定理** | $A = \int_{\sigma(A)} \lambda\,\mathrm{d}E_\lambda$（有界积分） | $A = \int_{\mathbb{R}} \lambda\,\mathrm{d}E_\lambda$，且定义域由 $\int_{\mathbb{R}}\lambda^2\,\mathrm{d}\lVert E_\lambda x\rVert^2 < \infty$ 刻画 |

### 第 7 章 4 道课后作业题考点快速索引

| 题号 | 教材页码 | 核心考点 | 所用关键定理 / 核心技巧 |
| :---: | :---: | :--- | :--- |
| **第 3 题** | p. 293 | Legendre 微分算子的全部特征值与特征向量 | 多项式空间 $\mathcal{P}_n$ 上三角矩阵 + Rodrigues 公式 + 分部积分对称性 + Weierstrass 多项式逼近 |
| **第 4 题** | p. 293 | $L^2(\mathbb{R})$ 上乘法（位置）算子的无界性与闭性 | 示性函数 $\chi_{[n,n+1]}$ 证无界 + 两次抽几乎处处收敛子列证闭性 |
| **第 6 题** | p. 293 | 对称算子的闭性 / 自伴性与 $A, A^{**}, A^*$ 的关系 | 前置验证 $\mathcal{D}(A^*)$ 稠密 + 定理 7.1（共轭必闭）+ 定理 7.2（稠定闭 $\implies A=A^{**}$） |
| **第 7 题** | p. 294 | 非负对称算子（$(Ax,x) \ge 0$）的范数估计与闭性充要条件 | 内积展开得双控不等式 + Cauchy 列与 $H$ 完备性 + 极限唯一性 |

---

## 第一部分：第 7 章核心知识零基础精讲

> **全章基础约定**：如无特别说明，本章中 $H$ 始终表示**复数域 $\mathbb{C}$ 上的 Hilbert 空间**（完备内积空间），内积 $(\cdot, \cdot)$ 对第一变元线性、对第二变元共轭线性（即 $(\alpha x, \beta y) = \alpha \bar{\beta}(x, y)$）。一个无界线性算子永远由两样东西共同决定：**线性子空间 $\mathcal{D}(A) \subset H$（定义域）** 与 **线性映射 $A: \mathcal{D}(A) \to H$**。

---

### 7.1.1 算子的定义域、扩张与稠定算子的共轭算子（p. 272–275）

#### 1. 带定义域的算子运算与“扩张”（定义 7.1，p. 272–273）

设 $A, B$ 是分别以线性子空间 $\mathcal{D}(A), \mathcal{D}(B) \subset H$ 为定义域的线性算子，$\alpha \in \mathbb{F}$（$\mathbb{R}$ 或 $\mathbb{C}$）：

- **(i) 数乘 $\alpha A$**：定义域为 $\mathcal{D}(A)$，作用 $x \mapsto \alpha Ax$。
  - **零基础必踩避坑点**：$0A$（数 $0$ 乘算子 $A$）的定义域仍是 $\mathcal{D}(A)$，而零算子 $0$ 的定义域是全空间 $H$！所以当 $\mathcal{D}(A) \subsetneq H$ 时，**$0A \neq 0$**，只能写 $0A \subset 0$！（当然作用在零向量上恒有 $A0 = 0$。）
- **(ii) 和 $A + B$**：定义域只能取公共部分 $\mathcal{D}(A+B) = \mathcal{D}(A) \cap \mathcal{D}(B)$，作用 $x \mapsto Ax + Bx$。
- **(iii) 复合（积）$BA$**：先算 $Ax$，再要求 $Ax$ 落进 $B$ 的定义域，故 $\mathcal{D}(BA) = \{x \in \mathcal{D}(A) \mid Ax \in \mathcal{D}(B)\}$，作用 $x \mapsto B(Ax)$。
  - 特别地，例 7.1 中乘法算子 $(Af)(t)=tf(t)$ 的 $n$ 次幂 $A^n$ 的定义域是 $\mathcal{D}(A^n) = \{f \in L^2(\mathbb{R}) \mid \int_{\mathbb{R}} |t^n f(t)|^2\,\mathrm{d}t < \infty\}$。
- **(iv) 延拓 / 扩张 $A \subset B$**：如果 **$\mathcal{D}(A) \subset \mathcal{D}(B)$**，并且对所有 $x \in \mathcal{D}(A)$ 都有 **$Ax = Bx$**，就称 $B$ 是 $A$ 的**延拓（或扩张）**，记作 **$A \subset B$**（或 $B \supset A$）。
  - *通俗直觉*：$B$ 干的活跟 $A$ 一模一样，只是 $B$ 的“管辖范围（定义域）”比 $A$ 更大（或相等）。
- **(v) 算子相等 $A = B$**：$A \subset B$ 且 $B \subset A$，即 **$\mathcal{D}(A) = \mathcal{D}(B)$ 且对一切 $x \in \mathcal{D}(A)$ 有 $Ax = Bx$**。

#### 2. 稠定算子与共轭算子 $A^*$（定义 7.2–7.3、引理 7.1，p. 273–274）

在有界算子里，共轭算子 $A^*$ 由公式 $(Ax, y) = (x, A^*y)$（$\forall x, y \in H$）定义。现在 $A$ 只定义在子空间 $\mathcal{D}(A)$ 上，我们还能用同一招来定义 $A^*$ 吗？

这里会遇到**两个必须解决的问题**：
1. **唯一性问题**：如果已知对一切 $x \in \mathcal{D}(A)$ 都有 $(Ax, y) = (x, y^*)$，这个 $y^*$ 会不会不唯一？
   - 设 $y_1^*, y_2^*$ 都满足上式，相减得 $(x, y_1^* - y_2^*) = 0$（$\forall x \in \mathcal{D}(A)$），即 $y_1^* - y_2^* \perp \mathcal{D}(A)$。
   - 要逼出 $y_1^* - y_2^* = 0$，当且仅当 $\mathcal{D}(A)^\perp = \{0\}$，也就是 **$\overline{\mathcal{D}(A)} = H$（$\mathcal{D}(A)$ 在 $H$ 中稠密）**！
2. **存在性问题**：即使 $\mathcal{D}(A)$ 稠密，也不是每个 $y \in H$ 都能找到对应的 $y^*$！
   - 固定 $y \in H$，看线性泛函 $\varphi_y(x) \triangleq (Ax, y)$（$x \in \mathcal{D}(A)$）。
   - 要存在 $y^* \in H$ 使 $\varphi_y(x) = (x, y^*)$，由 **Riesz 表示定理**（将连续线性泛函从稠密子空间唯一延拓到 $H$），充要条件是 $\varphi_y$ 在 $\mathcal{D}(A)$ 上**有界**，即存在常数 $M_y > 0$（依赖于 $y$）使得
     $$|(Ax, y)| \le M_y \|x\|, \qquad \forall x \in \mathcal{D}(A).$$
   - 只有那些能让上式成立的“好向量 $y$”，才有资格进入 $\mathcal{D}(A^*)$！

这就精准引出了教材的**定义 7.2、定义 7.3 与引理 7.1**：

- **定义 7.2（稠定算子，p. 273）**：若线性子空间 $\mathcal{D}(A)$ 在 $H$ 中稠密（$\overline{\mathcal{D}(A)} = H$），则称 $A: \mathcal{D}(A) \to H$ 是 $H$ 上的**稠定算子**。
- **定义 7.3（共轭算子，p. 273）**：设 $A$ 是 $H$ 上的稠定线性算子，令
  $$\mathcal{D}(A^*) \triangleq \left\{y \in H \;\middle|\; \text{存在 } y^* \in H \text{ 使得对任意 } x \in \mathcal{D}(A),\ (Ax, y) = (x, y^*)\right\},$$
  定义算子 $A^*: \mathcal{D}(A^*) \to H$ 为 $A^* y = y^*$，称 $A^*$ 为 $A$ 的**共轭算子**（又称伴随算子）。
- **引理 7.1（唯一性的充要条件，p. 273–274）**：对 $y \in H$，满足 $(Ax, y) = (x, y^*)$（$\forall x \in \mathcal{D}(A)$）的向量 $y^* \in H$ **最多只有一个**的充要条件是 **$A$ 稠定**。

> **零基础警示**：千万别以为“$A$ 稠定”就能推出“$\mathcal{D}(A^*)$ 稠密”！事实上存在稠定算子满足 $\mathcal{D}(A^*) = \{0\}$（习题 7 第 1 题就指出：$\mathcal{D}(A^*) = \{0\} \iff A$ 的图像在 $H \times H$ 中稠密）。所以写 $A^{**} = (A^*)^*$ 之前，**必须先验证 $\mathcal{D}(A^*)$ 在 $H$ 中稠密**！（这是习题 7 第 6 题的关键得分点！）

#### 3. 共轭算子的基本运算性质与“扩张反号铁律”（引理 7.2，p. 274–275）

- **(i) 线性性**：$A^*$ 是线性算子（即 $\mathcal{D}(A^*)$ 是 $H$ 的线性子空间，且 $A^*(\alpha y_1 + \beta y_2) = \alpha A^* y_1 + \beta A^* y_2$）。
- **(ii) 和的共轭**：若 $A, B$ 及 $A+B$ 均稠定（即 $\mathcal{D}(A)\cap\mathcal{D}(B)$ 稠密），则
  $$A^* + B^* \subset (A + B)^*.$$
- **(iii) 扩张反号铁律（极其高频！）**：设 $A, B$ 都是稠定线性算子，若 **$A \subset B$**，则 **$B^* \subset A^*$**（即 $A^* \supset B^*$）。
  - **一句话悟透“为什么会反号”**：
    想想进入 $\mathcal{D}(B^*)$ 的“入场考试”：你必须通过**所有 $x \in \mathcal{D}(B)$** 的考题 $(Bx, y) = (x, y^*)$；而进入 $\mathcal{D}(A^*)$ 只需通过**子集 $\mathcal{D}(A)$** 的考题 $(Ax, y) = (x, y^*)$！
    **考题越多（$A \subset B$，范围更大、约束越严），能通过考试的人反而越少（$\mathcal{D}(B^*) \subset \mathcal{D}(A^*)$）！** 记住这个直觉，这辈子都不会把包含号写反！
- **(iv) 积的共轭**：设 $A, B$ 和 $BA$ 均为稠定算子，则
  $$A^* B^* \subset (BA)^*.$$
  特别地，如果还有 $B \in \mathcal{L}(H)$（外层算子是全空间有界算子），则等号成立：**$A^* B^* = (BA)^*$**。

> **教材印刷勘误提示（p. 274 引理 7.2）**：教材 p. 274 引理 7.2 的 (ii) 陈述与末尾把包含符号印反成了 $(A+B)^* \subset A^* + B^*$，但查看该页正文证明第一行“如果 $y \in \mathcal{D}(A^*+B^*)$，则 $y \in \mathcal{D}((A+B)^*)$”，实际证明出的正确方向就是 $A^* + B^* \subset (A+B)^*$；(iv) 亦存在 $AB$ 与 $BA$ 混印，正确形式为上文给出的 $A^* B^* \subset (BA)^*$（当 $B \in \mathcal{L}(H)$ 时取等）。复习与答题请以正确结论为准！

---

### 7.1.2 对称算子 vs 自伴算子：一字之差，命运天壤之别（p. 275–277）

#### 1. 精确定义与等价判别（定义 7.4、引理 7.3–7.4，p. 275–277）

- **定义 7.4（p. 275）**：设 $A$ 是 $H$ 上的**稠定**线性算子：
  - 若 **$A \subset A^*$**，称 $A$ 为**对称算子**（Symmetric operator）；
  - 若 **$A = A^*$**，称 $A$ 为**自伴算子**（Self-adjoint operator）。
- **引理 7.3（对称算子的实操检验公式，p. 277）**：稠定线性算子 $A$ 是对称算子，当且仅当对任意 $x, y \in \mathcal{D}(A)$，都有
  $$(Ax, y) = (x, Ay).$$
  *证明直观*：只要对所有 $x, y \in \mathcal{D}(A)$ 满足 $(Ax, y) = (x, Ay)$，按定义 7.3 立即看出每个 $y \in \mathcal{D}(A)$ 自动满足 $y \in \mathcal{D}(A^*)$ 且 $A^* y = Ay$，即 $A \subset A^*$！但注意：这**完全无法保证反方向的 $\mathcal{D}(A^*) \subset \mathcal{D}(A)$**！
- **引理 7.4（自伴算子的极大性，p. 277）**：设 $A$ 是自伴算子，$B$ 是对称算子，如果 $A \subset B$，那么必有 **$A = B$**！
  - *一行神证明*：由 $A \subset B$ 及扩张反号律得 $B^* \subset A^*$，于是
    $$A = A^* \supset B^* \supset B \supset A \implies A = B!$$
  - *几何直觉*：对称算子一旦“长”到了自伴，就已经到达了对称性的天花板，再也不可能被真扩张成更大的对称算子！

#### 2. 醒目对比表：对称算子 vs 自伴算子

| 对比项 | 对称算子（Symmetric） | 自伴算子（Self-adjoint） |
| :--- | :--- | :--- |
| **核心关系式** | $A \subset A^*$ | $A = A^*$ |
| **定义域关系** | 允许严格更小：$\mathcal{D}(A) \subsetneq \mathcal{D}(A^*)$ | 必须完全相等：**$\mathcal{D}(A) = \mathcal{D}(A^*)$** |
| **验证难度** | **容易**：只需验证 $\mathcal{D}(A)$ 稠密且 $(Ax,y)=(x,Ay)$ | **困难得多**：除对称外，必须证明 $\mathcal{D}(A^*) \subset \mathcal{D}(A)$（或用定理 7.4 判别） |
| **闭性** | 未必是闭算子（但必**可闭**，因 $\mathcal{D}(A^*)$ 稠密） | **必然是闭算子**（定理 7.1） |
| **谱的性质** | 谱**未必**在实轴上！甚至可能是整个复平面 $\mathbb{C}$（见例 7.2 中的 $A_3$） | 谱必然全在实轴上 **$\sigma(A) \subset \mathbb{R}$**，且无剩余谱（推论 7.2） |
| **谱分解定理** | **不成立**（无法生成酉群 $\mathrm{e}^{\mathrm{i}tA}$，在量子力学中不能代表物理观测量） | **成立** $A = \int_{\mathbb{R}} \lambda\,\mathrm{d}E_\lambda$（定理 7.8） |

#### 3. 教科书级活例子：边界条件如何决定算子的命运（例 7.2，p. 275–277）

在 $H = L^2(0, 1)$ 上，考察**完全相同的微分运算公式** $A_k f = \mathrm{i} f'$，只改变边界条件（$AC[0,1]$ 表示 $[0,1]$ 上的绝对连续函数，可分部积分）：

$$\begin{aligned}
\mathcal{D}(A_1) &= \{f \in L^2(0,1) \mid f \in AC[0,1],\ f' \in L^2(0,1)\} && \text{（最松：无任何边界条件）}, \\
\mathcal{D}(A_2) &= \{f \in \mathcal{D}(A_1) \mid f(0) = f(1)\} && \text{（适中：周期边界条件）}, \\
\mathcal{D}(A_3) &= \{f \in \mathcal{D}(A_1) \mid f(0) = f(1) = 0\} && \text{（最紧：两端固定为零）}.
\end{aligned}$$

显然 $A_3 \subsetneq A_2 \subsetneq A_1$，且三者都是稠定的。做一次分部积分看边界项：
$$(A_k f, g) = \int_0^1 \mathrm{i} f'(x) \overline{g(x)}\,\mathrm{d}x = \underbrace{\mathrm{i}\left[f(1)\overline{g(1)} - f(0)\overline{g(0)}\right]}_{\text{必须消掉的边界项！}} + \int_0^1 f(x) \overline{\mathrm{i} g'(x)}\,\mathrm{d}x.$$

为了消掉这个边界项，三个算子走上了截然不同的命运：
- **$A_3$（约束太紧）**：$f(0)=f(1)=0$ 直接把边界项杀死了，所以对 $g$ 的端点值**没有任何要求**！结果 **$A_3^* = A_1 \supsetneq A_3$** —— $A_3$ 是**对称算子，但不是自伴算子**！（它的谱甚至是整个复平面 $\sigma(A_3) = \mathbb{C}$！）
- **$A_1$（约束太松）**：$f(0), f(1)$ 随便取，为了消掉边界项，只好强迫 $g(0)=g(1)=0$！结果 **$A_1^* = A_3 \subsetneq A_1$** —— $A_1$ **连对称算子都不是**！
- **$A_2$（松紧恰到好处）**：$f(0)=f(1)$，边界项变成 $\mathrm{i} f(0)[\overline{g(1)} - \overline{g(0)}]$，为使其对一切 $f(0)$ 为零，恰好等价于要求 $g(0)=g(1)$！结果 **$A_2^* = A_2$** —— **$A_2$ 是自伴算子**！（其谱恰好是离散实点集 $\{2n\pi \mid n \in \mathbb{Z}\}$！）

> **灵魂总结**：结合马上要学的 **Hellinger-Toeplitz 定理（定理 7.3 (i)）**——“定义在全空间 $H$ 上的对称算子必定有界”，我们立刻明白：**任何无界对称算子都绝不可能定义在整个 $H$ 上**。为了成为自伴算子，定义域不能太小（太小则共轭太大，只是对称），也不能太大（太大则共轭太小，失去对称），必须精准调配边界条件达到 $A = A^*$ 的平衡！

#### 4. 两大核心自伴算子范例与酉等价（例 7.3、7.1.3 节引理 7.5、定义 7.5、例 7.4，p. 277–279）

- **例 7.3（位置/乘法算子自伴，p. 277）**：$L^2(\mathbb{R})$ 上的乘法算子 $(Af)(t) = tf(t)$（定义域 $\mathcal{D}(A) = \{f \in L^2(\mathbb{R}) \mid tf \in L^2(\mathbb{R})\}$）不仅是闭的无界算子（习题 7 第 4 题），而且是**自伴算子**（$A = A^*$）。
- **引理 7.5 & 定义 7.5（酉等价，p. 277–278）**：若 $U$ 是 $H$ 上的酉算子（$U^* = U^{-1}$），称 $B = UAU^{-1}$（定义域 $\mathcal{D}(B) = U\mathcal{D}(A)$）与 $A$ **酉等价**。**若 $A$ 自伴，则 $UAU^{-1}$ 依然自伴**（且谱完全相同！）。
- **例 7.4（Fourier 变换把微分变成乘法，p. 278–279）**：在 $L^2(\mathbb{R})$ 中，动量算子 $A_1 f = \frac{1}{\mathrm{i}} f'$（定义域为局部绝对连续且 $f, f' \in L^2(\mathbb{R})$）通过酉的 Fourier 变换 $U$ 与乘法算子 $A_2 f(t) = tf(t)$ 满足 **$A_1 = U A_2 U^{-1}$**！由此一行秒证全直线上的微分算子 $\frac{1}{\mathrm{i}}\frac{\mathrm{d}}{\mathrm{d}t}$ 是自伴算子！

---

### 7.1.4 算子的图像、闭算子与可闭算子（p. 279–280）

#### 1. 什么是闭算子？——“退一步海阔天空”的连续性替身

在 $H \times H = \{\{a, b\} \mid a, b \in H\}$ 上赋予内积 $(\{a,b\}, \{c,d\}) = (a,c) + (b,d)$ 和范数 $\|\{a,b\}\|^2 = \|a\|^2 + \|b\|^2$，它仍是 Hilbert 空间。

- **算子的图像**：$G(A) \triangleq \{\{x, Ax\} \mid x \in \mathcal{D}(A)\} \subset H \times H$。
- **闭算子的等价定义**：称 $A$ 是**闭算子**，如果它的图像 $G(A)$ 是 $H \times H$ 中的闭子空间；用序列语言翻译，就是：
  $$\text{若 } \{x_n\} \subset \mathcal{D}(A) \text{ 满足 } x_n \to x \text{ 且 } Ax_n \to y \quad (x, y \in H),\quad \text{则必有 } x \in \mathcal{D}(A) \text{ 且 } Ax = y.$$

**零基础秒懂：对比“连续算子”与“闭算子”**
- **连续算子（高配版）**：只要你喂给我 $x_n \to x$，我**打包票** $Ax_n$ 一定收敛，而且收敛到 $Ax$。
- **闭算子（务实版）**：我承认我是无界算子，你给 $x_n \to x$ 时，$Ax_n$ 完全可能发散到无穷大（比如 $u_n(x)=\frac{1}{n}\sin(n^2 x) \to 0$，但导数 $u_n'(x)=n\cos(n^2 x)$ 发散）；**但**如果 $x_n \to x$ **而且** $Ax_n$ 碰巧也收敛到某个 $y$，那我担保两件事：① $x$ 绝不会逃出定义域 $\mathcal{D}(A)$；② 极限值绝不撒谎，$y$ 必定等于 $Ax$！

#### 2. 旋转算子 $V$ 与图像正交补（引理 7.6、推论 7.1，p. 279–280）

在 $H \times H$ 上定义一个把坐标旋转 $90^\circ$ 的酉算子：
$$V\{a, b\} = \{-b, a\}, \qquad \text{显然 } V^2 = -I.$$

观察 $(Ax, y) = (x, z)$ 等价于 $(-Ax, y) + (x, z) = 0$，即 $(\{-Ax, x\}, \{y, z\}) = 0$！也就是 $\{y, z\} \perp V\{x, Ax\}$！这一神来之笔直接给出整个 7.1.4 节的底层几何引擎：

- **引理 7.6（p. 279–280）**：设 $A$ 是 $H$ 上的稠定线性算子，则
  $$G(A^*) = [V G(A)]^\perp.$$
- **推论 7.1（p. 280）**：如果 $A$ 是 $H$ 上的**稠定闭算子**（从而 $VG(A)$ 是闭子空间），则有正交直和分解
  $$H \times H = V G(A) \oplus G(A^*).$$

#### 3. 两大核心柱石定理（定理 7.1–7.2，p. 280）

由引理 7.6 和推论 7.1，几乎不费吹灰之力就能摘下本章最重要的两颗果实（**习题 7 第 6 题直接考这两条！**）：

- **定理 7.1（共轭算子恒闭，p. 280）**：设 $A$ 是 $H$ 上的稠定线性算子，则 **$A^*$ 是闭算子**。特别地，**自伴算子必是闭算子**。
  - *一句话秒记证明*：因为 $G(A^*) = [VG(A)]^\perp$，而任何集合的正交补 $M^\perp$ 永远是闭子空间！
- **定理 7.2（稠定闭算子的二阶共轭还原，p. 280）**：设 $A$ 是 $H$ 上的**稠定闭算子**，则 **$\mathcal{D}(A^*)$ 在 $H$ 中稠密，且 $A = A^{**}$**。
  - *直观理解*：在 $H \times H = VG(A) \oplus G(A^*)$ 两边再乘一次 $V$（用 $V^2 = -I$），得 $H \times H = G(A) \oplus VG(A^*)$，这就迫使 $G(A) = [VG(A^*)]^\perp = G(A^{**})$，即 $A = A^{**}$！

#### 4. 进阶串联：可闭算子、闭包 $\bar{A} = A^{**}$ 与对称算子的四级阶梯

如果一个算子 $A$ 还不是闭算子，我们能不能“把缺的极限点补进去”，把它扩张成一个闭算子？

- **可闭算子**：若 $A$ 存在闭扩张，称 $A$ **可闭**。等价判别：图像闭包 $\overline{G(A)}$ 中不含形如 $\{0, y\}$（$y \neq 0$）的“一个输入对应两个输出”的坏点（即不存在 $x_n \in \mathcal{D}(A)$ 使 $x_n \to 0$ 而 $Ax_n \to y \neq 0$）。
- **最小闭扩张（闭包 $\bar{A}$）**：当 $A$ 可闭时，以 $\overline{G(A)}$ 为图像的算子 $\bar{A}$ 称为 $A$ 的闭包（最小闭扩张）。
- **稠定可闭算子的黄金公式**：稠定算子 $A$ 可闭 $\iff \mathcal{D}(A^*)$ 稠密，此时 **$\bar{A} = A^{**}$**（且 $\bar{A}^* = A^*$）！
- **特别地，对任意对称算子 $A$**：因 $A \subset A^*$ 自动保证 $\mathcal{D}(A) \subset \mathcal{D}(A^*)$ 稠密，故**所有对称算子都可闭**，且恒满足**黄金包含链**：
  $$A \subset \bar{A} = A^{**} \subset A^*.$$
  由此形成清晰的**四级阶梯**：
  1. **一般对称算子**：$A \subsetneq \bar{A} = A^{**} \subsetneq A^*$；
  2. **闭对称算子（非自伴）**：$A = A^{**} \subsetneq A^*$（习题 7 第 6 题 (i)）；
  3. **本质自伴算子**（闭包自伴）：$A \subsetneq \bar{A} = A^{**} = A^*$；
  4. **自伴算子**：$A = A^{**} = A^*$（习题 7 第 6 题 (ii)）。

---

### 7.1.5 对称算子成为自伴算子的判别准则（p. 280–282）

实际做题时，直接去算 $\mathcal{D}(A^*)$ 常常非常麻烦。本节给出判断一个对称算子到底是不是自伴算子的几把“快刀”：

#### 1. 两个充分条件与 Hellinger-Toeplitz 定理（定理 7.3，p. 280–281）

设 $A$ 是 $H$ 上的**对称算子**：
- **(i) Hellinger-Toeplitz 定理**：如果 **$\mathcal{D}(A) = H$**（到处定义），则 **$A$ 是自伴的，且 $A \in \mathcal{L}(H)$（必为有界算子）**！
  - *两行证明*：由 $H = \mathcal{D}(A) \subset \mathcal{D}(A^*) \subset H$ 得 $A = A^*$（自伴）；由定理 7.1 知 $A$ 是定义在全空间 $H$ 上的闭算子，再由第 3 章**闭图像定理**立得 $A$ 有界！
- **(ii) 满射判别法**：如果 **$\mathcal{R}(A) = H$**（值域为全空间），则 **$A$ 是自伴的，且 $A^{-1} \in \mathcal{L}(H)$**！

#### 2. 虚数移位 $A \pm \mathrm{i}I$ 的魔力与自伴基本判别准则（引理 7.7–7.9、定理 7.4，p. 281–282）

对于任意**对称算子 $A$** 和任意 $x \in \mathcal{D}(A)$，把 $\|(A \pm \mathrm{i}I)x\|^2$ 展开，交叉项因对称性精准抵消（$(Ax, \pm\mathrm{i}x) + (\pm\mathrm{i}x, Ax) = \mp\mathrm{i}(Ax,x) \pm \mathrm{i}(x,Ax) = 0$），得到贯穿后半章的**勾股恒等式（式 7.9，p. 281）**：
$$\|(A \pm \mathrm{i}I)x\|^2 = \|Ax\|^2 + \|x\|^2, \qquad \forall x \in \mathcal{D}(A).$$

由式 (7.9) 立刻推倒一串多米诺骨牌：
- **引理 7.7（p. 281）**：若 $A$ 对称，则 **$\operatorname{Ker}(A \pm \mathrm{i}I) = \{0\}$**（$A \pm \mathrm{i}I$ 必单射）。
- **引理 7.8（p. 281）**：若 $A$ 是**闭对称算子**，则 **$\mathcal{R}(A \pm \mathrm{i}I)$ 是闭子空间**。（*对比习题 7 第 7 题：当对称算子额外非负 $(Ax,x)\ge 0$ 时，实数移位 $A+I$ 扮演了完全相同的角色！*）
- **引理 7.9（p. 281）**：对任意稠定算子 $A$，恒有正交关系
  $$\operatorname{Ker} A^* = (\mathcal{R}(A))^\perp, \qquad \overline{\mathcal{R}(A)} = (\operatorname{Ker} A^*)^\perp.$$

将三块骨牌拼在一起，就得到了**自伴算子基本判别准则（定理 7.4，p. 281–282）**：

> **定理 7.4**：设 $A$ 是 Hilbert 空间 $H$ 上的**对称算子**，则以下三个命题**等价**：
> - **(i)** $A$ 是**自伴算子**（$A = A^*$）；
> - **(ii)** $A$ 是**闭算子**，且 **$\operatorname{Ker}(A^* + \mathrm{i}I) = \{0\}$** 与 **$\operatorname{Ker}(A^* - \mathrm{i}I) = \{0\}$**（合记为 $\operatorname{Ker}(A^* \pm \mathrm{i}I) = \{0\}$，即 $A^*$ 没有特征值 $+\mathrm{i}$ 和 $-\mathrm{i}$）；
> - **(iii)** **$\mathcal{R}(A + \mathrm{i}I) = H$** 且 **$\mathcal{R}(A - \mathrm{i}I) = H$**（合记为 $\mathcal{R}(A \mp \mathrm{i}I) = H$）。

**零基础易错点**：(ii) 中**千万不能漏掉“$A$ 是闭算子”**！如果 $A$ 不闭，即使 $\operatorname{Ker}(A^* \pm \mathrm{i}I) = \{0\}$，也只能推出 $\overline{\mathcal{R}(A \pm \mathrm{i}I)} = H$（$A$ 本质自伴，即 $\bar{A}$ 自伴），不能推出 $A$ 本身自伴！同时，(iii) 中**正负两个符号 $+\mathrm{i}$ 和 $-\mathrm{i}$ 缺一不可**！

---

### 7.2 无界自伴算子的谱、Cayley 变换与谱分解（p. 282–293）

#### 1. 无界算子的预解集与谱（定义 7.6–7.7，p. 282–283）

设 $X$ 是复 Banach 空间，$A: \mathcal{D}(A) \subset X \to X$ 是线性算子：

- **定义 7.6（预解集与谱，p. 282）**：
  $$\rho(A) = \left\{\lambda \in \mathbb{C} \;\middle|\; \lambda I - A \text{ 是 } \mathcal{D}(A) \to X \text{ 的双射，且 } (\lambda I - A)^{-1} \in \mathcal{L}(X)\right\}.$$
  称 $\rho(A)$ 为 $A$ 的**预解集**（其中点称**正则点**），$R(\lambda, A) \triangleq (\lambda I - A)^{-1}$ 为**预解式（预解算子）**，$\sigma(A) \triangleq \mathbb{C} \setminus \rho(A)$ 为 $A$ 的**谱集**。
  > **为什么这里必须手动补写 $(\lambda I - A)^{-1} \in \mathcal{L}(X)$？** 第 6 章里 $A$ 有界，Banach 逆算子定理自动保证逆算子有界；现在 $A$ 可能既无界又不闭！如果 $A$ 是**闭算子**，则 $\lambda I - A$ 是闭算子，当它是双射时，逆算子 $(\lambda I - A)^{-1}$ 是定义在**全空间 $X$** 上的闭算子，由**闭图像定理**仍自动有界（属于 $\mathcal{L}(X)$）！
- **定义 7.7（闭线性算子谱的分类，p. 282–283）**：对闭线性算子 $A$，每个 $\lambda \in \sigma(A)$ 恰好落入以下三种互不相交的情形之一：
  - **(i) 点谱（特征值）$\sigma_p(A)$**：$\lambda I - A$ 不是单射（即存在非零 $x \in \mathcal{D}(A)$ 使 $Ax = \lambda x$）；
  - **(ii) 连续谱 $\sigma_c(A)$**：$\lambda I - A$ 是单射，$\mathcal{R}(\lambda I - A) \neq X$ 但 $\overline{\mathcal{R}(\lambda I - A)} = X$；
  - **(iii) 剩余谱 $\sigma_r(A)$**：$\lambda I - A$ 是单射，且 $\overline{\mathcal{R}(\lambda I - A)} \neq X$。

#### 2. 自伴算子的谱必在实轴上（定理 7.5、推论 7.2，p. 283–284）

- **定理 7.5（下有界判别法 / Weyl 准则，p. 283）**：设 $A$ 是 $H$ 上的**自伴算子**，则
  $$\lambda \in \rho(A) \iff \text{存在常数 } C > 0 \text{ 使得 } \|(\lambda I - A)x\| \ge C \|x\|,\quad \forall x \in \mathcal{D}(A).$$
  - *逆否命题（求谱利器！）*：$\lambda \in \sigma(A) \iff$ 存在单位向量列 $\{x_n\} \subset \mathcal{D}(A)$（$\|x_n\|=1$）使得 $\|(\lambda I - A)x_n\| \to 0$！（这就是证明命题 7.1 $\sigma(A)=\mathbb{R}$ 的核心武器！）
- **推论 7.2（谱的实性，p. 283–284）**：设 $A$ 是 $H$ 上的自伴算子，则
  $$\sigma(A) \subset \mathbb{R}, \qquad \text{且 } \sigma_r(A) = \varnothing.$$
  - *证明直觉*：把 $\lambda = \alpha + \mathrm{i}\beta$（$\beta \neq 0$）代入展开，对称性使交叉项再次消去，得
    $$\|(\lambda I - A)x\|^2 = \|(\alpha I - A)x\|^2 + \beta^2 \|x\|^2 \ge |\beta|^2 \|x\|^2!$$
    由定理 7.5 立得一切非实复数 $\lambda$ 都属于 $\rho(A)$，从而 $\sigma(A) \subset \mathbb{R}$（且预解式满足估计 $\|R(\lambda, A)\| \le \frac{1}{|\operatorname{Im}\lambda|}$）！

#### 3. Cayley 变换：把无界自伴算子“压”成有界酉算子（定理 7.6–7.7、定义 7.9、推论 7.3，p. 285–286）

想把第 6 章的谱分解定理推广到无界自伴算子 $A$，最大障碍就是 $A$ **无界**（谱跑到无穷远，做不了有界多项式逼近）。怎么化无界为有界？

复变函数里有一个神奇的映射：
$$u(\lambda) = \frac{\lambda - \mathrm{i}}{\lambda + \mathrm{i}}.$$
当实数 $\lambda$ 从 $-\infty$ 跑到 $+\infty$ 时，$|u(\lambda)| = \frac{|\lambda - \mathrm{i}|}{|\lambda + \mathrm{i}|} = 1$ 始终落在**单位圆周**上，只是永远达不到点 $1$（$\lambda \to \pm\infty$ 时趋近于 $1$）。

把这个标量公式直接换成算子，就是 **Cayley 变换（定义 7.9，p. 286）**：
$$U = (A - \mathrm{i}I)(A + \mathrm{i}I)^{-1}.$$

- **定理 7.6（正向 Cayley 变换，p. 285–286）**：设 $A$ 是对称算子，令 $U = (A - \mathrm{i}I)(A + \mathrm{i}I)^{-1}$，则：
  - (i) $U$ 是 $\mathcal{R}(A+\mathrm{i}I) \to \mathcal{R}(A-\mathrm{i}I)$ 的**等距线性算子**（因式 7.9 保证 $\|(A-\mathrm{i}I)x\| = \|(A+\mathrm{i}I)x\|$）；
  - (ii) **$1 \notin \sigma_p(U)$**（$I - U$ 是单射）；
  - (iii) 可以精确反解回 $A$：**$A = \mathrm{i}(I+U)(I-U)^{-1}$**，且 $\mathcal{D}(A) = \mathcal{R}(I - U)$；
  - (iv) $A$ 闭 $\implies U$ 闭；**特别地，当 $A$ 自伴时**，由定理 7.4 知 $\mathcal{R}(A\pm\mathrm{i}I) = H$，从而 **$U$ 是全空间 $H$ 上的（有界）酉算子**！
- **推论 7.3（p. 286）**：$1 \in \rho(U) \iff A$ 是**有界**自伴算子。（因此无界自伴算子恰好对应 $1 \in \sigma_c(U)$！）
- **定理 7.7（逆向 Cayley 变换，p. 286）**：反之，若 $U$ 是等距算子且 $\mathcal{R}(I-U)$ 稠密，则 $1 \notin \sigma_p(U)$，且 $A = \mathrm{i}(I+U)(I-U)^{-1}$ 是定义在 $\mathcal{R}(I-U)$ 上的对称算子；**特别地，若 $U$ 是酉算子，则 $A$ 是自伴算子**！

#### 4. 无界自伴算子的谱分解定理与乘法算子范例（定义 7.10、定理 7.8、命题 7.1–7.2，p. 290–293）

有了 Cayley 变换，无界自伴算子的谱分解就水到渠成：
- **路线图**：无界自伴算子 $A$ $\xrightarrow{\text{Cayley变换}}$ 有界酉算子 $U = \int_0^{2\pi} \mathrm{e}^{\mathrm{i}\theta}\,\mathrm{d}F_\theta$（第 6 章推论 6.11，因 $1 \notin \sigma_p(U)$ 故 $F_0 = 0, \lim_{\theta\to 2\pi^-} F_\theta = I$）$\xrightarrow{\text{换元 }\lambda = -\cot(\theta/2),\ E_\lambda = F_\theta}$ 还原回实轴上的谱分解！
- **定理 7.8（自伴算子谱分解定理，p. 290–292）**：设 $A$ 是 $H$ 上的自伴算子，则必有 $H$ 上唯一的实轴谱系 $\{E_\lambda \mid \lambda \in \mathbb{R}\}$（单调递增、右连续、强极限 $E_{-\infty}=0, E_{+\infty}=I$ 的正交投影族）使得
  $$A = \int_{\mathbb{R}} \lambda\,\mathrm{d}E_\lambda,$$
  其定义域精确刻画为（按定义 7.10 无界函数谱积分）：
  $$\mathcal{D}(A) = \left\{x \in H \;\middle|\; \int_{\mathbb{R}} \lambda^2\,\mathrm{d}\|E_\lambda x\|^2 < +\infty\right\},$$
  且对 $x \in \mathcal{D}(A), y \in H$ 有 $\|Ax\|^2 = \int_{\mathbb{R}} \lambda^2\,\mathrm{d}\|E_\lambda x\|^2$，$(Ax, y) = \int_{\mathbb{R}} \lambda\,\mathrm{d}(E_\lambda x, y)$。
- **乘法算子的完整画像（命题 7.1–7.2，p. 292–293）**：对于例 7.1（即习题 7 第 4 题）中的 $L^2(\mathbb{R})$ 位置算子 $(Ax)(t) = tx(t)$：
  - **命题 7.1**：$\sigma_p(A) = \sigma_r(A) = \varnothing$，$\sigma_c(A) = \sigma(A) = \mathbb{R}$。（一个特征值都没有，整条实轴都是连续谱！）
  - **命题 7.2**：它的谱系就是把函数在 $(-\infty, \lambda]$ 上截断的乘法投影：
    $$(E_\lambda f)(t) = \chi_{(-\infty, \lambda]}(t) f(t), \qquad A = \int_{\mathbb{R}} \lambda\,\mathrm{d}E_\lambda.$$

---

### 第 7 章考前速记总表

| 考点模块 | 必背公式 / 核心结论 | 教材编号 | 考试最高频踩坑点 |
| :--- | :--- | :--- | :--- |
| **稠定与共轭** | $(Ax, y) = (x, A^*y)$；$A^*$ 唯一 $\iff \overline{\mathcal{D}(A)} = H$ | 定义 7.2–7.3、引理 7.1 | 没证 $\mathcal{D}(A^*)$ 稠密就直接写 $A^{**}$ |
| **扩张反号律** | $A \subset B \implies B^* \subset A^*$ | 引理 7.2 (iii) | 把包含方向写正了（记住：约束越多，满足者越少！） |
| **对称 vs 自伴** | 对称：$A \subset A^*$；自伴：$A = A^*$；$A$ 自伴、$B$ 对称、$A \subset B \implies A = B$ | 定义 7.4、引理 7.3–7.4 | 只验 $(Ax,y)=(x,Ay)$ 就敢宣称 $A$ 自伴，忘了核对定义域！ |
| **闭算子与 $A^{**}$** | $G(A^*) = [VG(A)]^\perp$；$A^*$ 恒闭；$A$ 稠定闭 $\implies \mathcal{D}(A^*)$ 稠密且 $A=A^{**}$ | 引理 7.6、定理 7.1–7.2 | 证闭性时只证了 $Ax=y$ 却漏证 $x \in \mathcal{D}(A)$ |
| **Hellinger-Toeplitz** | $A$ 对称且 $\mathcal{D}(A) = H \implies A$ 自伴且 $A \in \mathcal{L}(H)$ | 定理 7.3 (i) | 误以为无界自伴算子可以定义在整个 $H$ 上 |
| **自伴判别准则** | 对称算子 $A$ 自伴 $\iff A$ 闭且 $\operatorname{Ker}(A^*\pm\mathrm{i}I)=\{0\} \iff \mathcal{R}(A\mp\mathrm{i}I)=H$ | 定理 7.4 | 漏掉“闭算子”条件，或只验了 $+\mathrm{i}$ 漏了 $-\mathrm{i}$ |
| **谱与下有界性** | $A$ 自伴：$\lambda \in \rho(A) \iff \exists C>0,\ \lVert (\lambda I-A)x\rVert \ge C\lVert x\rVert$；$\sigma(A)\subset\mathbb{R},\ \sigma_r(A)=\varnothing$ | 定义 7.6–7.7、定理 7.5、推论 7.2 | 无界自伴算子的谱可以**无界**（如 $\sigma(A)=\mathbb{R}$） |
| **Cayley 变换** | $U = (A-\mathrm{i}I)(A+\mathrm{i}I)^{-1}$，$A = \mathrm{i}(I+U)(I-U)^{-1}$；$A$ 自伴 $\iff U$ 酉且 $1 \notin \sigma_p(U)$ | 定理 7.6–7.7、定义 7.9 | 正负号记混（记住分子是 $-\mathrm{i}$、分母是 $+\mathrm{i}$） |
| **谱分解定理** | $A = \int_{\mathbb{R}} \lambda\,\mathrm{d}E_\lambda$，$\mathcal{D}(A) = \{x \colon \int_{\mathbb{R}} \lambda^2\,\mathrm{d}\lVert E_\lambda x\rVert^2 < \infty\}$ | 定理 7.8、命题 7.1–7.2 | 乘法算子谱系 $E_\lambda f = \chi_{(-\infty,\lambda]} f$ |

---

## 第二部分：第 7 章课后作业逐题详解（习题 7 第 3、4、6、7 题）

> **全题前置约定**：按教材习题 7（p. 293）题首说明，以下各题如无特别说明，均假设 $H$ 是复 Hilbert 空间。

---

### 习题 7 第 3 题（教材 p. 293）

#### 【原题】

定义算子 $A : C^2([-1, 1]) \to C([-1, 1])$ 如下：
$$\mathcal{L}u(x) = -\frac{\mathrm{d}}{\mathrm{d}x}\left[(1 - x^2)\frac{\mathrm{d}u}{\mathrm{d}x}\right], \quad x \in [-1, 1].$$
求出 $A$ 的特征值和特征向量。

> **原题记号说明**：原题首句将算子记为 $A$，表达式中又沿用微分算子习惯记作 $\mathcal{L}u(x)$，二者指同一算子（即著名的 **Legendre 微分算子**）：
> $$(Au)(x) = \mathcal{L}u(x) = -\frac{\mathrm{d}}{\mathrm{d}x}\left[(1 - x^2)u'(x)\right] = -(1 - x^2)u''(x) + 2x u'(x), \qquad x \in [-1, 1].$$

#### 【零基础破题思路】

这道题要找所有的标量 $\lambda$ 和**非零**函数 $u \in C^2([-1,1])$，使得 $Au = \lambda u$，也就是解二阶变系数微分方程
$$-(1-x^2)u''(x) + 2x u'(x) = \lambda u(x), \qquad x \in [-1,1].$$

很多同学一看到变系数微分方程就懵了，其实这道题用**两步极其干净的泛函 + 线代组合拳**，连一丁点儿高深奇点理论都不需要，就能 100% 严密拿下：

1. **第一步（先找出“一批”特征值和特征函数）**：试试把最简单的单项式 $x^k$ 喂给 $A$！算一下立刻惊喜地发现：
   $$A(x^k) = k(k+1)x^k - k(k-1)x^{k-2}.$$
   次数不但没有变高，还只往低次掉！这说明对每个 $n \ge 0$，$n$ 次及以下多项式空间 $\mathcal{P}_n = \operatorname{span}\{1, x, \dots, x^n\}$ 是 $A$ 的**有限维不变子空间**，且在基 $\{1, x, \dots, x^n\}$ 下矩阵是**上三角矩阵**，对角线上明晃晃地写着 $0, 2, 6, \dots, n(n+1)$！
   于是每个 $\lambda_n = n(n+1)$（$n = 0, 1, 2, \dots$）铁定是特征值，而且对应一个唯一的 $n$ 次多项式特征函数 $P_n(x)$（就是大名鼎鼎的 **Legendre 多项式**，既可以用递推系数写出来，也可以用 Rodrigues 公式写出来并 3 行验完！）。
2. **第二步（证明“全找齐了”——再无其他特征值，且每个特征空间只有 1 维）**：万一还有某个非多项式函数 $u \in C^2([-1,1])$ 也是特征函数怎么办？
   这时**泛函分析的主角登场**：把 $C^2([-1,1])$ 放进 Hilbert 空间 $L^2(-1,1)$ 里！
   - 做两次分部积分，发现由于权函数 $(1-x^2)$ 在端点 $x = \pm 1$ **自己变成 $0$**，边界项自动消失，所以 $A$ 在 $L^2(-1,1)$ 上竟然**天然对称**：$(Au, v) = (u, Av)$！
   - 于是任何满足 $Au = \lambda u$ 的函数 $u$ 都自动满足 $[\lambda - n(n+1)](u, P_n) = 0$（$\forall n \ge 0$）。
   - 如果 $\lambda \neq n(n+1)$，那 $u$ 就与**所有多项式 $P_n$ 都正交**！而由 **Weierstrass 多项式逼近定理**，全体多项式在 $C([-1,1])$（从而在 $L^2(-1,1)$）里稠密——与稠密集正交的连续函数只能恒等于 $0$！
   - 同理，如果 $\lambda = n_0(n_0+1)$，把 $u$ 减去它在 $P_{n_0}$ 方向上的投影 $c P_{n_0}$，剩下的余项 $w = u - c P_{n_0}$ 与所有 $P_n$ 正交，立刻逼出 $w \equiv 0$，即 $u = c P_{n_0}$！

#### 【规范解答】

**结论**：算子 $A$ 的全部特征值为
$$\lambda_n = n(n + 1), \qquad n = 0, 1, 2, \dots,$$
对应于特征值 $\lambda_n = n(n+1)$ 的全部特征向量（特征函数）为
$$u(x) = c P_n(x), \qquad c \neq 0,$$
其中 $P_n(x)$ 为 $n$ 次 **Legendre 多项式**，其显式表达式（Rodrigues 公式）为
$$P_n(x) = \frac{1}{2^n n!} \frac{\mathrm{d}^n}{\mathrm{d}x^n}\left[(x^2 - 1)^n\right], \qquad n = 0, 1, 2, \dots$$
（前四个为：$P_0(x) = 1$，$P_1(x) = x$，$P_2(x) = \dfrac{1}{2}(3x^2 - 1)$，$P_3(x) = \dfrac{1}{2}(5x^3 - 3x)$。）

下面分两步严格证明上述结论（空间默认在复数域 $\mathbb{C}$ 上考虑；若限定实数域 $\mathbb{R}$，结论及证明完全相同）。

---

**第一步：证明每个 $\lambda_n = n(n+1)$（$n = 0, 1, 2, \dots$）都是 $A$ 的特征值，且 $P_n(x)$ 是对应的特征向量。**

先计算算子 $A$ 作用在单项式 $x^k \in C^2([-1,1])$（$k = 0, 1, 2, \dots$）上的结果：
- 当 $k = 0$ 时，$A(1) = 0 = 0 \cdot (0+1) \cdot 1$；
- 当 $k = 1$ 时，$A(x) = -\dfrac{\mathrm{d}}{\mathrm{d}x}(1 - x^2) = 2x = 1 \cdot (1+1) x$；
- 当 $k \ge 2$ 时，
  $$A(x^k) = -\frac{\mathrm{d}}{\mathrm{d}x}\left[(1 - x^2) k x^{k-1}\right] = -\frac{\mathrm{d}}{\mathrm{d}x}\left(k x^{k-1} - k x^{k+1}\right) = k(k + 1)x^k - k(k - 1)x^{k-2}.$$

对任意非负整数 $n$，记次数不超过 $n$ 的多项式子空间为 $\mathcal{P}_n = \operatorname{span}\{1, x, \dots, x^n\} \subset C^2([-1,1])$。由上式可知 $A(\mathcal{P}_n) \subset \mathcal{P}_n$，即 $\mathcal{P}_n$ 是 $A$ 的 $(n+1)$ 维不变子空间。在基 $\{1, x, x^2, \dots, x^n\}$ 下，限制算子 $A|_{\mathcal{P}_n}$ 的表示矩阵为**上三角矩阵**（主对角线元素依次为 $0, 2, 6, \dots, n(n+1)$）：
$$M_n = \begin{pmatrix} 0 & 0 & -2 & 0 & \cdots & \ast \\ 0 & 2 & 0 & -6 & \cdots & \ast \\ 0 & 0 & 6 & 0 & \cdots & \ast \\ 0 & 0 & 0 & 12 & \cdots & \ast \\ \vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & 0 & 0 & \cdots & n(n+1) \end{pmatrix}.$$

由于主对角线上的 $(n+1)$ 个数 $k(k+1)$（$k = 0, 1, \dots, n$）两两互异，故对每个 $n \ge 0$，$\lambda_n = n(n+1)$ 都是 $A|_{\mathcal{P}_n}$（从而也是 $A$）的特征值，且在 $\mathcal{P}_n$ 中存在唯一（相差非零常数倍）的 $n$ 次多项式特征向量。

下面直接验证 Rodrigues 公式给出的 $n$ 次多项式 $P_n(x) = \dfrac{1}{2^n n!} \dfrac{\mathrm{d}^n}{\mathrm{d}x^n}\left[(x^2-1)^n\right] \in C^2([-1,1])$ 确为对应于 $\lambda_n = n(n+1)$ 的非零特征向量：

1. **非零性与次数**：$(x^2 - 1)^n = x^{2n} + \text{低次项}$，求 $n$ 阶导数后首项为 $\dfrac{(2n)!}{n!} x^n$，故 $P_n(x)$ 是精确的 $n$ 次多项式（特别地 $P_n \not\equiv 0$；按 Leibniz 公式在 $x=1$ 处计算还可得归一化 $P_n(1)=1$）。
2. **满足特征方程**：记 $f(x) = (x^2 - 1)^n$，直接求导得 $f'(x) = 2nx(x^2-1)^{n-1}$，两边同乘 $(x^2 - 1)$ 得恒等式
   $$(x^2 - 1)f'(x) = 2n x f(x).$$
   对上式两边应用 Leibniz 乘积求导公式求 $(n+1)$ 阶导数：
   - 左边：
     $$(x^2 - 1)f^{(n+2)}(x) + \binom{n+1}{1}(2x)f^{(n+1)}(x) + \binom{n+1}{2}(2)f^{(n)}(x) = (x^2 - 1)f^{(n+2)}(x) + 2(n+1)x f^{(n+1)}(x) + n(n+1)f^{(n)}(x);$$
   - 右边：
     $$2n\left[x f^{(n+1)}(x) + \binom{n+1}{1} \cdot 1 \cdot f^{(n)}(x)\right] = 2nx f^{(n+1)}(x) + 2n(n+1)f^{(n)}(x).$$
   两边相等并移项整理，得
   $$(x^2 - 1)f^{(n+2)}(x) + 2x f^{(n+1)}(x) - n(n+1)f^{(n)}(x) = 0.$$
   两边同乘 $\dfrac{1}{2^n n!}$ 并注意到 $P_n(x) = \dfrac{1}{2^n n!} f^{(n)}(x)$，立得
   $$-(1 - x^2)P_n''(x) + 2x P_n'(x) = n(n+1)P_n(x), \quad \text{即 } A P_n = n(n+1)P_n.$$
   （注：亦可设 $P_n(x) = \sum_{k=0}^n a_k x^k$ 代入 $AP_n = n(n+1)P_n$ 比较系数，得到 $a_{n-1}=0$ 及递推公式 $a_k = -\dfrac{(k+2)(k+1)}{(n-k)(n+k+1)} a_{k+2}$ 唯一确定各项系数。）

---

**第二步：证明 $A$ 再无其他特征值，且每个 $\lambda_n$ 对应的特征空间恰为一维空间 $\operatorname{span}\{P_n\}$。**

将 $C^2([-1,1])$ 与 $C([-1,1])$ 嵌入复 Hilbert 空间 $L^2(-1, 1)$，内积与范数定义为
$$(u, v) = \int_{-1}^{1} u(x)\overline{v(x)}\,\mathrm{d}x, \qquad \|u\|_{L^2} = \left(\int_{-1}^{1} |u(x)|^2\,\mathrm{d}x\right)^{1/2}.$$

1. **验证 $A$ 的对称性及 $\{P_n\}$ 的正交性**：
   对任意 $u, v \in C^2([-1,1])$，由于 $(1 - x^2)\big|_{x=-1}^{x=1} = 0$，连续分部积分两次得边界项全部为零：
   $$\begin{aligned}
   (Au, v) &= -\int_{-1}^{1} \left[(1 - x^2)u'(x)\right]' \overline{v(x)}\,\mathrm{d}x \\
   &= -\Big[(1 - x^2)u'(x)\overline{v(x)}\Big]_{-1}^{1} + \int_{-1}^{1} (1 - x^2)u'(x)\overline{v'(x)}\,\mathrm{d}x \\
   &= \int_{-1}^{1} (1 - x^2)u'(x)\overline{v'(x)}\,\mathrm{d}x \\
   &= \Big[(1 - x^2)u(x)\overline{v'(x)}\Big]_{-1}^{1} - \int_{-1}^{1} u(x)\overline{\left[(1 - x^2)v'(x)\right]'}\,\mathrm{d}x = (u, Av).
   \end{aligned}$$
   特别地，对任意 $m \neq n$，由 $m(m+1)(P_m, P_n) = (AP_m, P_n) = (P_m, AP_n) = n(n+1)(P_m, P_n)$ 及 $m(m+1) \neq n(n+1)$，得
   $$(P_m, P_n) = 0 \qquad (m \neq n),$$
   且因为 $P_n$ 连续非零，有 $\|P_n\|_{L^2} > 0$。

2. **证明多项式系 $\{P_n\}_{n=0}^\infty$ 的完备性（唯一性引理）**：
   设连续函数 $w \in C([-1,1])$ 满足 $(w, P_n) = 0$ 对一切 $n = 0, 1, 2, \dots$ 成立。
   因为每个 $P_n$ 恰为 $n$ 次多项式，线性张成 $\operatorname{span}\{P_0, P_1, \dots\}$ 恰好是 $[-1,1]$ 上的全体多项式，故对一切多项式 $p(x)$ 均有 $(w, p) = 0$。
   由 **Weierstrass 多项式逼近定理**（对复值函数 $w$，分别对其连续的实部 $\operatorname{Re} w$ 与虚部 $\operatorname{Im} w$ 作实多项式一致逼近再合并），存在多项式列 $\{p_k\}$ 在 $[-1,1]$ 上一致收敛于 $w$。于是
   $$\|w - p_k\|_{L^2}^2 = \int_{-1}^{1} |w(x) - p_k(x)|^2\,\mathrm{d}x \le 2 \max_{x \in [-1,1]} |w(x) - p_k(x)|^2 \to 0 \qquad (k \to \infty).$$
   由内积的连续性得
   $$\|w\|_{L^2}^2 = (w, w) = \lim_{k\to\infty} (w, p_k) = 0.$$
   因 $w \in C([-1,1])$ 连续且 $\int_{-1}^1 |w(x)|^2\,\mathrm{d}x = 0$，必有 $w(x) \equiv 0$（对一切 $x \in [-1,1]$）。

3. **确定全部特征值与特征空间**：
   现设 $\lambda \in \mathbb{C}$ 与 $u \in C^2([-1,1])$ 满足 $Au = \lambda u$。对一切 $n = 0, 1, 2, \dots$，利用 $A$ 的对称性以及 $n(n+1) \in \mathbb{R}$，有
   $$\lambda (u, P_n) = (\lambda u, P_n) = (Au, P_n) = (u, A P_n) = \left(u, n(n+1)P_n\right) = n(n+1)(u, P_n),$$
   整理即得
   $$\left[\lambda - n(n+1)\right](u, P_n) = 0, \qquad \forall n = 0, 1, 2, \dots$$
   分两种情形讨论：
   - **情形 A**：$\lambda \notin \{n(n+1) \mid n = 0, 1, 2, \dots\}$。
     此时对一切 $n \ge 0$ 都有 $\lambda - n(n+1) \neq 0$，从而 $(u, P_n) = 0$（$\forall n \ge 0$）。由上述完备性引理立得 $u(x) \equiv 0$。故 $\lambda$ 不是特征值。
   - **情形 B**：$\lambda = n_0(n_0+1)$（其中 $n_0$ 为某个非负整数）。
     此时对一切 $n \neq n_0$ 都有 $(u, P_n) = 0$。令 $c = \dfrac{(u, P_{n_0})}{\|P_{n_0}\|_{L^2}^2}$，构造 $w = u - c P_{n_0} \in C^2([-1,1])$。利用正交性 $(P_{n_0}, P_n) = 0$（$n \neq n_0$）：
     - 当 $n \neq n_0$ 时，$(w, P_n) = (u, P_n) - c(P_{n_0}, P_n) = 0 - 0 = 0$；
     - 当 $n = n_0$ 时，$(w, P_{n_0}) = (u, P_{n_0}) - \dfrac{(u, P_{n_0})}{\|P_{n_0}\|_{L^2}^2}(P_{n_0}, P_{n_0}) = 0$。
     故 $(w, P_n) = 0$ 对一切 $n \ge 0$ 成立。再由上述完备性引理得 $w(x) \equiv 0$，即 $u(x) = c P_{n_0}(x)$。当 $u$ 为特征向量时要求 $u \not\equiv 0$，即 $c \neq 0$。

综上，$A$ 的全部特征值为 $\lambda_n = n(n+1)$（$n = 0, 1, 2, \dots$），对应的全部特征向量为 $c P_n(x)$（$c \neq 0$）。$\blacksquare$

#### 【考点复盘与避坑提示】

- **为什么题目没有给任何边界条件（如 $u(\pm 1)=0$），$A$ 却自动对称且选出了离散的特征值？**
  因为系数 $(1-x^2)$ 在端点 $x = \pm 1$ 自动退化为零（数理方程里称为**奇异 Sturm-Liouville 问题**），只要函数要求落在 $C^2([-1,1])$（端点有限且二阶连续可微），边界项就被 $(1-x^2)$ 自动抹平了！
- **与数理方程中“幂级数解法”的融会贯通**：
  学过《数学物理方法》的同学可能记得另一种做法：把 $u(x) = \sum_{k=0}^\infty a_k x^k$ 代入方程得两步递推 $a_{k+2} = \frac{k(k+1)-\lambda}{(k+2)(k+1)} a_k$；当 $\lambda \neq n(n+1)$ 时，级数永不截断，由比值判别法可证其发散于端点 $x = \pm 1$，从而不在 $C^2([-1,1])$ 中（实际上连 $C([-1,1])$ 都不在）。但幂级数法在考试中需要细致证明任意 $C^2([-1,1])$ 解在 $(-1,1)$ 上解析且在端点处发散，很容易漏掉严格性细节；相比之下，上面的**“多项式上三角 + $L^2$ 对称性 + Weierstrass 稠密性”**纯用泛函工具，无任何死角，是最稳妥的满分写法！
- **避坑提醒**：千万别只写完第一步“找到 $P_n(x)$”就收工！题干问的是“求出 $A$ 的特征值和特征向量”，必须证明**没有遗漏其他特征值**，且**每个特征空间都没有多项式以外的其他特征函数**（即第二步不可或缺！）。

---

### 习题 7 第 4 题（教材 p. 293）

#### 【原题】

设 $H = L^2(\mathbb{R})$，算子 $A$ 定义如下：
$$\begin{cases}
\mathcal{D}(A) = \left\{x \in H \;\middle|\; \displaystyle\int_{\mathbb{R}} t^2 |x(t)|^2\,\mathrm{d}t < +\infty\right\}, \\
(Ax)(t) = t x(t), \quad \forall x \in \mathcal{D}(A).
\end{cases}$$
证明：$A$ 是一个闭的无界算子。

#### 【零基础破题思路】

这就是贯穿全章的明星算子——量子力学里的**位置算子**（教材例 7.1、例 7.3、命题 7.1–7.2 的主角）。要证明三件事（其实第一件是隐含前提）：

1. **线性与良定义**（快速补全，显得严谨）：确认 $\mathcal{D}(A)$ 是线性子空间且 $Ax \in L^2(\mathbb{R})$。
2. **无界性**：要推翻“存在常数 $M$ 使 $\|Ax\| \le M\|x\|$”，只需造出一列“自己大小固定为 $1$，但被 $A$ 作用后大小爆到 $+\infty$”的测试函数！
   - $A$ 的本质是乘以坐标 $t$，哪里的 $t$ 最大？当然是远处 $[n, n+1]$！
   - 取示性函数 $x_n = \chi_{[n, n+1]}$（宽度 $1$、高度 $1$ 的小方块），$\|x_n\|_2 = 1$，而乘上 $t$ 后至少放大 $n$ 倍（$\|Ax_n\|_2 \ge n \to +\infty$），瞬间秒杀！
3. **闭性**：回忆闭算子的定义：
   $$\text{若 } x_n \in \mathcal{D}(A),\ x_n \xrightarrow{L^2} x,\ Ax_n \xrightarrow{L^2} y,\quad \text{需证 } x \in \mathcal{D}(A) \text{ 且 } Ax = y.$$
   - 怎么把“积分意义下的 $L^2$ 收敛”变成“逐点等式 $y(t) = tx(t)$”？用实变函数最经典的利器：**$L^2$ 收敛必含几乎处处（a.e.）逐点收敛的子列**！
   - **细节杀手锏**：不要一次糊弄说“取子列同时 a.e. 收敛”，而要**显式分两次抽子列**——先从 $\{x_n\}$ 抽子列 $\{x_{n_k}\}$ 使 $x_{n_k} \to x$ a.e.；再从对应的 $\{Ax_{n_k}\}$ 中**进一步抽子列** $\{Ax_{n_{k_j}}\}$ 使 $t x_{n_{k_j}} \to y$ a.e.！同一子列同时满足 $t x_{n_{k_j}} \to tx$ a.e. 和 $t x_{n_{k_j}} \to y$ a.e.，立得 $y(t) = t x(t)$ a.e.！
   - **临门一脚千万别忘**：利用已知的 $y \in L^2(\mathbb{R})$，得 $\int_{\mathbb{R}} t^2|x(t)|^2\,\mathrm{d}t = \int_{\mathbb{R}} |y(t)|^2\,\mathrm{d}t < +\infty$，这才真正证完 **$x \in \mathcal{D}(A)$** 且 $Ax = y$！

#### 【规范解答】

**第一步：验证 $A$ 是良定义的线性算子。**

显然 $0 \in \mathcal{D}(A)$。任取 $x, y \in \mathcal{D}(A)$ 及 $\alpha, \beta \in \mathbb{C}$，利用基本不等式 $|a + b|^2 \le 2(|a|^2 + |b|^2)$，有
$$\int_{\mathbb{R}} t^2 |\alpha x(t) + \beta y(t)|^2\,\mathrm{d}t \le 2|\alpha|^2 \int_{\mathbb{R}} t^2 |x(t)|^2\,\mathrm{d}t + 2|\beta|^2 \int_{\mathbb{R}} t^2 |y(t)|^2\,\mathrm{d}t < +\infty,$$
故 $\alpha x + \beta y \in \mathcal{D}(A)$，即 $\mathcal{D}(A)$ 是 $H = L^2(\mathbb{R})$ 的线性子空间（事实上 $\mathcal{D}(A)$ 包含全体紧支集连续函数 $C_c(\mathbb{R})$，故还是稠密子空间）。又对任意 $x \in \mathcal{D}(A)$，
$$\|Ax\|^2 = \int_{\mathbb{R}} |t x(t)|^2\,\mathrm{d}t = \int_{\mathbb{R}} t^2 |x(t)|^2\,\mathrm{d}t < +\infty,$$
故 $Ax \in L^2(\mathbb{R})$，且显然 $A(\alpha x + \beta y) = \alpha Ax + \beta Ay$，因此 $A: \mathcal{D}(A) \to H$ 是良定义的线性算子。

---

**第二步：证明 $A$ 是无界算子。**

对每个正整数 $n = 1, 2, \dots$，构造区间 $[n, n+1]$ 上的示性函数：
$$x_n(t) = \chi_{[n, n+1]}(t) = \begin{cases} 1, & t \in [n, n+1], \\ 0, & t \notin [n, n+1]. \end{cases}$$
直接计算它们的 $L^2$ 范数：
$$\|x_n\|^2 = \int_{\mathbb{R}} |x_n(t)|^2\,\mathrm{d}t = \int_n^{n+1} 1\,\mathrm{d}t = 1 \implies \|x_n\| = 1,$$
$$\int_{\mathbb{R}} t^2 |x_n(t)|^2\,\mathrm{d}t = \int_n^{n+1} t^2\,\mathrm{d}t = \frac{(n+1)^3 - n^3}{3} = n^2 + n + \frac{1}{3} < +\infty.$$
由此可知每个 $x_n \in \mathcal{D}(A)$ 且 $\|x_n\| = 1$，但
$$\|Ax_n\| = \sqrt{n^2 + n + \frac{1}{3}} \ge n \to +\infty \qquad (n \to \infty).$$
因此不存在常数 $M > 0$ 使得 $\|Ax\| \le M\|x\|$ 对一切 $x \in \mathcal{D}(A)$ 成立，故 $A$ 是**无界算子**。

---

**第三步：证明 $A$ 是闭算子。**

任取序列 $\{x_n\}_{n=1}^\infty \subset \mathcal{D}(A)$ 以及元素 $x, y \in L^2(\mathbb{R})$，使得在 $L^2(\mathbb{R})$ 范数意义下
$$\lim_{n\to\infty} \|x_n - x\| = 0, \qquad \lim_{n\to\infty} \|Ax_n - y\| = 0.$$
需要证明：$x \in \mathcal{D}(A)$ 且 $Ax = y$。

由 $L^2(\mathbb{R})$ 空间的基本性质（Riesz 定理：$L^2$ 强收敛序列必存在几乎处处逐点收敛的子列），分两步抽取子列：

1. 因为 $x_n \xrightarrow{L^2} x$，存在子列 $\{x_{n_k}\}_{k=1}^\infty$ 使得
   $$\lim_{k\to\infty} x_{n_k}(t) = x(t), \qquad \text{a.e. } t \in \mathbb{R}.$$
2. 因为原序列满足 $Ax_n \xrightarrow{L^2} y$，其子列 $\{Ax_{n_k}\}_{k=1}^\infty$ 依然在 $L^2(\mathbb{R})$ 中强收敛于 $y$，故可从 $\{n_k\}$ 中**进一步抽取子列** $\{n_{k_j}\}_{j=1}^\infty$ 使得
   $$\lim_{j\to\infty} (Ax_{n_{k_j}})(t) = \lim_{j\to\infty} t x_{n_{k_j}}(t) = y(t), \qquad \text{a.e. } t \in \mathbb{R}.$$

另一方面，作为 $\{x_{n_k}\}$ 的子列，$\{x_{n_{k_j}}\}$ 仍满足 $\lim_{j\to\infty} x_{n_{k_j}}(t) = x(t)$（a.e. $t \in \mathbb{R}$）。对每个使该式成立的有限实数 $t \in \mathbb{R}$，两端同乘 $t$ 即得
$$\lim_{j\to\infty} t x_{n_{k_j}}(t) = t x(t), \qquad \text{a.e. } t \in \mathbb{R}.$$

两个零测集之并仍为零测集。由逐点极限的唯一性，得
$$y(t) = t x(t), \qquad \text{a.e. } t \in \mathbb{R}.$$

最后验证 $x \in \mathcal{D}(A)$ 且 $Ax = y$：由于已知 $y \in L^2(\mathbb{R})$，代入上式有
$$\int_{\mathbb{R}} t^2 |x(t)|^2\,\mathrm{d}t = \int_{\mathbb{R}} |t x(t)|^2\,\mathrm{d}t = \int_{\mathbb{R}} |y(t)|^2\,\mathrm{d}t = \|y\|^2 < +\infty.$$
这就证明了 $x \in \mathcal{D}(A)$，进而 $(Ax)(t) = t x(t) = y(t)$（a.e. $t \in \mathbb{R}$），即在 $L^2(\mathbb{R})$ 中 $Ax = y$。

综上所述，$A$ 是 $L^2(\mathbb{R})$ 上的一个闭的无界算子。$\blacksquare$

#### 【考点复盘与避坑提示】

- **扣分高发区 1（抽子列不严谨）**：切忌写“$x_n \to x$ a.e. 且 $Ax_n \to y$ a.e.”——$L^2$ 收敛本身**不能**保证原序列几乎处处收敛！也不能各自独立抽两个不同的子列 $\{x_{n_k}\}$ 和 $\{Ax_{m_k}\}$！必须**嵌套两次抽子列**（从第一个子列里再抽子列），才能保证同一个下标序列 $\{n_{k_j}\}$ 上两条几乎处处极限同时成立。
- **扣分高发区 2（漏证 $x \in \mathcal{D}(A)$）**：证明闭算子有两个硬指标：① **$x \in \mathcal{D}(A)$**；② **$Ax = y$**。本题中恰恰因为 $y \in L^2(\mathbb{R})$ 且 $t x(t) = y(t)$ a.e.，才把 $\int t^2|x(t)|^2\mathrm{d}t < +\infty$ 给锁死，千万不能跳过这一行！
- **另解拓展（与教材例 7.3 + 定理 7.1 联动）**：也可以先证明 $A$ 稠定且自伴（$A = A^*$，见教材例 7.3，p. 277），再由**定理 7.1**（自伴算子必是闭算子）直接推出 $A$ 是闭算子；但直接用序列定义最简短自包含，考场首选直接法！

---

### 习题 7 第 6 题（教材 p. 293）

#### 【原题】

设 $A$ 是 $H$ 上的对称算子。证明：
- (i) $A$ 是闭算子当且仅当 $A = A^{**} \subset A^*$；
- (ii) $A$ 是自伴算子当且仅当 $A = A^{**} = A^*$。

#### 【零基础破题思路】

这道题篇幅很短，却是检验你是否真正搞懂第 7 章概念逻辑链的“照妖镜”！

看到题目里出现了二阶共轭 $A^{**} = (A^*)^*$，零基础最容易直接埋头套定理，却漏掉一个最要命的**前提检查**：
- **前置暗坑**：凭什么 $A^{**}$ 一定存在？按定义 7.3，只有**稠定算子**才有共轭算子！所以必须先花两行证明：$A$ 对称 $\implies A$ 稠定且 $A \subset A^* \implies \mathcal{D}(A) \subset \mathcal{D}(A^*) \implies \mathcal{D}(A^*)$ 也在 $H$ 中稠密 $\implies A^{**} = (A^*)^*$ 良定义存在！

迈过这个坑之后，(i) 和 (ii) 就像搭乐高一样，精准对接教材 **定理 7.1**（稠定算子的共轭必闭）和 **定理 7.2**（稠定闭算子满足 $A = A^{**}$）：
- **(i)**
  - $\Rightarrow$：$A$ 是稠定闭算子 $\xrightarrow{\text{定理 7.2}} A = A^{**}$；又 $A$ 对称 $\implies A \subset A^*$；拼起来就是 $A = A^{**} \subset A^*$！
  - $\Leftarrow$：已知 $A = A^{**}$；由于 $A^*$ 稠定，把**定理 7.1** 用在 $A^*$ 身上，立刻得出 $(A^*)^* = A^{**}$ 是闭算子，所以 $A$ 就是闭算子！
- **(ii)**
  - $\Rightarrow$：$A$ 自伴即 $A = A^*$ $\xrightarrow{\text{定理 7.1}} A$ 是闭算子 $\xrightarrow{\text{(i)}} A = A^{**}$，合并即 $A = A^{**} = A^*$！
  - $\Leftarrow$：题设直接送了你 $A = A^*$，按定义 7.4 这就是自伴算子！

#### 【规范解答】

**前置验证（证明 $A^{**}$ 良定义存在）：**

因为 $A$ 是 $H$ 上的对称算子，由定义 7.4（p. 275）知，$A$ 是 $H$ 上的稠定线性算子（即 $\overline{\mathcal{D}(A)} = H$），且满足 $A \subset A^*$。由算子扩张的定义（定义 7.1 (iv)）可知 $\mathcal{D}(A) \subset \mathcal{D}(A^*)$，于是
$$H = \overline{\mathcal{D}(A)} \subset \overline{\mathcal{D}(A^*)} \subset H,$$
故 $\overline{\mathcal{D}(A^*)} = H$。结合引理 7.2 (i) 知 $A^*$ 是 $H$ 上的稠定线性算子。因此，按定义 7.3，$A^*$ 的共轭算子 $A^{**} \triangleq (A^*)^*$ 良定义且唯一存在。

---

**证明 (i)：$A$ 是闭算子 $\iff A = A^{**} \subset A^*$。**

- **必要性（$\Rightarrow$）**：
  设 $A$ 是闭算子。因为 $A$ 是对称算子，所以 $A$ 是 $H$ 上的稠定闭线性算子。由**定理 7.2**（p. 280）立得 $A = A^{**}$。又因为 $A$ 是对称算子，由定义 7.4 知 $A \subset A^*$。两式合并即得
  $$A = A^{**} \subset A^*.$$
- **充分性（$\Leftarrow$）**：
  设 $A = A^{**} \subset A^*$。由前置验证知 $A^*$ 是 $H$ 上的稠定线性算子，对 $A^*$ 应用**定理 7.1**（p. 280，稠定线性算子的共轭算子必为闭算子），可知 $(A^*)^* = A^{**}$ 是闭算子。因为 $A = A^{**}$，所以 $A$ 是闭算子。

---

**证明 (ii)：$A$ 是自伴算子 $\iff A = A^{**} = A^*$。**

- **必要性（$\Rightarrow$）**：
  设 $A$ 是自伴算子，即按定义 7.4 有 $A = A^*$。由**定理 7.1**（p. 280）知，自伴算子必为闭算子。再由刚刚证明的 (i)（或直接对 $A = A^*$ 两边取共轭得 $A^{**} = (A^*)^* = A^* = A$）立得
  $$A = A^{**} = A^*.$$
- **充分性（$\Leftarrow$）**：
  设 $A = A^{**} = A^*$，则直接有 $A = A^*$。又已知 $A$ 是稠定线性算子，按定义 7.4 即知 $A$ 是自伴算子。$\blacksquare$

#### 【考点复盘与避坑提示】

- **加分亮点（一定要写的前置步）**：很多同学一上来就用“对 $A^*$ 用定理 7.1 知 $A^{**}$ 闭”，却忘了先说明 **为什么 $A^*$ 是稠定的**！对一般稠定算子 $T$，$T^*$ 并不一定稠定；本题之所以能无障碍地写 $A^{**}$，全靠题设里 **$A$ 是对称算子** 给出的 $\mathcal{D}(A) \subset \mathcal{D}(A^*)$！
- **彻底打通对称算子的“黄金包含链”与四级阶梯**：
  实际上，对**任意**对称算子 $A$（哪怕不闭），都自动成立
  $$A \subset \bar{A} = A^{**} \subset A^*.$$
  （速记理由：① 对任意 $x \in \mathcal{D}(A), y \in \mathcal{D}(A^*)$，有 $(A^* y, x) = \overline{(x, A^* y)} = \overline{(Ax, y)} = (y, Ax)$，直接推出 $A \subset A^{**}$；② 对 $A \subset A^*$ 两边取共轭并应用扩张反号律引理 7.2 (iii)，立得 $A^{**} \subset A^*$！）
  于是对称算子的全部类型一目了然：

| 算子层级 | 包含链状态 | 典型教材例子 |
| :--- | :--- | :--- |
| **1. 一般（非闭）对称算子** | $A \subsetneq \bar{A} = A^{**} \subsetneq A^*$ | 把例 7.2 的 $A_3$ 再缩小到 $C_0^\infty(0,1)$ |
| **2. 闭对称算子（非自伴）** | $A = A^{**} \subsetneq A^*$ | 例 7.2 中的 $A_3$（$A_3 = A_3^{**} \subsetneq A_1 = A_3^*$） |
| **3. 本质自伴算子（闭包自伴）** | $A \subsetneq \bar{A} = A^{**} = A^*$ | 把乘法算子 $A$ 限制在紧支光滑函数 $C_c^\infty(\mathbb{R})$ 上 |
| **4. 自伴算子** | $A = A^{**} = A^*$ | 例 7.2 的 $A_2$、例 7.3 的位置算子、例 7.4 的动量算子 |

---

### 习题 7 第 7 题（教材 p. 294）

#### 【原题】

设 $A$ 是 $H$ 上的对称算子，而且对任意 $x \in \mathcal{D}(A)$ 均有 $(Ax, x) \geqslant 0$。证明：
- (i) $\| (A + I)x \|^2 \geqslant \|x\|^2 + \|Ax\|^2$；
- (ii) $A$ 是闭算子的充要条件是 $\mathcal{R}(A + I)$ 是闭集。

#### 【零基础破题思路】

这道题是对教材 **式 (7.9) 与引理 7.8（p. 281）** 的超级经典“双胞胎变形题”！

1. **先看第 (i) 问**：直接把左边的范数平方按内积展开：
   $$\|(A+I)x\|^2 = (Ax + x, Ax + x) = \|Ax\|^2 + (Ax, x) + (x, Ax) + \|x\|^2.$$
   因为 $A$ 对称，所以 $(x, Ax) = (Ax, x)$；又题目给了**非负条件** $(Ax, x) \ge 0$，于是中间两项合起来 $2(Ax, x) \ge 0$，直接扔掉就得到
   $$\|(A+I)x\|^2 \ge \|x\|^2 + \|Ax\|^2!$$
   立刻收获两把万能扳手（**双控不等式**）：
   $$\|x\| \le \|(A+I)x\|, \qquad \|Ax\| \le \|(A+I)x\|.$$
   这两把扳手的威力是：**只要 $(A+I)x_n$ 乖乖收敛（是 Cauchy 列），$x_n$ 和 $Ax_n$ 就都被死死按住，同时变成 Cauchy 列！**

2. **再看第 (ii) 问（充要条件两个方向）**：
   - **必要性（$A$ 闭 $\implies \mathcal{R}(A+I)$ 闭）**：套路极稳！取 $y_n = (A+I)x_n \to y$，目标证明 $y \in \mathcal{R}(A+I)$。
     - $\{y_n\}$ 收敛 $\implies \{y_n\}$ 是 Cauchy 列；
     - 把 (i) 用在差向量 $x_m - x_n$ 上：$\|x_m - x_n\|^2 + \|Ax_m - Ax_n\|^2 \le \|y_m - y_n\|^2 \to 0$；
     - 立刻推出 $\{x_n\}$ 和 $\{Ax_n\}$ 都是 $H$ 中的 Cauchy 列；
     - 由 Hilbert 空间 $H$ 完备，$x_n \to x$ 且 $Ax_n \to z$；
     - 由 $A$ 是闭算子，$x \in \mathcal{D}(A)$ 且 $Ax = z$；
     - 于是 $y = \lim(Ax_n + x_n) = Ax + x = (A+I)x \in \mathcal{R}(A+I)$！收工！
   - **充分性（$\mathcal{R}(A+I)$ 闭 $\implies A$ 闭）**：设 $x_n \in \mathcal{D}(A)$，$x_n \to x$ 且 $Ax_n \to z$，目标证明 $x \in \mathcal{D}(A)$ 且 $Ax = z$。
     - 相加立得 $(A+I)x_n = Ax_n + x_n \to z + x$；
     - 因为 $\mathcal{R}(A+I)$ 是闭集，极限 $z+x$ 必落回 $\mathcal{R}(A+I)$，也就是**存在某个 $x_0 \in \mathcal{D}(A)$ 使 $(A+I)x_0 = z + x$**！
     - 现在怎么把陌生的 $x_0$ 和目标 $x$ 画上等号？再次祭出扳手 $\|u\| \le \|(A+I)u\|$，用在 $x_n - x_0$ 上：
       $$\|x_n - x_0\| \le \|(A+I)(x_n - x_0)\| = \|(A+I)x_n - (z+x)\| \to 0!$$
     - 这说明 $x_n \to x_0$！而题目本来就有 $x_n \to x$，由极限唯一性立得 **$x = x_0 \in \mathcal{D}(A)$**！
     - 最后顺手算 $Ax = Ax_0 = (A+I)x_0 - x_0 = (z+x) - x = z$，大功告成！

#### 【规范解答】

注意 $\mathcal{D}(A+I) = \mathcal{D}(A) \cap H = \mathcal{D}(A)$。

---

**证明 (i)：**

对任意 $x \in \mathcal{D}(A)$，将范数平方按内积展开：
$$\|(A + I)x\|^2 = (Ax + x, Ax + x) = (Ax, Ax) + (Ax, x) + (x, Ax) + (x, x) = \|Ax\|^2 + (Ax, x) + (x, Ax) + \|x\|^2.$$
因为 $A$ 是对称算子，由引理 7.3（p. 277）得 $(x, Ax) = (Ax, x)$；又由题设对任意 $x \in \mathcal{D}(A)$ 有 $(Ax, x) \ge 0$，因此
$$\|(A + I)x\|^2 = \|Ax\|^2 + 2(Ax, x) + \|x\|^2 \ge \|x\|^2 + \|Ax\|^2.$$

特别地，由 (i) 立即得到对任意 $x \in \mathcal{D}(A)$ 均成立两条估计：
$$\|x\| \le \|(A + I)x\|, \qquad \|Ax\| \le \|(A + I)x\|. \tag{$\ast$}$$

---

**证明 (ii)：$A$ 是闭算子 $\iff \mathcal{R}(A + I)$ 是闭集。**

- **必要性（$\Rightarrow$）：设 $A$ 是闭算子，证明 $\mathcal{R}(A+I)$ 是闭集。**

  任取序列 $\{y_n\}_{n=1}^\infty \subset \mathcal{R}(A+I)$ 及 $y \in H$，满足 $y_n \to y$（$n \to \infty$）。我们需要证明 $y \in \mathcal{R}(A+I)$。

  由 $y_n \in \mathcal{R}(A+I)$，对每个 $n \ge 1$ 存在 $x_n \in \mathcal{D}(A)$ 使得 $y_n = (A+I)x_n$。因为 $\{y_n\}$ 在 $H$ 中收敛，故 $\{y_n\}$ 是 Cauchy 列。对任意 $m, n \ge 1$，由于 $x_m - x_n \in \mathcal{D}(A)$，将 (i) 应用于 $x_m - x_n$ 得
  $$\|x_m - x_n\|^2 + \|Ax_m - Ax_n\|^2 \le \|(A+I)(x_m - x_n)\|^2 = \|y_m - y_n\|^2 \to 0 \qquad (m, n \to \infty).$$
  这表明 $\{x_n\}$ 与 $\{Ax_n\}$ 均为 $H$ 中的 Cauchy 列。由 Hilbert 空间 $H$ 的完备性，存在 $x, z \in H$ 使得
  $$\lim_{n\to\infty} x_n = x, \qquad \lim_{n\to\infty} Ax_n = z.$$
  因为 $\{x_n\} \subset \mathcal{D}(A)$ 且 $A$ 是闭算子，由闭算子定义立得 $x \in \mathcal{D}(A)$ 且 $Ax = z$。于是
  $$y = \lim_{n\to\infty} y_n = \lim_{n\to\infty} (Ax_n + x_n) = z + x = Ax + x = (A + I)x \in \mathcal{R}(A + I).$$
  故 $\mathcal{R}(A + I)$ 包含其一切极限点，是 $H$ 中的闭集。

- **充分性（$\Leftarrow$）：设 $\mathcal{R}(A+I)$ 是闭集，证明 $A$ 是闭算子。**

  任取序列 $\{x_n\}_{n=1}^\infty \subset \mathcal{D}(A)$ 及 $x, z \in H$，满足
  $$\lim_{n\to\infty} x_n = x, \qquad \lim_{n\to\infty} Ax_n = z.$$
  我们需要证明 $x \in \mathcal{D}(A)$ 且 $Ax = z$。

  考察序列 $(A+I)x_n = Ax_n + x_n \in \mathcal{R}(A+I)$，显然有
  $$\lim_{n\to\infty} (A + I)x_n = \lim_{n\to\infty} (Ax_n + x_n) = z + x.$$
  因为 $\mathcal{R}(A+I)$ 是闭集，所以极限点 $z + x \in \mathcal{R}(A+I)$，即存在 $x_0 \in \mathcal{D}(A)$ 使得
  $$(A + I)x_0 = z + x.$$
  对任意 $n \ge 1$，由于 $x_n - x_0 \in \mathcal{D}(A)$，应用 ($\ast$) 式的第一条不等式得
  $$\|x_n - x_0\| \le \|(A + I)(x_n - x_0)\| = \|(A + I)x_n - (A + I)x_0\| = \|(A + I)x_n - (z + x)\| \to 0 \qquad (n \to \infty).$$
  因此 $\lim_{n\to\infty} x_n = x_0$。由极限的唯一性（已有 $\lim_{n\to\infty} x_n = x$），立得
  $$x = x_0 \in \mathcal{D}(A).$$
  进而由 $(A+I)x_0 = z + x$ 计算得
  $$Ax = Ax_0 = (A + I)x_0 - x_0 = (z + x) - x = z.$$
  故 $A$ 是闭算子。$\blacksquare$

#### 【考点复盘与避坑提示】

- **为什么教材引理 7.8 用 $A \pm \mathrm{i}I$，而这里能用 $A + I$？**
  对一般对称算子，实数内积 $(Ax, x)$ 可正可负，因此 $(A+I)x$ 展开里的 $2(Ax,x)$ 可能变成负数导致下界失效（甚至出现特征值 $-1$ 使得 $A+I$ 连单射都不是！）；所以一般对称算子必须靠虚数单位 $\pm\mathrm{i}$ 让交叉项严格抵消（式 7.9：$\|(A\pm\mathrm{i}I)x\|^2 = \|Ax\|^2 + \|x\|^2$）。而本题多了**非负条件 $(Ax,x) \ge 0$**，相当于把谱挡在了非负半轴 $[0, +\infty)$ 上，于是点 $-1$ 自动远离谱集，实数移位 $A+I$ 就拥有了同样强大的双控不等式 $\|(A+I)x\|^2 \ge \|x\|^2 + \|Ax\|^2$！
- **充分性证明的避坑提示**：证明“$\mathcal{R}(A+I)$ 闭 $\implies A$ 闭”时，得到存在 $x_0 \in \mathcal{D}(A)$ 使 $(A+I)x_0 = z + x$ 后，**绝不能直接假设 $x_0 = x$**（因为此时还不知道 $x \in \mathcal{D}(A)$，更不能写 $(A+I)x$！），必须用不等式 $\|x_n - x_0\| \le \|(A+I)(x_n - x_0)\| \to 0$ 证出 $x_n \to x_0$，再借极限唯一性把 $x_0$ 和 $x$ 焊接在一起！
- **考场拓展（和定理 7.3 的强强联合）**：若非负对称算子 $A$ 进一步满足 $\mathcal{R}(A+I) = H$，则因为 $A+I$ 也是对称算子且满射，由**定理 7.3 (ii)** 立得 $A+I$ 自伴（从而 $A$ 自伴！）。这是判断偏微分方程中正椭圆算子（如 $-\Delta$）是否自伴的超级实用判据！
