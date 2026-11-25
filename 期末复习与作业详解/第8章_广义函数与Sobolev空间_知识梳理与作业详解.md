# 第 8 章 广义函数与 Sobolev 空间：零基础知识梳理与作业详解

> **定位与使用指南**：本文档面向**此前没有学过泛函分析或偏微分方程**的本科生，严格对齐教材第 8 章（p.303–343）的记号、定义与定理编号，将「知识零基础精讲」「高频公式速查」与「全部 8 道必做作业逐题详解」合为一体。复习时建议先花 15 分钟读完**两大核心约定**与第 1 部分的核心直觉，再对照第 2 部分逐题刷完并背熟「规范解答」的得分动作。

---

## ⚠️ 必读前置：本书第 8 章两大核心约定（防网课/参考书混坑！）

不同泛函分析与偏微分方程教材在这一章有**两处最容易打架的记号差异**。做作业和期末考试前，请先锁死本书约定：

1. **常规广义函数的配对与乘法带复共轭（教材 p.304 例 8.2，p.305）**
   - 对局部可积函数 $g \in L^1_{\mathrm{loc}}(\Omega)$，其诱导的广义函数 $\lambda_g \in \mathcal{D}'(\Omega)$ 定义为
     $$\lambda_g(f) \triangleq \int_{\mathbb{R}^n} f(x)\overline{g(x)}\,\mathrm{d}x, \quad \forall f \in \mathcal{D}(\Omega).$$
     （之所以在 $g$ 上取共轭，是为了与复 Hilbert 空间 $L^2(\Omega)$ 的内积 $\langle f, g\rangle_{L^2} = \int f\bar{g}\,\mathrm{d}x$ 保持一致。当 $g$ 是**实值函数**时，$\overline{g} = g$，直接退化为 $\lambda_g(f) = \int g(x)f(x)\,\mathrm{d}x$。）
   - 对光滑函数 $h \in C^\infty(\Omega)$ 与广义函数 $\lambda \in \mathcal{D}'(\Omega)$，乘积 $h\lambda$ 定义为
     $$(h\lambda)(f) \triangleq \lambda(\bar{h}f), \quad \forall f \in \mathcal{D}(\Omega).$$
     这样对常规分布 $\lambda_g$ 刚好满足 $(h\lambda_g)(f) = \int \bar{h}f\bar{g}\,\mathrm{d}x = \int f\overline{hg}\,\mathrm{d}x = \lambda_{hg}(f)$，与普通函数乘法完全相容！

2. **Fourier 变换在指数上带 $2\pi$ 因子（教材 p.312 式 (8.9)）**
   - 对 $f \in L^1(\mathbb{R}^n)$，本书定义 Fourier 变换为
     $$\hat{f}(\xi) = (\mathcal{F}f)(\xi) \triangleq \int_{\mathbb{R}^n} f(x)\mathrm{e}^{-2\pi\mathrm{i}\xi\cdot x}\,\mathrm{d}x.$$
   - 由此联动产生**全部标准系数**（考试写错 $2\pi$ 会直接扣结果分！）：
     - 求导变乘法：$\widehat{\partial_j f}(\xi) = 2\pi\mathrm{i}\xi_j \hat{f}(\xi)$（教材 p.313 命题 8.15）；
     - 乘坐标变求导：$\partial_{\xi_j}\hat{f}(\xi) = -2\pi\mathrm{i}\widehat{x_j f}(\xi)$（教材 p.312 命题 8.14）；
     - **反演公式与 Plancherel 等式均不带 $(2\pi)^{-n}$**：
       $$f(x) = \int_{\mathbb{R}^n} \hat{f}(\xi)\mathrm{e}^{2\pi\mathrm{i}\xi\cdot x}\,\mathrm{d}\xi, \qquad \|\hat{f}\|_{L^2(\mathbb{R}^n)} = \|f\|_{L^2(\mathbb{R}^n)}.$$
   - **换算秘籍**：如果网课或别的参考书采用不带 $2\pi$ 的版本 $\tilde{f}(\omega) = \int f(x)\mathrm{e}^{-\mathrm{i}\omega\cdot x}\,\mathrm{d}x$，二者关系极简单——只需代入 $\omega = 2\pi\xi$，即 $\tilde{f}(\omega) = \hat{f}\!\left(\frac{\omega}{2\pi}\right)$。

---

## 第一部分：第 8 章核心知识零基础精讲

### 一、为什么要学广义函数与 Sobolev 空间？（一句话看懂全章主线）

传统微积分有两个致命痛点：
1. **很多物理与工程对象「不能求导」**：比如阶跃函数 $\operatorname{sgn}(x)$ 在 $0$ 点跳变、绝对值 $|x|$ 在 $0$ 点有尖角、点电荷密度 $\delta$ 甚至不是普通函数。
2. **偏微分方程（PDE）的经典解空间不完备**：用逼近法解方程时，一列光滑解的极限往往失去经典导数，掉出了 $C^k$ 空间。

全章就是用一套**分步升级的巧妙设计**解决这两个痛点：
- **第 1 步（8.2 节）**：先挑出一批「极其乖巧」的探针——紧支光滑函数 $\mathcal{D}(\Omega)$；
- **第 2 步（8.3 节）**：不再问函数“在某一点等于几”，而改为问它“和每个探针做积分得到几”，于是诞生了**广义函数（分布）$\mathcal{D}'(\Omega)$**，任何广义函数都能无穷次求导！
- **第 3 步（8.4 节）**：为了把求导变成代数乘法（Fourier 变换），把探针稍微放大成**速降函数 $\mathcal{S}(\mathbb{R}^n)$**，对偶得到**缓增分布 $\mathcal{S}'(\mathbb{R}^n)$**；
- **第 4 步（8.6–8.9 节）**：从广义函数里挑出“本身和各阶弱导数都在 $L^p$（尤其 $L^2$）里”的优等生，建成完备的 **Sobolev 空间 $W^{m,p}(\Omega)$ 与实指数空间 $H^s(\mathbb{R}^n)$**，成为现代偏微分方程的绝对舞台。

---

### 二、8.2 试验函数空间 $\mathcal{D}(\Omega) = C_c^\infty(\Omega)$（教材 p.303–304）

#### 1. 零基础直觉：什么是「超级探针」？
试验函数 $\varphi \in C_c^\infty(\Omega)$（也记作 $C_0^\infty(\Omega)$）同时满足两大神级性质：
- **无穷次光滑（$C^\infty$）**：想对它求多少阶偏导都行，永远不会出错；
- **紧支集（subscript $c$ / $0$）**：$\operatorname{supp}\varphi \triangleq \overline{\{x \in \Omega : \varphi(x) \neq 0\}}$ 是包含在开集 $\Omega$ 内部的紧集（在 $\mathbb{R}^n$ 中即有界闭集）。换句话说，**只要离开某个有限区域，$\varphi$ 及其一切导数统统严格等于 $0$！**

> **为什么这个设计如此天才？** 因为做分部积分时，$\int_\Omega u\,\varphi'\,\mathrm{d}x = \underbrace{[u\varphi]_{\partial\Omega}}_{=\,0} - \int_\Omega u'\,\varphi\,\mathrm{d}x$，**边界项永远自动归零**！

#### 2. 严格收敛定义与核心性质
- **定义 8.6（$\mathcal{D}(\Omega)$ 中的收敛，教材 p.303）**：设 $\{f_k\}_{k=1}^\infty \subset C_c^\infty(\Omega)$，$f \in C_c^\infty(\Omega)$。称 $f_k \to f$ 于 $C_c^\infty(\Omega)$（记为 $\mathcal{D}(\Omega)$），当且仅当同时满足两条：
  1. **共同紧支集**：存在**同一个**紧集 $K \subset \Omega$，使得所有 $f_k$ 与 $f$ 的支集都包含在 $K$ 中；
  2. **任意阶导数一致收敛**：对每个非负整数 $j = 0, 1, 2, \dots$，$\nabla^j f_k$ 在 $K$ 上一致收敛于 $\nabla^j f$（即对一切多重指标 $\alpha$，$\sup_{x \in K} |\partial^\alpha f_k(x) - \partial^\alpha f(x)| \to 0$）。
- **命题 8.2（p.303）**：对紧集 $K \subset \mathbb{R}^n$ 与赋范空间 $X$，线性映射 $T: C_c^\infty(K) \to X$ 连续 $\iff \exists k \ge 0, C > 0$ 使得 $\|Tf\|_X \le C\|f\|_{C^k(K)}$。
- **命题 8.4（p.303）**：对任意 $1 \le p \le \infty$，嵌入 $C_c^\infty(\Omega) \hookrightarrow L^p(\Omega)$ 连续单射。
- **命题 8.6（p.304）**：光滑系数的线性微分算子是 $\mathcal{D}(\Omega)$ 上的连续线性算子。

---

### 三、8.3 广义函数（分布）$\mathcal{D}'(\Omega)$（教材 p.304–312）

#### 1. 基本定义与连续性判定
- **定义 8.7（广义函数，p.304）**：称连续线性泛函 $\lambda: \mathcal{D}(\Omega) \to \mathbb{C}$ 为 $\Omega$ 上的**广义函数**（也称分布），全体记为 $\mathcal{D}'(\Omega)$。
- **连续性实用判定（式 (8.3)，p.304）**：线性泛函 $\lambda: \mathcal{D}(\mathbb{R}^n) \to \mathbb{C}$ 是广义函数 $\iff$ 对每个紧集 $K \subset \mathbb{R}^n$，存在非负整数 $k$（可依赖 $K$）与常数 $C > 0$，使得
  $$|\lambda(\varphi)| \le C\|\varphi\|_{C^k(K)}, \quad \forall \varphi \in C_c^\infty(K).$$
  - 若 $k$ 可不依赖紧集 $K$ 统一选取，则称 $\lambda$ 至多 $k$ 阶（定义 8.12，p.311）：如 $\delta$ 为 $0$ 阶，$\delta'$ 与 $\operatorname{p.v.}\frac{1}{x}$ 为 $1$ 阶。
- **定义 8.8（广义函数列的收敛，p.304）**：称 $\lambda_k \to \lambda$ 于 $\mathcal{D}'(\Omega)$，当且仅当**对每个固定的试验函数** $\varphi \in \mathcal{D}(\Omega)$，数列满足
  $$\lim_{k\to\infty}\lambda_k(\varphi) = \lambda(\varphi).$$
- **收敛强弱关系（p.305，必考结论！）**：
  $$L^p \text{ 收敛或一致收敛} \implies L^1_{\mathrm{loc}}(\Omega) \text{ 收敛} \implies \mathcal{D}'(\Omega) \text{ 收敛}.$$
  只要 $g_k \to g$ 于 $L^1_{\mathrm{loc}}(\Omega)$（即在每个紧集 $K \subset \Omega$ 上 $\int_K |g_k - g|\,\mathrm{d}x \to 0$），就自动有 $\lambda_{g_k} \to \lambda_g$ 于 $\mathcal{D}'(\Omega)$！

#### 2. 五大必背经典广义函数（教材例 8.2–8.6，p.304–305）

| 名称 | 记号 | 对 $\varphi \in \mathcal{D}(\Omega)$ 的作用公式 | 教材出处与备注 |
| :-------- | :------------ | :------------------------------- | :------------------ |
| 常规广义函数 | $\lambda_g$（$g \in L^1_{\mathrm{loc}}(\Omega)$） | $\lambda_g(\varphi) = \int_{\Omega} \varphi(x)\overline{g(x)}\,\mathrm{d}x$ | **例 8.2（p.304）**；实值 $g$ 时为 $\int g\varphi\,\mathrm{d}x$ |
| Dirac 脉冲函数 | $\delta$（Radon 测度特例） | $\delta(\varphi) = \varphi(0)$ | **例 8.3（p.304）**；$0$ 阶，非函数 |
| $\delta$ 的导数 | $\delta'$（$n=1$） | $\delta'(\varphi) = -\varphi'(0)$ | **例 8.4（p.305）**；$1$ 阶，非 Radon 测度 |
| Cauchy 主值 | $\operatorname{p.v.}\dfrac{1}{x}$（$n=1$） | $\left(\operatorname{p.v.}\dfrac{1}{x}\right)(\varphi) = \lim\limits_{\varepsilon\to 0^+}\displaystyle\int_{\lvert x\rvert>\varepsilon}\dfrac{\varphi(x)}{x}\,\mathrm{d}x$ | **例 8.5（p.305）**；因 $\frac{1}{x} \notin L^1_{\mathrm{loc}}$ 须挖洞对称取极限 |
| 例 8.6 广义函数 | $\lambda_r$（$n=1, r>0$） | $\lambda_r(\varphi) = \displaystyle\int_{\lvert x\rvert<r}\dfrac{\varphi(x)-\varphi(0)}{\lvert x\rvert}\,\mathrm{d}x + \displaystyle\int_{\lvert x\rvert\ge r}\dfrac{\varphi(x)}{\lvert x\rvert}\,\mathrm{d}x$ | **例 8.6（p.305）**；**习题 8 第 11(3) 题直接考 $r=1$ 即 $\lambda_1$！** |

> **零基础秒懂例 8.5（$\operatorname{p.v.}\frac{1}{x}$）与例 8.6（$\lambda_r$）**：
> - **例 8.5（Cauchy 主值 $\operatorname{p.v.}\frac{1}{x}$，p.305）**：$\frac{1}{x}$ 在原点附近不可积（$\int_0^1 \frac{\mathrm{d}x}{x} = +\infty$），不能直接按普通积分定义。但 $\frac{1}{x}$ 是**奇函数**，左右两侧分别冲向 $-\infty$ 和 $+\infty$；只要**对称地挖掉小区间 $(-\varepsilon, \varepsilon)$**（只积 $|x|>\varepsilon$），两边的无穷大刚好精准抵消！设 $\operatorname{supp}\varphi \subset [-R, R]$，由奇函数对称积分为 $0$ 得 $\int_{\varepsilon<|x|\le R}\frac{\varphi(0)}{x}\,\mathrm{d}x = 0$，于是
>   $$\int_{|x|>\varepsilon}\frac{\varphi(x)}{x}\,\mathrm{d}x = \int_{\varepsilon<|x|\le R}\frac{\varphi(x)-\varphi(0)}{x}\,\mathrm{d}x \xrightarrow{\varepsilon\to 0^+} \int_{-R}^{R}\frac{\varphi(x)-\varphi(0)}{x}\,\mathrm{d}x,$$
>   由中值定理 $\left|\frac{\varphi(x)-\varphi(0)}{x}\right| \le \|\varphi'\|_\infty$，积分收敛且 $\left|\left(\operatorname{p.v.}\frac{1}{x}\right)(\varphi)\right| \le 2R\|\varphi'\|_\infty$（$1$ 阶广义函数，既非局部可积函数也非 Radon 测度）。
> - **例 8.6（$\lambda_r$，p.305）**：函数 $\frac{1}{|x|}$ 同样在原点不可积（且是偶函数，挖洞也不能抵消！）。因此在 $|x|<r$ 内直接把 $\varphi(x)$ 减去常数 $\varphi(0)$，由于 $|\varphi(x)-\varphi(0)| \le |x|\|\varphi'\|_\infty$，分子的 $|x|$ 恰好消掉分母的 $|x|$，积分立刻收敛；外围 $|x|\ge r$ 本来无奇点，直接积分即可。

#### 3. 广义函数的核心运算与雷区
1. **广义导数（p.307–308）**：
   $$\frac{\partial}{\partial x_j}\lambda(\varphi) \triangleq -\lambda\!\left(\frac{\partial}{\partial x_j}\varphi\right), \qquad \partial^\alpha \lambda(\varphi) = (-1)^{|\alpha|}\lambda(\partial^\alpha \varphi).$$
   - **万能记忆口诀**：「**求导甩给试验函数，每甩一次加一个负号！**」
   - **例 8.8（p.308）**：若 $f \in C^m(\mathbb{R}^n)$，其广义导数与经典导数完全一致。
   - **命题 8.9（p.310，式 (8.6)）**：任意广义函数的混合偏导数**无条件可交换**次序：$\partial_{x_i}\partial_{x_j}\lambda = \partial_{x_j}\partial_{x_i}\lambda$。
   - **🔥 超级大杀器——求导与分布极限无条件交换（p.307）**：
     $$\lambda_k \to \lambda \text{ 于 } \mathcal{D}'(\Omega) \implies \partial^\alpha \lambda_k \to \partial^\alpha \lambda \text{ 于 } \mathcal{D}'(\Omega).$$
     *证明只有一行*：$(\partial^\alpha\lambda_k)(\varphi) = (-1)^{|\alpha|}\lambda_k(\partial^\alpha\varphi) \to (-1)^{|\alpha|}\lambda(\partial^\alpha\varphi) = (\partial^\alpha\lambda)(\varphi)$！
     结合「$L^1_{\mathrm{loc}}$ 收敛 $\implies \mathcal{D}'$ 收敛」，遇到难算极限的 $f_\varepsilon$ 时，**先找它的原函数 $F_\varepsilon$，证 $F_\varepsilon \to F$ 于 $L^1_{\mathrm{loc}}$，再两边求广义导数**，一招秒杀（见习题 8 第 14 题）！
   - **定理 8.7（p.310）**：一维区间上任意广义函数必有原函数；导数为 $0$ 的广义函数必为常数。
   - **定理 8.8（p.310）**：一维区间上若 $\lambda' \in L^1_{\mathrm{loc}}((a,b))$，则 $\lambda$ 必为绝对连续函数。

2. **光滑函数乘法（p.305，p.309 命题 8.8）**：
   - 对 $h \in C^\infty(\Omega), \lambda \in \mathcal{D}'(\Omega)$，$(h\lambda)(\varphi) = \lambda(\bar{h}\varphi)$；特别有 $x\delta = 0$，以及 Leibniz 乘法法则 $\partial_{x_j}(h\lambda) = h(\partial_{x_j}\lambda) + (\partial_{x_j}h)\lambda$。
   - **⚠️ 两大禁区（p.311–312，注记 8.7）**：
     1. 不能对一般广义函数施加非线性算子（如取绝对值 $|\lambda|$）；
     2. **不能把两个一般广义函数相乘**（如 $\delta^2$ 无意义），也**不能用不光滑的函数乘一般广义函数**（如 $f \in L^1(\mathbb{R})$ 时 $f\delta$ 无法良定义，见习题 8 第 7 题）。注记 8.7 指出：即使 $\varphi_{1/n} \to \delta$，平方后的 $\varphi_{1/n}^2$ 却发散到无穷！

3. **支集与稠密性（p.306–307）**：
   - 定义 8.11 与注记 8.5：$\operatorname{supp}\lambda \subset K$ 要求 $\varphi$ 在 $K$ 的某个**开邻域**内恒为 $0$ 时 $\lambda(\varphi)=0$（不能只在 $K$ 本身为 $0$，例如 $\varphi(0)=0$ 时 $\delta'(\varphi)=-\varphi'(0)$ 仍可非零）。
   - 引理 8.7：$\mathcal{D}(\mathbb{R}^n)$ 在 $\mathcal{D}'(\mathbb{R}^n)$ 中稠密。

---

### 四、8.4 Schwartz 空间 $\mathcal{S}(\mathbb{R}^n)$、缓增分布 $\mathcal{S}'(\mathbb{R}^n)$ 与 Fourier 变换（教材 p.312–319）

#### 1. 为什么不能直接在 $\mathcal{D}'(\mathbb{R}^n)$ 上做 Fourier 变换？
一个非零紧支光滑函数 $\varphi \in \mathcal{D}(\mathbb{R}^n)$，它的 Fourier 变换 $\hat{\varphi}$ 是解析函数，**不可能再有紧支集**（除非 $\varphi \equiv 0$）！也就是说，$\mathcal{F}$ 把 $\mathcal{D}(\mathbb{R}^n)$ 映出了自己，无法通过 $(\mathcal{F}\lambda)(\varphi) = \lambda(\mathcal{F}\varphi)$ 给一般分布定义 Fourier 变换。

**解决方案**：把探针空间从「紧支光滑 $\mathcal{D}(\mathbb{R}^n)$」适度放宽为「无穷远比任何多项式衰减得都快的光滑函数」，即 **Schwartz 空间 $\mathcal{S}(\mathbb{R}^n)$**！
- **定义 8.14（Schwartz 空间 $\mathcal{S}(\mathbb{R}^n)$，p.313）**：$f \in C^\infty(\mathbb{R}^n; \mathbb{C})$ 称为速降函数，如果对任意多重指标 $\alpha, \beta$ 都有
  $$\lim_{|x|\to\infty} x^\alpha \partial^\beta f(x) = 0.$$
  最典型的非紧支速降函数是 Gauss 函数 $\mathrm{e}^{-a|x|^2}$（$a > 0$）。
- **命题 8.16–8.17（p.314–316）**：$\mathcal{S}(\mathbb{R}^n)$ 对求导、乘多项式、乘法、卷积封闭，在 $L^p(\mathbb{R}^n)$（$1\le p<\infty$）中稠密；最重要的是，**$\mathcal{F}: \mathcal{S}(\mathbb{R}^n) \to \mathcal{S}(\mathbb{R}^n)$ 是双射同胚，且逆变换就是伴随变换 $\mathcal{F}^{-1} = \mathcal{F}^*$**！

#### 2. 缓增分布 $\mathcal{S}'(\mathbb{R}^n)$（教材 p.317–318）
- **定义 8.15（p.317）**：$\mathcal{S}(\mathbb{R}^n)$ 上的连续线性泛函全体称为**缓增分布空间** $\mathcal{S}'(\mathbb{R}^n)$。
- **包含链**：因为 $\mathcal{D}(\mathbb{R}^n) \hookrightarrow \mathcal{S}(\mathbb{R}^n)$ 连续且稠密，取对偶后方向反转：
  $$\mathcal{D}(\mathbb{R}^n) \subset \mathcal{S}(\mathbb{R}^n) \subset L^2(\mathbb{R}^n) \subset \mathcal{S}'(\mathbb{R}^n) \subset \mathcal{D}'(\mathbb{R}^n).$$
- **直观辨析（例 8.14 & 注记 8.9，p.317）**：“缓增”粗略指在无穷远处**至多呈多项式增长**（如多项式、有界函数、$L^p$ 函数、$\delta$ 都是缓增分布）；但像 $\mathrm{e}^x$ 那样指数爆炸的函数属于 $\mathcal{D}'(\mathbb{R})$ 却**不属于** $\mathcal{S}'(\mathbb{R})$！
- **定义 8.16（缓增分布的 Fourier 变换，p.317）**：对 $\lambda \in \mathcal{S}'(\mathbb{R}^n)$，定义 $(\mathcal{F}\lambda)(f) \triangleq \lambda(\mathcal{F}f)$（$\forall f \in \mathcal{S}(\mathbb{R}^n)$）。
  > *注记*：在实值或双线性配对 $\langle \lambda, f\rangle = \int gf$ 下，上式与普通函数的 Fourier 变换完全一致；在复共轭配对 $\lambda_g(f) = \int f\bar{g}$ 下，由 Parseval 等式 $\langle f, \mathcal{F}g\rangle_{L^2} = \langle \mathcal{F}^* f, g\rangle_{L^2}$ 可知应对偶使用 $\mathcal{F}^*$。本章所有实值/偶函数及命题 8.21 的运算律均直接按教材公式使用，不会产生任何符号冲突。

#### 3. ⚡ 本书 $2\pi$ 约定下 Fourier 变换高频公式速查表（p.312–318）

| 运算 / 性质 | 物理空间（$x$） | 频域空间（$\xi$） | 教材编号 |
| :--- | :--- | :--- | :--- |
| 定义（正变换） | $f(x)$ | $\hat{f}(\xi) = \displaystyle\int_{\mathbb{R}^n} f(x)\mathrm{e}^{-2\pi\mathrm{i}\xi\cdot x}\,\mathrm{d}x$ | 式 (8.9)（p.312） |
| 反演公式 | $f(x) = \displaystyle\int_{\mathbb{R}^n} \hat{f}(\xi)\mathrm{e}^{2\pi\mathrm{i}\xi\cdot x}\,\mathrm{d}\xi$ | $\hat{f}(\xi)$ | 命题 8.17、式 (8.18) |
| Riemann-Lebesgue | $f \in L^1(\mathbb{R}^n)$ | $\lim\limits_{\lvert\xi\rvert\to\infty}\hat{f}(\xi) = 0$（且连续） | 命题 8.13（p.312） |
| 导数 $\leftrightarrow$ 乘多项式 | $\partial_{x_j} f$（或 $\partial^\alpha f$） | $2\pi\mathrm{i}\xi_j \hat{f}(\xi)$（或 $(2\pi\mathrm{i}\xi)^\alpha \hat{f}(\xi)$） | 命题 8.15、8.21(6) |
| 乘坐标 $\leftrightarrow$ 频域求导 | $-2\pi\mathrm{i}x_j f(x)$ | $\partial_{\xi_j}\hat{f}(\xi)$ | 命题 8.14、8.21(6) |
| 卷积 $\leftrightarrow$ 逐点乘积 | $(f * g)(x)$ | $\hat{f}(\xi)\hat{g}(\xi)$ | 命题 8.21(2)（p.317） |
| 乘积 $\leftrightarrow$ 卷积 | $f(x)g(x)$ | $(\hat{f} * \hat{g})(\xi)$ | 命题 8.21(2)（p.317） |
| 平移 $\leftrightarrow$ 相位调制 | $(\tau_{x_0}f)(x) = f(x - x_0)$ | $\mathrm{e}^{-2\pi\mathrm{i}\xi\cdot x_0}\hat{f}(\xi)$ | 命题 8.21(3)（p.317） |
| 线性伸缩 | $f_L(x) = f(Lx)$（$L$ 可逆） | $\dfrac{1}{\lvert\det L\rvert}\hat{f}\big((L^*)^{-1}\xi\big)$ | 命题 8.18、8.21(5) |
| Parseval / Plancherel | $\langle f, g\rangle_{L^2},\ \lVert f\rVert_{L^2}$ | $\langle \hat{f}, \hat{g}\rangle_{L^2},\ \lVert\hat{f}\rVert_{L^2}$（酉算子！） | 定理 8.9（p.316） |
| $\delta$ 与常数 $1$ | $\delta$ ；常数 $1$ | $\hat{\delta} = 1$ ；$\hat{1} = \delta$ | 命题 8.21(7)（p.318） |
| Gauss 核（$a>0$） | $\mathrm{e}^{-ax^2}$（$1$ 维） | $\sqrt{\dfrac{\pi}{a}}\,\mathrm{e}^{-\frac{\pi^2\xi^2}{a}}$ | 习题 8 第 17(3) 题 |
| Poisson 核（$a>0$） | $\dfrac{1}{a^2+x^2}$（$1$ 维） | $\dfrac{\pi}{a}\,\mathrm{e}^{-2\pi a\lvert\xi\rvert}$ | 习题 8 第 17(2) 题 |

---

### 五、8.6–8.7 整数阶 Sobolev 空间 $W^{m,p}(\Omega)$ 与嵌入定理（教材 p.321–335）

#### 1. 弱导数与 $W^{m,p}(\Omega)$
- **定义 8.17（弱导数，p.321）**：设 $f \in L^p(\Omega)$，$\alpha \in \mathbb{N}^n$。若存在函数 $g \in L^1(\Omega)$（更一般地 $L^1_{\mathrm{loc}}(\Omega)$，见 p.309）满足对所有试验函数 $\varphi \in C_0^\infty(\Omega)$ 均有
  $$\int_\Omega f\,\partial^\alpha\varphi\,\mathrm{d}x = (-1)^{|\alpha|}\int_\Omega g\,\varphi\,\mathrm{d}x,$$
  则称 $g$ 是 $f$ 的 $\alpha$ 阶**弱导数**，记作 $\partial^\alpha f = g$。
- **⚠️ 弱导数 vs. 广义导数（注记 8.15，p.321，概念常考！）**：
  - 任何局部可积函数都有**广义导数**（落在 $\mathcal{D}'(\Omega)$ 里）；
  - 只有当这个广义导数**恰好还是一个局部可积函数**时，才叫**弱导数**！
  - *秒懂对比*：$|x|$ 的广义导数是 $\operatorname{sgn}(x) \in L^1_{\mathrm{loc}}(\mathbb{R})$，所以 $|x|$ **有弱导数**；而 $\operatorname{sgn}(x)$ 的广义导数是 $2\delta$，不是函数，所以 $\operatorname{sgn}(x)$ **没有弱导数**！
- **定义 8.18（Sobolev 空间，p.321–322）**：设 $1 \le p \le \infty, m \in \mathbb{N}$，定义
  $$W^{m,p}(\Omega) \triangleq \{f \in L^p(\Omega) : \partial^\alpha f \in L^p(\Omega),\ \forall |\alpha| \le m\}, \qquad \|f\|_{W^{m,p}(\Omega)} \triangleq \sum_{|\alpha|\le m}\|\partial^\alpha f\|_{L^p(\Omega)}.$$
  - 特别记 $H^m(\Omega) \triangleq W^{m,2}(\Omega)$（是 Hilbert 空间；注：教材 p.338 行首偶将 $W^{m,2}$ 误印为 $W^{2,m}$，以定义 8.18 及 p.336 命题 8.28 为准）。
  - **注记 8.16**：$W^{0,p}(\Omega) = L^p(\Omega)$，且 $W^{m+1,p}(\Omega) \subset W^{m,p}(\Omega)$。
  - **引理 8.8（p.322）**：当 $1 \le p < \infty$ 时，$C_c^\infty(\mathbb{R}^n)$ 在 $W^{k,p}(\mathbb{R}^n)$ 中稠密。

#### 2. Sobolev 嵌入定理（p.327–335，直观口诀记忆！）
> **核心直觉**：「**导数阶数可以兑换可积指数（甚至连续性），维数 $n$ 越高兑换成本越贵！**」记住核心标度指标 $\frac{1}{p} - \frac{k}{n}$：

- **情形 A：低光滑 $kp < n$（可积性提升，定理 8.11 & 8.13）**
  - 一阶（定理 8.11，p.327）：$1 \le p < n \implies W^{1,p}(\Omega) \hookrightarrow L^{p^*}(\Omega)$ 连续嵌入，其中 Sobolev 共轭指数满足
    $$\frac{1}{p^*} = \frac{1}{p} - \frac{1}{n} \iff p^* = \frac{np}{n - p}.$$
  - 高阶（定理 8.13，p.331）：$kp < n \implies W^{k,p}(\Omega) \hookrightarrow L^{\frac{np}{n-kp}}(\Omega)$。
- **情形 B：高光滑 $kp > n$（跃升为连续/Hölder 光滑函数，定理 8.12 & 8.13，p.331）**
  - $W^{k,p}(\Omega) \hookrightarrow C^{m,\alpha}(\Omega)$，其中 $m$ 取小于 $k - \frac{n}{p}$ 的非负整数（特别一阶 $p > n$ 时嵌入到 $C^{0, 1 - \frac{n}{p}}(\Omega)$）。
- **紧嵌入（Rellich-Kondrachov 定理 8.14 & 推论 8.3，p.333–335）**：
  - 当 $\Omega \subset \mathbb{R}^n$ 是**有界光滑区域**时，把连续嵌入的临界指数**严格往下降一点**（$q < \frac{np}{n-p}$；或 $n \le p$ 时任意 $q \in [1, \infty)$），嵌入 $W^{1,p}(\Omega) \hookrightarrow L^q(\Omega)$ 不仅是连续的，还是**紧算子**！推论 8.3：$m_1 > m_2 \implies W^{m_1,p}(\Omega) \hookrightarrow W^{m_2,p}(\Omega)$ 紧。
  - **注记 8.18（p.335）**：区域 $\Omega$ 的**有界性绝对不能丢**，在无界区域（如 $\mathbb{R}^n$）上平移一列紧支函数即可破坏紧性。

---

### 六、8.8–8.9 实指数 Sobolev 空间 $H^s(\mathbb{R}^n)$、迹定理、$H^{-m}$ 与椭圆算子（教材 p.335–340）

#### 1. 用 Plancherel 定理把「求导」翻译到频域，解锁任意实数阶 $s \in \mathbb{R}$！
对于整数 $k \ge 0$ 与 $f \in \mathcal{S}(\mathbb{R}^n)$，由 Plancherel 定理和 $\widehat{\partial_j f} = 2\pi\mathrm{i}\xi_j\hat{f}$，有（教材 p.335–336 式 (8.47)–(8.48)）：
$$\|f\|_{W^{k,2}(\mathbb{R}^n)}^2 = \int_{\mathbb{R}^n}\sum_{j=0}^k (2\pi|\xi|)^{2j}|\hat{f}(\xi)|^2\,\mathrm{d}\xi \asymp \int_{\mathbb{R}^n} (1 + |\xi|^2)^k |\hat{f}(\xi)|^2\,\mathrm{d}\xi.$$
引入极其方便的**日本括号（Japanese bracket）**记号：
$$\langle\xi\rangle \triangleq (1 + |\xi|^2)^{1/2}.$$
于是整数 $k$ 瞬间解放为**任意实数 $s \in \mathbb{R}$（正数、分数、负数统统适用！）**：

- **实指数 Sobolev 空间定义（教材 p.336）**：对任意 $s \in \mathbb{R}$，
  $$H^s(\mathbb{R}^n) \triangleq \big\{f \in \mathcal{S}'(\mathbb{R}^n) : \langle\xi\rangle^s \hat{f}(\xi) \in L^2(\mathbb{R}^n)\big\}, \qquad \|f\|_{H^s(\mathbb{R}^n)} \triangleq \|\langle\xi\rangle^s \hat{f}\|_{L^2(\mathbb{R}^n)}.$$
  - **怎么凭直觉理解正负 $s$？**
    - $s > 0$（如 $H^2$）：要求 $\hat{f}(\xi)$ 在高频 $|\xi|\to\infty$ 衰减得更快 $\implies$ 物理空间函数 $f$ 更光滑；
    - $s = 0$：$H^0(\mathbb{R}^n) = L^2(\mathbb{R}^n)$；
    - $s < 0$（如 $H^{-m}$）：允许 $\hat{f}(\xi)$ 在高频增长（只需乘上衰减权重 $\langle\xi\rangle^{-m}$ 后落回 $L^2$）$\implies$ 容纳了带奇点的分布（如 $\delta \in H^{-s}(\mathbb{R}^n)$ 当且仅当 $s > \frac{n}{2}$！）。
  - **基本性质（p.336 命题 8.28）**：
    - $s > s' \implies H^s(\mathbb{R}^n) \subset H^{s'}(\mathbb{R}^n)$；
    - $\mathcal{S}(\mathbb{R}^n)$ 在每个 $H^s(\mathbb{R}^n)$ 中稠密；非负整数 $k$ 时 $W^{k,2}(\mathbb{R}^n) = H^k(\mathbb{R}^n)$ 且范数等价；
    - $H^s(\mathbb{R}^n)$ 是复 Hilbert 空间，内积为 $\langle f, g\rangle_{H^s} = \int_{\mathbb{R}^n}\langle\xi\rangle^{2s}\hat{f}(\xi)\overline{\hat{g}(\xi)}\,\mathrm{d}\xi$（*注：教材 p.336 内积公式漏印了帽子，应为 $\hat{f},\hat{g}$*）。
- **对偶刻画（命题 8.29，p.336）**：$H^s(\mathbb{R}^n)$ 上的连续线性泛函由唯一的 $g \in H^{-s}(\mathbb{R}^n)$ 通过 $\lambda(f) = \int_{\mathbb{R}^n}\hat{f}(\xi)\overline{\hat{g}(\xi)}\,\mathrm{d}\xi$ 表示，即 $(H^s(\mathbb{R}^n))^* \cong H^{-s}(\mathbb{R}^n)$，且
  $$\|f\|_{H^s(\mathbb{R}^n)} = \sup\big\{|\langle f, g\rangle_{L^2(\mathbb{R}^n)}| : g \in \mathcal{S}(\mathbb{R}^n),\ \|g\|_{H^{-s}(\mathbb{R}^n)} \le 1\big\}.$$
- **$H^s(\mathbb{R}^n)$ 的嵌入定理（定理 8.15–8.16，p.336–337）**：
  - **定理 8.15**：若 $0 < s < \frac{n}{2}$ 且 $\frac{n}{2} - s \le \frac{n}{q} \le \frac{n}{2}$，则 $H^s(\mathbb{R}^n) \hookrightarrow L^q(\mathbb{R}^n)$ 连续嵌入；
  - **定理 8.16**：若 $s > \frac{n}{2}$，则对 $0 < \alpha \le \min\!\left(s - \frac{n}{2}, 1\right)$，$H^s(\mathbb{R}^n) \hookrightarrow C^{0,\alpha}(\mathbb{R}^n)$ 连续嵌入（口诀：**「光滑度过半维数 $n/2$，立刻连续！」**）。

#### 2. 常系数 $m$ 阶椭圆算子（教材 p.337）
- 设 $m \ge 0$，$L = \sum_{|\alpha|=m} c_\alpha \partial^\alpha$ 是 $m$ 阶常系数齐次微分算子，定义其**象征（symbol）**为齐次多项式 $l(\xi) \triangleq \sum_{|\alpha|=m} c_\alpha \xi^\alpha$。
- 若存在常数 $c > 0$，使得对所有 $\xi \in \mathbb{R}^n$ 均有
  $$|l(\xi)| \ge c|\xi|^m,$$
  则称 $L$ 为 **$m$ 阶椭圆算子**。
- **例 8.17（p.337）**：
  - **是椭圆算子**：Laplace 算子 $\Delta = \sum_{j=1}^n \partial_j^2$（$2$ 阶，象征 $l(\xi)=|\xi|^2$）；$\mathbb{R}^2$ 上的 Cauchy-Riemann 算子 $\partial_{x_1} - \mathrm{i}\partial_{x_2}$（$1$ 阶，象征 $l(\xi)=\xi_1-\mathrm{i}\xi_2$，模长 $|l(\xi)|=|\xi|$）。
  - **非椭圆算子**：热算子 $\partial_t - \Delta$、Schrödinger 算子 $\mathrm{i}\partial_t + \Delta$、波算子 $-\partial_t^2 + \Delta$（在 $\mathbb{R}^{1+n}$ 上象征会在非零频率点处消失或阶数不齐）。

#### 3. 8.9 节迹定理、$H^{-m}(\Omega)$ 与 Poincaré 不等式（教材 p.338–340）
- **迹定理（定理 8.17–8.19，p.338–339）**：$L^2$ 函数在边界 $\partial\Omega$（零测集）上本来不能谈取值，但 $H^1(\Omega)$ 的函数多了一阶弱导数后，限制映射 $\gamma_0(\varphi) = \varphi|_{\partial\Omega}$ 可唯一延拓为有界线性算子 $\gamma_0: H^1(\Omega) \to L^2(\partial\Omega)$，从而严格赋予了 PDE 的边界值意义！将直到 $m-1$ 阶法向迹全为 $0$ 的子空间记为 $H_0^m(\Omega)$（**命题 8.30**：$C_0^\infty(\Omega)$ 在 $H_0^m(\Omega)$ 中稠密）。
- **定理 8.20（$H^{-m}(\Omega)$ 的结构，p.339）**：$H^{-m}(\Omega) \triangleq (H_0^m(\Omega))'$ 恰好等于 $\{\partial^\alpha f : |\alpha| \le m,\ f \in L^2(\Omega)\}$ 张成的广义函数子空间。
- **定理 8.21（Poincaré 不等式，p.340）**：若 $\Omega \subset \mathbb{R}^n$ 沿 $x_1$ 方向有界（$\sup_{x\in\Omega}|x_1| = K < \infty$），则 $\|\varphi\|_{L^2(\Omega)} \le 2K\|\partial_1\varphi\|_{L^2(\Omega)}$（$\forall \varphi \in H_0^1(\Omega)$）。

---

## 第二部分：第 8 章全部 8 道作业逐题详解

> **题目索引**：以下按顺序完整解答 **习题 8 第 7、9、11、14、17(2)(3)、20、28 题**（教材 p.341–343）以及 **教材 p.337 练习 8.8.2**，共 8 题。

---

### 题目 1：习题 8 第 7 题（教材 p.341）

#### 【原题】
设 $f \in L^1(\mathbb{R})$，能否定义 $f$ 与 $\delta$ 的乘积？

#### 【零基础破题思路】
回想光滑函数 $h \in C^\infty(\mathbb{R})$ 乘以 $\delta$ 发生了什么：
$$(h\delta)(\varphi) = \delta(\bar{h}\varphi) = \overline{h(0)}\varphi(0) = \big(\overline{h(0)}\delta\big)(\varphi)$$
（实值函数就是 $h(0)\delta$）。本质上，**乘以 $\delta$ 会把 $h$ 在单点 $x = 0$ 处的精确取值“拎出来”**！

现在把光滑函数 $h$ 换成 $L^1(\mathbb{R})$ 里的函数 $f$，立刻撞上三堵死墙：
1. **$L^1$ 里没有“单点取值”**：$L^1(\mathbb{R})$ 的元素是几乎处处相等的函数等价类，单点 $\{0\}$ 测度为 $0$，随便改 $f(0)$ 都不改变它在 $L^1$ 中的身份；
2. **奇点可能发散**：$L^1(\mathbb{R})$ 里有像 $|x|^{-1/2}\chi_{(-1,1)}(x)$ 这样在原点趋于 $+\infty$ 的函数，根本没法给 $f(0)$ 赋有限值；
3. **哪怕想“先用光滑函数列在 $L^1$ 里逼近 $f$，再对乘积取 $\mathcal{D}'$ 极限”也不行**：我们可以构造出**都以 $0$ 为 $L^1$ 极限**的两列光滑函数（如 $0$ 和越来越窄但顶点永远高度为 $1$ 的高斯尖峰 $\mathrm{e}^{-k^2 x^2}$），它们乘上 $\delta$ 后一个永远是 $0$，另一个永远是 $\delta$！极限不唯一，甚至还能构造出直接发散的序列！

#### 【规范解答】
**结论：不能**（无法在 $L^1(\mathbb{R})$ 上对所有 $f$ 定义与光滑函数乘法相容且满足连续性的乘积 $f\delta \in \mathcal{D}'(\mathbb{R})$）。下面从三个层面严密说明：

1. **直接沿用点值公式不良定义（等价类与无界奇点障碍）**
   按教材 p.305 光滑函数与广义函数的乘积定义，对任意 $h \in C^\infty(\mathbb{R})$ 及试验函数 $\varphi \in \mathcal{D}(\mathbb{R})$，有
   $$(h\delta)(\varphi) = \delta(\bar{h}\varphi) = \overline{h(0)}\varphi(0) = \big(\overline{h(0)}\delta\big)(\varphi).$$
   若试图将该公式推广到 $f \in L^1(\mathbb{R})$，存在两个本质矛盾：
   - **等价类歧义**：$L^1(\mathbb{R})$ 中的元素是几乎处处相等的可测函数等价类，单点集 $\{0\}$ 是零测集。例如函数 $f_1(x) \equiv 0$ 与 $f_2(x) = \begin{cases} 1, & x = 0, \\ 0, & x \neq 0 \end{cases}$ 在 $L^1(\mathbb{R})$ 中是**同一个元素**，但它们在 $x = 0$ 处的值分别为 $0$ 与 $1$，无法唯一确定 $f(0)$。
   - **原点无有限值**：存在函数 $f(x) = |x|^{-1/2}\chi_{(-1,1)}(x)$，因 $\int_{-1}^{1} |x|^{-1/2}\,\mathrm{d}x = 4 < \infty$ 故 $f \in L^1(\mathbb{R})$，但 $\lim_{x\to 0} f(x) = +\infty$，在原点处无法赋予任何有限值。

2. **通过 $L^1(\mathbb{R})$ 稠密逼近作连续延拓也不可能**
   由于 $C^\infty(\mathbb{R}) \cap L^1(\mathbb{R})$（乃至 $C_c^\infty(\mathbb{R})$）在 $L^1(\mathbb{R})$ 中稠密，考虑能否对 $f \in L^1(\mathbb{R})$ 取光滑函数列 $f_k \to f$（于 $L^1(\mathbb{R})$），定义 $f\delta \triangleq \lim_{k\to\infty} f_k\delta$（于 $\mathcal{D}'(\mathbb{R})$）。以下反例证明该极限**既不唯一，也不一定存在**：
   - 取 $g_k(x) \equiv 0 \in C^\infty(\mathbb{R}) \cap L^1(\mathbb{R})$，显然 $g_k \to 0$ 于 $L^1(\mathbb{R})$，且对所有 $k$ 有 $g_k\delta = 0 \to 0$（于 $\mathcal{D}'(\mathbb{R})$）。
   - 取 $h_k(x) = \mathrm{e}^{-k^2 x^2} \in C^\infty(\mathbb{R}) \cap L^1(\mathbb{R})$，计算其 $L^1$ 范数：
     $$\|h_k - 0\|_{L^1(\mathbb{R})} = \int_{-\infty}^{\infty}\mathrm{e}^{-k^2 x^2}\,\mathrm{d}x = \frac{\sqrt{\pi}}{k} \xrightarrow{k\to\infty} 0,$$
     即 $h_k \to 0$ 于 $L^1(\mathbb{R})$。然而对任意 $\varphi \in \mathcal{D}(\mathbb{R})$，
     $$(h_k\delta)(\varphi) = \overline{h_k(0)}\varphi(0) = 1 \cdot \varphi(0) = \delta(\varphi) \implies \lim_{k\to\infty} h_k\delta = \delta \neq 0.$$
   - 进一步取 $u_k(x) = \sqrt{k}\,\mathrm{e}^{-k^2 x^2} \in C^\infty(\mathbb{R}) \cap L^1(\mathbb{R})$，则 $\|u_k\|_{L^1(\mathbb{R})} = \sqrt{\frac{\pi}{k}} \to 0$，但对满足 $\varphi(0) \neq 0$ 的试验函数，$(u_k\delta)(\varphi) = \sqrt{k}\varphi(0)$ 当 $k \to \infty$ 时发散！

综上所述，映射 $h \mapsto h\delta$ 不能从光滑函数空间合理延拓到 $L^1(\mathbb{R})$，因此**不能定义**一般 $f \in L^1(\mathbb{R})$ 与 $\delta$ 的乘积。

#### 【考点复盘与避坑提示】
- **核心考点**：理解广义函数乘法的正则性门槛（教材 p.305、p.311–312 注记 8.7）以及 $L^1$ 拓扑不控制单点值的本质。
- **延伸思考**：如果把条件加强为 $f \in C(\mathbb{R})$（连续函数），那么**可以**定义 $f\delta \triangleq \overline{f(0)}\delta$！因为 $\delta$ 是 $0$ 阶分布（Radon 测度，命题 8.10），只需要乘子在原点连续（甚至只要求在原点逐点连续）就够了，不需要 $C^\infty$。但若是 $f\delta'$，则至少需要 $f \in C^1(\mathbb{R})$。

---

### 题目 2：习题 8 第 9 题（教材 p.341）

#### 【原题】
在广义函数意义下求下列极限：
(1) $\displaystyle\lim_{\varepsilon\to 0}\frac{\varepsilon}{(\varepsilon^2 + x^2)}$；
(2) $\displaystyle\lim_{\varepsilon\to 0}\frac{1}{\sqrt{\pi\varepsilon}}\mathrm{e}^{-\frac{x^2}{\varepsilon}}$。

#### 【零基础破题思路】
按定义 8.8，求一列常规函数 $f_\varepsilon$ 在 $\mathcal{D}'(\mathbb{R})$ 中的极限，**固定模板只有三步**：
1. 任取一个试验函数 $\varphi \in \mathcal{D}(\mathbb{R})$，写出积分 $\int_{-\infty}^{\infty} f_\varepsilon(x)\varphi(x)\,\mathrm{d}x$；
2. **变量替换“把 $\varepsilon$ 从分母/指数里赶到 $\varphi$ 的括号里”**：第 (1) 题令 $x = \varepsilon y$，第 (2) 题令 $x = \sqrt{\varepsilon}\,y$；
3. **用 Lebesgue 控制收敛定理把 $\lim_{\varepsilon\to 0^+}$ 放进积分号里**（因为 $\varphi$ 紧支连续有界，$|\varphi(\cdot)| \le M$，外层剩下的 $\frac{1}{1+y^2}$ 和 $\mathrm{e}^{-y^2}$ 都是经典可积函数），直接算出常数倍的 $\varphi(0)$！

> **特别提醒**：第 (1) 题的表达式关于 $\varepsilon$ 是**奇函数**，当 $\varepsilon \to 0^+$ 时极限是 $\pi\delta$，而当 $\varepsilon \to 0^-$ 时极限是 $-\pi\delta$！作业与教材通常默认参数 $\varepsilon > 0$（即 $\varepsilon \to 0^+$），我们在解答中先给出 $\varepsilon \to 0^+$ 的完整证明，再顺手补一句 $\varepsilon \to 0^-$ 的结论，滴水不漏！

#### 【规范解答】
任取试验函数 $\varphi \in \mathcal{D}(\mathbb{R}) = C_c^\infty(\mathbb{R})$。因为 $\varphi$ 具有紧支集且连续，故有界，记 $M \triangleq \sup_{x\in\mathbb{R}} |\varphi(x)| < \infty$。两题中的函数均为实值光滑可积函数，按例 8.2 均定义 $\mathcal{D}'(\mathbb{R})$ 中的常规广义函数。

**(1) 求解 $\displaystyle\lim_{\varepsilon\to 0}\frac{\varepsilon}{\varepsilon^2 + x^2}$**

- **情形一：$\varepsilon \to 0^+$（默认情形）**
  对任意 $\varepsilon > 0$，作变量替换 $x = \varepsilon y$（则 $\mathrm{d}x = \varepsilon\,\mathrm{d}y$，积分限仍为 $(-\infty, +\infty)$），得
  $$\int_{-\infty}^{\infty}\frac{\varepsilon}{\varepsilon^2 + x^2}\varphi(x)\,\mathrm{d}x = \int_{-\infty}^{\infty}\frac{\varepsilon}{\varepsilon^2(1 + y^2)}\varphi(\varepsilon y)\cdot\varepsilon\,\mathrm{d}y = \int_{-\infty}^{\infty}\frac{\varphi(\varepsilon y)}{1 + y^2}\,\mathrm{d}y.$$
  验证 Lebesgue 控制收敛定理的条件：
  1. 对每个固定的 $y \in \mathbb{R}$，由 $\varphi$ 的连续性，有 $\lim_{\varepsilon\to 0^+} \frac{\varphi(\varepsilon y)}{1 + y^2} = \frac{\varphi(0)}{1 + y^2}$；
  2. 对所有 $\varepsilon > 0$ 及 $y \in \mathbb{R}$，有统一控制估计
     $$\left|\frac{\varphi(\varepsilon y)}{1 + y^2}\right| \le \frac{M}{1 + y^2},$$
     且控制函数可积：$\int_{-\infty}^{\infty}\frac{M}{1 + y^2}\,\mathrm{d}y = M\big[\arctan y\big]_{-\infty}^{+\infty} = \pi M < \infty$。

  由 Lebesgue 控制收敛定理，得
  $$\lim_{\varepsilon\to 0^+}\int_{-\infty}^{\infty}\frac{\varepsilon}{\varepsilon^2 + x^2}\varphi(x)\,\mathrm{d}x = \int_{-\infty}^{\infty}\frac{\varphi(0)}{1 + y^2}\,\mathrm{d}y = \pi\varphi(0) = (\pi\delta)(\varphi).$$
  故当 $\varepsilon \to 0^+$ 时，在 $\mathcal{D}'(\mathbb{R})$ 意义下，
  $$\lim_{\varepsilon\to 0^+}\frac{\varepsilon}{\varepsilon^2 + x^2} = \pi\delta.$$

- **情形二：$\varepsilon \to 0^-$（单侧极限补充）**
  令 $\eta = -\varepsilon > 0$，则 $\frac{\varepsilon}{\varepsilon^2 + x^2} = -\frac{\eta}{\eta^2 + x^2}$。由情形一立刻得到，当 $\varepsilon \to 0^-$ 时在 $\mathcal{D}'(\mathbb{R})$ 中极限为 $-\pi\delta$。
  （因此严格的双侧极限 $\varepsilon \to 0$ 不存在；按教材默认 $\varepsilon > 0$ 理解，答案为 $\pi\delta$。）

---

**(2) 求解 $\displaystyle\lim_{\varepsilon\to 0}\frac{1}{\sqrt{\pi\varepsilon}}\mathrm{e}^{-\frac{x^2}{\varepsilon}}$**

由表达式中的 $\sqrt{\pi\varepsilon}$ 可知必有 $\varepsilon > 0$，即考虑 $\varepsilon \to 0^+$。作变量替换 $x = \sqrt{\varepsilon}\,y$（则 $\mathrm{d}x = \sqrt{\varepsilon}\,\mathrm{d}y$，积分限仍为 $(-\infty, +\infty)$），得
$$\int_{-\infty}^{\infty}\frac{1}{\sqrt{\pi\varepsilon}}\mathrm{e}^{-\frac{x^2}{\varepsilon}}\varphi(x)\,\mathrm{d}x = \frac{1}{\sqrt{\pi}}\int_{-\infty}^{\infty}\mathrm{e}^{-y^2}\varphi(\sqrt{\varepsilon}\,y)\,\mathrm{d}y.$$
验证 Lebesgue 控制收敛定理的条件：
1. 对每个固定的 $y \in \mathbb{R}$，由 $\varphi$ 的连续性，有 $\lim_{\varepsilon\to 0^+}\mathrm{e}^{-y^2}\varphi(\sqrt{\varepsilon}\,y) = \mathrm{e}^{-y^2}\varphi(0)$；
2. 对所有 $\varepsilon > 0$ 及 $y \in \mathbb{R}$，有统一控制估计
   $$\left|\mathrm{e}^{-y^2}\varphi(\sqrt{\varepsilon}\,y)\right| \le M\mathrm{e}^{-y^2},$$
   由熟知的 Gauss 积分 $\int_{-\infty}^{\infty}\mathrm{e}^{-y^2}\,\mathrm{d}y = \sqrt{\pi} < \infty$，控制函数属于 $L^1(\mathbb{R})$。

由 Lebesgue 控制收敛定理，得
$$\lim_{\varepsilon\to 0^+}\int_{-\infty}^{\infty}\frac{1}{\sqrt{\pi\varepsilon}}\mathrm{e}^{-\frac{x^2}{\varepsilon}}\varphi(x)\,\mathrm{d}x = \frac{\varphi(0)}{\sqrt{\pi}}\int_{-\infty}^{\infty}\mathrm{e}^{-y^2}\,\mathrm{d}y = \varphi(0) = \delta(\varphi).$$
故在 $\mathcal{D}'(\mathbb{R})$ 意义下，
$$\lim_{\varepsilon\to 0^+}\frac{1}{\sqrt{\pi\varepsilon}}\mathrm{e}^{-\frac{x^2}{\varepsilon}} = \delta.$$

#### 【考点复盘与避坑提示】
- **核心考点**：近似恒等元（Approximate Identity / $\delta$ 序列）与 Lebesgue 控制收敛定理的规范用法。其实这两题藏着 PDE 最著名的两大基本核：
  - 第 (1) 题除以 $\pi$ 就是**上半平面 Poisson 核** $P_\varepsilon(x) = \frac{1}{\pi}\frac{\varepsilon}{\varepsilon^2+x^2}$（直接预告**第 20 题**！）；
  - 第 (2) 题令 $\varepsilon = 4t$ 就是**一维热核** $K_t(x) = \frac{1}{\sqrt{4\pi t}}\mathrm{e}^{-x^2/(4t)}$（教材 p.319 式 (8.23)）。
- **避坑提示**：千万不要只写“把 $\lim$ 放进积分号”而不写**可积控制函数**（$M/(1+y^2)$ 和 $M\mathrm{e}^{-y^2}$）；也别忘了检查第 (1) 题 $\varepsilon$ 正负号导致积分换元是否反转上下限！

---

### 题目 3：习题 8 第 11 题（教材 p.341）

#### 【原题】
设 $n = 1$。证明下列结论：
(1) 符号函数 $\operatorname{sgn}(x)$ 的导数等于 $2\delta$。
(2) 局部可积函数 $\ln|x|$ 的导数等于 $\operatorname{p.v.}\dfrac{1}{x}$。
(3) 局部可积函数 $\ln|x|\operatorname{sgn}(x)$ 的导数等于例 8.6 中的广义函数 $\lambda_1$。
(4) $|x|$ 的导数是 $\operatorname{sgn}(x)$。

#### 【零基础破题思路】
四道小题全部围绕同一个口诀展开——**「求导甩给 $\varphi$ 加负号，挖掉奇点分部积，边界消零看中值！」**
- 第 (1)、(4) 题没有发散奇点，只在 $x=0$ 处有跳变或折点：直接把 $(-\infty, \infty)$ 拆成 $(-\infty, 0)$ 和 $(0, \infty)$ 两段，逐段分部积分即可。
- 第 (2) 题 $\ln|x|$ 在原点发散到 $-\infty$（虽然它是局部可积的！）：先用控制收敛定理挖掉 $(-\varepsilon, \varepsilon)$，在 $(-\infty, -\varepsilon)$ 和 $(\varepsilon, \infty)$ 上分别分部积分。神奇的事发生了：边界项凑成了 $\ln\varepsilon(\varphi(\varepsilon)-\varphi(-\varepsilon))$，靠 $\varphi$ 的光滑性（$|\varphi(\varepsilon)-\varphi(-\varepsilon)| = O(\varepsilon)$）压死对数奇性（$\varepsilon\ln\varepsilon \to 0$）！
- 第 (3) 题是**全题最精彩的一小问**！先抄出教材 p.305 例 8.6（取 $r=1$）：
  $$\lambda_1(\varphi) = \int_{|x|<1}\frac{\varphi(x)-\varphi(0)}{|x|}\,\mathrm{d}x + \int_{|x|\ge 1}\frac{\varphi(x)}{|x|}\,\mathrm{d}x.$$
  观察到 $\lambda_1$ 是在 $|x|=1$ 处分段的，而且在内段 $|x|<1$ 里藏着 $\varphi(x)-\varphi(0)$！这就直接告诉了我们**必胜绝招**：
  - 把积分拆成四段 $(-R, -1), (-1, -\varepsilon), (\varepsilon, 1), (1, R)$；
  - 在内侧两段 $(-1, -\varepsilon)$ 和 $(\varepsilon, 1)$ 上，把 $\varphi'(x)$ 改写成 $[\varphi(x)-\varphi(0)]'$ 再分部积分——这样原点端 $\pm\varepsilon$ 处出现 $\ln\varepsilon(\varphi(\pm\varepsilon)-\varphi(0)) = O(\varepsilon\ln\varepsilon) \to 0$，而接缝点 $\pm 1$ 处恰好因为 $\ln 1 = 0$，**所有边界项全部清零**！

#### 【规范解答】
本题中四个函数均为实值局部可积函数（对 (2)(3)，因 $\int_0^R |\ln x|\,\mathrm{d}x = [x\ln x - x]_0^R < \infty$，故 $\ln|x|$ 与 $\ln|x|\operatorname{sgn}(x)$ 均属于 $L^1_{\mathrm{loc}}(\mathbb{R})$），由例 8.2 知其诱导广义函数满足 $\lambda_g(\varphi) = \int_{-\infty}^{\infty} g(x)\varphi(x)\,\mathrm{d}x$。

以下统一任取 $\varphi \in \mathcal{D}(\mathbb{R})$，设 $\operatorname{supp}\varphi \subset [-R, R]$（取足够大的 $R > 1$），从而当 $|x| \ge R$ 时 $\varphi(x) = 0$。

---

**(1) 证明 $(\operatorname{sgn} x)' = 2\delta$**

符号函数定义为 $\operatorname{sgn}(x) = \begin{cases} 1, & x > 0, \\ 0, & x = 0, \\ -1, & x < 0. \end{cases}$ 按广义导数定义与 Newton-Leibniz 公式，有
$$\begin{aligned}
(\operatorname{sgn})'(\varphi) &= -\int_{-\infty}^{\infty}\operatorname{sgn}(x)\varphi'(x)\,\mathrm{d}x \\
&= -\left(-\int_{-R}^{0}\varphi'(x)\,\mathrm{d}x + \int_{0}^{R}\varphi'(x)\,\mathrm{d}x\right) \\
&= \big[\varphi(0) - \varphi(-R)\big] - \big[\varphi(R) - \varphi(0)\big] \\
&= 2\varphi(0) = (2\delta)(\varphi).
\end{aligned}$$
由 $\varphi \in \mathcal{D}(\mathbb{R})$ 的任意性，得 $(\operatorname{sgn} x)' = 2\delta$。

---

**(2) 证明 $(\ln|x|)' = \operatorname{p.v.}\dfrac{1}{x}$**

因为 $\ln|x| \in L^1_{\mathrm{loc}}(\mathbb{R})$ 且 $\varphi' \in C_c^\infty(\mathbb{R})$ 有界紧支，故 $\ln|x|\varphi'(x) \in L^1(\mathbb{R})$。利用指示函数 $\chi_{|x|>\varepsilon}(x) \to 1$（a.e.）及控制函数 $|\ln|x|\varphi'(x)| \in L^1(\mathbb{R})$，由 Lebesgue 控制收敛定理得
$$(\ln|x|)'(\varphi) = -\int_{-R}^{R}\ln|x|\,\varphi'(x)\,\mathrm{d}x = -\lim_{\varepsilon\to 0^+}\left(\int_{-R}^{-\varepsilon}\ln(-x)\,\varphi'(x)\,\mathrm{d}x + \int_{\varepsilon}^{R}\ln x\,\varphi'(x)\,\mathrm{d}x\right).$$
注意到当 $x < 0$ 时 $\frac{\mathrm{d}}{\mathrm{d}x}\ln(-x) = \frac{-1}{-x} = \frac{1}{x}$。对两段分别分部积分（利用 $\varphi(\pm R) = 0$）：
$$\begin{aligned}
\int_{-R}^{-\varepsilon}\ln(-x)\,\varphi'(x)\,\mathrm{d}x &= \Big[\ln(-x)\varphi(x)\Big]_{-R}^{-\varepsilon} - \int_{-R}^{-\varepsilon}\frac{\varphi(x)}{x}\,\mathrm{d}x = \ln\varepsilon\,\varphi(-\varepsilon) - \int_{-R}^{-\varepsilon}\frac{\varphi(x)}{x}\,\mathrm{d}x, \\
\int_{\varepsilon}^{R}\ln x\,\varphi'(x)\,\mathrm{d}x &= \Big[\ln x\,\varphi(x)\Big]_{\varepsilon}^{R} - \int_{\varepsilon}^{R}\frac{\varphi(x)}{x}\,\mathrm{d}x = -\ln\varepsilon\,\varphi(\varepsilon) - \int_{\varepsilon}^{R}\frac{\varphi(x)}{x}\,\mathrm{d}x.
\end{aligned}$$
两式相加并取负号，得
$$-\left(\int_{-R}^{-\varepsilon} + \int_{\varepsilon}^{R}\right)\ln|x|\,\varphi'(x)\,\mathrm{d}x = \ln\varepsilon\big(\varphi(\varepsilon) - \varphi(-\varepsilon)\big) + \int_{|x|>\varepsilon}\frac{\varphi(x)}{x}\,\mathrm{d}x.$$
由 Lagrange 中值定理，$|\varphi(\varepsilon) - \varphi(-\varepsilon)| \le 2\varepsilon\|\varphi'\|_\infty$，从而
$$\big|\ln\varepsilon\,\big(\varphi(\varepsilon) - \varphi(-\varepsilon)\big)\big| \le 2\|\varphi'\|_\infty \cdot \varepsilon|\ln\varepsilon| \xrightarrow{\varepsilon\to 0^+} 0.$$
令 $\varepsilon \to 0^+$，由例 8.5 得
$$(\ln|x|)'(\varphi) = \lim_{\varepsilon\to 0^+}\int_{|x|>\varepsilon}\frac{\varphi(x)}{x}\,\mathrm{d}x = \left(\operatorname{p.v.}\frac{1}{x}\right)(\varphi).$$
故 $(\ln|x|)' = \operatorname{p.v.}\dfrac{1}{x}$。

---

**(3) 证明 $(\ln|x|\operatorname{sgn}(x))' = \lambda_1$**

回忆教材 p.305 例 8.6（取 $r = 1$）：
$$\lambda_1(\varphi) = \int_{|x|<1}\frac{\varphi(x) - \varphi(0)}{|x|}\,\mathrm{d}x + \int_{|x|\ge 1}\frac{\varphi(x)}{|x|}\,\mathrm{d}x.$$
记 $g(x) = \ln|x|\operatorname{sgn}(x)$，则当 $x > 0$ 时 $g(x) = \ln x$，当 $x < 0$ 时 $g(x) = -\ln(-x)$；特别在 $x \neq 0$ 处恒有
$$g'(x) = \frac{1}{|x|}, \qquad \text{且 } g(\pm 1) = \pm\ln 1 = 0.$$
同 (2)，因 $g(x)\varphi'(x) \in L^1(\mathbb{R})$，由控制收敛定理挖去 $(-\varepsilon, \varepsilon)$（$0 < \varepsilon < 1 < R$），并拆为四段积分：
$$g'(\varphi) = -\lim_{\varepsilon\to 0^+}\big(I_1 + I_2 + I_3 + I_4\big),$$
其中四段分别作如下分部积分（内侧两段 $I_2, I_3$ 改写 $\varphi'(x) = \big(\varphi(x)-\varphi(0)\big)'$，且所有端点利用 $g(\pm 1)=0$ 及 $\varphi(\pm R)=0$）：
1. $I_1 = \displaystyle\int_{-R}^{-1} \big(-\ln(-x)\big)\varphi'(x)\,\mathrm{d}x = \Big[-\ln(-x)\varphi(x)\Big]_{-R}^{-1} - \int_{-R}^{-1}\frac{\varphi(x)}{|x|}\,\mathrm{d}x = -\int_{-R}^{-1}\frac{\varphi(x)}{|x|}\,\mathrm{d}x$；
2. $I_2 = \displaystyle\int_{-1}^{-\varepsilon} \big(-\ln(-x)\big)\big(\varphi(x)-\varphi(0)\big)'\,\mathrm{d}x = \Big[-\ln(-x)\big(\varphi(x)-\varphi(0)\big)\Big]_{-1}^{-\varepsilon} - \int_{-1}^{-\varepsilon}\frac{\varphi(x)-\varphi(0)}{|x|}\,\mathrm{d}x$
   $= -\ln\varepsilon\big(\varphi(-\varepsilon)-\varphi(0)\big) - \displaystyle\int_{-1}^{-\varepsilon}\frac{\varphi(x)-\varphi(0)}{|x|}\,\mathrm{d}x$；
3. $I_3 = \displaystyle\int_{\varepsilon}^{1} \ln x\,\big(\varphi(x)-\varphi(0)\big)'\,\mathrm{d}x = \Big[\ln x\,\big(\varphi(x)-\varphi(0)\big)\Big]_{\varepsilon}^{1} - \int_{\varepsilon}^{1}\frac{\varphi(x)-\varphi(0)}{|x|}\,\mathrm{d}x$
   $= -\ln\varepsilon\big(\varphi(\varepsilon)-\varphi(0)\big) - \displaystyle\int_{\varepsilon}^{1}\frac{\varphi(x)-\varphi(0)}{|x|}\,\mathrm{d}x$；
4. $I_4 = \displaystyle\int_{1}^{R} \ln x\,\varphi'(x)\,\mathrm{d}x = \Big[\ln x\,\varphi(x)\Big]_{1}^{R} - \int_{1}^{R}\frac{\varphi(x)}{|x|}\,\mathrm{d}x = -\int_{1}^{R}\frac{\varphi(x)}{|x|}\,\mathrm{d}x$。

将四式相加并取负号，得
$$-(I_1 + I_2 + I_3 + I_4) = \ln\varepsilon\big(\varphi(-\varepsilon)-\varphi(0)\big) + \ln\varepsilon\big(\varphi(\varepsilon)-\varphi(0)\big) + \int_{\varepsilon < |x| < 1}\frac{\varphi(x)-\varphi(0)}{|x|}\,\mathrm{d}x + \int_{|x|\ge 1}\frac{\varphi(x)}{|x|}\,\mathrm{d}x.$$
由中值定理，$|\varphi(\pm\varepsilon) - \varphi(0)| \le \varepsilon\|\varphi'\|_\infty$，故边界项之和满足
$$\Big|\ln\varepsilon\big(\varphi(-\varepsilon)-\varphi(0)\big) + \ln\varepsilon\big(\varphi(\varepsilon)-\varphi(0)\big)\Big| \le 2\|\varphi'\|_\infty \cdot \varepsilon|\ln\varepsilon| \xrightarrow{\varepsilon\to 0^+} 0.$$
又因在 $(-1, 1)\setminus\{0\}$ 上 $\left|\frac{\varphi(x)-\varphi(0)}{|x|}\right| \le \|\varphi'\|_\infty \in L^1((-1,1))$，令 $\varepsilon \to 0^+$ 并由控制收敛定理得
$$g'(\varphi) = \int_{|x|<1}\frac{\varphi(x)-\varphi(0)}{|x|}\,\mathrm{d}x + \int_{|x|\ge 1}\frac{\varphi(x)}{|x|}\,\mathrm{d}x = \lambda_1(\varphi).$$
故 $(\ln|x|\operatorname{sgn}(x))' = \lambda_1$。

---

**(4) 证明 $(|x|)' = \operatorname{sgn}(x)$**

因 $|x|$ 连续（从而局部可积），拆分为 $[-R, 0]$ 与 $[0, R]$ 两段分部积分（注意 $\varphi(\pm R)=0$ 且原点处 $0\cdot\varphi(0)=0$）：
$$\begin{aligned}
(|x|)'(\varphi) &= -\int_{-R}^{R}|x|\varphi'(x)\,\mathrm{d}x = -\left(\int_{-R}^{0}(-x)\varphi'(x)\,\mathrm{d}x + \int_{0}^{R}x\varphi'(x)\,\mathrm{d}x\right) \\
&= -\left(\Big[-x\varphi(x)\Big]_{-R}^{0} + \int_{-R}^{0}\varphi(x)\,\mathrm{d}x + \Big[x\varphi(x)\Big]_{0}^{R} - \int_{0}^{R}\varphi(x)\,\mathrm{d}x\right) \\
&= -\int_{-R}^{0}\varphi(x)\,\mathrm{d}x + \int_{0}^{R}\varphi(x)\,\mathrm{d}x = \int_{-\infty}^{\infty}\operatorname{sgn}(x)\varphi(x)\,\mathrm{d}x = \operatorname{sgn}(\varphi).
\end{aligned}$$
故 $(|x|)' = \operatorname{sgn}(x)$。

#### 【考点复盘与避坑提示】
- **核心考点**：广义导数的定义、挖洞分部积分法、利用 $\varphi(x)-\varphi(0)$ 消去边界发散项的经典技巧。
- **第 (3) 题为什么必须是 $\lambda_1$ 而不是 $\lambda_r$（$r \neq 1$）？** 因为只有切在 $|x|=1$ 处，$\ln 1 = 0$ 才刚好把内外交界处的边界项消成 $0$；如果硬切在 $\pm r$ 处，会多出一个 $2(\ln r)\varphi(0)$ 的边界项，实际上 $\lambda_1 = \lambda_r + 2(\ln r)\delta$！
- **弱导数对照**：顺手记牢，由 (4) 知 $|x|$ 有弱导数 $\operatorname{sgn}(x)$，由 (1) 知 $\operatorname{sgn}(x)$ 无弱导数（因 $2\delta$ 不是局部可积函数）。

---

### 题目 4：习题 8 第 14 题（教材 p.341）

#### 【原题】
设 $f_\varepsilon(x) = (x+\mathrm{i}\varepsilon)^{-1} + (x-\mathrm{i}\varepsilon)^{-1}$。证明：极限 $\lim_{\varepsilon\to 0} f_\varepsilon$ 属于 $\mathcal{D}'(\mathbb{R}^n)$。

> *教材印刷勘误提示*：原题公式中的变量 $x$ 是一维实变量，题末的 $\mathcal{D}'(\mathbb{R}^n)$（原书印作 $D'(\mathbb{R}^n)$）系 $\mathcal{D}'(\mathbb{R})$ 的笔误。

#### 【零基础破题思路】
第一步先做初中通分，立刻看穿纸老虎：
$$f_\varepsilon(x) = \frac{1}{x+\mathrm{i}\varepsilon} + \frac{1}{x-\mathrm{i}\varepsilon} = \frac{(x-\mathrm{i}\varepsilon) + (x+\mathrm{i}\varepsilon)}{x^2 + \varepsilon^2} = \frac{2x}{x^2 + \varepsilon^2}.$$
注意两个惊喜：① 它是**实值光滑函数**；② 它只含 $\varepsilon^2$，所以关于 $\varepsilon$ 是偶函数，**$\varepsilon \to 0$ 的双侧极限天然相同！**

接下来怎么证明它在 $\mathcal{D}'(\mathbb{R})$ 中收敛（并且明确算出极限是 $2\operatorname{p.v.}\frac{1}{x}$）？提供**两种满分证法**：
- **方法一（最优雅，3 步秒杀）**：一眼认出 $\frac{2x}{x^2+\varepsilon^2} = \frac{\mathrm{d}}{\mathrm{d}x}\ln(x^2 + \varepsilon^2)$！只要用控制收敛定理证原函数 $\ln(x^2+\varepsilon^2) \to 2\ln|x|$ 于 $L^1_{\mathrm{loc}}(\mathbb{R})$，立刻推得在 $\mathcal{D}'(\mathbb{R})$ 中收敛，再用「广义导数与分布极限可交换」及第 11(2) 题，一行写出极限 $= (2\ln|x|)' = 2\operatorname{p.v.}\frac{1}{x} \in \mathcal{D}'(\mathbb{R})$！
- **方法二（完全自包含，不依赖第 11 题）**：利用奇函数积分 $\int_{-R}^R \frac{2x}{x^2+\varepsilon^2}\varphi(0)\,\mathrm{d}x = 0$，把 $\varphi(x)$ 换成 $\varphi(x)-\varphi(0)$，瞬间用中值定理控住 $\frac{2x}{x^2+\varepsilon^2}(\varphi(x)-\varphi(0))$，直接控制收敛取极限并按式 (8.3) 验证极限泛函属于 $\mathcal{D}'(\mathbb{R})$！

#### 【规范解答】
先化简 $f_\varepsilon(x)$（$\varepsilon \in \mathbb{R}\setminus\{0\}$）：
$$f_\varepsilon(x) = \frac{1}{x+\mathrm{i}\varepsilon} + \frac{1}{x-\mathrm{i}\varepsilon} = \frac{2x}{x^2 + \varepsilon^2}.$$
显然对任意 $\varepsilon \neq 0$，$f_\varepsilon$ 为实值光滑函数，故定义 $\mathcal{D}'(\mathbb{R})$ 中的常规广义函数；又 $f_{-\varepsilon} = f_\varepsilon$，故 $\varepsilon \to 0$ 的双侧极限等价于 $\varepsilon \to 0^+$。下面给出两种严密证明（考试任选其一即可）。

---

##### 方法一：原函数 + 广义导数的连续性（最简洁）

1. **构造原函数**：对 $\varepsilon \neq 0$，令 $F_\varepsilon(x) \triangleq \ln(x^2 + \varepsilon^2) \in C^\infty(\mathbb{R})$。由经典求导与例 8.8，在 $\mathcal{D}'(\mathbb{R})$ 中有
   $$F_\varepsilon' = \frac{2x}{x^2 + \varepsilon^2} = f_\varepsilon.$$

2. **证明 $F_\varepsilon \to 2\ln|x|$ 于 $L^1_{\mathrm{loc}}(\mathbb{R})$**：
   任取紧区间 $[-R, R]$（$R > 0$），考虑 $0 < |\varepsilon| \le 1$。
   - 对每个 $x \in [-R, R]\setminus\{0\}$，有 $\lim_{\varepsilon\to 0}\ln(x^2 + \varepsilon^2) = \ln(x^2) = 2\ln|x|$（几乎处处收敛）。
   - 对任意 $x \in [-R, R]\setminus\{0\}$ 及 $0 < |\varepsilon| \le 1$，有 $x^2 \le x^2 + \varepsilon^2 \le R^2 + 1$。因 $\ln t$ 单调递增，故 $\ln(x^2) \le \ln(x^2+\varepsilon^2) \le \ln(R^2+1)$，从而
     $$\big|\ln(x^2 + \varepsilon^2)\big| \le 2|\ln|x|| + \ln(R^2 + 1) \in L^1([-R, R]).$$
   由 Lebesgue 控制收敛定理，得
   $$\lim_{\varepsilon\to 0}\int_{-R}^{R}\big|F_\varepsilon(x) - 2\ln|x|\big|\,\mathrm{d}x = 0,$$
   即 $F_\varepsilon \to 2\ln|x|$ 于 $L^1_{\mathrm{loc}}(\mathbb{R})$。按教材 p.305，$L^1_{\mathrm{loc}}(\mathbb{R})$ 收敛蕴含 $\mathcal{D}'(\mathbb{R})$ 收敛，因此 $F_\varepsilon \to 2\ln|x|$ 于 $\mathcal{D}'(\mathbb{R})$。

3. **两边取广义导数**：由广义导数在 $\mathcal{D}'(\mathbb{R})$ 上的连续性（教材 p.307）及第 11(2) 题结论 $(\ln|x|)' = \operatorname{p.v.}\frac{1}{x}$，对任意 $\varphi \in \mathcal{D}(\mathbb{R})$ 有
   $$\lim_{\varepsilon\to 0} f_\varepsilon(\varphi) = \lim_{\varepsilon\to 0} F_\varepsilon'(\varphi) = -\lim_{\varepsilon\to 0} F_\varepsilon(\varphi') = -(2\ln|x|)(\varphi') = (2\ln|x|)'(\varphi) = \left(2\operatorname{p.v.}\frac{1}{x}\right)(\varphi).$$
   由例 8.5 知 $2\operatorname{p.v.}\frac{1}{x} \in \mathcal{D}'(\mathbb{R})$，故 $\lim_{\varepsilon\to 0} f_\varepsilon = 2\operatorname{p.v.}\frac{1}{x} \in \mathcal{D}'(\mathbb{R})$。

---

##### 方法二：奇函数减常值 + 直接控制收敛（完全自包含）

任取 $\varphi \in \mathcal{D}(\mathbb{R})$，设 $\operatorname{supp}\varphi \subset [-R, R]$（$R > 0$）。因 $\frac{2x}{x^2+\varepsilon^2}$ 是奇函数，有 $\int_{-R}^{R}\frac{2x}{x^2+\varepsilon^2}\varphi(0)\,\mathrm{d}x = 0$，于是
$$f_\varepsilon(\varphi) = \int_{-R}^{R}\frac{2x}{x^2+\varepsilon^2}\varphi(x)\,\mathrm{d}x = \int_{-R}^{R}\frac{2x}{x^2+\varepsilon^2}\big(\varphi(x) - \varphi(0)\big)\,\mathrm{d}x.$$
对所有 $x \in [-R, R]\setminus\{0\}$ 及 $\varepsilon \neq 0$，由 $\frac{|x|}{x^2+\varepsilon^2} \le \frac{1}{|x|}$ 与中值定理 $|\varphi(x)-\varphi(0)| \le |x|\|\varphi'\|_\infty$，得统一可积控制：
$$\left|\frac{2x}{x^2+\varepsilon^2}\big(\varphi(x) - \varphi(0)\big)\right| \le \frac{2}{|x|}\cdot |x|\|\varphi'\|_\infty = 2\|\varphi'\|_\infty \in L^1([-R, R]).$$
且对每个 $x \neq 0$，被积函数逐点收敛于 $2\frac{\varphi(x)-\varphi(0)}{x}$。由 Lebesgue 控制收敛定理，得
$$\lim_{\varepsilon\to 0} f_\varepsilon(\varphi) = 2\int_{-R}^{R}\frac{\varphi(x) - \varphi(0)}{x}\,\mathrm{d}x \triangleq \lambda(\varphi).$$
最后验证 $\lambda \in \mathcal{D}'(\mathbb{R})$：
- $\lambda$ 显然是线性泛函；
- 对任意紧集 $K \subset [-R, R]$ 及 $\varphi \in C_c^\infty(K)$，有
  $$|\lambda(\varphi)| \le 2\int_{-R}^{R}\frac{|\varphi(x)-\varphi(0)|}{|x|}\,\mathrm{d}x \le 2\int_{-R}^{R}\|\varphi'\|_\infty\,\mathrm{d}x = 4R\|\varphi'\|_\infty \le 4R\|\varphi\|_{C^1(K)},$$
  由教材 p.304 式 (8.3) 立得 $\lambda \in \mathcal{D}'(\mathbb{R})$（且由对称挖洞即知 $\lambda = 2\operatorname{p.v.}\frac{1}{x}$）。证毕。

#### 【考点复盘与避坑提示】
- **知识串联（Sokhotski-Plemelj 公式）**：把第 9(1) 题和本题拼在一起看，你会看到数学物理中最著名的公式之一：
  $$\frac{1}{x \pm \mathrm{i}\varepsilon} = \frac{x}{x^2+\varepsilon^2} \mp \mathrm{i}\frac{\varepsilon}{x^2+\varepsilon^2} \xrightarrow{\varepsilon\to 0^+} \operatorname{p.v.}\frac{1}{x} \mp \mathrm{i}\pi\delta.$$
  本题把两者相加，虚部的 $\mp\mathrm{i}\pi\delta$ 恰好抵消，只留下实部的 $2\operatorname{p.v.}\frac{1}{x}$！
- **避坑提示**：用方法一控制 $|\ln(x^2+\varepsilon^2)|$ 时，必须把**正负两端都控住**（下界用 $\ln(x^2)$ 控原点附近，上界用 $\ln(R^2+1)$ 控远离原点处），千万别直接错写成 $|\ln(x^2+\varepsilon^2)| \le |\ln(x^2)|$（当 $x=1, \varepsilon=1$ 时左边是 $\ln 2$，右边是 $0$，不等式反了！）。

---

### 题目 5：习题 8 第 17 题 (2)(3)（教材 p.342）

#### 【原题】
计算下列函数的 Fourier 变换：
(2) $\dfrac{1}{a^2 + x^2}$；
(3) $\mathrm{e}^{-ax^2}$；

#### 【零基础破题思路】
这两题都是**不依赖复变函数留数定理**、纯靠大一微积分就能优雅搞定的经典计算！全程牢记本书式 (8.9) 约定：$\hat{f}(\xi) = \int_{-\infty}^{\infty} f(x)\mathrm{e}^{-2\pi\mathrm{i}x\xi}\,\mathrm{d}x$。

- **(2) 借花献佛（先算指数衰减再反演）**：直接积分 $\int \frac{\mathrm{e}^{-2\pi\mathrm{i}x\xi}}{a^2+x^2}\,\mathrm{d}x$ 需要留数法，但反过来算指数函数 $g(t) = \mathrm{e}^{-2\pi |a||t|}$ 的 Fourier 变换却只是大一微积分的 $\int \mathrm{e}^{\alpha t}\,\mathrm{d}t$！一算出来居然恰好等于 $\frac{|a|}{\pi(a^2+x^2)}$，再用 Fourier 反演公式一翻身，立刻秒得答案！
- **(3) 求导建立一阶常微分方程（ODE 大法）**：令 $f(x) = \mathrm{e}^{-ax^2}$（$a>0$），它满足最简单的微分关系 $f'(x) = -2ax f(x)$。两边做 Fourier 变换（命题 8.14 与 8.15），“乘 $x$”变成“对 $\xi$ 求导”，“对 $x$ 求导”变成“乘 $2\pi\mathrm{i}\xi$”，直接蹦出关于 $\hat{f}(\xi)$ 的一阶齐次 ODE，10 秒解完再代初值 $\hat{f}(0) = \sqrt{\frac{\pi}{a}}$！

#### 【规范解答】
均采用教材 p.312 式 (8.9) 的 Fourier 变换约定：$\hat{f}(\xi) = \int_{-\infty}^{\infty} f(x)\mathrm{e}^{-2\pi\mathrm{i}x\xi}\,\mathrm{d}x$。

---

**(2) 计算 $f(x) = \dfrac{1}{a^2 + x^2}$ 的 Fourier 变换（设实数 $a \neq 0$；注：$a = 0$ 时 $\frac{1}{x^2} \notin L^1_{\mathrm{loc}}(\mathbb{R})$）**

构造辅助函数 $g(t) \triangleq \mathrm{e}^{-2\pi |a||t|}$。因 $|a| > 0$，显然 $g \in L^1(\mathbb{R}) \cap L^2(\mathbb{R})$ 且连续。按定义分段直接计算 $\hat{g}(x)$：
$$\begin{aligned}
\hat{g}(x) &= \int_{-\infty}^{\infty}\mathrm{e}^{-2\pi |a||t|}\mathrm{e}^{-2\pi\mathrm{i}tx}\,\mathrm{d}t \\
&= \int_{-\infty}^{0}\mathrm{e}^{2\pi(|a| - \mathrm{i}x)t}\,\mathrm{d}t + \int_{0}^{\infty}\mathrm{e}^{-2\pi(|a| + \mathrm{i}x)t}\,\mathrm{d}t \\
&= \left.\frac{\mathrm{e}^{2\pi(|a| - \mathrm{i}x)t}}{2\pi(|a| - \mathrm{i}x)}\right|_{-\infty}^{0} + \left.\frac{\mathrm{e}^{-2\pi(|a| + \mathrm{i}x)t}}{-2\pi(|a| + \mathrm{i}x)}\right|_{0}^{+\infty} \\
&= \frac{1}{2\pi(|a| - \mathrm{i}x)} + \frac{1}{2\pi(|a| + \mathrm{i}x)} \\
&= \frac{1}{2\pi}\cdot\frac{2|a|}{a^2 + x^2} = \frac{|a|}{\pi(a^2 + x^2)}.
\end{aligned}$$
注意到 $\hat{g} \in L^1(\mathbb{R}) \cap L^2(\mathbb{R})$。由 Plancherel 定理（定理 8.9）及 $\mathcal{S}(\mathbb{R})$ 在 $L^2(\mathbb{R})$ 中的稠密性（命题 8.16(6)），命题 8.17 中的反演关系 $\mathcal{F}^*\mathcal{F} = I$ 延拓到整个 $L^2(\mathbb{R})$；又由注记 8.8，对 $L^1 \cap L^2$ 函数 $L^2$-Fourier 变换与积分定义一致，故几乎处处成立
$$g(\xi) = \int_{-\infty}^{\infty}\hat{g}(x)\mathrm{e}^{2\pi\mathrm{i}x\xi}\,\mathrm{d}x = \frac{|a|}{\pi}\int_{-\infty}^{\infty}\frac{1}{a^2 + x^2}\mathrm{e}^{2\pi\mathrm{i}x\xi}\,\mathrm{d}x.$$
因左右两端都是 $\mathbb{R}$ 上的连续函数（右端由控制收敛定理连续），故上式对一切 $\xi \in \mathbb{R}$ 逐点成立。将 $\xi$ 替换为 $-\xi$ 并利用 $g(-\xi) = g(\xi)$，即得
$$\widehat{\frac{1}{a^2 + x^2}}(\xi) = \int_{-\infty}^{\infty}\frac{1}{a^2 + x^2}\mathrm{e}^{-2\pi\mathrm{i}x\xi}\,\mathrm{d}x = \frac{\pi}{|a|}\mathrm{e}^{-2\pi |a||\xi|}.$$
特别地，当 $a > 0$ 时，结果为 $\dfrac{\pi}{a}\mathrm{e}^{-2\pi a|\xi|}$。

---

**(3) 计算 $f(x) = \mathrm{e}^{-ax^2}$ 的 Fourier 变换（设 $a > 0$）**

*参数说明*：若 $a < 0$，函数超指数增长，既不属于 $L^1(\mathbb{R})$ 也不属于 $\mathcal{S}'(\mathbb{R})$；若 $a = 0$，常数函数 $1$ 的缓增分布 Fourier 变换为 $\delta$（命题 8.21(7)）。以下求解标准情形 $a > 0$。

当 $a > 0$ 时，$f(x) = \mathrm{e}^{-ax^2} \in \mathcal{S}(\mathbb{R})$，故 $f, xf, f' \in L^1(\mathbb{R})$ 且 $f'(x) = -2ax f(x)$，即 $x f(x) = -\frac{1}{2a}f'(x)$。
依次应用教材 p.312 命题 8.14 与 p.313 命题 8.15，得
$$\hat{f}'(\xi) = -2\pi\mathrm{i}\,\widehat{xf}(\xi) = -2\pi\mathrm{i}\left(-\frac{1}{2a}\right)\widehat{f'}(\xi) = \frac{\pi\mathrm{i}}{a}\cdot (2\pi\mathrm{i}\xi)\hat{f}(\xi) = -\frac{2\pi^2}{a}\xi\hat{f}(\xi).$$
（若直接分部积分，同样有 $\hat{f}'(\xi) = \frac{\pi\mathrm{i}}{a}\int_{-\infty}^{\infty}(\mathrm{e}^{-ax^2})'\mathrm{e}^{-2\pi\mathrm{i}x\xi}\,\mathrm{d}x = -\frac{\pi\mathrm{i}}{a}(-2\pi\mathrm{i}\xi)\hat{f}(\xi) = -\frac{2\pi^2}{a}\xi\hat{f}(\xi)$。）

由此得到一阶常微分方程 $\hat{f}'(\xi) + \frac{2\pi^2}{a}\xi\hat{f}(\xi) = 0$。两边同乘积分因子 $\mathrm{e}^{\frac{\pi^2\xi^2}{a}}$，得
$$\frac{\mathrm{d}}{\mathrm{d}\xi}\left(\mathrm{e}^{\frac{\pi^2\xi^2}{a}}\hat{f}(\xi)\right) = 0 \implies \hat{f}(\xi) = \hat{f}(0)\mathrm{e}^{-\frac{\pi^2\xi^2}{a}}.$$
令 $\xi = 0$ 计算初值（作换元 $t = \sqrt{a}\,x$）：
$$\hat{f}(0) = \int_{-\infty}^{\infty}\mathrm{e}^{-ax^2}\,\mathrm{d}x = \frac{1}{\sqrt{a}}\int_{-\infty}^{\infty}\mathrm{e}^{-t^2}\,\mathrm{d}t = \sqrt{\frac{\pi}{a}}.$$
代入即得最终结果：
$$\widehat{\mathrm{e}^{-ax^2}}(\xi) = \sqrt{\frac{\pi}{a}}\,\mathrm{e}^{-\frac{\pi^2\xi^2}{a}} \quad (a > 0).$$

#### 【考点复盘与避坑提示】
- **快速自检两绝招（考场防算错常数！）**：
  1. **代 $\xi = 0$ 检验**：(2) 左边 $\int_{-\infty}^\infty \frac{\mathrm{d}x}{a^2+x^2} = \frac{\pi}{|a|}$，与右边 $\frac{\pi}{|a|}\mathrm{e}^0$ 吻合；(3) 取 $a = \pi$ 得 $\widehat{\mathrm{e}^{-\pi x^2}}(\xi) = \mathrm{e}^{-\pi\xi^2}$（Gauss 函数是 Fourier 变换的不动点，和教材 p.314 证明完全一致！）。
  2. **若对比其他不带 $2\pi$ 的教材**（$\tilde{f}(\omega) = \int f(x)\mathrm{e}^{-\mathrm{i}x\omega}\,\mathrm{d}x$）：把结果中的 $\xi$ 换成 $\frac{\omega}{2\pi}$，两题分别对应 $\frac{\pi}{|a|}\mathrm{e}^{-|a||\omega|}$ 和 $\sqrt{\frac{\pi}{a}}\mathrm{e}^{-\frac{\omega^2}{4a}}$。**切忌把别的教材公式直接抄到本题试卷上！**

---

### 题目 6：习题 8 第 20 题（教材 p.342）

#### 【原题】
利用 Fourier 变换求解以下方程：
$$\begin{cases}
-\Delta u(x,y) = 0, & (x,y) \in \mathbb{R}_+^2 = (-\infty,\infty) \times (0,\infty), \\
u(x,0) = g(x), & x \in \mathbb{R}, \\
u(x,y) \to 0, & \text{若 } |(x,y)| \to \infty.
\end{cases} \tag{8.50}$$

#### 【零基础破题思路】
这是偏微分方程里最经典的**「半空间偏 Fourier 变换法」**！记住一个原则：**哪个变量跑满整个 $(-\infty, +\infty)$，就对哪个变量做 Fourier 变换！**
- 这里 $x \in (-\infty, \infty)$ 跑满整条实轴，而 $y \in (0, \infty)$ 只有半条轴。
- 所以我们对 $x$ 做一维 Fourier 变换（把 $y$ 当普通参数看）：
  1. $u_{xx}$ 瞬间变成代数项 $-4\pi^2\xi^2 \hat{u}$，原本的偏微分方程 $u_{xx}+u_{yy}=0$ 立刻退化成关于 $y$ 的大一常微分方程 $\hat{u}_{yy} - (2\pi\xi)^2 \hat{u} = 0$！
  2. 写出通解 $C_1\mathrm{e}^{-2\pi|\xi|y} + C_2\mathrm{e}^{2\pi|\xi|y}$，靠 $y\to+\infty$ 不爆炸砍掉 $C_2$，靠 $y=0$ 边值定出 $C_1 = \hat{g}(\xi)$，得到 $\hat{u}(\xi, y) = \hat{g}(\xi)\mathrm{e}^{-2\pi y|\xi|}$；
  3. **梦幻联动第 17(2) 题**：因子 $\mathrm{e}^{-2\pi y|\xi|}$ 刚好是 Poisson 核 $P_y(x) = \frac{1}{\pi}\frac{y}{x^2+y^2}$ 的 Fourier 变换！频域相乘 = 物理空间卷积，直接写出上半平面 Poisson 积分公式！
  4. **再联动第 9(1) 题**：当 $y \to 0^+$ 时 $P_y \to \delta$，所以 $P_y * g \to g$ 自动还原边值！第 9(1)、17(2)、20 三题首尾呼应，形成完美闭环！

#### 【规范解答】
假定边值函数 $g$ 满足适当的正则与衰减条件（例如 $g \in \mathcal{S}(\mathbb{R})$，或 $g \in L^1(\mathbb{R}) \cap C_0(\mathbb{R})$），且待求经典解 $u(x,y)$ 及其关于 $x$ 的偏导数充分衰减以保证 Fourier 变换与积分号下求导合法。

1. **对水平变量 $x$ 作一维 Fourier 变换（按本书式 (8.9) 约定）**
   记
   $$\hat{u}(\xi, y) \triangleq \int_{-\infty}^{\infty} u(x, y)\mathrm{e}^{-2\pi\mathrm{i}x\xi}\,\mathrm{d}x, \qquad \xi \in \mathbb{R},\ y \ge 0.$$
   方程 $-\Delta u = 0$ 即 $u_{xx}(x, y) + u_{yy}(x, y) = 0$。两边关于 $x$ 取 Fourier 变换，由命题 8.15 得 $\widehat{u_{xx}}(\xi, y) = (2\pi\mathrm{i}\xi)^2 \hat{u}(\xi, y) = -4\pi^2\xi^2 \hat{u}(\xi, y)$，且 $\widehat{u_{yy}}(\xi, y) = \frac{\partial^2 \hat{u}}{\partial y^2}(\xi, y)$，从而化为以 $\xi$ 为参数的关于 $y$ 的二阶常微分方程：
   $$\frac{\partial^2 \hat{u}}{\partial y^2}(\xi, y) - (2\pi\xi)^2 \hat{u}(\xi, y) = 0, \qquad y > 0.$$

2. **求解常微分方程并利用边界与无穷远条件定常数**
   对固定的 $\xi \neq 0$，特征方程 $r^2 - (2\pi\xi)^2 = 0$ 的两根为 $r = \pm 2\pi|\xi|$，故通解为
   $$\hat{u}(\xi, y) = C_1(\xi)\mathrm{e}^{-2\pi|\xi|y} + C_2(\xi)\mathrm{e}^{2\pi|\xi|y}.$$
   - 当 $y \to +\infty$ 时，$\mathrm{e}^{2\pi|\xi|y} \to +\infty$。由题设无穷远条件 $u(x,y) \to 0$（要求 $\hat{u}(\xi, y)$ 在 $y \to +\infty$ 时有界），必有 $C_2(\xi) = 0$。
   - 再由 $y = 0$ 处的边值条件 $u(x, 0) = g(x)$，两边关于 $x$ 取 Fourier 变换得 $\hat{u}(\xi, 0) = \hat{g}(\xi)$，从而 $C_1(\xi) = \hat{g}(\xi)$。

   因此（由连续性包含 $\xi = 0$），
   $$\hat{u}(\xi, y) = \hat{g}(\xi)\mathrm{e}^{-2\pi y|\xi|}, \qquad \xi \in \mathbb{R},\ y > 0.$$

3. **Fourier 反演与卷积定理求显式解**
   对任意固定的 $y > 0$，引入上半平面 **Poisson 核**：
   $$P_y(x) \triangleq \frac{1}{\pi}\frac{y}{x^2 + y^2}.$$
   由习题 8 第 17(2) 题（令 $a = y > 0$），$\widehat{\frac{1}{x^2+y^2}}(\xi) = \frac{\pi}{y}\mathrm{e}^{-2\pi y|\xi|}$，因此恰好有
   $$\widehat{P_y}(\xi) = \frac{y}{\pi}\cdot\frac{\pi}{y}\mathrm{e}^{-2\pi y|\xi|} = \mathrm{e}^{-2\pi y|\xi|}.$$
   于是
   $$\hat{u}(\xi, y) = \hat{g}(\xi)\widehat{P_y}(\xi) = \widehat{(P_y * g)}(\xi).$$
   由 Fourier 反演公式，得边值问题 (8.50) 的解（**上半平面 Poisson 积分公式**）：
   $$u(x, y) = (P_y * g)(x) = \frac{y}{\pi}\int_{-\infty}^{\infty}\frac{g(t)}{(x - t)^2 + y^2}\,\mathrm{d}t, \qquad (x, y) \in \mathbb{R}_+^2.$$

4. **解的合理性简明核验（对 $g \in C_0(\mathbb{R})$，即连续且 $\lim_{|t|\to\infty} g(t) = 0$）**
   - **调和性**：当 $y > 0$ 时，直接求偏导易验证 $(\partial_{xx} + \partial_{yy})P_y(x-t) = 0$，且各阶偏导在 $y \ge y_0 > 0$ 上可积，积分号下求导即得 $-\Delta u = 0$。
   - **边值还原**：当 $y \to 0^+$ 时，如第 9(1) 题所示 $P_y \to \delta$（且 $\int P_y = 1, P_y > 0$），由 $g$ 的连续有界性易证 $u(x,y) \to g(x)$。
   - **无穷远趋于 $0$**：任给 $\eta > 0$，取 $M > 0$ 使 $|t| > M$ 时 $|g(t)| < \eta$。将积分拆成 $|t| > M$ 与 $|t| \le M$ 两部分：外段绝对值不超过 $\eta\int_{\mathbb{R}} P_y(x-t)\,\mathrm{d}t = \eta$；内段利用初等不等式 $\frac{y}{(x-t)^2+y^2} \le \min\!\left\{\frac{1}{y},\ \frac{1}{2|x-t|}\right\}$ 可知当 $|(x,y)| \to \infty$ 时内段一致趋于 $0$。故 $\lim_{|(x,y)|\to\infty} u(x,y) = 0$。

#### 【考点复盘与避坑提示】
- **核心考点**：用 Fourier 变换求解常系数偏微分方程的标准四步（变换降维 $\to$ 解 ODE $\to$ 边界与衰减定常数 $\to$ 卷积反演）。
- **避坑提示**：写 $\mathrm{e}^{-2\pi|\xi|y}$ 时**绝对不能漏掉 $|\xi|$ 上的绝对值**！如果错写成 $\mathrm{e}^{-2\pi\xi y}$，当 $\xi < 0$ 时它反而随 $y \to +\infty$ 指数爆炸，全盘皆错！

---

### 题目 7：习题 8 第 28 题（教材 p.343）

#### 【原题】
设 $k \in \mathbb{N}$。求使得 $\dfrac{1}{|x|^k} \in H^{-m}(\mathbb{R}^n)$ 的最小自然数 $m$。

#### 【零基础破题思路】
这道题是全章最能体现**实指数 Sobolev 空间威力**的压轴好题，也是最容易踩「漏查一半条件」大坑的题！请务必分清**两件事**：

1. **先查入门资格**：函数 $f_k(x) = \frac{1}{|x|^k}$ 要成为 $\mathbb{R}^n$ 上合法的常规缓增分布（即局部可积 $L^1_{\mathrm{loc}}(\mathbb{R}^n)$），在原点用球坐标算一下：
   $$\int_{|x|<1}\frac{\mathrm{d}x}{|x|^k} = \omega_{n-1}\int_0^1 r^{n-1-k}\,\mathrm{d}r < \infty \iff k < n.$$
   （若 $k = 0$ 则 $f_0 \equiv 1$，$\hat{1} = \delta$ 不是函数，不在任何 $H^{-m}$ 里；若 $k \ge n$，它在原点不可积，作为普通函数根本不能定义常规广义函数，后面我们还会单独讨论其正则化情形！）

2. **用 $H^{-m}(\mathbb{R}^n)$ 定义直击频域（教材 p.336）**：
   $$f_k \in H^{-m}(\mathbb{R}^n) \iff \int_{\mathbb{R}^n} (1 + |\xi|^2)^{-m} |\widehat{f_k}(\xi)|^2\,\mathrm{d}\xi < \infty.$$
   关键是要算出 $\widehat{|x|^{-k}}(\xi)$！下面我们用第 17(3) 题的 Gauss 变换配合 $\Gamma$ 函数**4 行严密推导出** $\widehat{|x|^{-k}}(\xi) = c_{n,k}|\xi|^{k-n}$（零黑箱！）。代入后转球坐标，积分变成：
   $$\int_0^\infty (1 + r^2)^{-m} r^{2k - n - 1}\,\mathrm{d}r < \infty.$$
   **⚠️ 超级防坑关键**：这个积分在 $r \to 0^+$ 和 $r \to +\infty$ **两端都有可能发散**，必须两边同时检查！
   - **低频端 $r \to 0^+$**（对应物理空间无穷远 $|x|\to\infty$ 的衰减）：需要 $2k - n - 1 > -1 \iff \mathbf{k > \frac{n}{2}}$！这一关和 $m$ 无关——也就是说，**当 $k \le \frac{n}{2}$ 时，不管 $m$ 取多大都无解！**
   - **高频端 $r \to +\infty$**（对应原点 $x = 0$ 的奇性）：需要 $2k - n - 1 - 2m < -1 \iff \mathbf{m > k - \frac{n}{2}}$！由此立刻定出最小自然数 $m = \left\lfloor k - \frac{n}{2}\right\rfloor + 1$！

#### 【规范解答】
按教材惯例记自然数集 $\mathbb{N} = \{0, 1, 2, \dots\}$（下文可见，凡有解情形求得的最小 $m$ 均满足 $m \ge 1$，故无论 $\mathbb{N}$ 是否含 $0$，有解时的答案完全一致）。记 $S^{n-1}$ 的表面积为 $\omega_{n-1}$。

##### 第一步：先排查 $k = 0$ 与常规缓增分布的存在范围（$1 \le k < n$）
- **若 $k = 0$**：此时 $\frac{1}{|x|^0} \equiv 1$。由命题 8.21(7) 知 $\hat{1} = \delta$，而对任意 $m \in \mathbb{N}$，$\langle\xi\rangle^{-m}\delta = \delta$ 不是局部可积函数，更不属于 $L^2(\mathbb{R}^n)$。故 $k = 0$ 时对任何 $m \in \mathbb{N}$ 均有 $1 \notin H^{-m}(\mathbb{R}^n)$（无解）。
- **检查局部可积性**：对正整数 $k \ge 1$，在单位球上用球坐标得
  $$\int_{|x|<1}\frac{\mathrm{d}x}{|x|^k} = \omega_{n-1}\int_0^1 r^{n-1-k}\,\mathrm{d}r < \infty \iff n - 1 - k > -1 \iff k < n.$$
  而在 $|x| \ge 1$ 上有 $0 < |x|^{-k} \le 1$ 有界。因此当 $1 \le k < n$ 时，$f_k(x) \triangleq \frac{1}{|x|^k} \in L^1_{\mathrm{loc}}(\mathbb{R}^n) \cap \mathcal{S}'(\mathbb{R}^n)$ 是良定义的实值常规缓增分布。

##### 第二步：严格计算 $\widehat{|x|^{-k}}(\xi)$（$1 \le k < n$）
利用 Gamma 函数定义 $\Gamma(\alpha) = \int_0^\infty u^{\alpha-1}\mathrm{e}^{-u}\,\mathrm{d}u$（$\alpha > 0$），取 $\alpha = \frac{k}{2}$ 并作换元 $u = \pi |x|^2 t$（$x \neq 0$），得 Gauss 积分表示：
$$\frac{1}{|x|^k} = \frac{\pi^{k/2}}{\Gamma(k/2)}\int_0^\infty t^{\frac{k}{2}-1}\mathrm{e}^{-\pi t|x|^2}\,\mathrm{d}t.$$
由第 17(3) 题（取 $a = \pi t > 0$ 并对 $n$ 个坐标连乘），$n$ 维 Gauss 函数的 Fourier 变换为
$$\int_{\mathbb{R}^n}\mathrm{e}^{-\pi t|x|^2}\mathrm{e}^{-2\pi\mathrm{i}x\cdot\xi}\,\mathrm{d}x = t^{-n/2}\mathrm{e}^{-\frac{\pi}{t}|\xi|^2}.$$
任取 $\varphi \in \mathcal{S}(\mathbb{R}^n)$。因 $1 \le k < n$，有 $\int_{\mathbb{R}^n}|x|^{-k}|\hat{\varphi}(x)|\,\mathrm{d}x < \infty$ 且（对下述换元 $s = 1/t$ 后）$\int_{\mathbb{R}^n}|\xi|^{k-n}|\varphi(\xi)|\,\mathrm{d}\xi < \infty$，两次应用 Fubini 定理得
$$\begin{aligned}
(\mathcal{F}f_k)(\varphi) &= \int_{\mathbb{R}^n}\frac{1}{|x|^k}\hat{\varphi}(x)\,\mathrm{d}x \\
&= \frac{\pi^{k/2}}{\Gamma(k/2)}\int_0^\infty t^{\frac{k}{2}-1}\left(\int_{\mathbb{R}^n}\mathrm{e}^{-\pi t|x|^2}\hat{\varphi}(x)\,\mathrm{d}x\right)\mathrm{d}t \\
&= \frac{\pi^{k/2}}{\Gamma(k/2)}\int_0^\infty t^{\frac{k-n}{2}-1}\left(\int_{\mathbb{R}^n}\mathrm{e}^{-\frac{\pi}{t}|\xi|^2}\varphi(\xi)\,\mathrm{d}\xi\right)\mathrm{d}t \\
&= \frac{\pi^{k/2}}{\Gamma(k/2)}\int_{\mathbb{R}^n}\left(\int_0^\infty s^{\frac{n-k}{2}-1}\mathrm{e}^{-\pi s|\xi|^2}\,\mathrm{d}s\right)\varphi(\xi)\,\mathrm{d}\xi \\
&= c_{n,k}\int_{\mathbb{R}^n}|\xi|^{k-n}\varphi(\xi)\,\mathrm{d}\xi,
\end{aligned}$$
其中 $s = \frac{1}{t}$，正常数
$$c_{n,k} = \pi^{k - \frac{n}{2}}\frac{\Gamma\!\left(\frac{n-k}{2}\right)}{\Gamma\!\left(\frac{k}{2}\right)} > 0.$$
于是，在 $\mathcal{S}'(\mathbb{R}^n)$ 中有
$$\widehat{|x|^{-k}}(\xi) = c_{n,k}|\xi|^{k-n}.$$

##### 第三步：按 $H^{-m}(\mathbb{R}^n)$ 定义检查两端的积分收敛性（$1 \le k < n$）
按教材 p.336 定义，$f_k \in H^{-m}(\mathbb{R}^n) \iff \langle\xi\rangle^{-m}\widehat{f_k}(\xi) \in L^2(\mathbb{R}^n)$。因 $c_{n,k} > 0$，这等价于如下球坐标积分收敛：
$$I \triangleq \int_{\mathbb{R}^n}(1 + |\xi|^2)^{-m}|\xi|^{2(k-n)}\,\mathrm{d}\xi = \omega_{n-1}\int_0^\infty (1 + r^2)^{-m} r^{2k - n - 1}\,\mathrm{d}r < \infty.$$
将 $(0, \infty)$ 拆分为 $(0, 1]$ 与 $[1, \infty)$ 两段分别检验：

1. **低频端（$r \to 0^+$）**：在 $(0, 1]$ 上，$(1 + r^2)^{-m} \in [2^{-m}, 1]$，故 $\int_0^1 (1+r^2)^{-m} r^{2k-n-1}\,\mathrm{d}r < \infty$ 的充要条件为
   $$2k - n - 1 > -1 \iff k > \frac{n}{2}.$$
   - 特别地，当 $1 \le k \le \frac{n}{2}$ 时，无论自然数 $m$ 取多大，$I = +\infty$，故**不存在任何 $m$ 使得 $\frac{1}{|x|^k} \in H^{-m}(\mathbb{R}^n)$**！

2. **高频端（$r \to +\infty$）**：在 $[1, \infty)$ 上，$(1 + r^2)^{-m} \asymp r^{-2m}$，故被积函数等价于 $r^{2k - n - 1 - 2m}$，积分 $\int_1^\infty (1+r^2)^{-m} r^{2k-n-1}\,\mathrm{d}r < \infty$ 的充要条件为
   $$2k - n - 1 - 2m < -1 \iff 2m > 2k - n \iff m > k - \frac{n}{2}.$$

综合两端条件：当 **$\dfrac{n}{2} < k < n$** 时，$m \in \mathbb{N}$ 满足 $\frac{1}{|x|^k} \in H^{-m}(\mathbb{R}^n)$ 当且仅当 $m > k - \frac{n}{2}$，因此**最小自然数 $m$** 为：
$$m = \left\lfloor k - \frac{n}{2}\right\rfloor + 1 = \begin{cases}
k - \dfrac{n}{2} + 1, & n \text{ 为偶数}, \\[8pt]
k - \dfrac{n - 1}{2}, & n \text{ 为奇数}.
\end{cases}$$

##### 第四步：关于 $k \ge n$ 的完整补全说明
- **作为常规函数**：当 $k \ge n$ 时，$\frac{1}{|x|^k} \notin L^1_{\mathrm{loc}}(\mathbb{R}^n)$，在原点附近不可积，故不能按例 8.2 直接视为 $\mathcal{D}'(\mathbb{R}^n)$ 或 $\mathcal{S}'(\mathbb{R}^n)$ 中的常规广义函数。
- **若理解为正则化缓增分布（Hadamard 有限部 $\operatorname{Pf}|x|^{-k} \in \mathcal{S}'(\mathbb{R}^n)$）**：其 Fourier 变换在 $k - n \notin 2\mathbb{N}$ 时为 $c|\xi|^{k-n}$，在 $k - n \in 2\mathbb{N}$ 时为 $c_1|\xi|^{k-n}\ln|\xi| + c_2|\xi|^{k-n}$。此时因 $k \ge n \implies 2k - n - 1 \ge n - 1 \ge 0 > -1$，低频端自动收敛；高频端 $\int_1^\infty (1+r^2)^{-m} r^{2k-n-1}(\ln r)^{2\text{ 或 }0}\,\mathrm{d}r < \infty$ 的充要条件依然是严格不等式 $m > k - \frac{n}{2}$（对数因子不影响临界严格不等式），故最小自然数仍为
  $$m = \left\lfloor k - \frac{n}{2}\right\rfloor + 1.$$

##### 最终结论汇总
- 当 $k \le \dfrac{n}{2}$（含 $k = 0$）时：**不存在**满足条件的自然数 $m$（因无穷远处衰减不足，低频积分发散）；
- 当 $\dfrac{n}{2} < k < n$ 时（常规缓增分布）：最小自然数为
  $$m = \left\lfloor k - \frac{n}{2}\right\rfloor + 1 = \begin{cases} k - \dfrac{n}{2} + 1, & n \text{ 为偶数}, \\ k - \dfrac{n-1}{2}, & n \text{ 为奇数}; \end{cases}$$
- 当 $k \ge n$ 时：作为普通函数不在 $L^1_{\mathrm{loc}}(\mathbb{R}^n)$ 中；若理解为 Hadamard 有限部正则化缓增分布，最小自然数仍为 $m = \left\lfloor k - \frac{n}{2}\right\rfloor + 1$。

#### 【考点复盘与避坑提示】
- **用具体数字秒验结论**：
  - $n = 3, k = 1$（三维牛顿势 $\frac{1}{|x|}$）：$1 \le \frac{3}{2} \implies$ **无解**（$\widehat{|x|^{-1}} = \frac{1}{\pi|\xi|^2}$，在 $\xi=0$ 处 $|\xi|^{-4} \cdot r^2 = r^{-2}$ 发散！）；
  - $n = 3, k = 2$（$\frac{1}{|x|^2}$）：$\frac{3}{2} < 2 < 3 \implies m > 2 - \frac{3}{2} = 0.5 \implies$ 最小自然数 $m = 1$；
  - $n = 4, k = 3$（$\frac{1}{|x|^3}$）：$2 < 3 < 4 \implies m > 3 - 2 = 1 \implies$ 最小自然数 $m = 2$。
- **为什么 $k \le \frac{n}{2}$ 无解？（物理空间直观验证）**：若取非负截断 $\psi \in C_c^\infty(\{1 \le |x| \le 2\})$ 并放大为 $\psi_R(x) = \psi(x/R)$，则配对 $\int |x|^{-k}\psi_R\,\mathrm{d}x \sim R^{n-k}$，而 $\|\psi_R\|_{H^m} \lesssim R^{n/2}$；当 $k < \frac{n}{2}$ 时比值 $\sim R^{\frac{n}{2}-k} \to \infty$，根本不是 $H^m(\mathbb{R}^n)$ 上的有界线性泛函！
- **Sobolev 嵌入对偶印证高频条件**：为什么原点奇性对应的恰好是 $m > k - \frac{n}{2}$？由定理 8.15，当 $m < \frac{n}{2}$ 时 $H^m \hookrightarrow L^q$（$\frac{1}{q} = \frac{1}{2} - \frac{m}{n}$）；由 Hölder 不等式，要在局部与 $L^q$ 配对，需要 $|x|^{-k} \in L^{q'}(B_1)$，即 $kq' < n \iff \frac{k}{n} < \frac{1}{q'} = 1 - \frac{1}{q} = \frac{1}{2} + \frac{m}{n} \iff m > k - \frac{n}{2}$！频域和物理空间完美印证！

---

### 题目 8：教材第 337 页 练习 8.8.2

#### 【原题】
（教材 p.337 **练习 8.8.2**）如果 $L$ 是一个 $m$ 阶椭圆算子，$f$ 是满足 $f, Lf \in H^s(\mathbb{R}^n)$ 的缓增分布，则 $f \in H^{s+m}(\mathbb{R}^n)$，且有如下估计：
$$\|f\|_{H^{s+m}(\mathbb{R}^n)} \le C\big(\|f\|_{H^s(\mathbb{R}^n)} + \|Lf\|_{H^s(\mathbb{R}^n)}\big), \tag{8.49}$$
其中 $C$ 依赖于 $s, m, n, L$。

> *教材 p.337 上方定义的两处印刷笔误提示*：
> 1. 原书将 $m$ 阶常系数齐次微分算子误印为 $\sum_{j_1+\dots+j_n=m} c_{j_1,\dots,j_n}\frac{\partial^n}{\partial x_{j_1}\cdots\partial x_{j_n}}$，按标准多重指标记法应为 $L = \sum_{|\alpha|=m} c_\alpha \partial^\alpha = \sum_{\alpha_1+\dots+\alpha_n=m} c_\alpha \frac{\partial^m}{\partial x_1^{\alpha_1}\cdots\partial x_n^{\alpha_n}}$，其象征为 $l(\xi) = \sum_{|\alpha|=m} c_\alpha \xi^\alpha$；
> 2. 原书将椭圆性下界印作 $l(\xi) \ge c|\xi|^m$，但结合紧随其后的**例 8.17**（Cauchy-Riemann 算子 $\partial_{x_1} - \mathrm{i}\partial_{x_2}$ 的象征 $l(\xi) = \xi_1 - \mathrm{i}\xi_2$ 为复值函数），正确定义应取绝对值（模长）：存在常数 $c > 0$ 使得对所有 $\xi \in \mathbb{R}^n$ 有 $|l(\xi)| \ge c|\xi|^m$。

#### 【零基础破题思路】
千万不要被「$m$ 阶椭圆算子」「实指数 $H^s$」「缓增分布」这堆吓人的名字唬住——**一转到频域，这道题其实是初中代数不等式！**

看清它的翻译词典：
- 已知 $f \in H^s$ $\iff$ $\langle\xi\rangle^s |\hat{f}(\xi)| \in L^2$（已经控住了**低频** $1 \cdot \langle\xi\rangle^s |\hat{f}|$）；
- 已知 $Lf \in H^s$ $\iff$ $\langle\xi\rangle^s |\widehat{Lf}(\xi)| = \langle\xi\rangle^s (2\pi)^m |l(\xi)|\,|\hat{f}(\xi)| \in L^2$；
- 椭圆性条件 $|l(\xi)| \ge c|\xi|^m$ $\iff$ 用 $|\widehat{Lf}(\xi)|$ 就能反过来压住**高频** $|\xi|^m|\hat{f}(\xi)|$！
- 目标 $f \in H^{s+m}$ $\iff$ 要控住 $\langle\xi\rangle^{s+m}|\hat{f}| = \langle\xi\rangle^s (1 + |\xi|^2)^{m/2} |\hat{f}|$！
- 而初等不等式 $(1 + |\xi|^2)^{m/2} \le 2^{m/2}(1 + |\xi|^m)$ 刚好把它一拆为二：**「1」交给 $f$ 的 $H^s$ 范数，「$|\xi|^m$」交给 $Lf$ 的 $H^s$ 范数**，两边取 $L^2$ 范数用三角不等式，瞬间收官！

#### 【规范解答】
设 $m$ 阶常系数齐次椭圆算子为
$$L = \sum_{|\alpha|=m} c_\alpha \partial^\alpha \quad (m \ge 0),$$
其象征 $l(\xi) = \sum_{|\alpha|=m} c_\alpha \xi^\alpha$ 满足椭圆性条件：存在常数 $c > 0$，使得
$$|l(\xi)| \ge c|\xi|^m, \qquad \forall \xi \in \mathbb{R}^n.$$

1. **将假设翻译为频域可测函数**
   由 $H^s(\mathbb{R}^n)$ 的定义（教材 p.336），$f \in H^s(\mathbb{R}^n)$ 意味着 $f \in \mathcal{S}'(\mathbb{R}^n)$ 且 $\langle\xi\rangle^s \hat{f}(\xi) \in L^2(\mathbb{R}^n)$，其中 $\langle\xi\rangle = (1 + |\xi|^2)^{1/2} \ge 1$。由此 $\hat{f}(\xi) = \langle\xi\rangle^{-s}\big(\langle\xi\rangle^s \hat{f}(\xi)\big)$ 是 $\mathbb{R}^n$ 上的局部 $L^2$（从而局部可积）可测函数；同理，$Lf \in H^s(\mathbb{R}^n)$ 意味着 $\langle\xi\rangle^s \widehat{Lf}(\xi) \in L^2(\mathbb{R}^n)$ 且 $\widehat{Lf}$ 是局部可积可测函数。

2. **计算 $\widehat{Lf}(\xi)$ 并应用椭圆性下界**
   在 $\mathcal{S}'(\mathbb{R}^n)$ 中，反复应用教材 p.318 命题 8.21(6)（$\mathcal{F}(\partial_{x_j}\lambda) = 2\pi\mathrm{i}\xi_j \mathcal{F}\lambda$），对任意 $|\alpha| = m$ 有 $\widehat{\partial^\alpha f}(\xi) = (2\pi\mathrm{i})^m \xi^\alpha \hat{f}(\xi)$。线性叠加得
   $$\widehat{Lf}(\xi) = \sum_{|\alpha|=m} c_\alpha (2\pi\mathrm{i})^m \xi^\alpha \hat{f}(\xi) = (2\pi\mathrm{i})^m l(\xi)\hat{f}(\xi).$$
   因 $l(\xi)$ 为多项式且 $\hat{f}, \widehat{Lf} \in L^1_{\mathrm{loc}}(\mathbb{R}^n)$，上式作为可测函数在 $\mathbb{R}^n$ 上几乎处处成立。两边取模长并利用椭圆性条件 $|l(\xi)| \ge c|\xi|^m$，得
   $$|\widehat{Lf}(\xi)| = (2\pi)^m |l(\xi)|\,|\hat{f}(\xi)| \ge c(2\pi)^m |\xi|^m |\hat{f}(\xi)| \quad \text{a.e. } \xi \in \mathbb{R}^n,$$
   即
   $$|\xi|^m |\hat{f}(\xi)| \le \frac{1}{c(2\pi)^m}|\widehat{Lf}(\xi)| \quad \text{a.e. } \xi \in \mathbb{R}^n.$$

3. **利用初等不等式估计 $\langle\xi\rangle^{s+m}|\hat{f}(\xi)|$**
   对一切 $\xi \in \mathbb{R}^n$，显然 $1 + |\xi|^2 \le 2\max\{1, |\xi|^2\}$。因 $m \ge 0$，两边取 $\frac{m}{2}$ 次幂，得
   $$\langle\xi\rangle^m = (1 + |\xi|^2)^{m/2} \le 2^{m/2}\max\{1, |\xi|^m\} \le 2^{m/2}\big(1 + |\xi|^m\big).$$
   （当 $m = 0$ 时，左端 $= 1$，右端 $= 2$，不等式亦显然成立。）

   将上式与步骤 2 的估计相乘结合，对几乎处处的 $\xi \in \mathbb{R}^n$ 有
   $$\begin{aligned}
   \langle\xi\rangle^{s+m}|\hat{f}(\xi)| &= \langle\xi\rangle^s \langle\xi\rangle^m |\hat{f}(\xi)| \\
   &\le 2^{m/2}\langle\xi\rangle^s |\hat{f}(\xi)| + 2^{m/2}\langle\xi\rangle^s |\xi|^m |\hat{f}(\xi)| \\
   &\le 2^{m/2}\langle\xi\rangle^s |\hat{f}(\xi)| + \frac{2^{m/2}}{c(2\pi)^m}\langle\xi\rangle^s |\widehat{Lf}(\xi)|.
   \end{aligned}$$

4. **取 $L^2(\mathbb{R}^n)$ 范数得出结论**
   由已知 $f, Lf \in H^s(\mathbb{R}^n)$，右端两项均属于 $L^2(\mathbb{R}^n)$。由 $L^2(\mathbb{R}^n)$ 范数的单调性与三角不等式（Minkowski 不等式），立得 $\langle\xi\rangle^{s+m}\hat{f} \in L^2(\mathbb{R}^n)$，即 **$f \in H^{s+m}(\mathbb{R}^n)$**，且满足估计
   $$\begin{aligned}
   \|f\|_{H^{s+m}(\mathbb{R}^n)} &= \|\langle\xi\rangle^{s+m}\hat{f}\|_{L^2(\mathbb{R}^n)} \\
   &\le 2^{m/2}\|\langle\xi\rangle^s \hat{f}\|_{L^2(\mathbb{R}^n)} + \frac{2^{m/2}}{c(2\pi)^m}\|\langle\xi\rangle^s \widehat{Lf}\|_{L^2(\mathbb{R}^n)} \\
   &= 2^{m/2}\|f\|_{H^s(\mathbb{R}^n)} + \frac{2^{m/2}}{c(2\pi)^m}\|Lf\|_{H^s(\mathbb{R}^n)} \\
   &\le C\big(\|f\|_{H^s(\mathbb{R}^n)} + \|Lf\|_{H^s(\mathbb{R}^n)}\big),
   \end{aligned}$$
   其中常数可显式取为
   $$C = 2^{m/2}\max\left\{1,\ \frac{1}{c(2\pi)^m}\right\}.$$
   （该常数只依赖于阶数 $m$ 与 $L$ 的椭圆常数 $c$，自然满足题设要求。）证毕。

#### 【考点复盘与避坑提示】
- **核心考点**：实指数 Sobolev 空间 $H^s(\mathbb{R}^n)$ 的频域范数定义（p.336）、椭圆算子的象征下界与**全局椭圆正则性先验估计（式 (8.49)）**。
- **为什么需要式中的 $\|f\|_{H^s}$ 项，不能只用 $\|Lf\|_{H^s}$ 控制？**
  因为当 $\xi \to 0$ 时，$|\xi|^m \to 0$，$Lf$ 根本管不住低频附近的 $\hat{f}(\xi)$（例如常函数满足 $\Delta f = 0$，但 $f \neq 0$）！所以必须由 $\|f\|_{H^s}$ 兜底低频的「$1$」，由 $\|Lf\|_{H^s}$ 压住高频的「$|\xi|^m$」，二者缺一不可！
- **为什么非椭圆算子不成立（教材 p.337 注记 8.19）？**
  以波算子 $L = -\partial_t^2 + \Delta$ 为例，其象征 $l(\tau, \xi) = \tau^2 - |\xi|^2$ 在特征锥 $|\tau| = |\xi| \neq 0$ 上恒等于 $0$！即使 $Lf = 0$，解 $\hat{f}$ 也可以把任意高频率的能量集中堆在特征锥 $|\tau|=|\xi|$ 上，完全无法由 $Lf$ 控制，因此估计 (8.49) 失效。
- **往深一步（衔接练习 8.8.3）**：若 $Lf = 0$ 且 $f \in H^s(\mathbb{R}^n)$，反复套用本题结论可得 $f \in H^{s+m} \implies H^{s+2m} \implies \cdots \implies \bigcap_{k\ge 1} H^{s+km}$，再由 Sobolev 嵌入定理 8.16 立得 $f \in C^\infty(\mathbb{R}^n)$——这就是 PDE 中大名鼎鼎的**椭圆正则性自举（Bootstrap）思想**！

---

## 第三部分：期末考前 3 分钟冲刺速查清单

1. **做广义函数求导题**：先写「任取 $\varphi \in \mathcal{D}(\mathbb{R})$，设 $\operatorname{supp}\varphi \subset [-R, R]$」$\to$ 甩锅求导加负号 $\to$ 奇点处挖去 $(-\varepsilon, \varepsilon)$ 用控制收敛定理 $\to$ 分部积分用中值定理消掉 $\varepsilon\ln\varepsilon$ 边界项。
2. **做分布极限题**：优先试两种大招——① 换元 + 控制收敛定理（第 9 题）；② 找原函数证 $L^1_{\mathrm{loc}}$ 收敛再求广义导数（第 14 题）。
3. **做 Fourier 变换与 PDE 题**：时刻检查本书 **$2\pi$ 约定**——求导乘 $2\pi\mathrm{i}\xi$，反演与 Plancherel **不带** $(2\pi)^{-n}$；背熟 Gauss 核与 Poisson 核两对变换；半空间方程对无界变量做变换，写衰减通解别忘 $|\xi|$ 的绝对值。
4. **做 $H^s(\mathbb{R}^n)$ 题**：第一时间切到频域看 $\langle\xi\rangle^s \hat{f}(\xi) \in L^2(\mathbb{R}^n)$；查积分收敛时**低频 $r \to 0^+$ 和高频 $r \to +\infty$ 两端都必须查**（第 28 题）；证椭圆估计用 $\langle\xi\rangle^m \le 2^{m/2}(1 + |\xi|^m)$ 一分为二（练习 8.8.2）。
