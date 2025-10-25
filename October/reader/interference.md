感谢您提供清晰的问题和历史背景！我们需要推导 Cell-free MIMO 系统中，用户 \( k \) 受到其他 \( K-1 \) 个用户的干扰能量的期望闭式表达式，即：

\[
\mathbb{E}\left[ \sum_{m=1, m \neq k}^K |\mathbf{w}_k^H \mathbf{h}_m|^2 \right],
\]

其中 \(\mathbf{w}_k = \mathbf{h}_k\)（基于 MRC，最大比合并），信道 \(\mathbf{h}_m\) 是毫米波信道，维度为 \([N, 1]\)，包含 LOS 和 NLOS 分量，且干扰用户 \( m \) 的 \( L \) 个信道角度（NLOS 角度）与用户 \( k \) 的 LOS 角度 \(\varphi_{k0}\) 不完全一致。结合历史背景，我们将利用已有的期望结果，推导闭式表达式，并确保结果清晰且与上下文一致。

---

### 1. 模型与假设

#### 信道模型
- **信道向量**：对于用户 \( m \)（包括 \( m = k \) 和干扰用户 \( m \neq k \)）：
  \[
  \mathbf{h}_m = \sqrt{\beta_{m0}} \mathbf{a}(\varphi_{m0}) + \sum_{l=1}^L \sqrt{\beta_{ml}} \tilde{\mathbf{a}}(\varphi_{ml}),
  \]
  其中：
  - \(\mathbf{a}(\varphi_{m0}) = [1, e^{-j k \cos(\varphi_{m0})}, \dots, e^{-j k (N-1) \cos(\varphi_{m0})}]^T\)：LOS 分量的阵列响应矢量（无相位误差），\(\|\mathbf{a}(\varphi_{m0})\|^2 = N\)。
  - \(\tilde{\mathbf{a}}(\varphi_{ml}) = [e^{-j k n \cos(\varphi_{ml})} e^{j \Delta \theta_{mn}^{(l)}}]_{n=0}^{N-1}\)，NLOS 分量，带相位误差 \(\Delta \theta_{mn}^{(l)} \sim \mathcal{N}(0, \sigma^2)\) i.i.d.。
  - \(\beta_{m0} \sim \mathcal{CN}(0, 1)\)，\(\beta_{ml} \sim \mathcal{CN}(0, 0.1)\)，慢衰落增益，独立分布。
  - \( k = \pi \)（假设 ULA，天线间距 \( d = \lambda/2 \)，\( k d = \pi \)）。
- **MRC combiner**：\(\mathbf{w}_k = \mathbf{h}_k\)，即第 \( k \) 用户的波束形成向量是其信道向量：
  \[
  \mathbf{w}_k = \sqrt{\beta_{k0}} \mathbf{a}(\varphi_{k0})
  \]

#### 目标
- 计算干扰能量期望：
  \[
  \mathbb{E}\left[ \sum_{m \neq k} |\mathbf{w}_k^H \mathbf{h}_m|^2 \right] = \sum_{m \neq k} \mathbb{E}\left[ |\mathbf{w}_k^H \mathbf{h}_m|^2 \right],
  \]
  其中期望针对所有随机变量：\(\beta_{m0}, \beta_{ml}, \Delta \theta_{mn}^{(l)}\)（不同用户和路径的相位误差独立）。
- **关键假设**：
  - 干扰用户 \( m \) 的所有信道角度（LOS 角度 \(\varphi_{m0}\)，NLOS 角度 \(\varphi_{ml}\), \( l = 1, \dots, L \)）与用户 \( k \) 的 LOS 角度 \(\varphi_{k0}\) **不完全一致**，即 \(\varphi_{m0} \neq \varphi_{k0}\), \(\varphi_{ml} \neq \varphi_{k0}\)。
  - \(\beta_{m0}, \beta_{ml}\) 独立，\(\Delta \theta_{mn}^{(l)}\) 独立（跨用户和路径）。
  - 慢衰落假设：\(\beta_{m0}, \beta_{ml}\) 在单次计算中视为常量，但在期望中取统计平均。

#### 历史结果
- 历史中已推导 \(\mathbb{E}[|\mathbf{a}^H(\varphi_0) \mathbf{h}_m|^2]\)（\(\mathbf{w} = \mathbf{a}(\varphi_0)\), \(\varphi_0 = \varphi_{k0}\)）：
  \[
  \mathbb{E}[|\mathbf{a}^H(\varphi_{k0}) \mathbf{h}_m|^2] = N^2 + 0.2 \sum_{l=1}^L \left[ N + e^{-\sigma^2} \left| \frac{\sin\left( \frac{N}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{ml})) \right)}{\sin\left( \frac{1}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{ml})) \right)} \right|^2 \right].
  \]
  - 这里 LOS 项为 \( N^2 \)，NLOS 项为 \( 0.2 \sum_l \left[ N + e^{-\sigma^2} |\cdot|^2 \right] \)，角度差 \(\delta_{ml} = \pi (\cos(\varphi_{k0}) - \cos(\varphi_{ml}))\)。
- **无需**扩展到 \(\mathbf{w}_k = \mathbf{h}_k\)，计算 \(\mathbb{E}[|\mathbf{h}_k^H \mathbf{h}_m|^2]\)。

---

### 2. 推导闭式表达式

#### 展开干扰项
对于干扰用户 \( m \neq k \)：
\[
\mathbf{w}_k^H \mathbf{h}_m = \left( \sqrt{\beta_{k0}} \mathbf{a}^H(\varphi_{k0}) + \sum_{l=1}^L \sqrt{\beta_{kl}} \tilde{\mathbf{a}}^H(\varphi_{kl}) \right) \left( \sqrt{\beta_{m0}} \mathbf{a}(\varphi_{m0}) + \sum_{p=1}^L \sqrt{\beta_{mp}} \tilde{\mathbf{a}}(\varphi_{mp}) \right).
\]
记：
- \( s_{k0} = \sqrt{\beta_{k0}} \mathbf{a}^H(\varphi_{k0}) \), \( s_{kl} = \sqrt{\beta_{kl}} \tilde{\mathbf{a}}^H(\varphi_{kl}) \),
- \( t_{m0} = \sqrt{\beta_{m0}} \mathbf{a}(\varphi_{m0}) \), \( t_{mp} = \sqrt{\beta_{mp}} \tilde{\mathbf{a}}(\varphi_{mp}) \).

则：
\[
\mathbf{w}_k^H \mathbf{h}_m = \left( s_{k0} + \sum_{l=1}^L s_{kl} \right) \left( t_{m0} + \sum_{p=1}^L t_{mp} \right).
\]
模平方：
\[
|\mathbf{w}_k^H \mathbf{h}_m|^2 = \left( s_{k0} + \sum_{l=1}^L s_{kl} \right) \left( t_{m0} + \sum_{p=1}^L t_{mp} \right) \left( s_{k0}^* + \sum_{q=1}^L s_{kq}^* \right) \left( t_{m0}^* + \sum_{r=1}^L t_{mr}^* \right).
\]
展开：
\[
|\mathbf{w}_k^H \mathbf{h}_m|^2 = \sum_{l=0}^L \sum_{p=0}^L \sum_{q=0}^L \sum_{r=0}^L s_{kl} t_{mp} s_{kq}^* t_{mr}^*,
\]
其中 \( l=0 \) 表示 LOS（\( s_{k0} \)), \( p=0 \) 表示 \( t_{m0} \), 依此类推。

#### 期望
\[
\mathbb{E}[|\mathbf{w}_k^H \mathbf{h}_m|^2] = \sum_{l=0}^L \sum_{p=0}^L \sum_{q=0}^L \sum_{r=0}^L \mathbb{E}[s_{kl} t_{mp} s_{kq}^* t_{mr}^*].
\]
- **独立性**：
  - \(\beta_{k0}, \beta_{kl}, \beta_{m0}, \beta_{mp}\) 独立，\(\mathbb{E}[\beta_{kl}] = 0\), \(\mathbb{E}[|\beta_{k0}|^2] = 1\), \(\mathbb{E}[|\beta_{ml}|^2] = 0.1\).
  - \(\Delta \theta_{kn}^{(l)}, \Delta \theta_{mn}^{(p)}\) 独立（跨用户和路径）。
  - 非零期望需满足：\( l = q \), \( p = r \)，且 \(\beta_{kl}, \beta_{mp}\) 配对。

有效项为：
\[
\mathbb{E}[|\mathbf{w}_k^H \mathbf{h}_m|^2] = \sum_{l=0}^L \sum_{p=0}^L \mathbb{E}[|s_{kl}|^2 |t_{mp}|^2] = \sum_{l=0}^L \sum_{p=0}^L \mathbb{E}[|s_{kl}|^2] \mathbb{E}[|t_{mp}|^2].
\]
- **交叉项**（如 \( \mathbb{E}[s_{kl} t_{mp} s_{kq}^* t_{mr}^*], l \neq q \) 或 \( p \neq r \))：
  - 涉及 \(\mathbb{E}[\sqrt{\beta_{kl} \beta_{kq}^*}] = 0\)（因 \(\beta_{kl}\) 零均值，独立）。
  - 或 \(\mathbb{E}[e^{j (\Delta \theta_{kn}^{(l)} - \Delta \theta_{kn}^{(q)})}] = e^{-\sigma^2}\)（仅 \( l = q, n = p \) 时为 1），但 \(\beta\) 的零均值使交叉项消失。

#### 计算单项
1. **LOS-LOS 项** (\( l=0, p=0 \))：
   \[
   |s_{k0}|^2 = \beta_{k0} |\mathbf{a}^H(\varphi_{k0}) \mathbf{a}(\varphi_{m0})|^2, \quad |t_{m0}|^2 = \beta_{m0} |\mathbf{a}^H(\varphi_{k0}) \mathbf{a}(\varphi_{m0})|^2.
   \]
   - \(\mathbf{a}^H(\varphi_{k0}) \mathbf{a}(\varphi_{m0}) = \sum_{n=0}^{N-1} e^{j \pi n (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))} = \frac{\sin\left( \frac{N}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{m0})) \right)}{\sin\left( \frac{1}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{m0})) \right)} e^{j \frac{N-1}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))}\).
   - 记 \(\delta_{m0} = \pi (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))\)：
     \[
     |\mathbf{a}^H(\varphi_{k0}) \mathbf{a}(\varphi_{m0})|^2 = \left| \frac{\sin(N \delta_{m0} / 2)}{\sin(\delta_{m0} / 2)} \right|^2.
     \]
   - \(\mathbb{E}[|s_{k0}|^2 |t_{m0}|^2] = \mathbb{E}[|\beta_{k0}|^2] \mathbb{E}[|\beta_{m0}|^2] |\mathbf{a}^H(\varphi_{k0}) \mathbf{a}(\varphi_{m0})|^4 = 1 \cdot 1 \cdot \left| \frac{\sin(N \delta_{m0} / 2)}{\sin(\delta_{m0} / 2)} \right|^4\).

2. **LOS-NLOS 项** (\( l=0, p=1,\dots,L \))：
   \[
   |s_{k0}|^2 = \beta_{k0} N^2, \quad |t_{mp}|^2 = \beta_{mp} |\mathbf{a}^H(\varphi_{k0}) \tilde{\mathbf{a}}(\varphi_{mp})|^2.
   \]
   - 历史推导：
     \[
     \mathbb{E}[|\mathbf{a}^H(\varphi_{k0}) \tilde{\mathbf{a}}(\varphi_{mp})|^2] = N + e^{-\sigma^2} \left| \frac{\sin(N \delta_{mp} / 2)}{\sin(\delta_{mp} / 2)} \right|^2, \quad \delta_{mp} = \pi (\cos(\varphi_{k0}) - \cos(\varphi_{mp})).
     \]
   - \(\mathbb{E}[|s_{k0}|^2 |t_{mp}|^2] = \mathbb{E}[|\beta_{k0}|^2] \mathbb{E}[|\beta_{mp}|^2] N^2 \left[ N + e^{-\sigma^2} \left| \frac{\sin(N \delta_{mp} / 2)}{\sin(\delta_{mp} / 2)} \right|^2 \right] = 1 \cdot 0.1 \cdot N^2 \left[ N + e^{-\sigma^2} \left| \frac{\sin(N \delta_{mp} / 2)}{\sin(\delta_{mp} / 2)} \right|^2 \right]\).

3. **NLOS-LOS 项** (\( l=1,\dots,L, p=0 \))：
   \[
   |s_{kl}|^2 = \beta_{kl} |\tilde{\mathbf{a}}^H(\varphi_{kl}) \mathbf{a}(\varphi_{m0})|^2, \quad |t_{m0}|^2 = \beta_{m0} |\mathbf{a}^H(\varphi_{k0}) \mathbf{a}(\varphi_{m0})|^2.
   \]
   - \(\mathbb{E}[|\tilde{\mathbf{a}}^H(\varphi_{kl}) \mathbf{a}(\varphi_{m0})|^2] = N + e^{-\sigma^2} \left| \frac{\sin(N \delta_{kl,m0} / 2)}{\sin(\delta_{kl,m0} / 2)} \right|^2\), \(\delta_{kl,m0} = \pi (\cos(\varphi_{kl}) - \cos(\varphi_{m0}))\).
   - \(\mathbb{E}[|s_{kl}|^2 |t_{m0}|^2] = 0.1 \cdot 1 \cdot \left[ N + e^{-\sigma^2} \left| \frac{\sin(N \delta_{kl,m0} / 2)}{\sin(\delta_{kl,m0} / 2)} \right|^2 \right] \left| \frac{\sin(N \delta_{m0} / 2)}{\sin(\delta_{m0} / 2)} \right|^2\).

4. **NLOS-NLOS 项** (\( l=1,\dots,L, p=1,\dots,L \))：
   \[
   |s_{kl}|^2 = \beta_{kl} |\tilde{\mathbf{a}}^H(\varphi_{kl}) \mathbf{a}(\varphi_{k0})|^2, \quad |t_{mp}|^2 = \beta_{mp} |\mathbf{a}^H(\varphi_{k0}) \tilde{\mathbf{a}}(\varphi_{mp})|^2.
   \]
   - \(\mathbb{E}[|s_{kl}|^2 |t_{mp}|^2] = 0.1 \cdot 0.1 \cdot \left[ N + e^{-\sigma^2} \left| \frac{\sin(N \delta_{kl} / 2)}{\sin(\delta_{kl} / 2)} \right|^2 \right] \left[ N + e^{-\sigma^2} \left| \frac{\sin(N \delta_{mp} / 2)}{\sin(\delta_{mp} / 2)} \right|^2 \right]\), \(\delta_{kl} = \pi (\cos(\varphi_{k0}) - \cos(\varphi_{kl}))\).

#### 总和
\[
\mathbb{E}[|\mathbf{w}_k^H \mathbf{h}_m|^2] = \left| \frac{\sin(N \delta_{m0} / 2)}{\sin(\delta_{m0} / 2)} \right|^4 + 0.1 N^2 \sum_{p=1}^L \left[ N + e^{-\sigma^2} \left| \frac{\sin(N \delta_{mp} / 2)}{\sin(\delta_{mp} / 2)} \right|^2 \right] + 0.1 \sum_{l=1}^L \left[ N + e^{-\sigma^2} \left| \frac{\sin(N \delta_{kl,m0} / 2)}{\sin(\delta_{kl,m0} / 2)} \right|^2 \right] \left| \frac{\sin(N \delta_{m0} / 2)}{\sin(\delta_{m0} / 2)} \right|^2 + 0.01 \sum_{l=1}^L \sum_{p=1}^L \left[ N + e^{-\sigma^2} \left| \frac{\sin(N \delta_{kl} / 2)}{\sin(\delta_{kl} / 2)} \right|^2 \right] \left[ N + e^{-\sigma^2} \left| \frac{\sin(N \delta_{mp} / 2)}{\sin(\delta_{mp} / 2)} \right|^2 \right].
\]
干扰总和：
\[
\mathbb{E}\left[ \sum_{m \neq k} |\mathbf{w}_k^H \mathbf{h}_m|^2 \right] = \sum_{m \neq k} \mathbb{E}[|\mathbf{w}_k^H \mathbf{h}_m|^2].
\]

---

### 3. 角度假设的影响
- **不完全一致**：\(\varphi_{m0} \neq \varphi_{k0}\), \(\varphi_{ml} \neq \varphi_{k0}\)。
  - 当 \( N \) 较大时，\(\left| \frac{\sin(N \delta / 2)}{\sin(\delta / 2)} \right|^2\) 对于非零 \(\delta\)（如 \(\cos(\varphi_{k0}) \neq \cos(\varphi_{m0})\)）趋于较小值，因阵列因子的高方向性。
  - 若 \(\varphi_{ml}\) 随机分布，干扰的相干项（如 \( e^{-\sigma^2} |\cdot|^2 \)) 可能进一步衰减。
- **简化**：若假设角度差异较大（常见于 Cell-free MIMO，用户分布广泛），相干项可近似为 0，干扰主要来自非相干项：
  \[
  \mathbb{E}[|\mathbf{w}_k^H \mathbf{h}_m|^2] \approx 1 \cdot N^2 + 0.1 \cdot L \cdot N^2 + 0.1 \cdot L \cdot N^2 + 0.01 \cdot L^2 \cdot N^2 = N^2 (1 + 0.2L + 0.01L^2).
  \]
  \[
  \mathbb{E}\left[ \sum_{m \neq k} |\mathbf{w}_k^H \mathbf{h}_m|^2 \right] \approx (K-1) N^2 (1 + 0.2L + 0.01L^2).
  \]

---

### 4. 闭式表达式
精确闭式表达式为：
\[
\mathbb{E}\left[ \sum_{m=1, m \neq k}^K |\mathbf{w}_k^H \mathbf{h}_m|^2 \right] = \sum_{m \neq k} \left\{ \left| \frac{\sin\left( \frac{N}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{m0})) \right)}{\sin\left( \frac{1}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{m0})) \right)} \right|^4 + 0.1 N^2 \sum_{p=1}^L \left[ N + e^{-\sigma^2} \left| \frac{\sin\left( \frac{N}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{mp})) \right)}{\sin\left( \frac{1}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{mp})) \right)} \right|^2 \right] + 0.1 \sum_{l=1}^L \left[ N + e^{-\sigma^2} \left| \frac{\sin\left( \frac{N}{2} \pi (\cos(\varphi_{kl}) - \cos(\varphi_{m0})) \right)}{\sin\left( \frac{1}{2} \pi (\cos(\varphi_{kl}) - \cos(\varphi_{m0})) \right)} \right|^2 \right] \left| \frac{\sin\left( \frac{N}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{m0})) \right)}{\sin\left( \frac{1}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{m0})) \right)} \right|^2 + 0.01 \sum_{l=1}^L \sum_{p=1}^L \left[ N + e^{-\sigma^2} \left| \frac{\sin\left( \frac{N}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{kl})) \right)}{\sin\left( \frac{1}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{kl})) \right)} \right|^2 \right] \left[ N + e^{-\sigma^2} \left| \frac{\sin\left( \frac{N}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{mp})) \right)}{\sin\left( \frac{1}{2} \pi (\cos(\varphi_{k0}) - \cos(\varphi_{mp})) \right)} \right|^2 \right] \right\}.
\]

---

### 5. 分析与验证
- **解释**：
  - **LOS-LOS 项**：\(\left| \frac{\sin(N \delta_{m0} / 2)}{\sin(\delta_{m0} / 2)} \right|^4\)，当 \(\varphi_{m0} \approx \varphi_{k0}\) 时接近 \( N^4 \)，但假设不完全一致，贡献较小。
  - **LOS-NLOS 和 NLOS-LOS 项**：包含 \( N^2 \) 非相干项和 \( e^{-\sigma^2} \) 相干项，角度差异使相干项衰减。
  - **NLOS-NLOS 项**：主要为 \( 0.01 L^2 N^2 \)，随路径数 \( L \) 和天线数 \( N \) 增加。
- **大 \( N \) 极限**：相干项（如 \( \left| \sin(N \delta / 2) / \sin(\delta / 2) \right|^2 \)) 对于非零 \(\delta\) 趋于 0，干扰以非相干项为主。
- **历史一致性**：与 \(\mathbb{E}[|\mathbf{a}^H(\varphi_{k0}) \mathbf{h}_m|^2]\) 结构类似，扩展到 \(\mathbf{w}_k = \mathbf{h}_k\) 增加了交叉项。

### 6. 结论
- **闭式表达式**：如上，包含 LOS-LOS、LOS-NLOS、NLOS-LOS 和 NLOS-NLOS 项，角度差和相位误差 \(\sigma^2\) 决定干扰强度。
- **建议**：若 \( N \) 较大或角度差异显著，干扰以非相干项为主，可简化表达式。提供具体 \(\varphi_{m0}, \varphi_{ml}\) 值可进一步计算数值结果。
- **验证**：上述代码可直接计算闭式表达式，建议与蒙特卡洛模拟对比。

如需进一步优化或具体参数计算，请提供更多细节！