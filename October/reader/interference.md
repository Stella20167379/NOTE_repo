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
  \mathbf{h}_m = \beta_{m0} \tilde{\mathbf{a}}(\varphi_{m0}) + \sum_{l=1}^L \beta_{ml} \tilde{\mathbf{a}}(\varphi_{ml}),
  \]
  其中：
  - \(\mathbf{a}(\varphi_{m0}) = [1, e^{-j k \cos(\varphi_{m0})}, \dots, e^{-j k (N-1) \cos(\varphi_{m0})}]^T\)：LOS 分量的阵列响应矢量（无相位误差），\(\|\mathbf{a}(\varphi_{m0})\|^2 = N\)。
  - \(\tilde{\mathbf{a}}(\varphi_{ml}) = [e^{-j k n \cos(\varphi_{ml})} e^{j \Delta \theta_{mn}^{(l)}}]_{n=0}^{N-1}\)，NLOS 分量，带相位误差 \(\Delta \theta_{mn}^{(l)} \sim \mathcal{N}(0, \sigma^2)\) i.i.d.。
  - \(\beta_{m0} \sim \mathcal{CN}(0, 1)\)，\(\beta_{ml} \sim \mathcal{CN}(0, 0.1)\)，慢衰落增益，独立分布。
  - \( k = \pi \)（假设 ULA，天线间距 \( d = \lambda/2 \)，\( k d = \pi \)）。
- **MRC combiner**：\(\mathbf{w}_k = \mathbf{h}_k\)，即第 \( k \) 用户的波束形成向量是其信道向量：
  \[
  \mathbf{w}_k = \beta_{k0} \mathbf{a}(\varphi_{k0})
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

---

### 2. 推导闭式表达式

#### 展开干扰项
对于干扰用户 \( m \neq k \)：
\[
\mathbf{w}_k^H \mathbf{h}_m = \beta_{k0}^* \mathbf{a}^H(\varphi_{k0})\left( \beta_{m0} \tilde{\mathbf{a}}(\varphi_{m0}) + \sum_{p=1}^L \beta_{mp} \tilde{\mathbf{a}}(\varphi_{mp}) \right).
\]
记：
- \( s_{k0} = \beta_{k0}^* \mathbf{a}^H(\varphi_{k0}) \). 
- \( t_{m0} = \beta_{m0} \tilde{\mathbf{a}}(\varphi_{m0}) \), \( t_{mp} = \beta_{mp} \tilde{\mathbf{a}}(\varphi_{mp}) \).

则：
\[
\mathbf{w}_k^H \mathbf{h}_m = s_{k0} \left( t_{m0} + \sum_{p=1}^L t_{mp} \right).
\]
模平方：
\[
|\mathbf{w}_k^H \mathbf{h}_m|^2 = s_{k0} \left( t_{m0} + \sum_{p=1}^L t_{mp} \right) s_{k0}^* \left( t_{m0}^* + \sum_{r=1}^L t_{mr}^* \right).
\]

We need to compute:

\[
\mathbb{E}\left[ |\mathbf{w}_k^H \mathbf{h}_m|^2 \right] = \mathbb{E} \left[ \left| s_{k0} \left( t_{m0} + \sum_{p=1}^L t_{mp} \right) \right|^2 \right],
\]

Let’s denote:

\[
z_{m0} = \mathbf{a}^H(\varphi_{k0}) \tilde{\mathbf{a}}(\varphi_{m0}) = \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))} e^{j \Delta \theta_{mn}^{(0)}},
\]

\[
z_{mp} = \mathbf{a}^H(\varphi_{k0}) \tilde{\mathbf{a}}(\varphi_{mp}) = \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_{k0}) - \cos(\varphi_{mp}))} e^{j \Delta \theta_{mn}^{(p)}},
\]

\[
z_{k0} = \mathbf{a}^H(\varphi_{k0}) \tilde{\mathbf{a}}(\varphi_{k0}) = \sum_{n=0}^{N-1} e^{j \Delta \theta_{kn}^{(0)}}.
\]

Thus:

\[
s_{k0}s_{k0}^* = |\beta_{k0}|^2 z_{k0}, \quad t_{m0} = \beta_{m0} \tilde{\mathbf{a}}(\varphi_{m0}), \quad t_{mp} = \beta_{mp} \tilde{\mathbf{a}}(\varphi_{mp}),
\]

However, we need the inner product:

\[
\mathbf{w}_k^H \mathbf{h}_m = \beta_{k0}^* \left( \beta_{m0} \mathbf{a}^H(\varphi_{k0}) \tilde{\mathbf{a}}(\varphi_{m0}) + \sum_{p=1}^L \beta_{mp} \mathbf{a}^H(\varphi_{k0}) \tilde{\mathbf{a}}(\varphi_{mp}) \right) = \beta_{k0}^* \left( \beta_{m0} z_{m0} + \sum_{p=1}^L \beta_{mp} z_{mp} \right).
\]

So:

\[
|\mathbf{w}_k^H \mathbf{h}_m|^2 = |\beta_{k0}|^2 \left| \beta_{m0} z_{m0} + \sum_{p=1}^L \beta_{mp} z_{mp} \right|^2.
\]

### Step 1: Expand the Squared Term

Let:

\[
u_m = \beta_{m0} z_{m0} + \sum_{p=1}^L \beta_{mp} z_{mp},
\]

so:

\[
|\mathbf{w}_k^H \mathbf{h}_m|^2 = |\beta_{k0}|^2 |u_m|^2.
\]

Compute:

\[
|u_m|^2 = \left( \beta_{m0} z_{m0} + \sum_{p=1}^L \beta_{mp} z_{mp} \right) \left( \beta_{m0}^* z_{m0}^* + \sum_{r=1}^L \beta_{mr}^* z_{mr}^* \right).
\]

Expand:

\[
|u_m|^2 = |\beta_{m0}|^2 |z_{m0}|^2 + \sum_{p=1}^L |\beta_{mp}|^2 |z_{mp}|^2 + \sum_{p=1}^L \beta_{m0} \beta_{mp}^* z_{m0} z_{mp}^* + \sum_{p=1}^L \beta_{mp} \beta_{m0}^* z_{mp} z_{m0}^* + \sum_{p=1}^L \sum_{r=1, r \neq p}^L \beta_{mp} \beta_{mr}^* z_{mp} z_{mr}^*.
\]

### Step 2: Compute the Expectation

\[
\mathbb{E}[|\mathbf{w}_k^H \mathbf{h}_m|^2] = \mathbb{E}[|\beta_{k0}|^2 |u_m|^2] = \mathbb{E}[|\beta_{k0}|^2] \mathbb{E}[|u_m|^2],
\]

since \(\beta_{k0}\) is independent of \(\beta_{m0}, \beta_{mp}, \Delta \theta_{mn}^{(l)}\). Given \(\beta_{k0} \sim \mathcal{CN}(0, 1)\), \(\mathbb{E}[|\beta_{k0}|^2] = 1\), so:

\[
\mathbb{E}[|\mathbf{w}_k^H \mathbf{h}_m|^2] = \mathbb{E}[|u_m|^2].
\]

Now compute:

\[
\mathbb{E}[|u_m|^2] = \mathbb{E}[|\beta_{m0}|^2 |z_{m0}|^2] + \sum_{p=1}^L \mathbb{E}[|\beta_{mp}|^2 |z_{mp}|^2] + \sum_{p=1}^L \mathbb{E}[\beta_{m0} \beta_{mp}^* z_{m0} z_{mp}^*] + \sum_{p=1}^L \mathbb{E}[\beta_{mp} \beta_{m0}^* z_{mp} z_{m0}^*] + \sum_{p=1}^L \sum_{r=1, r \neq p}^L \mathbb{E}[\beta_{mp} \beta_{mr}^* z_{mp} z_{mr}^*].
\]

Since \(\beta_{m0}, \beta_{mp}\) are independent and zero-mean (\(\mathbb{E}[\beta_{m0}] = 0\), \(\mathbb{E}[\beta_{mp}] = 0\)), the cross terms vanish:

\[
\mathbb{E}[\beta_{m0} \beta_{mp}^* z_{m0} z_{mp}^*] = \mathbb{E}[\beta_{m0} \beta_{mp}^*] \mathbb{E}[z_{m0} z_{mp}^*] = 0,
\]

\[
\mathbb{E}[\beta_{mp} \beta_{m0}^* z_{mp} z_{m0}^*] = \mathbb{E}[\beta_{mp} \beta_{m0}^*] \mathbb{E}[z_{mp} z_{m0}^*] = 0,
\]

\[
\mathbb{E}[\beta_{mp} \beta_{mr}^* z_{mp} z_{mr}^*] = \mathbb{E}[\beta_{mp} \beta_{mr}^*] \mathbb{E}[z_{mp} z_{mr}^*] = 0, \quad (p \neq r),
\]

because \(\mathbb{E}[\beta_{mp} \beta_{mr}^*] = 0\) for all \( p, r \) (including \( p = 0 \)) when \( p \neq r \). Thus:

\[
\mathbb{E}[|u_m|^2] = \mathbb{E}[|\beta_{m0}|^2] \mathbb{E}[|z_{m0}|^2] + \sum_{p=1}^L \mathbb{E}[|\beta_{mp}|^2] \mathbb{E}[|z_{mp}|^2].
\]

Given:

- \(\mathbb{E}[|\beta_{m0}|^2] = 1\),
- \(\mathbb{E}[|\beta_{mp}|^2] = 0.1\),

we need \(\mathbb{E}[|z_{m0}|^2]\) and \(\mathbb{E}[|z_{mp}|^2]\).

### Step 3: Compute \(\mathbb{E}[|z_{m0}|^2]\)

\[
z_{m0} = \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))} e^{j \Delta \theta_{mn}^{(0)}},
\]

\[
|z_{m0}|^2 = \sum_{n=0}^{N-1} \sum_{q=0}^{N-1} e^{j k (n - q) (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))} e^{j (\Delta \theta_{mn}^{(0)} - \Delta \theta_{mq}^{(0)})}.
\]

\[
\mathbb{E}[|z_{m0}|^2] = \sum_{n=0}^{N-1} \sum_{q=0}^{N-1} e^{j k (n - q) (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))} \mathbb{E}[e^{j (\Delta \theta_{mn}^{(0)} - \Delta \theta_{mq}^{(0)})}].
\]

For \( n = q \):

\[
\mathbb{E}[e^{j (\Delta \theta_{mn}^{(0)} - \Delta \theta_{mn}^{(0)})}] = 1,
\]

contributing \( N \) terms. For \( n \neq q \), \(\Delta \theta_{mn}^{(0)} - \Delta \theta_{mq}^{(0)} \sim \mathcal{N}(0, 2\sigma^2)\), so:

\[
\mathbb{E}[e^{j (\Delta \theta_{mn}^{(0)} - \Delta \theta_{mq}^{(0)})}] = e^{-\sigma^2}.
\]

Let \(\Delta \cos_{m0} = \cos(\varphi_{k0}) - \cos(\varphi_{m0})\):

\[
\mathbb{E}[|z_{m0}|^2] = N + e^{-\sigma^2} \sum_{n \neq q} e^{j k (n - q) \Delta \cos_{m0}}.
\]

\[
\sum_{n \neq q} e^{j k (n - q) \Delta \cos_{m0}} = \sum_{n=0}^{N-1} \sum_{q=0}^{N-1} e^{j k (n - q) \Delta \cos_{m0}} - N = \left| \sum_{n=0}^{N-1} e^{j k n \Delta \cos_{m0}} \right|^2 - N.
\]

\[
\sum_{n=0}^{N-1} e^{j k n \Delta \cos_{m0}} = \frac{1 - e^{j k N \Delta \cos_{m0}}}{1 - e^{j k \Delta \cos_{m0}}},
\]

\[
\left| \sum_{n=0}^{N-1} e^{j k n \Delta \cos_{m0}} \right|^2 = \frac{\sin^2 \left( \frac{k N \Delta \cos_{m0}}{2} \right)}{\sin^2 \left( \frac{k \Delta \cos_{m0}}{2} \right)}.
\]

So:

\[
\mathbb{E}[|z_{m0}|^2] = N + e^{-\sigma^2} \left( \frac{\sin^2 \left( \frac{k N \Delta \cos_{m0}}{2} \right)}{\sin^2 \left( \frac{k \Delta \cos_{m0}}{2} \right)} - N \right).
\]

Since \(\mathbb{E}[|\beta_{m0}|^2] = 1\):

\[
\mathbb{E}[|\beta_{m0}|^2 |z_{m0}|^2] = N + e^{-\sigma^2} \left( \frac{\sin^2 \left( \frac{k N (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))}{2} \right)}{\sin^2 \left( \frac{k (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))}{2} \right)} - N \right).
\]

### Step 4: Compute \(\mathbb{E}[|z_{mp}|^2]\)

Similarly, for \( p = 1, \ldots, L \):

\[
z_{mp} = \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_{k0}) - \cos(\varphi_{mp}))} e^{j \Delta \theta_{mn}^{(p)}},
\]

\[
\mathbb{E}[|z_{mp}|^2] = N + e^{-\sigma^2} \sum_{n \neq q} e^{j k (n - q) (\cos(\varphi_{k0}) - \cos(\varphi_{mp}))}.
\]

Let \(\Delta \cos_{mp} = \cos(\varphi_{k0}) - \cos(\varphi_{mp})\):

\[
\mathbb{E}[|z_{mp}|^2] = N + e^{-\sigma^2} \left( \frac{\sin^2 \left( \frac{k N \Delta \cos_{mp}}{2} \right)}{\sin^2 \left( \frac{k \Delta \cos_{mp}}{2} \right)} - N \right).
\]

Since \(\mathbb{E}[|\beta_{mp}|^2] = 0.1\):

\[
\mathbb{E}[|\beta_{mp}|^2 |z_{mp}|^2] = 0.1 \left[ N + e^{-\sigma^2} \left( \frac{\sin^2 \left( \frac{k N (\cos(\varphi_{k0}) - \cos(\varphi_{mp}))}{2} \right)}{\sin^2 \left( \frac{k (\cos(\varphi_{k0}) - \cos(\varphi_{mp}))}{2} \right)} - N \right) \right].
\]

### Step 5: Total Interference Expectation

\[
\mathbb{E}[|u_m|^2] = \mathbb{E}[|\beta_{m0}|^2 |z_{m0}|^2] + \sum_{p=1}^L \mathbb{E}[|\beta_{mp}|^2 |z_{mp}|^2]
= \left[ N + e^{-\sigma^2} \left( \frac{\sin^2 \left( \frac{k N (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))}{2} \right)}{\sin^2 \left( \frac{k (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))}{2} \right)} - N \right) \right]+ \sum_{p=1}^L 0.1 \left[ N + e^{-\sigma^2} \left( \frac{\sin^2 \left( \frac{k N (\cos(\varphi_{k0}) - \cos(\varphi_{mp}))}{2} \right)}{\sin^2 \left( \frac{k (\cos(\varphi_{k0}) - \cos(\varphi_{mp}))}{2} \right)} - N \right) \right].
\]

Since \(\mathbb{E}[|\mathbf{w}_k^H \mathbf{h}_m|^2] = \mathbb{E}[|u_m|^2]\), the total interference energy expectation is:

\[
\mathbb{E}\left[ \sum_{m \neq k} |\mathbf{w}_k^H \mathbf{h}_m|^2 \right] = \sum_{m \neq k} \mathbb{E}[|u_m|^2].
\]

### Final Answer

\[
\mathbb{E}[|\mathbf{w}_k^H \mathbf{h}_m|^2] = \left[ N + e^{-\sigma^2} \left( \frac{\sin^2 \left( \frac{\pi N (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))}{2} \right)}{\sin^2 \left( \frac{\pi (\cos(\varphi_{k0}) - \cos(\varphi_{m0}))}{2} \right)} - N \right) \right] + \sum_{p=1}^L 0.1 \left[ N + e^{-\sigma^2} \left( \frac{\sin^2 \left( \frac{\pi N (\cos(\varphi_{k0}) - \cos(\varphi_{mp}))}{2} \right)}{\sin^2 \left( \frac{\pi (\cos(\varphi_{k0}) - \cos(\varphi_{mp}))}{2} \right)} - N \right) \right] 
\]

令\({\delta _{mp}} = \pi (\cos ({\varphi _{k0}}) - \cos ({\varphi _{mp}}))\)，则

\[
\mathbb{E}\left[ \sum_{m \neq k} |\mathbf{w}_k^H \mathbf{h}_m|^2 \right] = \sum_{m \neq k} \left\{ \left[ N + e^{-\sigma^2} \left( \left( \frac{\sin(N \delta_{m0} / 2)}{\sin(\delta_{m0} / 2)} \right)^2 - N \right) \right] + \sum_{p=1}^L 0.1 \left[ N + e^{-\sigma^2} \left( \left| \frac{\sin(N \delta_{ml} / 2)}{\sin(\delta_{ml} / 2)} \right|^2 - N \right) \right] \right\},
\]

where \( k = \pi \) (since \( d = \lambda/2 \), \( k = 2\pi/\lambda \), so \( k d = \pi \)). The expectation accounts for the LOS and NLOS components of the interfering users, with the phase errors reducing the interference power by \( e^{-\sigma^2} \), and the geometric series terms reflecting the angular separation between \(\varphi_{k0}\) and \(\varphi_{m0}, \varphi_{mp}\).