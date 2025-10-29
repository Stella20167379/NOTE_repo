
### 1. 模型与假设

#### 信道模型
- **信道向量**：
  \[
  \mathbf{h} = \sqrt{\beta_0} \tilde{\mathbf{a}}(\varphi_0) + \sum_{l=1}^L \sqrt{\beta_l} \tilde{\mathbf{a}}(\varphi_l),
  \]
  其中：
  - \(\tilde{\mathbf{a}}(\varphi_l) = [e^{-j \pi n \cos(\varphi_l)} e^{j \Delta \theta_{n}^{(l)}}]_{n=0}^{N-1}\)，维度为 \( N \times 1 \)，是带相位误差的阵列响应矢量，\(\Delta \theta_{n}^{(l)} \sim \mathcal{N}(0, \sigma^2)\) i.i.d.（跨天线 \( n \) 和路径 \( l \) 独立）。
  - \(\beta_0 \sim \mathcal{CN}(0, 1)\)，\(\beta_l \sim \mathcal{CN}(0, 0.1)\)（\( l = 1, \dots, L \)，历史中提到的慢衰落增益）。
  - 假设 ULA（均匀线性阵列），天线间距 \( d = \lambda/2 \)，波数 \( k d = \pi \)。
  - \(\varphi_0\) 是 LOS 路径角度，\(\varphi_l\) 是 NLOS 路径角度。
- **矩阵定义**：
  - \(\mathbf{h} \mathbf{h}^H\) 是 \( N \times N \) 矩阵，元素为 \([\mathbf{h} \mathbf{h}^H]_{n,p} = h_n h_p^*\)。
  - \( p_u \mathbf{h} \mathbf{h}^H + \mathbf{I} \) 是 \( N \times N \) 矩阵，\(\mathbf{I}\) 是单位矩阵。
- **目标**：推导矩阵 \( p_u \mathbf{h} \mathbf{h}^H + \mathbf{I} \) 的元素 \([p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,p}\)。
- **假设**：
  - 相位误差 \(\Delta \theta_{n}^{(l)}\) 独立，\(\mathbb{E}[e^{j \Delta \theta_{n}^{(l)}}] = e^{-\sigma^2 / 2}\)。
  - \(\beta_0, \beta_l\) 独立，\(\mathbb{E}[|\beta_0|^2] = 1\), \(\mathbb{E}[|\beta_l|^2] = 0.1\)。
  - 不同路径 \( l \) 的 \(\tilde{\mathbf{a}}(\varphi_l)\) 独立（因 \(\Delta \theta_{n}^{(l)}\) 独立）。
  
---

### 2. 推导矩阵元素

矩阵 \( p_u \mathbf{h} \mathbf{h}^H + \mathbf{I} \) 的元素为：
\[
[p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,p} = p_u h_n h_p^* + \delta_{n,p},
\]
其中 \(\delta_{n,p} = 1\)（若 \( n = p \)），否则为 0（单位矩阵的定义）。我们需要计算 \( h_n h_p^* \)，然后乘以 \( p_u \) 并加上单位矩阵项。

#### 信道向量元素
\[
\mathbf{h} = \sqrt{\beta_0} \tilde{\mathbf{a}}(\varphi_0) + \sum_{l=1}^L \sqrt{\beta_l} \tilde{\mathbf{a}}(\varphi_l).
\]
第 \( n \) 个元素：
\[
h_n = \sqrt{\beta_0} \tilde{a}_n(\varphi_0) + \sum_{l=1}^L \sqrt{\beta_l} \tilde{a}_n(\varphi_l),
\]
其中：
\[
\tilde{a}_n(\varphi_l) = e^{-j \pi n \cos(\varphi_l)} e^{j \Delta \theta_{n}^{(l)}}.
\]
共轭：
\[
h_p^* = \sqrt{\beta_0^*} \tilde{a}_p^*(\varphi_0) + \sum_{l=1}^L \sqrt{\beta_l^*} \tilde{a}_p^*(\varphi_l),
\]
\[
\tilde{a}_p^*(\varphi_l) = e^{j \pi p \cos(\varphi_l)} e^{-j \Delta \theta_{p}^{(l)}}.
\]

#### 计算 \( h_n h_p^* \)
\[
h_n h_p^* = \left( \sqrt{\beta_0} e^{-j \pi n \cos(\varphi_0)} e^{j \Delta \theta_{n}^{(0)}} + \sum_{l=1}^L \sqrt{\beta_l} e^{-j \pi n \cos(\varphi_l)} e^{j \Delta \theta_{n}^{(l)}} \right) \left( \sqrt{\beta_0^*} e^{j \pi p \cos(\varphi_0)} e^{-j \Delta \theta_{p}^{(0)}} + \sum_{m=1}^L \sqrt{\beta_m^*} e^{j \pi p \cos(\varphi_m)} e^{-j \Delta \theta_{p}^{(m)}} \right).
\]
展开：
\[
h_n h_p^* = \sum_{l=0}^L \sum_{m=0}^L \sqrt{\beta_l \beta_m^*} e^{-j \pi n \cos(\varphi_l)} e^{j \pi p \cos(\varphi_m)} e^{j (\Delta \theta_{n}^{(l)} - \Delta \theta_{p}^{(m)})}.
\]
- 包含 \( (L+1) \times (L+1) \) 项，分为：
  - 对角项（\( l = m \)）：如 \( \beta_0 |\tilde{a}_n(\varphi_0)| |\tilde{a}_p(\varphi_0)| \)。
  - 非对角项（\( l \neq m \)）：如 \( \sqrt{\beta_l \beta_m^*} \tilde{a}_n(\varphi_l) \tilde{a}_p^*(\varphi_m) \)。

#### 矩阵元素
\[
[p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,p} = p_u \sum_{l=0}^L \sum_{m=0}^L \sqrt{\beta_l \beta_m^*} e^{-j \pi n \cos(\varphi_l)} e^{j \pi p \cos(\varphi_m)} e^{j (\Delta \theta_{n}^{(l)} - \Delta \theta_{p}^{(m)})} + \delta_{n,p}.
\]
- **对角元素** (\( n = p \))：
  \[
  h_n h_n^* = |h_n|^2 = \left| \sqrt{\beta_0} e^{-j \pi n \cos(\varphi_0)} e^{j \Delta \theta_{n}^{(0)}} + \sum_{l=1}^L \sqrt{\beta_l} e^{-j \pi n \cos(\varphi_l)} e^{j \Delta \theta_{n}^{(l)}} \right|^2.
  \]
  展开：
  \[
  |h_n|^2 = \sum_{l=0}^L \sum_{m=0}^L \sqrt{\beta_l \beta_m^*} .
  \]
  因此：
  \[
  [p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,n} = p_u \sum_{l=0}^L \sum_{m=0}^L \sqrt{\beta_l \beta_m^*}  + 1.
  \]
  - 当 \( l = m \), \( e^{j (\Delta \theta_{n}^{(l)} - \Delta \theta_{n}^{(l)})} = 1 \), \( e^{-j \pi n (\cos(\varphi_l) - \cos(\varphi_l))} = 1 \)，得 \( |\beta_l| |\tilde{a}_n(\varphi_l)|^2 = |\beta_l| \).
  - 非对角项（\( l \neq m \)) 包含相位误差和角度差。

- **非对角元素** (\( n \neq p \))：
  \[
  [p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,p} = p_u \sum_{l=0}^L \sum_{m=0}^L \sqrt{\beta_l \beta_m^*} e^{-j \pi n \cos(\varphi_l)} e^{j \pi p \cos(\varphi_m)} e^{j (\Delta \theta_{n}^{(l)} - \Delta \theta_{p}^{(m)})}.
  \]

---

### 3. 期望（可选分析）

如果需要矩阵元素的期望（历史中常计算期望，如 SINR 分析），我们可以计算：
\[
\mathbb{E}[p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,p} = p_u \mathbb{E}[h_n h_p^*] + \delta_{n,p}.
\]
- **计算 \(\mathbb{E}[h_n h_p^*]\)**：
  \[
  \mathbb{E}[h_n h_p^*] = \sum_{l=0}^L \sum_{m=0}^L \mathbb{E}[\sqrt{\beta_l \beta_m^*}] \mathbb{E}[e^{j (\Delta \theta_{n}^{(l)} - \Delta \theta_{p}^{(m)})}] e^{-j \pi n \cos(\varphi_l)} e^{j \pi p \cos(\varphi_m)}.
  \]
  - \(\mathbb{E}[\sqrt{\beta_l \beta_m^*}] = 0\)（若 \( l \neq m \)，因 \(\beta_l\) 独立，零均值）。
  - 仅 \( l = m \) 非零：
    \[
    \mathbb{E}[|\beta_0|^2] = 1, \quad \mathbb{E}[|\beta_l|^2] = 0.1 \quad (l=1,\dots,L).
    \]
  - 相位误差：
    - 若 \( n = p, l = m \), \(\mathbb{E}[e^{j (\Delta \theta_{n}^{(l)} - \Delta \theta_{n}^{(l)})}] = 1\).
    - 若 \( n \neq p, l = m \), \(\mathbb{E}[e^{j (\Delta \theta_{n}^{(l)} - \Delta \theta_{p}^{(l)})}] = e^{-\sigma^2}\).
    - 若 \( l \neq m \), \(\mathbb{E}[e^{j (\Delta \theta_{n}^{(l)} - \Delta \theta_{p}^{(m)})}] = e^{-\sigma^2 / 2} e^{-\sigma^2 / 2} = e^{-\sigma^2}\).

  因此：
  \[
  \mathbb{E}[h_n h_p^*] = \sum_{l=0}^L \mathbb{E}[|\beta_l|^2] e^{-j \pi (n-p) \cos(\varphi_l)} \mathbb{E}[e^{j (\Delta \theta_{n}^{(l)} - \Delta \theta_{p}^{(l)})}].
  \]
  - **对角** (\( n = p \))：
    \[
    \mathbb{E}[h_n h_n^*] = \mathbb{E}[|\beta_0|^2] \cdot 1 + \sum_{l=1}^L \mathbb{E}[|\beta_l|^2] \cdot 1 = 1 + 0.1L.
    \]
    \[
    \mathbb{E}[p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,n} = p_u (1 + 0.1L) + 1.
    \]
  - **非对角** (\( n \neq p \))：
    \[
    \mathbb{E}[h_n h_p^*] = \sum_{l=0}^L \mathbb{E}[|\beta_l|^2] e^{-j \pi (n-p) \cos(\varphi_l)} e^{-\sigma^2} = e^{-\sigma^2} \left[ e^{-j \pi (n-p) \cos(\varphi_0)} + 0.1 \sum_{l=1}^L e^{-j \pi (n-p) \cos(\varphi_l)} \right].
    \]
    \[
    \mathbb{E}[p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,p} = p_u e^{-\sigma^2} \left[ e^{-j \pi (n-p) \cos(\varphi_0)} + 0.1 \sum_{l=1}^L e^{-j \pi (n-p) \cos(\varphi_l)} \right].
    \]

---

### 6. 结论

- **矩阵元素**：
  - 对角 (\( n = p \))：
    \[
    [p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,n} = p_u \sum_{l=0}^L  \sqrt{\beta_l \beta_l^*} + 1.
    \]
  - 非对角 (\( n \neq p \))：
    \[
    [p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,p} = p_u \sum_{l=0}^L \sum_{m=0}^L \sqrt{\beta_l \beta_m^*} e^{-j \pi n \cos(\varphi_l)} e^{j \pi p \cos(\varphi_m)} e^{j (\Delta \theta_{n}^{(l)} - \Delta \theta_{p}^{(m)})}.
    \]
- **期望**（若需要）：
  - 对角：
    \[
    \mathbb{E}[p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,n} = p_u (1 + 0.1L) + 1.
    \]
  - 非对角：
    \[
    \mathbb{E}[p_u \mathbf{h} \mathbf{h}^H + \mathbf{I}]_{n,p} = p_u e^{-\sigma^2} \left[ e^{-j \pi (n-p) \cos(\varphi_0)} + 0.1 \sum_{l=1}^L e^{-j \pi (n-p) \cos(\varphi_l)} \right].
    \]
- **验证**：
  - 可通过蒙特卡洛模拟生成 \(\beta_l, \Delta \theta_{n}^{(l)}\)，计算 \( h_n h_p^* \)，验证解析表达式。
  - 历史中的 Python 代码可稍作修改以计算矩阵元素。

如果需要进一步简化（如大 \( N \) 近似）、数值验证代码，或具体参数（如 \( p_u, \varphi_l \)）的计算，请提供更多细节！

---

### 7. 量化噪声期望

量化噪声功率因子如下所示：
\[{\bf{w}}_k^H diag\left( {{{\bf{R}}_{{{\bf{n}}_{\bf{q}}}{{\bf{n}}_{\bf{q}}}}}} \right){{\bf{w}}_k}\]

展开为元素累计计算形式，则为：
\[\sum\limits_{n = 1}^N {{{\left[ {{{\bf{R}}_{{{\bf{n}}_{\bf{q}}}{{\bf{n}}_{\bf{q}}}}}} \right]}_{nn}}{{\left| {{w_{k,n}}} \right|}^2}} \]

又因为
  \[
  \mathbf{w} = \beta_0 \mathbf{a}(\varphi_0) =  \beta_0 \begin{pmatrix} 1 \\ e^{-j k \cos(\varphi_0)} \\ \vdots \\ e^{-j k (N-1) \cos(\varphi_0)} \end{pmatrix}, \quad a_n(\varphi_0) = e^{-j k n \cos(\varphi_0)}, \quad \|\mathbf{a}(\varphi_0)\|^2 = N,
  \]

所以
\[\begin{array}{l}
\mathbb{E} \left[ {{\bf{w}}_k^Hdiag\left( {{{\bf{R}}_{{{\bf{n}}_{\bf{q}}}{{\bf{n}}_{\bf{q}}}}}} \right){{\bf{w}}_k}} \right]\\
 = \sum\limits_{n = 0}^{N - 1} { |\beta_0 |^2 {{\left| {{e^{ - jkn\cos ({\varphi _0})}}} \right|}^2} \cdot \mathbb{E} {{[{p_u}{\bf{h}}{{\bf{h}}^H} + {\bf{I}}]}_{n,n}}} \\
 = {N\left[ {{p_u}(1 + 0.1L) + 1} \right]}.
\end{array}\]
