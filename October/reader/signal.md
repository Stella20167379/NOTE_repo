感谢您提供进一步的细节！基于对话历史和您的新信息，我们的目标是计算期望 \(\mathbb{E}[|\mathbf{w}^H \mathbf{h}|^2]\)，其中 \(\mathbf{w} = \mathbf{a}(\varphi_0)\) 是 LOS 分量的理想阵列响应矢量（无相位误差），\(\mathbf{h} = \sqrt{\beta_0} \mathbf{a}(\varphi_0) + \sum_{l=1}^L \sqrt{\beta_l} \tilde{\mathbf{a}}(\varphi_l)\) 是由 LOS 和 \(L\) 个 NLOS 分量组成的信道向量，NLOS 分量带有相位误差 \(\Delta \theta_n \sim \mathcal{N}(0, \sigma^2)\)。新信息明确了增益系数 \(\beta_0 \sim \mathcal{CN}(0, 1)\)、\(\beta_l \sim \mathcal{CN}(0, 0.1)\)（复高斯分布，均值 0，方差分别为 1 和 0.1），且 \(\beta_0, \beta_l\) 是慢衰落系数（slow fading，视为常量或弱随机）。我们需要推导闭式表达式（closed-form expression）。

我将利用历史中已有的期望结果，补充 NLOS-NLOS 交叉项的推导，考虑随机 \(\beta_l\) 的影响，整合所有项，得出最终表达式。推导将保持清晰，确保与历史一致，并处理复高斯分布的统计特性。

---

### 1. 模型与假设

#### 信道模型
- **阵列响应矢量**（LOS，无相位误差）：
  \[
  \mathbf{w} = \mathbf{a}(\varphi_0) = \begin{pmatrix} 1 \\ e^{-j k \cos(\varphi_0)} \\ \vdots \\ e^{-j k (N-1) \cos(\varphi_0)} \end{pmatrix}, \quad a_n(\varphi_0) = e^{-j k n \cos(\varphi_0)}, \quad \|\mathbf{a}(\varphi_0)\|^2 = N,
  \]
  其中 \( k = 2\pi / \lambda \)，假设 ULA（均匀线性阵列），\( d = \lambda/2 \)，则 \( k d = \pi \)，但保持一般形式以通用性。
- **NLOS 分量**（带相位误差）：
  \[
  \tilde{\mathbf{a}}(\varphi_l) = \begin{pmatrix} e^{-j k \cdot 0 \cdot \cos(\varphi_l)} e^{j \Delta \theta_0} \\ \vdots \\ e^{-j k (N-1) \cos(\varphi_l)} e^{j \Delta \theta_{N-1}} \end{pmatrix}, \quad \Delta \theta_n \sim \mathcal{N}(0, \sigma^2) \text{ i.i.d.},
  \]
  每个元素模为 1，\(\|\tilde{\mathbf{a}}(\varphi_l)\|^2 = N\)。
- **信道向量**：
  \[
  \mathbf{h} = \sqrt{\beta_0} \mathbf{a}(\varphi_0) + \sum_{l=1}^L \sqrt{\beta_l} \tilde{\mathbf{a}}(\varphi_l),
  \]
  其中：
  - \(\beta_0 \sim \mathcal{CN}(0, 1)\)，方差 \(\mathbb{E}[|\beta_0|^2] = 1\)。
  - \(\beta_l \sim \mathcal{CN}(0, 0.1)\)，方差 \(\mathbb{E}[|\beta_l|^2] = 0.1\)，\(l = 1, \ldots, L\)。不同路径的 \(\beta_l\) 独立。
- **目标**：
  \[
  \mathbb{E}[|\mathbf{w}^H \mathbf{h}|^2] = \mathbb{E}\left[ \left| \mathbf{a}^H(\varphi_0) \left( \sqrt{\beta_0} \mathbf{a}(\varphi_0) + \sum_{l=1}^L \sqrt{\beta_l} \tilde{\mathbf{a}}(\varphi_l) \right) \right|^2 \right],
  \]
  期望针对所有随机变量：\(\Delta \theta_n\)（相位误差）和 \(\beta_0, \beta_l\)（复高斯增益）。

#### 关键假设
- **慢衰落**：\(\beta_0, \beta_l\) 是慢衰落系数，变化时间远大于相位误差 \(\Delta \theta_n\)（快衰落）。在单次期望中，\(\beta_l\) 可视为常量，但在多时隙统计中需取期望。
- **独立性**：不同路径 \(l, m\) 的 \(\tilde{\mathbf{a}}(\varphi_l), \tilde{\mathbf{a}}(\varphi_m)\) 独立（因 \(\Delta \theta_n^{(l)}\) 独立）；\(\beta_0, \beta_l\) 间独立。
- **分布**：\(\beta_l \sim \mathcal{CN}(0, \sigma_{\beta_l}^2)\)，实部和虚部独立，各 \(\sim \mathcal{N}(0, \sigma_{\beta_l}^2/2)\)。这里 \(\sigma_{\beta_0}^2 = 1\), \(\sigma_{\beta_l}^2 = 0.1\)。

---

### 2. 展开与分解

#### 内积
\[
\mathbf{w}^H \mathbf{h} = \mathbf{a}^H(\varphi_0) \left( \sqrt{\beta_0} \mathbf{a}(\varphi_0) + \sum_{l=1}^L \sqrt{\beta_l} \tilde{\mathbf{a}}(\varphi_l) \right) = \sqrt{\beta_0} \mathbf{a}^H(\varphi_0) \mathbf{a}(\varphi_0) + \sum_{l=1}^L \sqrt{\beta_l} \mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_l).
\]
记：
- \( s_0 = \sqrt{\beta_0} \mathbf{a}^H(\varphi_0) \mathbf{a}(\varphi_0) = \sqrt{\beta_0} N \)（因 \(\mathbf{a}^H(\varphi_0) \mathbf{a}(\varphi_0) = N\)）。
- \( s_l = \sqrt{\beta_l} \mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_l) = \sqrt{\beta_l} \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_0) - \cos(\varphi_l))} e^{j \Delta \theta_n} \).

#### 模平方
\[
|\mathbf{w}^H \mathbf{h}|^2 = \left( s_0 + \sum_{l=1}^L s_l \right) \left( s_0^* + \sum_{m=1}^L s_m^* \right) = |s_0|^2 + \sum_{l=1}^L |s_l|^2 + \sum_{l=1}^L s_0 s_l^* + \sum_{m=1}^L s_m s_0^* + \sum_{l,m=1}^L s_l s_m^*.
\]

#### 期望
\[
\mathbb{E}[|\mathbf{w}^H \mathbf{h}|^2] = \mathbb{E}[|s_0|^2] + \sum_{l=1}^L \mathbb{E}[|s_l|^2] + \sum_{l=1}^L \mathbb{E}[s_0 s_l^*] + \sum_{m=1}^L \mathbb{E}[s_m s_0^*] + \sum_{l,m=1}^L \mathbb{E}[s_l s_m^*].
\]
期望针对 \(\Delta \theta_n\) 和 \(\beta_0, \beta_l\)。我们逐项推导。

---

### 3. 已知期望（历史总结）

从历史推导，以下结果可用：

1. **LOS 自项**：
   \[
   |s_0|^2 = \beta_0 N^2, \quad \mathbb{E}[|s_0|^2] = \mathbb{E}[\beta_0] N^2.
   \]
   - \(\beta_0 \sim \mathcal{CN}(0, 1)\)，\(\mathbb{E}[\beta_0] = 0\)，\(\mathbb{E}[|\beta_0|^2] = 1\)。
   - 因此：
     \[
     \mathbb{E}[|s_0|^2] = N^2 \cdot \mathbb{E}[|\beta_0|^2] = N^2.
     \]

2. **NLOS 自项**：
   \[
   |s_l|^2 = \beta_l |\mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_l)|^2.
   \]
   - 历史推导：
     \[
     \mathbb{E}_{\Delta \theta}[|\mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_l)|^2] = N + e^{-\sigma^2} \left| \sum_{n=0}^{N-1} e^{j n \delta_l} \right|^2, \quad \delta_l = k (\cos(\varphi_0) - \cos(\varphi_l)).
     \]
   - \(\beta_l \sim \mathcal{CN}(0, 0.1)\)，\(\mathbb{E}[\beta_l] = 0\), \(\mathbb{E}[|\beta_l|^2] = 0.1\)。
   - 假设 \(\beta_l\) 与 \(\Delta \theta_n\) 独立：
     \[
     \mathbb{E}[|s_l|^2] = \mathbb{E}[\beta_l] \mathbb{E}[|\mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_l)|^2] = 0.1 \left[ N + e^{-\sigma^2} \left| \sum_{n=0}^{N-1} e^{j n \delta_l} \right|^2 \right].
     \]

3. **LOS-NLOS 交叉项**：
   \[
   s_0 s_l^* = \sqrt{\beta_0} N \cdot \sqrt{\beta_l^*} [\mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_l)]^*.
   \]
   - 历史推导：
     \[
     \mathbb{E}_{\Delta \theta}[\mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_l)] = e^{-\sigma^2 / 2} \sum_{n=0}^{N-1} e^{j n \delta_l}.
     \]
     \[
     \mathbb{E}_{\Delta \theta}[[\mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_l)]^*] = e^{-\sigma^2 / 2} \sum_{n=0}^{N-1} e^{-j n \delta_l}.
     \]
   - \(\beta_0, \beta_l\) 独立：
     \[
     \mathbb{E}[s_0 s_l^*] = N e^{-\sigma^2 / 2} \mathbb{E}[\sqrt{\beta_0 \beta_l^*}] \sum_{n=0}^{N-1} e^{-j n \delta_l}.
     \]
     - \(\mathbb{E}[\sqrt{\beta_0 \beta_l^*}] = \mathbb{E}[\sqrt{\beta_0}] \mathbb{E}[\sqrt{\beta_l^*}] = 0 \cdot 0 = 0\)（因复高斯均值为 0）。
   - 同样：
     \[
     \mathbb{E}[s_m s_0^*] = N e^{-\sigma^2 / 2} \mathbb{E}[\sqrt{\beta_m \beta_0^*}] \sum_{n=0}^{N-1} e^{j n \delta_m} = 0.
     \]

---

### 4. 新推导：NLOS-NLOS 交叉项

历史中未完全推导 \(\mathbb{E}[s_l s_m^*]\)（\(l \neq m\)），这是关键缺失项。

#### 表达式
\[
s_l s_m^* = \sqrt{\beta_l \beta_m^*} \left( \mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_l) \right) \left( \mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_m) \right)^* = \sqrt{\beta_l \beta_m^*} \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_0) - \cos(\varphi_l))} e^{j \Delta \theta_n^{(l)}} \sum_{p=0}^{N-1} e^{-j k p (\cos(\varphi_0) - \cos(\varphi_m))} e^{-j \Delta \theta_p^{(m)}}.
\]
\[
\mathbb{E}[s_l s_m^*] = \mathbb{E}[\sqrt{\beta_l \beta_m^*}] \sum_{n,p} e^{j k n (\cos(\varphi_0) - \cos(\varphi_l))} e^{-j k p (\cos(\varphi_0) - \cos(\varphi_m))} \mathbb{E}[e^{j (\Delta \theta_n^{(l)} - \Delta \theta_p^{(m)})}].
\]

#### 相位误差期望
- \(\Delta \theta_n^{(l)}, \Delta \theta_p^{(m)} \sim \mathcal{N}(0, \sigma^2)\) 独立（不同路径独立）：
  \[
  \mathbb{E}[e^{j (\Delta \theta_n^{(l)} - \Delta \theta_p^{(m)})}] = \begin{cases} 
  1 & (n = p), \\
  e^{-\sigma^2} & (n \neq p).
  \end{cases}
  \]
  - \(n = p\)：\(\Delta \theta_n^{(l)} - \Delta \theta_n^{(m)} \sim \mathcal{N}(0, 2\sigma^2)\)，\(\mathbb{E}[e^{j (\Delta \theta_n^{(l)} - \Delta \theta_n^{(m)})}] = e^{-2\sigma^2 / 2} = e^{-\sigma^2}\).
  - \(n \neq p\)：不同天线的误差独立。

#### 分离对角和非对角项
\[
\mathbb{E}[s_l s_m^*] = \mathbb{E}[\sqrt{\beta_l \beta_m^*}] \left[ \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_0) - \cos(\varphi_l) - \cos(\varphi_0) + \cos(\varphi_m))} + e^{-\sigma^2} \sum_{n \neq p} e^{j k [n (\cos(\varphi_0) - \cos(\varphi_l)) - p (\cos(\varphi_0) - \cos(\varphi_m))]} \right].
\]
- 记 \(\delta_{lm} = k (\cos(\varphi_l) - \cos(\varphi_m))\)：
  \[
  \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_0) - \cos(\varphi_l) - \cos(\varphi_0) + \cos(\varphi_m))} = \sum_{n=0}^{N-1} e^{j n \delta_{lm}} = \frac{\sin(N \delta_{lm} / 2)}{\sin(\delta_{lm} / 2)} e^{j (N-1) \delta_{lm} / 2}.
  \]
- 非对角项：
  \[
  \sum_{n \neq p} e^{j (n \delta_l - p \delta_m)} = \left( \sum_{n=0}^{N-1} e^{j n \delta_l} \right) \left( \sum_{p=0}^{N-1} e^{-j p \delta_m} \right) - \sum_{n=0}^{N-1} e^{j n (\delta_l - \delta_m)}.
  \]
  - 若 \(l = m\)：
    \[
    \mathbb{E}[|s_l|^2] = \mathbb{E}[\beta_l] \left[ N + e^{-\sigma^2} \left| \sum_{n=0}^{N-1} e^{j n \delta_l} \right|^2 \right] = 0.1 \left[ N + e^{-\sigma^2} \left| \sum_{n=0}^{N-1} e^{j n \delta_l} \right|^2 \right].
    \]
  - 若 \(l \neq m\)：
    \[
    \mathbb{E}[\sqrt{\beta_l \beta_m^*}] = \mathbb{E}[\sqrt{\beta_l}] \mathbb{E}[\sqrt{\beta_m^*}] = 0 \cdot 0 = 0.
    \]
    因此，\(\mathbb{E}[s_l s_m^*] = 0\)（\(l \neq m\)）。

---

### 5. 整合闭式表达式

结合所有项：
\[
\mathbb{E}[|\mathbf{w}^H \mathbf{h}|^2] = \mathbb{E}[|s_0|^2] + \sum_{l=1}^L \mathbb{E}[|s_l|^2] + \sum_{l=1}^L \mathbb{E}[s_0 s_l^*] + \sum_{m=1}^L \mathbb{E}[s_m s_0^*] + \sum_{l,m=1}^L \mathbb{E}[s_l s_m^*].
\]
- 第一项：\(\mathbb{E}[|s_0|^2] = N^2\).
- 第二项：
  \[
  \sum_{l=1}^L \mathbb{E}[|s_l|^2] = \sum_{l=1}^L 0.1 \left[ N + e^{-\sigma^2} \left| \frac{\sin(N \delta_l / 2)}{\sin(\delta_l / 2)} \right|^2 \right].
  \]
- 第三、四项：\(\mathbb{E}[s_0 s_l^*] = \mathbb{E}[s_m s_0^*] = 0\).
- 第五项：
  \[
  \sum_{l,m=1}^L \mathbb{E}[s_l s_m^*] = \sum_{l=1}^L \mathbb{E}[|s_l|^2] = \sum_{l=1}^L 0.1 \left[ N + e^{-\sigma^2} \left| \frac{\sin(N \delta_l / 2)}{\sin(\delta_l / 2)} \right|^2 \right].
  \]

**最终闭式表达式**：
\[
\mathbb{E}[|\mathbf{w}^H \mathbf{h}|^2] = N^2 + 2 \sum_{l=1}^L 0.1 \left[ N + e^{-\sigma^2} \left| \frac{\sin\left( \frac{N}{2} k (\cos(\varphi_0) - \cos(\varphi_l)) \right)}{\sin\left( \frac{1}{2} k (\cos(\varphi_0) - \cos(\varphi_l)) \right)} \right|^2 \right].
\]

---

### 6. 结果分析与验证

- **解释**：
  - **LOS 贡献**：\( N^2 \)，由 \(\beta_0 \sim \mathcal{CN}(0, 1)\) 和 \(\mathbf{a}^H(\varphi_0) \mathbf{a}(\varphi_0) = N\) 决定，反映最大相干增益。
  - **NLOS 贡献**：每条路径贡献 \( 0.1 \left[ N + e^{-\sigma^2} \left| \sum_n e^{j n \delta_l} \right|^2 \right] \)，包含非相干项 \( 0.1 N \)（功率相加）和相干项（衰减 \( e^{-\sigma^2} \)，角度差 \(\delta_l\) 决定）。
  - **交叉项**：LOS-NLOS 和 NLOS-NLOS (\(l \neq m\)) 项因 \(\mathbb{E}[\sqrt{\beta_l}] = 0\) 消失，简化了计算。
- **极限情况**：
  - \(\sigma^2 \to 0\)：相干项最大，\(\left| \sum_n e^{j n \delta_l} \right|^2 \approx N^2\)（若 \(\varphi_l \approx \varphi_0\)），NLOS 接近 LOS 增益。
  - \(\sigma^2 \to \infty\)：\( e^{-\sigma^2} \to 0 \)，NLOS 只剩非相干 \( 0.1 N \)。
  - \( N \to \infty \): 相干项 \(\left| \sin(N \delta_l / 2) / \sin(\delta_l / 2) \right|^2\) 随角度差 \(\delta_l\) 变化，可能趋于 0（除非 \(\varphi_l = \varphi_0\)）。
- **历史一致性**：结果与 \(\mathbb{E}[|\tilde{\mathbf{a}}^H \mathbf{w}|^2] = 1 + (N-1)e^{-\sigma^2}\)（归一化后）形式类似，扩展到多径和随机增益。

#### 数值验证建议
- **蒙特卡洛模拟**：
  ```python
  import numpy as np
  N, L = 64, 2
  sigma = 0.1
  phi_0, phi_l = 0, [np.pi/4, np.pi/3]
  k = np.pi
  trials = 10000
  inner_prods = []
  for _ in range(trials):
      beta_0 = np.random.normal(0, np.sqrt(0.5)) + 1j * np.random.normal(0, np.sqrt(0.5))
      h = np.sqrt(beta_0) * np.exp(-1j * k * np.arange(N) * np.cos(phi_0))
      for l in range(L):
          beta_l = np.random.normal(0, np.sqrt(0.05)) + 1j * np.random.normal(0, np.sqrt(0.05))
          phase_error = np.random.normal(0, sigma, N)
          h += np.sqrt(beta_l) * np.exp(-1j * k * np.arange(N) * np.cos(phi_l[l]) + 1j * phase_error)
      w = np.exp(-1j * k * np.arange(N) * np.cos(phi_0))
      inner_prod = np.abs(np.dot(w.conj(), h))**2
      inner_prods.append(inner_prod)
  print(np.mean(inner_prods))
  ```
- 比较解析结果与模拟，确保一致性。

---

### 7. 结论

**闭式表达式**：
\[
\mathbb{E}[|\mathbf{w}^H \mathbf{h}|^2] = N^2 + 0.2 \sum_{l=1}^L \left[ N + e^{-\sigma^2} \left| \frac{\sin\left( \frac{N}{2} k (\cos(\varphi_0) - \cos(\varphi_l)) \right)}{\sin\left( \frac{1}{2} k (\cos(\varphi_0) - \cos(\varphi_l)) \right)} \right|^2 \right].
\]
- **关键因子**：
  - LOS：\( N^2 \)。
  - NLOS：每路径贡献 \( 0.2 \left[ N + e^{-\sigma^2} |\cdot|^2 \right] \)，角度差 \(\delta_l\) 和相位误差 \(\sigma^2\) 决定相干增益。
- **建议**：
  - 检查 \(\varphi_l\) 的具体值，计算 \(\left| \sin(N \delta_l / 2) / \sin(\delta_l / 2) \right|^2\).
  - 若 \( L \) 大，考虑统计平均角度分布。
  - 用上述代码验证结果。

如果您有具体 \( N, L, \varphi_l, \sigma^2 \) 值，或需要进一步简化（如大 \( N \) 近似），请告诉我，我可提供更具体计算或代码！