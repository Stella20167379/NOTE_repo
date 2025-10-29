求
  \[
  {{\left\| {{\bf{w}}_k^H} \right\|}^2}
  \]
的期望

而
  \[
  \mathbf{w} = \beta_0 \mathbf{a}(\varphi_0) =  \beta_0 \begin{pmatrix} 1 \\ e^{-j k \cos(\varphi_0)} \\ \vdots \\ e^{-j k (N-1) \cos(\varphi_0)} \end{pmatrix}, \quad a_n(\varphi_0) = e^{-j k n \cos(\varphi_0)}, \quad \|\mathbf{a}(\varphi_0)\|^2 = N,
  \]
  其中 \(\beta_0 \sim \mathcal{CN}(0, 1)\)

易得

  \[{\left\| {{\bf{w}}_k^H} \right\|^2} = {\bf{w}}_k^H{{\bf{w}}_k} = N{\left| {{\beta _0}} \right|^2}\]

所以
  \[{\mathbb{E}}[{\left\| {{\bf{w}}_k^H} \right\|^2}] = N\]