Let’s dive into deriving \( |s_0|^2 \), \( |s_l|^2 \), and their expectations for the given MU-MISO communication scenario. The goal is to compute these expressions and their expected values, considering the phase errors and channel characteristics.

---

### Step 1: Compute \( |s_0|^2 \)

The expression for \( s_0 \) is given as:

\[
s_0 = \sqrt{\beta_0} \mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_0) = \sqrt{\beta_0} \sum_{n=0}^{N-1} e^{j \Delta \theta_n},
\]

where:
- \( \beta_0 \sim \mathcal{CN}(0, 1) \), so \( \mathbb{E}[|\beta_0|^2] = 1 \).
- \( \Delta \theta_n \sim \mathcal{N}(0, \sigma^2) \), i.i.d. across \( n \).
- \( \mathbf{a}(\varphi_0) = \begin{pmatrix} 1 \\ e^{-j k \cos(\varphi_0)} \\ \vdots \\ e^{-j k (N-1) \cos(\varphi_0)} \end{pmatrix} \).
- \( \tilde{\mathbf{a}}(\varphi_0) = \begin{pmatrix} e^{-j k \cdot 0 \cdot \cos(\varphi_0)} e^{j \Delta \theta_0} \\ \vdots \\ e^{-j k (N-1) \cos(\varphi_0)} e^{j \Delta \theta_{N-1}} \end{pmatrix} \).

The Hermitian inner product \( \mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_0) \) is:

\[
\mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_0) = \sum_{n=0}^{N-1} \left( e^{-j k n \cos(\varphi_0)} \right)^* e^{-j k n \cos(\varphi_0)} e^{j \Delta \theta_n} = \sum_{n=0}^{N-1} e^{j k n \cos(\varphi_0)} e^{-j k n \cos(\varphi_0)} e^{j \Delta \theta_n} = \sum_{n=0}^{N-1} e^{j \Delta \theta_n}.
\]

Thus:

\[
s_0 = \sqrt{\beta_0} \sum_{n=0}^{N-1} e^{j \Delta \theta_n}.
\]

Now, compute \( |s_0|^2 \):

\[
|s_0|^2 = \left| \sqrt{\beta_0} \sum_{n=0}^{N-1} e^{j \Delta \theta_n} \right|^2 = |\beta_0| \left| \sum_{n=0}^{N-1} e^{j \Delta \theta_n} \right|^2.
\]

Let’s denote:

\[
z_0 = \sum_{n=0}^{N-1} e^{j \Delta \theta_n},
\]

so:

\[
|s_0|^2 = |\beta_0| |z_0|^2.
\]

To compute \( |z_0|^2 \):

\[
|z_0|^2 = \left( \sum_{n=0}^{N-1} e^{j \Delta \theta_n} \right) \left( \sum_{m=0}^{N-1} e^{-j \Delta \theta_m} \right) = \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} e^{j (\Delta \theta_n - \Delta \theta_m)}.
\]

Thus:

\[
|s_0|^2 = |\beta_0| \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} e^{j (\Delta \theta_n - \Delta \theta_m)}.
\]

---

### Step 2: Expectation of \( |s_0|^2 \)

To find \( \mathbb{E}[|s_0|^2] \), we need:

\[
\mathbb{E}[|s_0|^2] = \mathbb{E} \left[ |\beta_0| \left| \sum_{n=0}^{N-1} e^{j \Delta \theta_n} \right|^2 \right].
\]

Since \( \beta_0 \) and \( \Delta \theta_n \) are independent, and noting that \( |\beta_0| = |\beta_0|^2 \) (since \( |\beta_0| \) is the magnitude squared in the context of power), we have:

\[
\mathbb{E}[|s_0|^2] = \mathbb{E}[|\beta_0|] \mathbb{E} \left[ \left| \sum_{n=0}^{N-1} e^{j \Delta \theta_n} \right|^2 \right] = \mathbb{E}[|\beta_0|^2] \mathbb{E} \left[ \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} e^{j (\Delta \theta_n - \Delta \theta_m)} \right].
\]

Given \( \beta_0 \sim \mathcal{CN}(0, 1) \), we have \( \mathbb{E}[|\beta_0|^2] = 1 \). Now compute:

\[
\mathbb{E} \left[ |z_0|^2 \right] = \mathbb{E} \left[ \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} e^{j (\Delta \theta_n - \Delta \theta_m)} \right].
\]

Separate the terms where \( n = m \) and \( n \neq m \):

\[
|z_0|^2 = \sum_{n=0}^{N-1} e^{j \Delta \theta_n} e^{-j \Delta \theta_n} + \sum_{n \neq m} e^{j (\Delta \theta_n - \Delta \theta_m)} = N + \sum_{n \neq m} e^{j (\Delta \theta_n - \Delta \theta_m)}.
\]

Taking the expectation:

\[
\mathbb{E} \left[ |z_0|^2 \right] = N + \sum_{n \neq m} \mathbb{E} \left[ e^{j (\Delta \theta_n - \Delta \theta_m)} \right].
\]

Since \( \Delta \theta_n \sim \mathcal{N}(0, \sigma^2) \), the difference \( \Delta \theta_n - \Delta \theta_m \sim \mathcal{N}(0, 2\sigma^2) \) for \( n \neq m \), because the variances add for independent Gaussian variables. The characteristic function of a Gaussian random variable \( X \sim \mathcal{N}(0, \sigma_x^2) \) is:

\[
\mathbb{E} [e^{j t X}] = e^{- \frac{t^2 \sigma_x^2}{2}}.
\]

For \( \Delta \theta_n - \Delta \theta_m \):

\[
\mathbb{E} \left[ e^{j (\Delta \theta_n - \Delta \theta_m)} \right] = e^{- \frac{1^2 \cdot 2\sigma^2}{2}} = e^{-\sigma^2}.
\]

There are \( N(N-1) \) terms where \( n \neq m \), so:

\[
\mathbb{E} \left[ |z_0|^2 \right] = N + \sum_{n \neq m} e^{-\sigma^2} = N + N(N-1) e^{-\sigma^2}.
\]

Thus:

\[
\mathbb{E}[|s_0|^2] = 1 \cdot \left[ N + N(N-1) e^{-\sigma^2} \right] = N + N(N-1) e^{-\sigma^2}.
\]

---

### Step 3: Compute \( |s_l|^2 \)

For the non-LOS paths (\( l = 1, \ldots, L \)):

\[
s_l = \sqrt{\beta_l} \mathbf{a}^H(\varphi_0) \tilde{\mathbf{a}}(\varphi_l) = \sqrt{\beta_l} \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_0) - \cos(\varphi_l))} e^{j \Delta \theta_n},
\]

where \( \beta_l \sim \mathcal{CN}(0, 0.1) \), so \( \mathbb{E}[|\beta_l|^2] = 0.1 \). Compute:

\[
|s_l|^2 = |\beta_l| \left| \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_0) - \cos(\varphi_l))} e^{j \Delta \theta_n} \right|^2.
\]

Let:

\[
z_l = \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_0) - \cos(\varphi_l))} e^{j \Delta \theta_n},
\]

so:

\[
|s_l|^2 = |\beta_l| |z_l|^2.
\]

Compute \( |z_l|^2 \):

\[
|z_l|^2 = \left( \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_0) - \cos(\varphi_l))} e^{j \Delta \theta_n} \right) \left( \sum_{m=0}^{N-1} e^{-j k m (\cos(\varphi_0) - \cos(\varphi_l))} e^{-j \Delta \theta_m} \right).
\]

\[
= \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} e^{j k (n - m) (\cos(\varphi_0) - \cos(\varphi_l))} e^{j (\Delta \theta_n - \Delta \theta_m)}.
\]

Thus:

\[
|s_l|^2 = |\beta_l| \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} e^{j k (n - m) (\cos(\varphi_0) - \cos(\varphi_l))} e^{j (\Delta \theta_n - \Delta \theta_m)}.
\]

---

### Step 4: Expectation of \( |s_l|^2 \)

\[
\mathbb{E}[|s_l|^2] = \mathbb{E}[|\beta_l|] \mathbb{E} \left[ \left| \sum_{n=0}^{N-1} e^{j k n (\cos(\varphi_0) - \cos(\varphi_l))} e^{j \Delta \theta_n} \right|^2 \right] = \mathbb{E}[|\beta_l|^2] \mathbb{E} [|z_l|^2].
\]

Since \( \mathbb{E}[|\beta_l|^2] = 0.1 \), compute:

\[
\mathbb{E} [|z_l|^2] = \mathbb{E} \left[ \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} e^{j k (n - m) (\cos(\varphi_0) - \cos(\varphi_l))} e^{j (\Delta \theta_n - \Delta \theta_m)} \right].
\]

Again, separate terms:

\[
|z_l|^2 = \sum_{n=0}^{N-1} 1 + \sum_{n \neq m} e^{j k (n - m) (\cos(\varphi_0) - \cos(\varphi_l))} e^{j (\Delta \theta_n - \Delta \theta_m)}.
\]

\[
\mathbb{E} [|z_l|^2] = N + \sum_{n \neq m} e^{j k (n - m) (\cos(\varphi_0) - \cos(\varphi_l))} \mathbb{E} \left[ e^{j (\Delta \theta_n - \Delta \theta_m)} \right].
\]

Using \( \mathbb{E} \left[ e^{j (\Delta \theta_n - \Delta \theta_m)} \right] = e^{-\sigma^2} \), we get:

\[
\mathbb{E} [|z_l|^2] = N + e^{-\sigma^2} \sum_{n \neq m} e^{j k (n - m) (\cos(\varphi_0) - \cos(\varphi_l))}.
\]

The sum is:

\[
\sum_{n \neq m} e^{j k (n - m) \Delta \cos} = \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} e^{j k (n - m) \Delta \cos} - \sum_{n=m} 1,
\]

where \( \Delta \cos = \cos(\varphi_0) - \cos(\varphi_l) \). The double sum is:

\[
\sum_{n=0}^{N-1} \sum_{m=0}^{N-1} e^{j k (n - m) \Delta \cos} = \left( \sum_{n=0}^{N-1} e^{j k n \Delta \cos} \right) \left( \sum_{m=0}^{N-1} e^{-j k m \Delta \cos} \right) = \left| \sum_{n=0}^{N-1} e^{j k n \Delta \cos} \right|^2.
\]

For a geometric series:

\[
\sum_{n=0}^{N-1} e^{j k n \Delta \cos} = \frac{1 - e^{j k N \Delta \cos}}{1 - e^{j k \Delta \cos}},
\]

\[
\left| \sum_{n=0}^{N-1} e^{j k n \Delta \cos} \right|^2 = \frac{\sin^2 \left( \frac{k N \Delta \cos}{2} \right)}{\sin^2 \left( \frac{k \Delta \cos}{2} \right)}.
\]

So:

\[
\mathbb{E} [|z_l|^2] = N + e^{-\sigma^2} \left( \frac{\sin^2 \left( \frac{k N \Delta \cos}{2} \right)}{\sin^2 \left( \frac{k \Delta \cos}{2} \right)} - N \right).
\]

\[
\mathbb{E}[|s_l|^2] = 0.1 \left[ N + e^{-\sigma^2} \left( \frac{\sin^2 \left( \frac{k N (\cos(\varphi_0) - \cos(\varphi_l))}{2} \right)}{\sin^2 \left( \frac{k (\cos(\varphi_0) - \cos(\varphi_l))}{2} \right)} - N \right) \right].
\]

---

### Final Answer

\[
|s_0|^2 = |\beta_0| \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} e^{j (\Delta \theta_n - \Delta \theta_m)},
\]

\[
\mathbb{E}[|s_0|^2] = N + N(N-1) e^{-\sigma^2},
\]

\[
|s_l|^2 = |\beta_l| \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} e^{j k (n - m) (\cos(\varphi_0) - \cos(\varphi_l))} e^{j (\Delta \theta_n - \Delta \theta_m)},
\]

\[
\mathbb{E}[|s_l|^2] = 0.1 \left[ N + e^{-\sigma^2} \left( \frac{\sin^2 \left( \frac{k N (\cos(\varphi_0) - \cos(\varphi_l))}{2} \right)}{\sin^2 \left( \frac{k (\cos(\varphi_0) - \cos(\varphi_l))}{2} \right)} - N \right) \right].
\]
