# Error Models and Node Estimates

This page summarizes a practical asymptotic error model for tanh-sinh quadrature and turns it into point-count estimates that are useful when choosing `N` for fixed-grid routines such as `integrate1D` and `integrate1D_avx`.

The derivation follows the same decomposition used by Vanherck et al. (2020), with notation aligned to this package.

## Setup and Assumptions

Let

```math
I=\int_a^b f(x)\,dx,
\qquad
Q_n^{(h)}=\sum_{i=-n}^{n} h\,g(ih),
```

be the exact integral and an $N=2n+1$ point tanh-sinh approximation with spacing $h$ in the transformed variable $t$.

Vanherck et al. use the decomposition

```math
|I-Q_n^{(h)}| \le |\Delta I_h| + |\varepsilon_t|,
```

where:

- $\Delta I_h$ is the discretization error of the infinite trapezoidal rule,
- $\varepsilon_t$ is the truncation error from restricting the infinite sum to $|i|\le n$.

For integrands analytic in a strip $|\Im z|<d$, a standard asymptotic model is

```math
E(n,h)\equiv |I-Q_n^{(h)}|
\approx
A e^{-2\pi d/h}
+
B e^{-(\pi/2)e^{nh}}.
```

The transformed integrand decay term $e^{-(\pi/2)e^{|t|}}$ is the double-exponential mechanism behind tanh-sinh.

## Optimal Spacing from Error Balancing

Balancing the two exponents gives

```math
\frac{2\pi d}{h}\approx \frac{\pi}{2}e^{nh}.
```

With $N=2n+1$, this yields the Lambert-$W$ expression

```math
h_{\mathrm{opt}}(N)=\frac{2}{N}W(2dN),
```

and for large $N$,

```math
h_{\mathrm{opt}}(N)\approx \frac{2}{N}\ln(2dN).
```

This is consistent with the standard optimal-spacing asymptotic law:

```math
|I-Q_n^{(h_{\mathrm{opt}})}|
=
\mathcal{O}\!\left(
\exp\!\left(-\frac{\pi d N}{\ln(2dN)}\right)
\right).
```

## Maximal Spacing Asymptotic (Derived)

In finite precision, the practical window constraint is often written as

```math
nh=t_{\max},
\qquad
h_{\max}(n)=\frac{t_{\max}}{n}\approx\frac{2t_{\max}}{N}.
```

Substitute this into the two-term model:

```math
E_{\max}(N)
\approx
A\exp\!\left(-\frac{\pi d}{t_{\max}}N\right)
+
B\exp\!\left(-\frac{\pi}{2}e^{t_{\max}}\right).
```

Equivalently,

```math
E_{\max}(N)
=
\mathcal{O}\!\left(
e^{-\pi d N/t_{\max}}
+
e^{-(\pi/2)e^{t_{\max}}}
\right).
```

Interpretation:

- A first regime with exponential decay in $N$ at rate $\pi d/t_{\max}$.
- A second regime where error plateaus at a fixed floor set by $t_{\max}$ and machine precision.

This explains why maximal spacing is often more robust in practice, even if optimal spacing is asymptotically faster.

## Point-Count Estimates for Target Accuracy

Let $\varepsilon$ denote the requested target error level and define

```math
L=\ln(1/\varepsilon).
```

### Estimate with Optimal Spacing

A useful inversion of the optimal-spacing law is

```math
N_{\mathrm{opt}}
\approx
\frac{L}{\pi d}\ln\!\left(\frac{2L}{\pi}\right),
\qquad
n_{\mathrm{opt}}\approx \frac{N_{\mathrm{opt}}-1}{2}.
```

### Estimate with Maximal Spacing

Ignoring unknown prefactors $A,B$, the leading estimate is

```math
N_{\max}
\approx
\frac{t_{\max}}{\pi d}L.
```

The truncation floor is approximately

```math
\varepsilon_{\mathrm{floor}}
\approx
\exp\!\left(-\frac{\pi}{2}e^{t_{\max}}\right).
```

If $\varepsilon \le \varepsilon_{\mathrm{floor}}$, increasing $N$ alone will not improve accuracy; increase $t_{\max}$ (if possible) and/or precision.

## Practical Workflow in This Package

1. Prefer `quad` when possible: it adapts automatically and uses the stopping test

```math
|I_{\mathrm{new}}-I_{\mathrm{old}}|
\le
\max(\mathrm{atol},\mathrm{rtol}|I_{\mathrm{new}}|).
```

2. For fixed-grid kernels (`integrate1D`, `integrate1D_avx`, ...), use the formulas above to choose an initial `N`.
3. Validate with a doubling test (`N` then `2N`) and stop when successive estimates satisfy the same tolerance criterion.
4. For repeated calls on similar integrands, calibrate `N` once and reuse precomputed nodes.

## Limitations

The formulas above are asymptotic and hide problem-dependent constants. They are best used as engineering estimates, not strict guarantees. Unknown strip width $d$, non-ideal analyticity, and floating-point effects can shift the required `N`.
