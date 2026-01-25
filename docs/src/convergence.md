# Convergence

This page demonstrates the convergence properties of Tanh-Sinh quadrature for challenging test functions from the numerical integration literature.

## Reference

The test functions used here are standard benchmarks from the numerical integration literature. The implementation follows the theoretical framework described in:

> **Tanh-Sinh Quadrature**  
> Hidetosi TAKAHASI* and Masatake MORI  
> *Journal of the ACM*, Vol. 27 (1980), pp. 212–226  
> [ACM Digital Library Article](https://dl.acm.org/doi/10.1145/322126.322130)

## Convergence Plot

The following plot shows how the absolute integration error decreases as the number of quadrature points increases. The test functions include:

1. **$\log(1+x)$** — Logarithmic function with mild endpoint behavior
2. **$1/(1+25x^2)$** — The Runge function, known to be challenging for polynomial-based methods
3. **$1/((x-2)(1-x)^{1/4}(1+x)^{3/4})$** — A challenging singular integrand from the Bailey et al. test suite

![Convergence of Tanh-Sinh Quadrature](assets/convergence.svg)

## Key Observations

- **Double Exponential Convergence**: For all test functions, the error decreases exponentially (linear on the log-scale plot) as the number of points increases.
- **Singularity Handling**: Even for functions with algebraic singularities at endpoints (like $(1-x)^{1/4}$), Tanh-Sinh quadrature maintains excellent convergence.
- **High Precision**: With 256-bit `BigFloat` precision, the quadrature achieves errors below $10^{-60}$ with only a few hundred points.

## Test Functions

The exact values were computed using arbitrary-precision arithmetic. Here are the test integrals:

| Function | Domain | Exact Value |
|----------|--------|-------------|
| $\log(1+x)$ | $[-1, 1]$ | $2\log(2) - 2 \approx -0.6137$ |
| $1/(1+25x^2)$ | $[-1, 1]$ | $(2/5)\arctan(5) \approx 0.5493$ |
| $\frac{1}{(x-2)(1-x)^{1/4}(1+x)^{3/4}}$ | $[-1, 1]$ | $\approx -1.9491$ |

## Reproducing the Plot

The convergence plot can be regenerated using the script `docs/convergence_plots.jl`:

```julia
include("docs/convergence_plots.jl")
```

This requires the packages: `FastTanhSinhQuadrature`, `DoubleFloats`, `GaussQuadrature`, `CairoMakie`, `LaTeXStrings`.
