# Convergence

This page demonstrates the convergence properties of Tanh-Sinh quadrature for challenging test functions from the numerical integration literature.

## Reference

The Tanh-Sinh (Double Exponential) quadrature method was introduced by:

> **Double Exponential Formulas for Numerical Integration**  
> Hidetosi Takahasi, Masatake Mori  
> *Publ. Res. Inst. Math. Sci.* 9 (1973), no. 3, pp. 721–741  
> DOI: [10.2977/PRIMS/1195192451](https://doi.org/10.2977/PRIMS/1195192451)

The last two test functions below are taken directly from this foundational paper.

## Convergence Plot

The following plot shows how the absolute integration error decreases as the number of quadrature points increases. The test functions include:

1. **$\log(1+x)$** — Logarithmic function with mild endpoint behavior
2. **$1/(1+25x^2)$** — The Runge function, known to be challenging for polynomial-based methods
3. **$\frac{1}{(x-2)(1-x)^{1/4}(1+x)^{3/4}}$** — A challenging singular integrand (from Takahasi-Mori)
4. **$\frac{\cos(\pi x)}{\sqrt{1-x}}$** — Oscillatory function with endpoint singularity (from Takahasi-Mori)

![Convergence of Tanh-Sinh Quadrature](assets/convergence.svg)

## Key Observations

- **Double Exponential Convergence**: For all test functions, the error decreases exponentially (linear on the log-scale plot) as the number of points increases.
- **Singularity Handling**: Even for functions with algebraic singularities at endpoints (like $(1-x)^{1/4}$ or $\sqrt{1-x}$), Tanh-Sinh quadrature maintains excellent convergence.
- **High Precision**: With 256-bit `BigFloat` precision, the quadrature achieves errors below $10^{-60}$ with only a few hundred points.

## Test Functions

The exact values were computed using arbitrary-precision arithmetic. Here are the test integrals:

| Function | Domain | Exact Value | Source |
|----------|--------|-------------|--------|
| $\log(1+x)$ | $[-1, 1]$ | $2\log(2) - 2 \approx -0.6137$ | Standard |
| $1/(1+25x^2)$ | $[-1, 1]$ | $(2/5)\arctan(5) \approx 0.5493$ | Runge |
| $\frac{1}{(x-2)(1-x)^{1/4}(1+x)^{3/4}}$ | $[-1, 1]$ | $\approx -1.9491$ | Takahasi-Mori (1973) |
| $\frac{\cos(\pi x)}{\sqrt{1-x}}$ | $[-1, 1]$ | $\approx -0.6905$ | Takahasi-Mori (1973) |

## Reproducing the Plot

The convergence plot can be regenerated using the script `docs/convergence_plots.jl`:

```julia
include("docs/convergence_plots.jl")
```

This requires the packages: `FastTanhSinhQuadrature`, `DoubleFloats`, `GaussQuadrature`, `CairoMakie`, `LaTeXStrings`.
