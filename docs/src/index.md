```@meta
CurrentModule = FastTanhSinhQuadrature
```

# FastTanhSinhQuadrature.jl


![Stable](https://img.shields.io/badge/docs-stable-blue.svg)
![Dev](https://img.shields.io/badge/docs-dev-blue.svg)
![Build Status](https://github.com/svretina/FastTanhSinhQuadrature.jl/actions/workflows/CI.yml/badge.svg?branch=master)
![Coverage](https://codecov.io/gh/svretina/FastTanhSinhQuadrature.jl/branch/master/graph/badge.svg)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

**FastTanhSinhQuadrature.jl** is a high-performance Julia library for numerical integration using the [Tanh-Sinh (Double Exponential) quadrature](https://en.wikipedia.org/wiki/Tanh-sinh_quadrature) method.

It handles **singularities at endpoints** robustly, supports **arbitrary precision** arithmetic (e.g., `BigFloat`, `Double64`), and leverages **SIMD** for speed.

The implementation follows the method introduced by [Takahasi & Mori (1973)](https://doi.org/10.2977/PRIMS/1195192451).

![Convergence of Tanh-Sinh Quadrature](assets/convergence.svg)

## Quick Start

```julia
using FastTanhSinhQuadrature

# High-level adaptive integration
val = quad(exp, 0.0, 1.0)  # ≈ e - 1

# Match the quadrature type to your bounds
val32 = quad(exp, 0.0f0, 1.0f0)  # Float32

# Handle singularities
f(x) = 1 / sqrt(abs(x))
val = quad_split(f, 0.0, -1.0, 1.0)  # Split at singularity

# Pre-computed nodes
x, w, h = tanhsinh(Float64, 80)
val = integrate1D(x -> sin(x)^2, 0.0, π, x, w, h)  # bounds can include π

# SIMD-accelerated (2-3x faster)
val = integrate1D_avx(sin, x, w, h)

# Use Val{N} for maximum performance with small N (< 128)
x_static, w_static, h_static = tanhsinh(Float64, Val(80))
val_static = integrate1D_avx(sin, x_static, w_static, h_static)
```

## Choosing an API

- Use `quad` when you want adaptive refinement without manually choosing `N`.
- For repeated adaptive calls, build `adaptive_cache_1D/2D/3D` once and pass `cache=...`.
- Use `integrate1D`/`integrate2D`/`integrate3D` with pre-computed `(x, w, h)` when reusing the same nodes across many integrals.
- Use `_avx` variants for `Float32`/`Float64` if your integrand works with `LoopVectorization`.
- Use `quad_split` for interior singularities and `quad_cmpl` for endpoint-sensitive formulations.

Pre-computed and high-level interfaces accept mixed real bounds (`Int`, `Float64`, `π`, etc.) and convert them to the numerical type used by the quadrature nodes. For concrete floating-point inputs such as `Float32`, the quadrature stays in that same type and returns that same type.

## Key Features

- **Arbitrary Precision**: `Float32`, `Float64`, `BigFloat`, `Double64`
- **High Performance**: SIMD-accelerated `_avx` variants
- **Multidimensional**: 1D, 2D, and 3D integration
- **Adaptive Integration**: `quad` and `quad_split` functions
- **Singularity Handling**: Robust at endpoints and interior points
- **Double Exponential Convergence**: Machine precision with few points
- **Type Stability**: Rigorously tested with `JET.jl` for zero runtime dispatch

## Contents

- [Theory](theory.md): Understand the mathematics behind the method.
- [Convergence](convergence.md): Convergence analysis and test functions.
- [Basic Examples](examples/basics.md): Learn how to integrate simple 1D functions.
- [Advanced Examples](examples/advanced.md): Multidimensional integration and performance tips.
- [Benchmarks](benchmarks.md): Performance comparison against other libraries.
- [API Reference](api.md): Detailed function documentation.

## Other Julia quadrature packages

No single quadrature method dominates every integrand class.

Benchmark-based guidance from this repository (`rtol=1e-6`, `atol=1e-8`):

- `FastTanhSinhQuadrature.jl` is strongest for endpoint-dominated cases (for example `1/sqrt(1-x^2)` and `log(1-x)`), where `quad` is much faster than `QuadGK` in our benchmark set.
- For strongly oscillatory 1D integrals (for example `sin^2(1000x)`), `QuadGK.jl` is the better default in our tests.
- For smooth low-order 1D functions, calibrated fixed Gaussian rules from `FastGaussQuadrature.jl` can be the fastest option.
- In 2D/3D box integrals, this package is often very competitive and robust on endpoint-singular products in our suite.
- If values are only available on a predetermined grid (not callable at arbitrary points), use sampled-data integration packages such as `Trapz.jl`, `Romberg.jl`, or `NumericalIntegration.jl`.

Practical rule of thumb:

- Choose this package for endpoint singularities, endpoint-sensitive formulas, and repeated integrations with precomputed nodes/caches.
- Prefer `QuadGK.jl` first for many oscillatory 1D workloads.
- Compare against `HCubature.jl`, `Cubature.jl`, and `Cuba.jl` for general multidimensional cubature problems.
