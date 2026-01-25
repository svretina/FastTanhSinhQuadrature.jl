```@meta
CurrentModule = FastTanhSinhQuadrature
```

# FastTanhSinhQuadrature.jl


![Stable](https://img.shields.io/badge/docs-stable-blue.svg)
![Dev](https://img.shields.io/badge/docs-dev-blue.svg)
![Build Status](https://github.com/svretina/FastTanhSinhQuadrature.jl/actions/workflows/CI.yml/badge.svg?branch=master)
![Coverage](https://codecov.io/gh/svretina/FastTanhSinhQuadrature.jl/branch/master/graph/badge.svg)

**FastTanhSinhQuadrature.jl** is a high-performance Julia library for numerical integration using the [Tanh-Sinh (Double Exponential) quadrature](https://en.wikipedia.org/wiki/Tanh-sinh_quadrature) method.

It handles **singularities at endpoints** robustly, supports **arbitrary precision** arithmetic (e.g., `BigFloat`, `Double64`), and leverages **SIMD** for speed.

The implementation follows the theoretical framework described in [Bailey et al. (2005)](https://ems.press/content/serial-article-files/41766?nt=1).

![Convergence of Tanh-Sinh Quadrature](assets/convergence.svg)

## Quick Start

```julia
using FastTanhSinhQuadrature

# High-level adaptive integration
val = quad(exp, 0.0, 1.0)  # ≈ e - 1

# Handle singularities
f(x) = 1 / sqrt(abs(x))
val = quad_split(f, 0.0, -1.0, 1.0)  # Split at singularity

# Pre-computed nodes for maximum performance
x, w, h = tanhsinh(Float64, 80)
val = integrate1D(sin, x, w, h)

# SIMD-accelerated (2-3x faster)
val = integrate1D_avx(sin, x, w, h)
```

## Key Features

- **Arbitrary Precision**: `Float32`, `Float64`, `BigFloat`, `Double64`
- **High Performance**: SIMD-accelerated `_avx` variants
- **Multidimensional**: 1D, 2D, and 3D integration
- **Adaptive Integration**: `quad` and `quad_split` functions
- **Singularity Handling**: Robust at endpoints and interior points
- **Double Exponential Convergence**: Machine precision with few points

## Contents

- [Theory](theory.md): Understand the mathematics behind the method.
- [Convergence](convergence.md): Convergence analysis and test functions.
- [Basic Examples](examples/basics.md): Learn how to integrate simple 1D functions.
- [Advanced Examples](examples/advanced.md): Multidimensional integration and performance tips.
- [Benchmarks](benchmarks.md): Performance comparison against other libraries.
- [API Reference](api.md): Detailed function documentation.

```@index
```
