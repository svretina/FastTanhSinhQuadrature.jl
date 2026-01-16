```@meta
CurrentModule = FastTanhSinhQuadrature
```

# FastTanhSinhQuadrature.jl

<p align="center">
  <img src="assets/logo.svg" alt="FastTanhSinhQuadrature Logo" width="250">
</p>

![Stable](https://img.shields.io/badge/docs-stable-blue.svg)
![Dev](https://img.shields.io/badge/docs-dev-blue.svg)
![Build Status](https://github.com/svretina/FastTanhSinhQuadrature.jl/actions/workflows/CI.yml/badge.svg?branch=master)
![Coverage](https://codecov.io/gh/svretina/FastTanhSinhQuadrature.jl/branch/master/graph/badge.svg)

**FastTanhSinhQuadrature.jl** is a high-performance Julia library for numerical integration using the [Tanh-Sinh (Double Exponential) quadrature](https://en.wikipedia.org/wiki/Tanh-sinh_quadrature) method.

It handles **singularities at endpoints** robustly, supports **arbitrary precision** arithmetic (e.g., `BigFloat`, `Double64`), and leverages **SIMD** for speed.

## Usage at a Glance

```julia
using FastTanhSinhQuadrature

# 1. Define function
f(x) = x * exp(x)

# 2. Integrate on [-1, 1]
val = integrate(f, 10) # 10 levels
println(val)
```

## Contents

- [Theory](theory.md): Understand the mathematics behind the method.
- [Basic Examples](examples/basics.md): Learn how to integrate simple 1D functions.
- [Advanced Examples](examples/advanced.md): Multidimensional integration and performance tips.
- [Benchmarks](benchmarks.md): Performance comparison against other libraries.
- [API Reference](api.md): Detailed function documentation.

```@index
```


