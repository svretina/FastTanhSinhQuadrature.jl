<img src="docs/FastTanhSinh_Logo.svg" alt="FastTanhSinhQuadrature Logo" width="300">

# FastTanhSinhQuadrature.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://svretina.github.io/FastTanhSinhQuadrature.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://svretina.github.io/FastTanhSinhQuadrature.jl/dev/)
[![Build Status](https://github.com/svretina/FastTanhSinhQuadrature.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/svretina/FastTanhSinhQuadrature.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/svretina/FastTanhSinhQuadrature.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/svretina/FastTanhSinhQuadrature.jl)


Fast and high-precision numerical integration using **Tanh-Sinh (Double Exponential) quadrature** in Julia.

## Overview

`FastTanhSinhQuadrature.jl` is designed for high-performance numerical integration, particularly effective for functions with endpoint singularities. It combines the rigorous accuracy of Tanh-Sinh quadrature with modern Julia performance features like SIMD acceleration.

### Key Features
- **Arbitrary Precision Support**: Seamlessly works with `Float32`, `Float64`, `BigFloat`, and extended precision types like `Double64` (from `DoubleFloats.jl`).
- **High Performance**: specialized `integrate_avx` routines utilize `LoopVectorization.jl` for maximum speed on supported hardware.
- **Multidimensional Support**: Built-in support for **1D**, **2D**, and **3D** integration domains.
- **Efficiency**: Pre-compute quadrature nodes and weights once and reuse them for multiple integrations to save computation time.
- **Singularity Handling**: robust handling of functions with singularities at integration boundaries.

## Installation

Install directly from the Julia REPL:

```julia
import Pkg
Pkg.add("FastTanhSinhQuadrature")
```

## Usage

### 1. Basic 1D Integration

For simple integration of a function $f(x)$ over $[-1, 1]$ (default) or an arbitrary interval $[a, b]$:

```julia
using FastTanhSinhQuadrature

# Define a function to integrate
f(x) = 3x^2

# Quick integration with a specified number of levels (n=10)
val = integrate(f, 10) 
println(val) # Should be approximately 2.0 (integral of 3x^2 from -1 to 1)
```

### 2. High-Precision Integration

You can easily switch to higher precision types like `BigFloat` or `Double64`.

```julia
using FastTanhSinhQuadrature, DoubleFloats

f(x) = exp(x)

# Pre-compute nodes and weights for Double64 precision
# n defines the resolution (higher n = more points)
x, w, h = tanhsinh(Double64, 12) 

# Integrate over [0, 1]
val = integrate(f, 0.0, 1.0, x, w, h)
println(val) # Result in Double64 precision
```

### 3. Adaptive Integration

If you don't know the required resolution `n`, use `adaptive_integrate`. It doubles the number of points until the relative difference between iterations is below `tol`.

```julia
# Integrate exp(x) from 0 to 1 with tolerance 1e-10
val = adaptive_integrate(exp, 0.0, 1.0, tol=1e-10)
println(val)
```

### 4. Pre-computation for Efficiency

If you need to integrate many functions (or the same function over different intervals) using the same precision and resolution, it is highly recommended to pre-compute the quadrature nodes and weights.

```julia
# 1. Generate nodes (x), weights (w), and step size (h) once
x, w, h = tanhsinh(Float64, 40)

# 2. Integrate multiple functions efficiently
f1(x) = sin(x)^2
f2(x) = cos(x)^2

# Use cached x, w, h
res1 = integrate(f1, 0.0, π, x, w, h)
res2 = integrate(f2, 0.0, π, x, w, h)

println("Integrals: $res1, $res2")
```

### 4. Multidimensional Integration (2D & 3D)

The library supports defining integration bounds using `StaticArrays`.

```julia
using StaticArrays

# 2D Integration Example
# Integrate f(x, y) = x*y over a square [-1, 1] x [-1, 1]
f_2d(x, y) = x * y

# Define bounds as SVector
low = SVector(-1.0, -1.0)
up  = SVector(1.0, 1.0)

# Generate quadrature scaling
x, w, h = tanhsinh(Float64, 40)

val_2d = integrate(f_2d, low, up, x, w, h)
println("2D Integral: $val_2d") # Should be 0.0
```

## API Reference

### Core Functions

- **`tanhsinh(T::Type, n::Int)`**
  Generates the quadrature nodes `x`, weights `w`, and step size `h` for a given type `T` and level `n`.
  - Returns: `(x, w, h)`

- **`integrate(f, n::Int)`**
  Convenience function to integrate `f` over $[-1, 1]$ using `Float64` with level `n`.

- **`integrate(T::Type, f, n::Int)`**
  Convenience function to integrate `f` over $[-1, 1]$ using type `T` with level `n`.

- **`integrate(f, x, w, h)`**
  Integrate `f` over $[-1, 1]$ using pre-computed nodes and weights.

- **`integrate(f, a, b, x, w, h)`**
  Integrate `f` over the interval `[a, b]` using pre-computed nodes and weights.

### Advanced / Optimized

- **`integrate_avx(f, ...)`**
  SIMD-accelerated version of `integrate`. Use this for performance-critical standard float operations. Note that this requires `f` to be compatible with `LoopVectorization`.

- **`quad(f, a, b, x, w, h)`**
  A safer wrapper around `integrate` that handles domain checks (e.g., if `a > b` or `a == b`) and automatically calls the appropriate integration routine.

### Singularities

- **`remove_endpoints!(xmin, xmax, x, w)`**
- **`remove_left_endpoint!(xmin, xmax, x, w)`**
- **`remove_right_endpoint!(xmin, xmax, x, w)`**
  Modify the weights `w` in-place to ignore endpoints where a function might be singular or undefined.

---

## Benchmarks

### Benchmark Environment
- **CPU**: Intel(R) Core(TM) Ultra 7 155U

### Results

Comparison against `FastGaussQuadrature.jl` for various functions.
- `TS`: `integrate` (standard Tanh-Sinh)
- `TS SIMD`: `integrate_avx` (SIMD-optimized Tanh-Sinh)
- `GQ`: `FastGaussQuadrature.jl` (Gauss-Legendre)

Timings are in nanoseconds (ns).

| Function | Domain | Points | TS (ns) | TS SIMD (ns) | GQ (ns) | Ratio (TS/GQ) | Ratio (TS SIMD/GQ) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| exp(x) | [-1, 1] | 9 | 54.24 | 18.98 | 46.01 | 1.18 | 0.41 |
| exp(x) | [-1, 1] | 100 | 656.00 | 151.80 | 310.92 | 2.11 | 0.49 |
| exp(x) | [-1, 1] | 1000 | 6941.80 | 1530.50 | 3075.75 | 2.26 | 0.50 |
| sin(x)^2 | [-1, 1] | 9 | 74.93 | 41.13 | 46.04 | 1.63 | 0.89 |
| sin(x)^2 | [-1, 1] | 100 | 795.19 | 287.41 | 353.35 | 2.25 | 0.81 |
| sin(x)^2 | [-1, 1] | 1000 | 8615.33 | 2749.33 | 3958.00 | 2.18 | 0.69 |
| 1/(1+25x^2) | [-1, 1] | 9 | 8.47 | 5.42 | 21.22 | 0.40 | 0.26 |
| 1/(1+25x^2) | [-1, 1] | 100 | 83.82 | 44.54 | 77.23 | 1.09 | 0.58 |
| 1/(1+25x^2) | [-1, 1] | 1000 | 1027.28 | 466.15 | 612.57 | 1.68 | 0.76 |
| sqrt(1-x^2) | [-1, 1] | 9 | 12.04 | 8.02 | 26.30 | 0.46 | 0.31 |
| sqrt(1-x^2) | [-1, 1] | 100 | 127.94 | 68.33 | 144.65 | 0.88 | 0.47 |
| sqrt(1-x^2) | [-1, 1] | 1000 | 1396.40 | 625.68 | 1528.60 | 0.91 | 0.41 |
| x^2 | [-1, 1] | 9 | 5.12 | 3.98 | 20.59 | 0.25 | 0.19 |
| x^2 | [-1, 1] | 100 | 34.88 | 9.74 | 69.44 | 0.50 | 0.14 |
| x^2 | [-1, 1] | 1000 | 442.34 | 117.87 | 578.60 | 0.76 | 0.20 |
| x^3 | [-1, 1] | 9 | 6.26 | 3.96 | 23.39 | 0.27 | 0.17 |
| x^3 | [-1, 1] | 100 | 38.67 | 11.80 | 65.03 | 0.59 | 0.18 |
| x^3 | [-1, 1] | 1000 | 418.34 | 106.36 | 540.83 | 0.77 | 0.20 |
| x^3+x^2+x+1 | [-1, 1] | 9 | 8.66 | 4.62 | 23.41 | 0.37 | 0.20 |
| x^3+x^2+x+1 | [-1, 1] | 100 | 89.29 | 17.90 | 71.05 | 1.26 | 0.25 |
| x^3+x^2+x+1 | [-1, 1] | 1000 | 915.95 | 180.46 | 560.56 | 1.63 | 0.32 |
| log(1-x) | [-1, 1] | 9 | 75.40 | 53.18 | 56.30 | 1.34 | 0.94 |
| log(1-x) | [-1, 1] | 100 | 752.70 | 327.72 | 458.73 | 1.64 | 0.71 |
| log(1-x) | [-1, 1] | 1000 | 7470.00 | 3503.38 | 3730.50 | 2.00 | 0.94 |

*Note: `FastTanhSinhQuadrature`'s SIMD-optimized `integrate_avx` uses `LoopVectorization.jl`. It achieves massive speedups for polynomials (up to 30x faster than GQ) and maintains excellent performance for singular integrals like `log(1-x)` and Runge functions.*

---

*This library was built to be a robust tool for scientific computing in Julia where precision and speed are paramount.*
