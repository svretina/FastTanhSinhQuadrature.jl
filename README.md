<p align="center">
  <img src="docs/FastTanhSinh_Logo.svg" alt="FastTanhSinhQuadrature Logo" width="300">
</p>

# FastTanhSinhQuadrature.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://svretina.github.io/FastTanhSinhQuadrature.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://svretina.github.io/FastTanhSinhQuadrature.jl/dev/)
[![Build Status](https://github.com/svretina/FastTanhSinhQuadrature.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/svretina/FastTanhSinhQuadrature.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/svretina/FastTanhSinhQuadrature.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/svretina/FastTanhSinhQuadrature.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)
[![status](https://joss.theoj.org/papers/84c42780ad01b47c0c3b5f1ab5c26260/status.svg)](https://joss.theoj.org/papers/84c42780ad01b47c0c3b5f1ab5c26260)
[![DOI](https://zenodo.org/badge/708778433.svg)](https://doi.org/10.5281/zenodo.19420465)

Fast and high-precision numerical integration using **Tanh-Sinh (Double Exponential) quadrature** in Julia.

## Reproducibility

For exact, command-level reproduction of the JOSS paper artifacts (tests,
benchmarks, figures, and PDF workflow), see
[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).

## Overview

`FastTanhSinhQuadrature.jl` is designed for high-performance numerical integration, particularly effective for functions with endpoint singularities. It combines the rigorous accuracy of Tanh-Sinh quadrature with modern Julia performance features like SIMD acceleration.

### Key Features
- **Arbitrary Precision Support**: Seamlessly works with `Float32`, `Float64`, `BigFloat`, and extended precision types like `Double64` (from `DoubleFloats.jl`).
- **High Performance**: Specialized `integrate1D_avx`, `integrate2D_avx`, `integrate3D_avx` routines utilize `LoopVectorization.jl` for maximum speed.
- **Multidimensional Support**: Built-in support for **1D**, **2D**, and **3D** integration domains.
- **Flexible Bound Types**: Pre-computed and high-level interfaces accept mixed `Real` bounds (for example `Int`, `Float64`, and `π`) and convert them safely to the node type.
- **Memory Efficiency**: Pre-compute quadrature nodes and weights once and reuse them for multiple integrations.
- **Singularity Handling**: Robust handling of functions with singularities at integration boundaries via `quad_split` and the new **Complement Interface** `quad_cmpl`.
- **Double Exponential Convergence**: Achieve machine precision with few points even for singular integrands.
- **Adaptive Integration**: Highly optimized adaptive routines for 1D, 2D, and 3D.
- **Transcendental Caching**: 2D and 3D adaptive routines cache node and weight mappings, reducing transcendental overhead by $O(N^{d-1})$.
- **Optimal and maximal spacing**: Optimal and maximal spacing of nodes and weights.
- **Underflow/Overflow Handling**: Robust handling of underflow/overflow when generating nodes and weights.
- **Type Stable**: Rigorously tested with `JET.jl` to ensure type stability and zero runtime dispatch.

<p align="center">
  <img src="docs/src/assets/convergence.svg" alt="Convergence Plot" width="400">
</p>

## Installation

Install directly from the Julia REPL:

```julia
import Pkg
Pkg.add("FastTanhSinhQuadrature")
```

Supported Julia versions: `1.9` - `1.12`.

## Quick Start

### Choosing an Interface

- Use `quad` for one-off integrations and automatic refinement.
- For repeated adaptive calls, prebuild an adaptive cache (`adaptive_cache_1D/2D/3D`) and pass `cache=...`.
- Use `integrate1D`/`integrate2D`/`integrate3D` with pre-computed `(x, w, h)` when evaluating many integrals on the same geometry.
- Use `_avx` variants for `Float32`/`Float64` when your integrand is compatible with `LoopVectorization`.
- Use `quad_split` for interior singularities and `quad_cmpl` for endpoint-sensitive formulas involving `1-x` / `1+x`.

### High-Level API: `quad`

The simplest way to integrate is using the `quad` function, which provides adaptive integration:

```julia
using FastTanhSinhQuadrature

# Integrate exp(x) from 0 to 1
val = quad(exp, 0.0, 1.0)
println(val)  # ≈ e - 1 ≈ 1.7182818...

# The result type follows the floating-point bound type
val32 = quad(exp, 0.0f0, 1.0f0)
println(typeof(val32))  # Float32

# Integrate over default domain [-1, 1]
val = quad(x -> 3x^2)
println(val)  # ≈ 2.0

# Handle singularities with quad_split
f(x) = 1 / sqrt(abs(x))  # Singular at x=0
val = quad_split(f, 0.0, -1.0, 1.0)  # Split at singularity
println(val)  # ≈ 4.0
```

### Reusing Adaptive Caches

For repeated adaptive integrations, prebuild the cache outside your loop:

```julia
using FastTanhSinhQuadrature

cache = adaptive_cache_1D(Float64; max_levels=16)

for α in (0.5, 1.0, 2.0, 4.0)
    f(x) = exp(-α * x^2)
    val = quad(f, -1.0, 1.0; rtol=1e-10, atol=1e-12, cache=cache)
    println((α, val))
end
```

### High-Accuracy Interface: `quad_cmpl`

Use `quad_cmpl` when your integrand is naturally written in endpoint distances.
You call it exactly like `quad`:

```julia
val = quad_cmpl(f, a, b)
```

The difference is the callback signature. For interval `[a, b]`, `f` must accept:
`f(x, b_minus_x, x_minus_a)`.

At each quadrature node `x`, the package passes:

- `x`: evaluation point in `[a, b]`
- `b_minus_x = b - x`: distance to the right endpoint
- `x_minus_a = x - a`: distance to the left endpoint

For `[-1, 1]`, this is `f(x, 1-x, 1+x)`.

Why this interface exists:
Tanh-Sinh places many nodes extremely close to endpoints. In that regime, computing `b - x`
or `x - a` directly from rounded `x` can lose relative precision (subtractive cancellation).
`quad_cmpl` computes and passes these complements in a numerically stable way, improving robustness for
endpoint-sensitive formulas such as `log(b-x)`, `1/sqrt((b-x)(x-a))`, or `exp(-1/(b-x))`.

```julia
using FastTanhSinhQuadrature

# Integrate f(x) = 1/sqrt(1-x^2) using passed endpoint distances
# 1-x^2 = (1-x)(1+x)
f(x, b_minus_x, x_minus_a) = 1 / sqrt(b_minus_x * x_minus_a)
val = quad_cmpl(f, -1.0, 1.0)
println(val)  # ≈ π ≈ 3.14159...
```

### Pre-computed Nodes for Maximum Performance

For repeated integrations, pre-compute nodes and weights once:

```julia
using FastTanhSinhQuadrature

# Generate nodes (x), weights (w), and step size (h)
x, w, h = tanhsinh(Float64, 80)

# Integrate multiple functions efficiently
f1(x) = sin(x)^2
f2(x) = cos(x)^2

# Bounds can be any Real type (Float64, Int, Irrational such as π).
res1 = integrate1D(f1, 0.0, π, x, w, h)
res2 = integrate1D(f2, 0.0, π, x, w, h)
println("Integrals: $res1, $res2")  # Both ≈ π/2

# Use Val{N} for maximum performance with small N (< 128)
# This returns StaticArrays for nodes and weights, allowing for full SIMD acceleration and zero allocations on the heap.
x_static, w_static, h_static = tanhsinh(Float64, Val(80))
val_static = integrate1D_avx(f1, 0.0, π, x_static, w_static, h_static)
```

### SIMD-Accelerated Integration

For `Float32`/`Float64`, use the `_avx` variants for maximum speed:

```julia
x, w, h = tanhsinh(Float64, 100)

# Standard integration
val1 = integrate1D(exp, x, w, h)

# SIMD-accelerated (2-3x faster)
val2 = integrate1D_avx(exp, x, w, h)
```

### High-Precision Integration

Switch to higher precision types like `BigFloat` or `Double64`:

```julia
using FastTanhSinhQuadrature, DoubleFloats

# Double64 precision (~32 decimal digits)
val = quad(exp, Double64(0), Double64(1); rtol=1e-30)

# BigFloat precision (arbitrary)
setprecision(BigFloat, 256)
x, w, h = tanhsinh(BigFloat, 120)
val = integrate1D(exp, x, w, h)
```

If you want single precision rather than higher precision, pass `Float32` bounds to `quad` or generate `Float32` nodes with `tanhsinh(Float32, N)`. The result stays in `Float32` for those concrete `Float32` entry points.

### Multidimensional Integration (2D & 3D)

Use `StaticArrays` for defining integration bounds:

```julia
using FastTanhSinhQuadrature, StaticArrays

# 2D: Integrate f(x,y) = x*y over [-1,1] × [-1,1]
x, w, h = tanhsinh(Float64, 40)
low = SVector(-1.0, -1.0)
up  = SVector(1.0, 1.0)
val = integrate2D((x, y) -> x * y, low, up, x, w, h)
println(val)  # ≈ 0.0

# 3D: Integrate constant 1 over unit cube
val = quad((x, y, z) -> 1.0, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
println(val)  # ≈ 1.0
```

## API Reference

### High-Level Functions

| Function | Description |
|----------|-------------|
| `quad(f; rtol, atol, max_levels, cache)` | Adaptive 1D integration over `[-1, 1]` |
| `quad(f, low, up; rtol, atol, max_levels, cache)` | Adaptive 1D integration over `[low, up]` for `low, up <: Real` |
| `quad_cmpl(f, low, up; rtol, atol, max_levels, cache)` | High-accuracy 1D integration for `f(x, b-x, x-a)` |
| `quad(f, low, up; rtol, atol, max_levels, cache)` | Adaptive 2D/3D integration (accepts `SVector` or vectors of reals) |
| `quad_split(f, c; rtol, atol, max_levels, cache)` | Split domain `[-1, 1]` at singularity `c` and integrate |
| `quad_split(f, c, low, up; rtol, atol, max_levels, cache)` | Split domain `[low, up]` at singularity `c` and integrate |

### Pre-computed Integration Functions
| Function | Description |
|----------|-------------|
| `tanhsinh(N)` | Generate Float64 nodes/weights for `N` points |
| `tanhsinh(T, N)` | Generate nodes `x`, weights `w`, step `h` for type `T` |
| `tanhsinh(T, Val(N))` | Generate `SVector` nodes/weights for small `N` (SIMD-ready) |
| `integrate1D(f, N)` | Integrate `f` over `[-1, 1]` using `N` points |
| `integrate1D(T, f, N)` | Integrate `f` over `[-1, 1]` using `N` points in type `T` |
| `integrate1D(f, x, w, h)` | Integrate `f` over `[-1, 1]` using pre-computed nodes |
| `integrate1D(f, low, up, x, w, h)` | Integrate `f` over `[low, up]` (`low, up <: Real`) using pre-computed nodes |
| `integrate2D(f, x, w, h)` | 2D integration over `[-1, 1]^2` |
| `integrate2D(f, low, up, x, w, h)` | 2D integration over rectangle (`SVector` or vector bounds) |
| `integrate3D(f, x, w, h)` | 3D integration over `[-1, 1]^3` |
| `integrate3D(f, low, up, x, w, h)` | 3D integration over box (`SVector` or vector bounds) |

### SIMD-Accelerated Variants

| Function | Description |
|----------|-------------|
| `integrate1D_avx(f, x, w, h)` | SIMD 1D integration over `[-1, 1]` |
| `integrate1D_avx(f, low, up, x, w, h)` | SIMD 1D integration over `[low, up]` |
| `integrate2D_avx(f, low, up, x, w, h)` | SIMD 2D integration over rectangle |
| `integrate3D_avx(f, x, w, h)` | SIMD 3D integration over `[-1, 1]^3` |
| `integrate3D_avx(f, low, up, x, w, h)` | SIMD 3D integration over box |

### Adaptive Integration Functions

| Function | Description |
| :--- | :--- |
| `adaptive_integrate_1D(T, f, a, b; rtol, atol, max_levels, warn, cache)` | Adaptive 1D with explicit type |
| `adaptive_integrate_1D_cmpl(T, f, a, b; rtol, atol, max_levels, warn, cache)` | Adaptive 1D with complement interface |
| `adaptive_integrate_2D(T, f, low, up; rtol, atol, max_levels, warn, cache)` | Adaptive 2D integration |
| `adaptive_integrate_3D(T, f, low, up; rtol, atol, max_levels, warn, cache)` | Adaptive 3D integration |

### Adaptive Cache Utilities

| Function | Description |
| :--- | :--- |
| `adaptive_cache_1D(T; max_levels=16, complement=false)` | Reusable cache for 1D adaptive calls |
| `adaptive_cache_2D(T; max_levels=8)` | Reusable cache for 2D adaptive calls |
| `adaptive_cache_3D(T; max_levels=5)` | Reusable cache for 3D adaptive calls |

---

## Benchmarks

### Benchmark Environment
- **CPU**: Intel(R) Core(TM) Ultra 7 155U

### Results

Comparison across:
- `FastTanhSinhQuadrature.jl` (`adaptive_integrate_*`, `quad`, `_avx`)
- `QuadGK.jl`
- `HCubature.jl`
- `Cubature.jl` (`h` and `p`)
- `Cuba.jl` (`Vegas`, `Divonne`, `Cuhre`, for >1D)
- `FastGaussQuadrature.jl` (1D)

Methodology:
- Domain: case-dependent (most tests use `[-1,1]^d`; oscillatory 1D uses `[-π,π]`)
- Tolerances: `rtol = 1e-6`, `atol = 1e-8`
- Max evaluations (external adaptive solvers): `200000`
- Timing: `@belapsed` with interpolation (`samples=3`, `evals=1`)
- One warm call is executed before each timed benchmark.
- Adaptive cache construction is excluded from timed regions (caches are prebuilt).

The plot below summarizes the speedup of the `quad` and `_avx` interfaces relative to the fastest accurate competing method on each benchmark case.

<p align="center">
  <img src="JOSS_paper/benchmark_summary.svg" alt="Benchmark summary plot" width="850">
</p>

| Dim | Function | FTS adaptive | FTS quad | FTS avx | QuadGK | HCubature | Cubature h | Cubature p | Cuba Vegas | Cuba Divonne | Cuba Cuhre | FastGauss |
| :-- | :------- | ----------: | -------: | ------: | -----: | --------: | ---------: | ---------: | ---------: | -----------: | ---------: | --------: |
| 1D | 1/(1+25x^2) | 735.0000 | 628.0000 | 474.0000 | 1612.0000 | 1.620e+04 | 1.644e+04 | 2.007e+04 | n/a | n/a | n/a | 176.0000 |
| 1D | 1/sqrt(1-x^2) | 574.0000 | 528.0000 | 523.0000 | 1.154e+04 | 2.189e+05 | 2.389e+05 | 2834.0000 * | n/a | n/a | n/a | 1.146e+04 * |
| 1D | log(1-x) | 1013.0000 | 592.0000 | 442.0000 | 5237.0000 | 6.031e+04 | 7.402e+04 | 1356.0000 * | n/a | n/a | n/a | 1.208e+04 |
| 1D | sin^2(1000x) | 1.770e+05 | 2.542e+05 | 2.909e+04 | 1.142e+04 | 1.313e+05 | 1.041e+05 | 1301.0000 * | n/a | n/a | n/a | 3.748e+04 |
| 1D | x^6 - 2x^3 + 0.5 | 1306.0000 | 803.0000 | 335.0000 | 1898.0000 | 4813.0000 | 2249.0000 | 8949.0000 | n/a | n/a | n/a | 118.0000 |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | 3663.0000 | 3660.0000 | 1865.0000 | n/a | 5.218e+07 | 2.832e+07 | 2433.0000 * | 9.357e+07 * | 1.140e+08 * | 9.744e+07 * | n/a |
| 2D | exp(x+y) | 1.864e+04 | 2.438e+04 | 5720.0000 | n/a | 8.769e+04 | 4.420e+04 | 3.861e+04 | 9.402e+07 * | 8.362e+07 | 6.697e+04 | n/a |
| 2D | x^2 + y^2 | 2678.0000 | 2153.0000 | 1968.0000 | n/a | 9986.0000 | 5980.0000 | 9098.0000 | 8.497e+07 * | 7.713e+07 | 6.858e+04 | n/a |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | 1.168e+05 | 1.153e+05 | 7.367e+04 | n/a | 3.339e+07 * | 3.151e+07 * | 5400.0000 * | 1.320e+08 * | 1.391e+08 * | 1.111e+08 * | n/a |
| 3D | exp(x+y+z) | 1.049e+06 | 1.069e+06 | 3.401e+05 | n/a | 3.619e+05 | 3.170e+05 | 7.207e+05 | 1.178e+08 * | 1.283e+08 * | 1.795e+05 | n/a |
| 3D | x^2*y^2*z^2 | 7.994e+04 | 7.634e+04 | 2.207e+04 | n/a | 1.345e+07 | 8.682e+06 | 2.441e+04 | 1.055e+08 * | 1.318e+08 * | 4.454e+07 | n/a |

`*` indicates the method did not meet the requested tolerance.

Full raw results are written by the benchmark script to:
- `benchmark/results/timings.csv`
- `benchmark/results/timings_full.md`
- `benchmark/results/timings_summary.md`

## Other Julia quadrature packages

Use the package choice that matches your integrand class.

From this repository’s current benchmark suite (`rtol=1e-6`, `atol=1e-8`):

- This package is excellent for endpoint-dominated 1D integrands. Examples:
  - `1/sqrt(1-x^2)`: `quad` (`528 ns`) vs `QuadGK` (`1.154e+04 ns`)
  - `log(1-x)`: `quad` (`592 ns`) vs `QuadGK` (`5237 ns`)
- This package is not the best adaptive default for strongly oscillatory 1D integrals. Example:
  - `sin^2(1000x)`: `quad` (`2.542e+05 ns`) vs `QuadGK` (`1.142e+04 ns`)
  - Even precomputed `integrate1D_avx` (`2.909e+04 ns`) is still slower than `QuadGK` on this case.
- For smooth low-order 1D functions, all three approaches can be good, and fixed Gaussian rules can win when `N` is calibrated. Example:
  - `x^6 - 2x^3 + 0.5`: `FastGauss` (`118 ns`) vs `quad` (`803 ns`) vs `QuadGK` (`1898 ns`)
- In 2D/3D box integrals, this package is often very competitive and consistently meets tolerance in our endpoint-singular test cases where several alternatives exceed tolerance/eval budgets.

Practical selection guide:

- Use `FastTanhSinhQuadrature.jl` when endpoint behavior is the main difficulty, or when you can reuse precomputed nodes/caches across many calls.
- Use `QuadGK.jl` as the first choice for many oscillatory or locally difficult 1D integrals.
- Use `FastGaussQuadrature.jl` when a fixed high-order Gaussian rule (possibly weighted) matches your problem well.
- For multidimensional cubature, also compare `HCubature.jl`, `Cubature.jl`, and `Cuba.jl`.
- If function values are only available on a fixed grid (not callable at arbitrary points), use sampled-data packages like `Trapz.jl`, `Romberg.jl`, or `NumericalIntegration.jl`.

## Contributing

Contributions are welcome through GitHub issues and pull requests. Please include tests or benchmark updates when relevant, especially for changes that affect numerical accuracy or performance.

Particularly useful contributions include:
- new benchmark cases and reproducible performance reports
- additional tests for edge cases, precision-specific behavior, and multidimensional interfaces
- documentation improvements and extra usage examples
- performance work such as SIMD-oriented optimizations, threading, and architecture-specific tuning
- extensions aligned with the current roadmap, including additional exponential or double-exponential formulas and more general `n`-dimensional support
