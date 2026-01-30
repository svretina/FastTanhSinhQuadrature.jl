---
title: 'FastTanhSinhQuadrature.jl: High-performance Tanh-Sinh numerical integration in Julia'
tags:
  - Julia
  - numerical integration
  - quadrature
  - tanh-sinh
  - SIMD
authors:
  - name: Stamatis Vretinaris
    orcid: 0000-0001-7575-813X
    affiliation: 1
affiliations:
  - name: Institute for Mathematics, Astrophysics and Particle Physics, Radboud University, Heyendaalseweg 135, 6525 AJ Nijmegen, The Netherlands
    index: 1
  - name: Albert-Einstein-Institut, Max-Planck-Institut für Gravitationsphysik, Callinstraße 38, 30167 Hannover, Germany
    index: 2
date: 26 January 2026
bibliography: paper.bib
---

# Summary

Numerical integration is a cornerstone of scientific computing, essential for evaluating integrals that cannot be solved analytically. The Tanh-Sinh (or Double Exponential) quadrature, originally proposed by @Takahasi1973, is a powerful technique known for its high accuracy and efficiency, particularly for integrands with endpoint singularities. `FastTanhSinhQuadrature.jl` provides a high-performance, arbitrary-precision implementation of this method in Julia. It leverages modern compiler technologies to achieve significant speedups over traditional implementations while maintaining rigorous mathematical precision.

# Statement of Need

In many fields of physics and engineering, researchers encounter integrals with singularities at the boundaries. Standard Gaussian quadrature rules often fail or require an excessive number of points to converge in these cases. While other integration libraries exist, they don't offer the combination of:
1.  **Robustness**: Handling singularities automatically without manual coordinate transformations.
2.  **Performance**: Utilizing SIMD (Single Instruction, Multiple Data) instructions for rapid evaluation.
3.  **Flexibility**: Supporting arbitrary precision types (e.g., `BigFloat`) and multidimensional integration.

`FastTanhSinhQuadrature.jl` addresses these needs by implementing a rigorous Tanh-Sinh scheme with an optimized "window selection" strategy enabling SIMD-accelerated execution paths. This makes it an ideal tool for large-scale simulations where both speed and precision are critical.

# State of the field

Numerical integration is a widely solved problem, with mature implementations available in libraries such as **Boost** (C++), **SciPy** (Python), and **mpmath** (Python). Within the Julia ecosystem, packages like `QuadGK.jl` and `HCubature.jl` provide efficient adaptive Gauss-Kronrod and h-adaptive cubature methods, respectively. However, implementations specifically of the **Tanh-Sinh quadrature** are less common and often face a fundamental performance bottleneck.

Most existing Tanh-Sinh implementations, both in Julia and other languages, rely on **dynamic condition checks** within the summation loop. These algorithms typically iterate until the current term drops below a certain tolerance to handle the infinite domain truncation and detect underflow. While mathematically sound, this approach introduces conditional branching inside the innermost loop, which prevents modern compilers from applying **Single Instruction, Multiple Data (SIMD)** vectorization. Consequently, these solvers are restricted to scalar execution speeds.

`FastTanhSinhQuadrature.jl` was developed to bridge this gap. Rather than modifying existing scalar implementations, we built a solver from the ground up centered on the **"window selection"** strategy described by @Vanherck2020. This method analytically pre-calculates the optimal step size $h$ and the truncation index $n$ based strictly on the floating-point precision (e.g., `Float64`, `BigFloat`). By determining the bounds *a priori*, we eliminate the need for runtime convergence checks. This creates a branch-free inner loop that allows the compiler (specifically via `LoopVectorization.jl`) to generate highly optimized SIMD instructions.

This design distinguishes `FastTanhSinhQuadrature.jl` from existing tools. By ensuring the algorithm is branch-free, we preserve the high accuracy of the Tanh-Sinh method while allowing for significant SIMD acceleration, a combination rarely found in arbitrary-precision quadrature libraries. Furthermore, the package is integrated into the SciML ecosystem as a backend for `Integrals.jl`.

# Software Design

## API Architecture

The package exposes a two-tier API designed to balance ease of use and maximum performance:

1. **High-Level API** (`quad`, `quad_split`): Intended for general purpose use, these functions perform fully adaptive integration. They automatically handle edge cases such as zero-width domains (`low == up`) and inverted integration bounds (`low > up`), normalizing them before integration. They also transparently manage domain singularities through `quad_split`, which splits the domain into subdomains and integrates them separately.
2. **Low-Level API** (`tanhsinh`, `integrateND`, `integrateND_avx`): Targeted at performance-critical loops, these functions allow users to pre-compute quadrature nodes and weights once and reuse them across millions of integrals. This avoids the overhead of node generation and memory allocation in tight loops. `integrateND_avx` is a SIMD-optimized version of `integrateND`.

## Efficient Adaptive Integration

The adaptive integration routines (`adaptive_integrate_ND`) implement a refinement scheme that maximizes computational efficiency. At each level of refinement, the step size $h$ is halved. Crucially, the algorithm reuses all function evaluations from previous levels. Since $t_{new} = k \cdot (h_{old}/2)$, the odd-indexed points at the new level are unique, while the even-indexed points correspond exactly to the nodes from the previous level. The solver only evaluates functions at these new odd indices and combines them with the stored weighted sum from the previous step.

Additionally, the algorithm exploits the intrinsic symmetry of the Tanh-Sinh transformation ($w(t) = w(-t)$ and $x(-t) = -x(t)$ in the symmetric domain). For every evaluation point $t_k > 0$, the solver computes contributions for both the positive and negative coordinate pairs (e.g., $f(x_0 + \Delta x \cdot x_k) + f(x_0 - \Delta x \cdot x_k)$) using a single weight lookup. This effectively halves the number of expensive `tanh`/`cosh` calls and memory accesses for the weights.

## Memory Efficiency & Static Allocation

To further optimize performance, the library stores only the non-negative quadrature nodes and weights, reducing memory footprint by 50%. This compactness allows for a powerful optimization: for moderate node counts (typically $N < 128$), the pre-computed nodes and weights can be stored in **StaticArrays**. This enables the `tanhsinh(Float64, Val(N))` interface to return stack-allocated, fixed-size arrays. Using these static arrays eliminates heap allocations entirely during the integration phase, a critical feature for high-performance nested loops or multi-threaded workloads.

## Window Selection

A critical implementation detail in Tanh-Sinh quadrature is determining the truncation limits for the infinite sum $I_h$, effectively defining the finite window $[-n, n]$ used in the approximation $I_{h}^n$. If $n$ is too small, significant contributions to the integral are discarded; if $n$ is too large, the weights $\Psi'(t_i)$ decay below machine precision, leading to unnecessary computations.

Standard practice often involves conditional checks within the quadrature loop to detect underflow and terminate the summation. However, such branching logic prevents the compiler from utilizing Single Instruction, Multiple Data (SIMD) vectorization, severely limiting performance. `FastTanhSinhQuadrature.jl` avoids this issue by adopting the methodology described by @Vanherck2020. This method pre-calculates the optimal step size $h$ and the truncation index $n$ to minimize error for a given floating-point precision. This ensures consistent accuracy across `Float64`, `BigFloat`, and other numeric types while maintaining a loop structure that is amenable to SIMD optimization.

## SIMD Optimization

To maximize performance on modern CPUs, the library utilizes `LoopVectorization.jl`. The core integration loops are designed to be "SIMD-friendly," allowing the compiler to process multiple data points simultaneously. This is achieved by carefully structuring the arrays of nodes and weights to align with vector registers. Benchmarks indicate that this approach yields a 2-3x speedup compared to standard scalar implementations for common integrands.

## Type Stability

One of the key strengths of the Julia language is its support for generic programming. `FastTanhSinhQuadrature.jl` is designed to be fully type-stable, meaning that the Julia compiler can infer return types solely from input types.
This design extends to the support of arbitrary numeric types. The package does not rely on hardcoded constants for `Float64`; instead, the quadrature parameters $h$ and $n$ are derived dynamically based on the machine epsilon of the number type used. Consequently, the package supports any custom data type $T$, such as `BigFloat`, `Double64` from `DoubleFloats.jl`, or other extended-precision types, provided that the type implements standard arithmetic operations alongside `eps(T)` and `prevfloat(T)`. The `eps(T)` function is essential for determining the target error tolerance and the summation window size, while `prevfloat(T)` is utilized to safely handle integration bounds near singularities without triggering domain errors.

# Research Impact

The impact of `FastTanhSinhQuadrature.jl` is evidenced by its immediate adoption within the Julia scientific computing ecosystem and its validated performance improvements over existing standard libraries.

## Community Adoption

The package has been integrated as a backend for `Integrals.jl`, the primary numerical integration library in the **SciML** (Scientific Machine Learning) ecosystem. This integration positions `FastTanhSinhQuadrature.jl` as a core tool for users requiring robust handling of boundary singularities. It allows researchers to seamlessly switch to Tanh-Sinh quadrature via the unified interface of `Integrals.jl` without modifying their problem definitions, facilitating widespread use in differential equation solvers and physics simulations.

## Performance Benchmarks

We conducted comparative benchmarks against `FastGaussQuadrature.jl`, a widely used Julia library implementing Gauss-Legendre quadrature. While Gaussian quadrature is typically favored for smooth functions, our SIMD-optimized Tanh-Sinh implementation (`integrate1D_avx`) demonstrates competitive and often superior runtime performance on modern hardware.

For example, when integrating $e^x$ over $[-1, 1]$ with 500 points, our SIMD solver achieves a **2.1x speedup** compared to `FastGaussQuadrature.jl` (782 ns vs 1665 ns). For integrands with endpoint singularities, such as $\sqrt{1-x^2}$, the speedup is approximately **2.4x** (313 ns vs 743 ns). These results confirm that `FastTanhSinhQuadrature.jl` effectively eliminates the traditional performance overhead associated with Tanh-Sinh quadrature, validating it as a high-performance solution for both singular and smooth integrands.

Comparison against `FastGaussQuadrature.jl` for various functions.

- `TS`: `integrate1D` (standard Tanh-Sinh)
- `TS SIMD`: `integrate1D_avx` (SIMD-optimized Tanh-Sinh)
- `GQ`: `FastGaussQuadrature.jl` (Gauss-Legendre)

Timings are in nanoseconds (ns) (smaller is better).

| Function | Domain | Points | TS (ns) | TS SIMD (ns) | GQ (ns) | Ratio (TS/GQ) | Ratio (TS SIMD/GQ) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| exp(x) | [-1, 1] | 5 | 33.09 | 14.27 | 33.03 | 1.00 | 0.43 |
| exp(x) | [-1, 1] | 50 | 322.87 | 81.99 | 176.68 | 1.83 | 0.46 |
| exp(x) | [-1, 1] | 500 | 3374.75 | 782.08 | 1665.50 | 2.03 | 0.47 |
| sin(x)^2 | [-1, 1] | 5 | 40.39 | 23.38 | 32.63 | 1.24 | 0.72 |
| sin(x)^2 | [-1, 1] | 50 | 404.32 | 145.48 | 176.75 | 2.29 | 0.82 |
| sin(x)^2 | [-1, 1] | 500 | 4361.71 | 1409.10 | 1768.40 | 2.47 | 0.80 |
| 1/(1+25x^2) | [-1, 1] | 5 | 3.00 | 3.75 | 17.35 | 0.17 | 0.22 |
| 1/(1+25x^2) | [-1, 1] | 50 | 41.82 | 39.17 | 34.77 | 1.20 | 1.13 |
| 1/(1+25x^2) | [-1, 1] | 500 | 477.21 | 213.24 | 320.39 | 1.49 | 0.67 |
| sqrt(1-x^2) | [-1, 1] | 5 | 4.06 | 4.69 | 20.31 | 0.20 | 0.23 |
| sqrt(1-x^2) | [-1, 1] | 50 | 68.01 | 45.31 | 69.68 | 0.98 | 0.65 |
| sqrt(1-x^2) | [-1, 1] | 500 | 636.53 | 313.13 | 743.43 | 0.86 | 0.42 |
| x^2 | [-1, 1] | 5 | 2.01 | 3.04 | 19.86 | 0.10 | 0.15 |
| x^2 | [-1, 1] | 50 | 17.72 | 21.68 | 38.07 | 0.47 | 0.57 |
| x^2 | [-1, 1] | 500 | 213.31 | 54.97 | 273.19 | 0.78 | 0.20 |
| log(1-x) | [-1, 1] | 5 | 38.06 | 30.29 | 34.32 | 1.11 | 0.88 |
| log(1-x) | [-1, 1] | 50 | 371.88 | 187.29 | 214.70 | 1.73 | 0.87 |
| log(1-x) | [-1, 1] | 500 | 4117.38 | 1978.80 | 2117.80 | 1.94 | 0.93 |

# Mathematics

Tanh-sinh quadrature is designed to compute integrals of the form:

$$I = \int_{-1}^{1} f(x) \, \mathrm{d}x$$

The method relies on a variable transformation $x = \Psi(t)$ that maps the finite interval $x \in (-1, 1)$ onto the entire real axis $t \in (-\infty, +\infty)$. The integral is rewritten as:

$$I = \int_{-\infty}^{\infty} g(t) \, \mathrm{d}t, \quad \text{where } g(t) := f(\Psi(t))\Psi'(t)$$

To approximate this integral, the trapezoidal rule is applied over the infinite domain. Given a step size $h$, the approximation is:

$$I_h = h \sum_{n=-\infty}^{\infty} \Psi'(t_n) f(\Psi(t_n))$$

where the evaluation points are equidistant:

$$t_n := nh, \quad n = 0, \pm 1, \pm 2, \dots$$

For computational implementation, the infinite sum is truncated to a finite window $[-n, n]$. This yields the quadrature formula:

$$I_{h}^n = h \sum_{n=-n}^{n} \Psi'(t_n) f(\Psi(t_n))$$

The Tanh-sinh quadrature specifically employs the transformation proposed by @Takahasi1973:

$$\Psi(t) = \tanh(\lambda \sinh(t))$$

which has the derivative:

$$\Psi'(t) = \frac{\lambda \cosh(t)}{\cosh^2(\lambda \sinh(t))}$$

where $\lambda = \pi/2$.


# Usage

## Installation

To install `FastTanhSinhQuadrature.jl`, use the Julia Package Manager. You can add it via the REPL by entering the package manager mode (press `]`) and typing:

```
pkg> add FastTanhSinhQuadrature
```

Alternatively, you can use the `Pkg` API:

```julia
using Pkg
Pkg.add("FastTanhSinhQuadrature")
```

## Basic Integration
The simplest entry point is the `quad` function, which performs adaptive integration:

```julia
using FastTanhSinhQuadrature

# Integrate exp(x) from 0 to 1
val = quad(exp, 0.0, 1.0)
println(val)  # ≈ 1.7182818...

# Handle singularities: 1/sqrt(x) from 0 to 1
val = quad(x -> 1/sqrt(x), 0.0, 1.0)
println(val)  # ≈ 2.0
```

## High-Performance Pre-computation
For repeated integrations, users can pre-compute nodes and weights:

```julia
# Pre-compute nodes (x) and weights (w) for Float64
x, w, h = tanhsinh(Float64, 80)

# Efficiently integrate multiple functions reusing x and w
f1(t) = sin(t)^2
f2(t) = cos(t)^2
integral1 = integrate1D(f1, 0.0, π, x, w, h)
integral2 = integrate1D(f2, 0.0, π, x, w, h)

## For maximum performance, pre-compute nodes and weights with Val() and use '_avx' variants
x, w, h = tanhsinh(Float64, Val(80))
integral1 = integrate1D_avx(f1, 0.0, π, x, w, h)
integral2 = integrate1D_avx(f2, 0.0, π, x, w, h)
```

# Convergence
We test the software on a range of integrands over $[-1, 1]$ and compare the results to the exact values.

![Convergence of Tanh-Sinh Quadrature compared to other methods. The plot illustrates the rapid error decay.](convergence.svg)

# AI usage disclosure

During the development of this package, the author utilized Gemini (Google) for assistance with documentation, debugging and the generation of the first draft of this paper. The author has reviewed and edited all AI-generated content to ensure accuracy and adherence to the package's coding standards.

# Acknowledgements

The author would like to acknowledge the developers of `LoopVectorization.jl` for providing the tools that enabled the performance optimizations in this package.

# References
