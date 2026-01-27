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
  - name: Institute for Mathematics, Astrophysics and Particle Physics,
    Radboud University, Heyendaalseweg 135, 6525 AJ Nijmegen, The Netherlands
    index: 1
  - name: Albert-Einstein-Institut, Max-Planck-Institut für Gravitationsphysik, Callinstraße 38, 30167 Hannover, Germany
    index: 2
date: 26 January 2026
bibliography: paper.bib
---

# Summary

Numerical integration is a cornerstone of scientific computing, essential for evaluating integrals that cannot be solved analytically. The Tanh-Sinh (or Double Exponential) quadrature, originally proposed by Takahasi and Mori (@Takahasi1973), is a powerful technique known for its high accuracy and efficiency, particularly for integrands with endpoint singularities. `FastTanhSinhQuadrature.jl` provides a high-performance, arbitrary-precision implementation of this method in Julia. It leverages modern compiler technologies to achieve significant speedups over traditional implementations while maintaining rigorous mathematical precision.

# Statement of Need

In many fields of physics and engineering, researchers encounter integrals with singularities at the boundaries, such as those found in potential theory or quantum field calculations. Standard Gaussian quadrature rules often fail or require an excessive number of points to converge in these cases. While other integration libraries exist, they don't offer the combination of:
1.  **Robustness**: Handling singularities automatically without manual coordinate transformations.
2.  **Performance**: Utilizing SIMD (Single Instruction, Multiple Data) instructions for rapid evaluation.
3.  **Flexibility**: Supporting arbitrary precision types (e.g., `BigFloat`) and multidimensional integration.

`FastTanhSinhQuadrature.jl` addresses these needs by implementing a rigorous Tanh-Sinh scheme with an optimized "window selection" strategy enabling SIMD-accelerated execution paths. This makes it an ideal tool for large-scale simulations where both speed and precision are critical.

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

The Tanh-sinh quadrature specifically employs the transformation proposed by Takahasi and Mori (@Takahasi1973):

$$\Psi(t) = \tanh(\lambda \sinh(t))$$

which has the derivative:

$$\Psi'(t) = \frac{\lambda \cosh(t)}{\cosh^2(\lambda \sinh(t))}$$

where $\lambda = \pi/2$.

# Software Design

## Window Selection

A critical implementation detail in Tanh-Sinh quadrature is determining the truncation limits for the infinite sum $I_h$, effectively defining the finite window $[-n, n]$ used in the approximation $I_{h}^n$. If $n$ is too small, significant contributions to the integral are discarded; if $n$ is too large, the weights $\Psi'(t_i)$ decay below machine precision, leading to unnecessary computations. 

Standard practice often involves conditional checks within the quadrature loop to detect underflow and terminate the summation. However, such branching logic prevents the compiler from utilizing Single Instruction, Multiple Data (SIMD) vectorization, severely limiting performance. `FastTanhSinhQuadrature.jl` avoids this issue by adopting the methodology described by @Vanherck2020. This method pre-calculates the optimal step size $h$ and the truncation index $n$ to minimize error for a given floating-point precision. This ensures consistent accuracy across `Float64`, `BigFloat`, and other numeric types while maintaining a loop structure that is amenable to SIMD optimization.

## SIMD Optimization

To maximize performance on modern CPUs, the library utilizes `LoopVectorization.jl`. The core integration loops are designed to be "SIMD-friendly," allowing the compiler to process multiple data points simultaneously. This is achieved by carefully structuring the arrays of nodes and weights to align with vector registers. Benchmarks indicate that this approach yields a 2-3x speedup compared to standard scalar implementations for common integrands.

## Type Stability

One of the key strengths of the Julia language is its support for generic programming. `FastTanhSinhQuadrature.jl` is designed to be fully type-stable, meaning that the Julia compiler can infer return types solely from input types.
This design extends to the support of arbitrary numeric types. The package does not rely on hardcoded constants for `Float64`; instead, the quadrature parameters $h$ and $n$ are derived dynamically based on the machine epsilon of the number type used. Consequently, the package supports any custom datatype $T$—such as `BigFloat`, `Double64` from `DoubleFloats.jl`, or other extended-precision types, provided that the type implements standard arithmetic operations alongside `eps(T)` and `prevfloat(T)`. The `eps(T)` function is essential for determining the target error tolerance and the summation window size, while `prevfloat(T)` is utilized to safely handle integration bounds near singularities without triggering domain errors.

# Usage

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

# Acknowledgements

The author would like to acknowledge the developers of `LoopVectorization.jl` for providing the tools that enabled the performance optimizations in this package.

# References
