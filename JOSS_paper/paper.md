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

Numerical integration is a cornerstone of scientific computing, essential for evaluating integrals that cannot be solved analytically. The Tanh-Sinh (or Double Exponential) quadrature, originally proposed by @Takahasi1973, is a powerful technique known for its high accuracy and efficiency, particularly for integrands with endpoint singularities. `FastTanhSinhQuadrature.jl` provides a high-performance, arbitrary-precision implementation of this method in Julia. It leverages modern compiler technologies to achieve significant speedups over traditional implementations while maintaining rigorous mathematical precision.

# Statement of Need

In many fields of physics and engineering, researchers encounter integrals with singularities at the boundaries, such as those found in potential theory or quantum field calculations. Standard Gaussian quadrature rules often fail or require an excessive number of points to converge in these cases. While other integration libraries exist, they may not offer the combination of:
1.  **Robustness**: Handling singularities automatically without manual coordinate transformations.
2.  **Performance**: Utilizing SIMD (Single Instruction, Multiple Data) instructions for rapid evaluation.
3.  **Flexibility**: Supporting arbitrary precision types (e.g., `BigFloat`) and multidimensional integration.

`FastTanhSinhQuadrature.jl` addresses these needs by implementing a rigorous Tanh-Sinh scheme with an optimized "window selection" strategy and SIMD-accelerated execution paths. This makes it an ideal tool for large-scale simulations where both speed and precision are critical.

# Mathematics

The Tanh-Sinh quadrature is based on the variable transformation 
$$x = \tanh\left(\frac{\pi}{2} \sinh(t)\right)$$
which maps the finite interval $x \in [-1, 1]$ to the infinite interval $t \in (-\infty, \infty)$. The transformed integral becomes:
$$ \int_{-1}^{1} f(x) dx = \int_{-\infty}^{\infty} f\left(\tanh\left(\frac{\pi}{2} \sinh(t)\right)\right) \cosh\left(\frac{\pi}{2} \sinh(t)\right) \frac{\pi}{2} \cosh(t) dt $$
This transformation results in an integrand that decays double-exponentially as $|t| \to \infty$. Consequently, the trapezoidal rule applied to this transformed integral converges exceptionally fast, even for functions with endpoint singularities.

# Software Design

## Window Selection
A key implementation detail in Tanh-Sinh quadrature is determining the truncation of the infinite sum (the "window"). If the window is too small, accuracy is lost; if too large, numerical underflow occurs. `FastTanhSinhQuadrature.jl` adopts the methodology described by @Vanherck2020. This method dynamically calculates the optimal step size $h$ and the number of points $N$ to minimize error for a given floating-point precision, ensuring consistent accuracy across `Float64`, `BigFloat`, and other types.

## SIMD Optimization
To maximize performance on modern CPUs, the library utilizes `LoopVectorization.jl`. The core integration loops are designed to be "SIMD-friendly," allowing the compiler to process multiple data points simultaneously. This is achieved by carefully structuring the arrays of nodes and weights to align with vector registers. Benchmarks indicate that this approach yields a 2-3x speedup compared to standard scalar implementations for common integrands.

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
```

# Figures

![Convergence of Tanh-Sinh Quadrature compared to other methods. The plot illustrates the rapid error decay.](convergence.svg)

# Research Impact Statement

`FastTanhSinhQuadrature.jl` streamlines the workflow for researchers dealing with complex or singular integrals. It has already been adopted in high-precision physics simulations where standard quadrature methods were insufficient. By providing a reliable and fast integration engine, the library enables:
- More accurate modeling of quantum mechanical systems.
- Faster iteration times in numerical experiments.
- Verification of theoretical results requiring high-precision arithmetic.

The library is designed to be easily extensible and to integrate seamlessly with other Julia packages, fostering a collaborative ecosystem for numerical analysis.

# Acknowledgements

The author would like to acknowledge the developers of `LoopVectorization.jl` for providing the tools that enabled the performance optimizations in this package.

# References
