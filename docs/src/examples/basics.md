# Basic Usage Examples

This section provides practical examples for using `FastTanhSinhQuadrature.jl` in common scenarios.

## 1. Simple 1D Integration with `quad`

The easiest way to integrate is using the `quad` function, which automatically adapts to your desired tolerance:

```julia
using FastTanhSinhQuadrature

# Integrate exp(x) from 0 to 1
val = quad(exp, 0.0, 1.0)
println("Integral of exp(x) on [0, 1]: $val")  # ≈ e - 1 ≈ 1.7183

# Integrate over default domain [-1, 1]
val = quad(x -> 3x^2)
println("Integral of 3x^2 on [-1, 1]: $val")  # ≈ 2.0
```

## 2. Pre-computed Nodes for Maximum Performance

For repeated integrations, pre-compute nodes and weights once:

```julia
using FastTanhSinhQuadrature

# Generate nodes (x), weights (w), and step size (h)
x, w, h = tanhsinh(Float64, 80)

# Integrate multiple functions efficiently
f1(x) = sin(x)^2
f2(x) = cos(x)^2

res1 = integrate1D(f1, 0.0, π, x, w, h)
res2 = integrate1D(f2, 0.0, π, x, w, h)
println("Integrals: $res1, $res2")  # Both ≈ π/2

# Integration over [-1, 1] (default domain)
val = integrate1D(exp, x, w, h)
println("Integral of exp(x) on [-1, 1]: $val")
```

## 3. High Precision Integration (`Double64`, `BigFloat`)

One of the main strengths of Tanh-Sinh quadrature is its ability to handle high-precision arithmetic efficiently.

```julia
using FastTanhSinhQuadrature
using DoubleFloats

f(x) = exp(x)

# Use Double64 for extended precision (~32 decimal digits)
val = quad(f, Double64(0), Double64(1); tol=1e-30)
println("High precision result: $val")

# Or with pre-computed nodes
x, w, h = tanhsinh(Double64, 100)
val = integrate1D(f, Double64(0), Double64(1), x, w, h)
println("Pre-computed result: $val")
```

## 4. Dealing with Singularities

Tanh-Sinh quadrature excels at handling endpoint singularities automatically.

### Logarithmic Singularity `log(1-x)`

This function has a singularity at $x=1$:

```julia
f_sing(x) = log(1-x)

# The singularity at x=1 is automatically handled
val = quad(f_sing, -1.0, 1.0)
println("Integral of log(1-x) on [-1, 1]: $val")  # ≈ -0.6137
```

### Inverse Square Root `1/sqrt(x)` at `x=0`

```julia
# Integrate 1/sqrt(x) from 0 to 1
f_sqrt(x) = 1.0 / sqrt(x)

x, w, h = tanhsinh(Float64, 80)
val = integrate1D(f_sqrt, 0.0, 1.0, x, w, h)
println("Integral of 1/sqrt(x) on [0, 1]: $val")  # ≈ 2.0
```

### Internal Singularities with `quad_split`

For singularities inside the domain, use `quad_split`:

```julia
# 1/sqrt(|x|) has a singularity at x=0
f_abs(x) = 1 / sqrt(abs(x))

# Split at the singularity point
val = quad_split(f_abs, 0.0, -1.0, 1.0)
println("Integral of 1/sqrt(|x|) on [-1, 1]: $val")  # ≈ 4.0
```

## 5. SIMD-Accelerated Integration

For `Float32`/`Float64`, use the `_avx` variants for maximum speed:

```julia
x, w, h = tanhsinh(Float64, 100)

# Standard integration
val1 = integrate1D(exp, x, w, h)

# SIMD-accelerated (2-3x faster)
val2 = integrate1D_avx(exp, x, w, h)

println("Standard: $val1, SIMD: $val2")
```
