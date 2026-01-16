# Advanced & Multidimensional Integration

## 1. Multidimensional Integration (2D, 3D)

`FastTanhSinhQuadrature.jl` supports multidimensional integration natively.

```julia
using FastTanhSinhQuadrature
using StaticArrays

# Define a 2D function f(x, y) = x^2 + y^2
f_2d(x, y) = x^2 + y^2

# Integration bounds: [-1, 1] x [-1, 1]
low = SVector(-1.0, -1.0)
up  = SVector(1.0, 1.0)

# Generate quadrature scaling
x, w, h = tanhsinh(Float64, 10)

val_2d = integrate(f_2d, low, up, x, w, h)
println("2D Integral: $val_2d") 
```

## 2. Pre-computing Nodes for Performance

For performance-critical code where you integrate many functions or run loops, always pre-calculate the quadrature nodes `(x, w, h)`.

```julia
# Pre-calculate once (expensive operation)
x, w, h = tanhsinh(Float64, 15)

# Reuse many times (cheap operation)
for i in 1:100
    param = i / 100.0
    f(t) = exp(-param * t^2)
    val = integrate(f, 0.0, 10.0, x, w, h)
    # ... use val
end
```

## 3. SIMD Acceleration with `integrate_avx`
For functions compatible with `LoopVectorization.jl`, you can achieve significant speedups.

```julia
using FastTanhSinhQuadrature

f_poly(x) = x^12 + 3x^5 - 2x

x, w, h = tanhsinh(Float64, 50) # Heavy quadrature

# Standard
@time integrate(f_poly, x, w, h)

# AVX Optimized
@time integrate_avx(f_poly, x, w, h) 
```
