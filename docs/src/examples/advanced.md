# Advanced & Multidimensional Integration

## 1. Multidimensional Integration (2D, 3D)

`FastTanhSinhQuadrature.jl` supports multidimensional integration natively using `StaticArrays`.

### 2D Integration

```julia
using FastTanhSinhQuadrature
using StaticArrays

# Define a 2D function f(x, y) = x^2 + y^2
f_2d(x, y) = x^2 + y^2

# Integration bounds: [-1, 1] × [-1, 1]
low = SVector(-1.0, -1.0)
up  = SVector(1.0, 1.0)

# Pre-computed nodes
x, w, h = tanhsinh(Float64, 40)
val = integrate2D(f_2d, low, up, x, w, h)
println("2D Integral: $val")  # ≈ 8/3

# Or using the high-level quad interface
val = quad(f_2d, low, up)
println("2D Integral (quad): $val")
```

### 3D Integration

```julia
# Integrate constant 1 over unit cube [0,1]³
f_3d(x, y, z) = 1.0
val = quad(f_3d, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
println("3D Volume: $val")  # ≈ 1.0

# Pre-computed for performance
x, w, h = tanhsinh(Float64, 30)
low3 = SVector(0.0, 0.0, 0.0)
up3  = SVector(1.0, 1.0, 1.0)
val = integrate3D((x,y,z) -> x*y*z, low3, up3, x, w, h)
println("3D Integral of xyz: $val")  # ≈ 1/8 = 0.125
```

## 2. Pre-computing Nodes for Performance

For performance-critical code where you integrate many functions or run loops, always pre-calculate the quadrature nodes `(x, w, h)`.

```julia
# Pre-calculate once (more expensive operation)
x, w, h = tanhsinh(Float64, 100)

# Reuse many times (cheap operation)
for i in 1:100
    param = i / 100.0
    f(t) = exp(-param * t^2)
    val = integrate1D(f, 0.0, 10.0, x, w, h)
    # ... use val
end
```

## 3. SIMD Acceleration with `_avx` Variants

For functions compatible with `LoopVectorization.jl`, you can achieve significant speedups (2-3x):

```julia
using FastTanhSinhQuadrature

f_poly(x) = x^12 + 3x^5 - 2x

x, w, h = tanhsinh(Float64, 500)

# Standard
@time val1 = integrate1D(f_poly, x, w, h)

# SIMD-accelerated
@time val2 = integrate1D_avx(f_poly, x, w, h)

println("Standard: $val1, SIMD: $val2")
```

### 2D/3D SIMD

```julia
using StaticArrays

low = SVector(-1.0, -1.0)
up  = SVector(1.0, 1.0)

# SIMD-accelerated 2D integration
val = integrate2D_avx((x,y) -> x^2 + y^2, low, up, x, w, h)
```

## 4. Handling Internal Singularities

For functions with singularities inside the domain, use `quad_split`:

```julia
# 2D example: singularity at (0, 0)
f_sing(x, y) = 1 / sqrt(x^2 + y^2 + 0.01)

center = SVector(0.0, 0.0)
low = SVector(-1.0, -1.0)
up  = SVector(1.0, 1.0)

val = quad_split(f_sing, center, low, up)
println("2D integral with near-singularity: $val")
```

## 5. Arbitrary Precision with BigFloat

```julia
setprecision(BigFloat, 256)  # 256-bit precision

x, w, h = tanhsinh(BigFloat, 150)
f(x) = log(1 + x)

val = integrate1D(f, x, w, h)
println("High precision: $val")
# Achieves ~50+ correct decimal digits
```
