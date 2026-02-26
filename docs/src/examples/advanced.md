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
println("Error: $(val - 8//3)")
# Error: 4.440892098500626e-16

# Or using the high-level quad interface
val = quad(f_2d, low, up)
println("Error: $(val - 8//3)")
# Error: -1.3322676295501878e-15
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
println("Error: $(val - 1//8)")
# Error: 0.0
```

### 3D Anisotropic Box with Analytic Check

```julia
using StaticArrays

x, w, h = tanhsinh(Float64, 80)
low = SVector(-2.0, -1.0, 0.0)
up  = SVector(3.0, 2.0, 4.0)

f(x, y, z) = x^2 + 2y + 3z^2
val = integrate3D(f, low, up, x, w, h)

# Exact integral:
exact = 1160//1

println("Error: $(val - exact)")
# Error: -1.8189894035458565e-12
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

For functions compatible with `LoopVectorization.jl`, you can achieve significant speedups:

```julia
using FastTanhSinhQuadrature
using BenchmarkTools

f_poly(x) = x^12 + 3x^5 - 2x

x, w, h = tanhsinh(Float64, 500)

julia> @benchmark integrate1D_avx($f_poly, $x, $w, $h)
BenchmarkTools.Trial: 10000 samples with 980 evaluations per sample.
 Range (min … max):  65.323 ns … 224.431 ns  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     73.924 ns               ┊ GC (median):    0.00%
 Time  (mean ± σ):   76.228 ns ±   7.036 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

           █ ▄▇                                                 
  ▂▂▂▂▂▂▂▂██▅██▅▇▆▄▅▄▃▄▅▃▃▅▄▃▄▄▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂ ▃
  65.3 ns         Histogram: frequency by time          106 ns <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark integrate1D($f_poly, $x, $w, $h)
BenchmarkTools.Trial: 10000 samples with 8 evaluations per sample.
 Range (min … max):  3.586 μs …  16.761 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     4.526 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   4.680 μs ± 617.014 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

            ▁▁█▄                                               
  ▂▂▂▂▁▂▁▂▃▄████▆▅▄▄▃▄▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂ ▃
  3.59 μs         Histogram: frequency by time        7.84 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.
t1 = @belapsed integrate1D_avx($f_poly, $x, $w, $h)
t2 = @belapsed integrate1D($f_poly, $x, $w, $h)
println("Speedup: $(t2/t1)")
# Speedup: 60.6060266328906
```

### 2D/3D SIMD

```julia
using StaticArrays

low = SVector(-1.0, -1.0)
up  = SVector(1.0, 1.0)

# SIMD-accelerated 2D integration
val = integrate2D_avx((x,y) -> x^2 + y^2, low, up, x, w, h)
aval = 8//3
println("Error: $(val - aval)")
# Error: -3.162230512660876733012324488637462686944251599484320216324290926764451115142467e-57
```

## 4. Handling Internal Singularities

For functions with singularities inside the domain, use `quad_split`:

```julia
# 2D example: singularity at (0, 0)
f_sing(x, y) = 1 / sqrt(x^2 + y^2)

center = SVector(0.0, 0.0)
low = SVector(-1.0, -1.0)
up  = SVector(1.0, 1.0)

val = quad_split(f_sing, center, low, up)
aval = 8*asinh(one(val))
println("Error: $(val - aval)")
# Error: -9.14823772291129e-14
```

For 3D:

```julia
using StaticArrays

f_abs3(x, y, z) = 1 / sqrt(abs(x * y * z))
center3 = SVector(0.0, 0.0, 0.0)
low3    = SVector(-1.0, -1.0, -1.0)
up3     = SVector(1.0, 1.0, 1.0)

val3 = quad_split(f_abs3, center3, low3, up3)
println("3D split integral: $val3")  # ≈ 64
```

## 5. Arbitrary Precision with BigFloat

```julia
setprecision(BigFloat, 256)  # 256-bit precision

x, w, h = tanhsinh(BigFloat, 150)
f(x) = log(1 + x)

val = integrate1D(f, x, w, h)
println("High precision: $val")
println("Exact value: $-2 + log(4)")
println("Error: $(val - (-2+ log(BigFloat(4))))")
# 1.351669432265398138420644011346412477230033406653628853271033861612521133030749e-62
```
