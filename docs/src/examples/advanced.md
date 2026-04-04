# Advanced Usage and Performance Tuning

This page focuses on fixed-grid workflows, repeated integrations, SIMD, higher
dimensions, and arbitrary precision.

## 1. Multidimensional Integration

`FastTanhSinhQuadrature.jl` supports adaptive and fixed-grid 2D/3D integration.

### 2D Integration

```julia
using FastTanhSinhQuadrature
using StaticArrays

f_2d(x, y) = x^2 + y^2

low = SVector(-1.0, -1.0)
up  = SVector(1.0, 1.0)

# High-level adaptive interface
val_quad = quad(f_2d, low, up)

# Fixed-grid interface
x, w, h = tanhsinh(Float64, 40)
val_fixed = integrate2D(f_2d, low, up, x, w, h)

println("quad error: $(val_quad - 8/3)")
println("fixed-grid error: $(val_fixed - 8/3)")
```

### 3D Integration

```julia
using FastTanhSinhQuadrature
using StaticArrays

f_3d(x, y, z) = x * y * z

low = SVector(0.0, 0.0, 0.0)
up  = SVector(1.0, 1.0, 1.0)

val = quad(f_3d, low, up)
println("Error: $(val - 1/8)")
```

## 2. Reusing Pre-Computed Nodes

If you integrate many related functions, pre-compute nodes once and reuse them:

```julia
using FastTanhSinhQuadrature

x, w, h = tanhsinh(Float64, 100)

for i in 1:100
    α = i / 100
    f(t) = exp(-α * t^2)
    val = integrate1D(f, 0.0, 10.0, x, w, h)
    # ... use val
end
```

This is often the best pattern when the integration domain and precision are
fixed but the integrand changes repeatedly.

### Reusing Adaptive Caches

If you prefer adaptive refinement (`quad`) but still want low overhead across
many calls, prebuild and reuse the adaptive cache:

```julia
using FastTanhSinhQuadrature

cache = adaptive_cache_2D(Float64; max_levels=8)
low = [-1.0, -1.0]
up  = [1.0, 1.0]

for α in (0.5, 1.0, 2.0)
    f(x, y) = exp(-α * (x^2 + y^2))
    val = quad(f, low, up; rtol=1e-8, atol=1e-10, cache=cache)
    println((α, val))
end
```

For 1D and 3D, use `adaptive_cache_1D` and `adaptive_cache_3D` respectively.

## 3. Choosing `N` for `integrate1D` / `integrate1D_avx`

The fixed-grid interfaces do not take a tolerance and do not estimate their own
error. A practical approach is to double `N` until two successive fixed-grid
results satisfy the same stopping rule used by `quad`:

```julia
using FastTanhSinhQuadrature

function choose_nodes_1d_avx(f, low, up;
    T=Float64,
    candidates=(16, 32, 64, 128, 256, 512),
    rtol=sqrt(eps(T)),
    atol=zero(T),
)
    prev_val = nothing

    for N in candidates
        x, w, h = tanhsinh(T, N)
        val = integrate1D_avx(f, T(low), T(up), x, w, h)

        if prev_val !== nothing
            err_est = abs(val - prev_val)
            target = max(T(atol), T(rtol) * abs(val))
            if err_est <= target
                return (N=N, value=val, x=x, w=w, h=h)
            end
        end

        prev_val = val
    end

    error("No candidate N met the requested tolerance.")
end

cfg = choose_nodes_1d_avx(exp, 0.0, 1.0; rtol=1e-8, atol=1e-12)
println(cfg.N)
```

If you already trust `quad` on a representative integrand family, another
practical workflow is to use `quad(...; rtol, atol)` as a reference and choose
the smallest `N` for which `integrate1D` or `integrate1D_avx` matches that
reference to the desired accuracy.

## 4. SIMD Acceleration with `_avx`

For `Float32` and `Float64`, the `_avx` variants use `LoopVectorization.jl` and
can be substantially faster for repeated fixed-grid integrations:

```julia
using FastTanhSinhQuadrature
using BenchmarkTools

f_poly(x) = x^12 + 3x^5 - 2x
x, w, h = tanhsinh(Float64, 256)

t_scalar = @belapsed integrate1D($f_poly, $x, $w, $h)
t_simd   = @belapsed integrate1D_avx($f_poly, $x, $w, $h)

println("Speedup: $(t_scalar / t_simd)")
```

Use the `_avx` path after you have already selected a suitable `N`.

## 5. Arbitrary Precision

Tanh-Sinh quadrature is often attractive when you need more than `Float64`
precision:

```julia
using FastTanhSinhQuadrature
using DoubleFloats

val = quad(exp, Double64(0), Double64(1); rtol=1e-30)
println(val)
```

Or with `BigFloat`:

```julia
using FastTanhSinhQuadrature

setprecision(BigFloat, 256)

x, w, h = tanhsinh(BigFloat, 150)
val = integrate1D(log1p, BigFloat(-1), BigFloat(1), x, w, h)
println(val)
```

## 6. Splitting Multidimensional Singularities

For internal singularities in 2D or 3D, use `quad_split` to turn them into
boundary singularities on subdomains:

```julia
using FastTanhSinhQuadrature
using StaticArrays

f_sing(x, y) = 1 / sqrt(x^2 + y^2)

center = SVector(0.0, 0.0)
low = SVector(-1.0, -1.0)
up  = SVector(1.0, 1.0)

val = quad_split(f_sing, center, low, up)
println(val)
```
