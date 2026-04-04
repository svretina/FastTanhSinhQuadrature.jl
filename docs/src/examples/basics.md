# Basic Usage Examples

These examples focus on the high-level adaptive interfaces: `quad`, `quad_split`,
and `quad_cmpl`.

## 1. Adaptive 1D Integration with `quad`

For most scalar integrals, `quad` is the easiest entry point:

```julia
using FastTanhSinhQuadrature

# Integrate exp(x) from 0 to 1
val = quad(exp, 0.0, 1.0)
println("Integral of exp(x) on [0, 1]: $val")  # ≈ e - 1

# Integrate over the default interval [-1, 1]
val = quad(x -> 3x^2)
println("Integral of 3x^2 on [-1, 1]: $val")  # ≈ 2.0
```

## 2. Tolerance Control with `rtol` and `atol`

`quad` uses an adaptive successive-refinement estimate. You can control the
requested accuracy with relative and absolute tolerances:

```julia
using FastTanhSinhQuadrature

val = quad(exp, 0.0, 1.0; rtol=1e-10, atol=1e-12)
println(val)
```

If `rtol` is omitted and `atol == 0`, the default is `sqrt(eps(T))`, where `T`
is the floating-point type of the bounds.

## 3. Reusing Adaptive Caches with `quad`

If you evaluate many similar integrals adaptively, reuse a cache:

```julia
using FastTanhSinhQuadrature

cache = adaptive_cache_1D(Float64; max_levels=16)

for α in (0.5, 1.0, 2.0)
    f(x) = exp(-α * x^2)
    val = quad(f, -1.0, 1.0; rtol=1e-10, atol=1e-12, cache=cache)
    println((α, val))
end
```

The same pattern works for `quad_split` and `quad_cmpl`.

## 4. Preserving the Input Floating-Point Type

The adaptive API follows the floating-point type implied by the bounds:

```julia
using FastTanhSinhQuadrature

val32 = quad(exp, 0.0f0, 1.0f0)
println(val32)
println(typeof(val32))  # Float32
```

## 5. Endpoint Singularities

Tanh-Sinh quadrature is especially effective when the integrand has a
singularity at an endpoint.

### Logarithmic Singularity

```julia
using FastTanhSinhQuadrature

f_log(x) = log(1 - x)
val = quad(f_log, -1.0, 1.0)
exact = -2 + log(4)
println("Error: $(val - exact)")
```

### Endpoint-Aware Integrands with `quad_cmpl`

When the formula naturally contains distances to the endpoints, `quad_cmpl`
avoids cancellation by passing them directly to the callback:

```julia
using FastTanhSinhQuadrature

f_cmpl(x, b_minus_x, x_minus_a) = inv(sqrt(b_minus_x * x_minus_a))

val = quad_cmpl(f_cmpl, -1.0, 1.0)
println("Error: $(val - π)")
```

## 6. Internal Singularities with `quad_split`

If the singularity is inside the integration domain, split the domain so that
each subproblem has a boundary singularity instead:

```julia
using FastTanhSinhQuadrature

f_abs(x) = 1 / sqrt(abs(x))
val = quad_split(f_abs, 0.0, -1.0, 1.0)
println("Error: $(val - 4)")
```

## 7. Mixed Real Bounds

The high-level API accepts mixed real bound types and promotes them internally:

```julia
using FastTanhSinhQuadrature

val = quad(x -> x, 0, π)
println("Error: $(val - π^2 / 2)")
```
