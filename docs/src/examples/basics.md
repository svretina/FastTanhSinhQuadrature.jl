# Basic Usage Examples

This section provides practical examples for using `FastTanhSinhQuadrature.jl` in common scenarios.

## 1. Simple 1D Integration

Integrating standard mathematical functions is straightforward.

```julia
using FastTanhSinhQuadrature

# Function to integrate: f(x) = exp(-x^2)
f(x) = exp(-x^2)

# Integrate over [-1, 1]
result = integrate(f, 10) # 10 levels of recursion
println("Integral of exp(-x^2) on [-1, 1]: $result")
```

Integrating over an arbitrary interval $[a, b]$:

```julia
# Integrate sin(x) from 0 to pi
result_sin = integrate(sin, 0.0, π, 10)
println("Integral of sin(x) on [0, π]: $result_sin") # Should be 2.0
```

## 2. High Precision Integration (`Double64`, `BigFloat`)

One of the main strengths of Tanh-Sinh quadrature is its ability to handle high-precision arithmetic efficiently.

```julia
using FastTanhSinhQuadrature
using DoubleFloats

f(x) = exp(x)

# Use Double64 for extended precision
# N=12 typically gives ~32 digits of precision
x, w, h = tanhsinh(Double64, 12)

# Integrate exp(x) on [0, 1]
val = integrate(f, 0.0, 1.0, x, w, h)
println("High precision result: $val")
```

## 3. Dealing with Singularities

Tanh-Sinh quadrature excels at handling endpoint singularities effectively because the quadrature nodes approach the endpoints exponentially fast but never reach them.

### Logarithmic Singularity `log(1-x)`
This function has a singularity at $x=1$.

```julia
f_sing(x) = log(1-x)

# Integrate on [-1, 1]
# The singularity at x=1 is automatically handled
val = integrate(f_sing, 10) 
println("Integral of log(1-x) on [-1, 1]: $val")
```

### Inverse Square Root `1/sqrt(x)` at `x=0`

```julia
# Integrate 1/sqrt(x) from 0 to 1
f_sqrt(x) = 1.0 / sqrt(x)

x, w, h = tanhsinh(Float64, 10)
val = integrate(f_sqrt, 0.0, 1.0, x, w, h)
println("Integral of 1/sqrt(x) on [0, 1]: $val") # Should be 2.0
```

## 4. Adaptive Integration

If you require a specific tolerance rather than specifying a fixed number of points, use `adaptive_integrate`.

```julia
# Integrate generic function to 1e-12 tolerance
val = adaptive_integrate(x -> cos(x)^2, 0.0, 2π, tol=1e-12)
println(val)
```
