# Benchmarks

Comparison of `FastTanhSinhQuadrature.jl` vs `FastGaussQuadrature.jl`.

**System**:
- **CPU**: Intel(R) Core(TM) Ultra 7 155U
- **Threads**: 1 (Single-threaded execution)

## Results

**Legend**:
- **TS**: `FastTanhSinhQuadrature.integrate`
- **TS SIMD**: `FastTanhSinhQuadrature.integrate_avx` (using `LoopVectorization`)
- **GQ**: `FastGaussQuadrature.gausslegendre`

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
| x^3 | [-1, 1] | 5 | 2.23 | 3.14 | 18.65 | 0.12 | 0.17 |
| x^3 | [-1, 1] | 50 | 21.86 | 22.09 | 35.69 | 0.61 | 0.62 |
| x^3 | [-1, 1] | 500 | 239.28 | 60.84 | 267.62 | 0.89 | 0.23 |
| x^3+x^2+x+1 | [-1, 1] | 5 | 3.98 | 4.00 | 20.12 | 0.20 | 0.20 |
| x^3+x^2+x+1 | [-1, 1] | 50 | 41.75 | 26.31 | 36.58 | 1.14 | 0.72 |
| x^3+x^2+x+1 | [-1, 1] | 500 | 450.22 | 80.03 | 263.54 | 1.71 | 0.30 |

## Analysis

- **Polynomials**: `FastTanhSinhQuadrature` with SIMD acceleration (`x^2`, `x^3`) is significantly faster (up to ~3-5x) than Gauss-Legendre quadrature.
- **Singularities**: Functions like `sqrt(1-x^2)` and `log(1-x)` are handled handled very efficiently, often matching or outperforming Gaussian quadrature due to the double exponential clustering of nodes.
- **Runge Function**: `1/(1+25x^2)` also shows competitive performance, especially with SIMD.

In summary, for smooth analytic functions, Gaussian quadrature (standard non-SIMD `integrate` vs GQ) is faster due to fewer nodes required for exactness using polynomials. However, `FastTanhSinhQuadrature`'s **SIMD implementation** often bridges or exceeds this gap, and it is the superior choice for **singular integrands**.
