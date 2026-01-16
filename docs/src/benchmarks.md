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
| exp(x) | [-1, 1] | 9 | 54.24 | 18.98 | 46.01 | 1.18 | 0.41 |
| exp(x) | [-1, 1] | 100 | 656.00 | 151.80 | 310.92 | 2.11 | 0.49 |
| exp(x) | [-1, 1] | 1000 | 6941.80 | 1530.50 | 3075.75 | 2.26 | 0.50 |
| sin(x)^2 | [-1, 1] | 9 | 74.93 | 41.13 | 46.04 | 1.63 | 0.89 |
| sin(x)^2 | [-1, 1] | 100 | 795.19 | 287.41 | 353.35 | 2.25 | 0.81 |
| sin(x)^2 | [-1, 1] | 1000 | 8615.33 | 2749.33 | 3958.00 | 2.18 | 0.69 |
| 1/(1+25x^2) | [-1, 1] | 9 | 8.47 | 5.42 | 21.22 | 0.40 | 0.26 |
| 1/(1+25x^2) | [-1, 1] | 100 | 83.82 | 44.54 | 77.23 | 1.09 | 0.58 |
| 1/(1+25x^2) | [-1, 1] | 1000 | 1027.28 | 466.15 | 612.57 | 1.68 | 0.76 |
| sqrt(1-x^2) | [-1, 1] | 9 | 12.04 | 8.02 | 26.30 | 0.46 | 0.31 |
| sqrt(1-x^2) | [-1, 1] | 100 | 127.94 | 68.33 | 144.65 | 0.88 | 0.47 |
| sqrt(1-x^2) | [-1, 1] | 1000 | 1396.40 | 625.68 | 1528.60 | 0.91 | 0.41 |
| x^2 | [-1, 1] | 9 | 5.12 | 3.98 | 20.59 | 0.25 | 0.19 |
| x^2 | [-1, 1] | 100 | 34.88 | 9.74 | 69.44 | 0.50 | 0.14 |
| x^2 | [-1, 1] | 1000 | 442.34 | 117.87 | 578.60 | 0.76 | 0.20 |
| x^3 | [-1, 1] | 9 | 6.26 | 3.96 | 23.39 | 0.27 | 0.17 |
| x^3 | [-1, 1] | 100 | 38.67 | 11.80 | 65.03 | 0.59 | 0.18 |
| x^3 | [-1, 1] | 1000 | 418.34 | 106.36 | 540.83 | 0.77 | 0.20 |
| x^3+x^2+x+1 | [-1, 1] | 9 | 8.66 | 4.62 | 23.41 | 0.37 | 0.20 |
| x^3+x^2+x+1 | [-1, 1] | 100 | 89.29 | 17.90 | 71.05 | 1.26 | 0.25 |
| x^3+x^2+x+1 | [-1, 1] | 1000 | 915.95 | 180.46 | 560.56 | 1.63 | 0.32 |
| log(1-x) | [-1, 1] | 9 | 75.40 | 53.18 | 56.30 | 1.34 | 0.94 |
| log(1-x) | [-1, 1] | 100 | 752.70 | 327.72 | 458.73 | 1.64 | 0.71 |
| log(1-x) | [-1, 1] | 1000 | 7470.00 | 3503.38 | 3730.50 | 2.00 | 0.94 |

## Analysis

- **Polynomials**: `FastTanhSinhQuadrature` with SIMD acceleration (`x^2`, `x^3`) is significantly faster (up to ~3-5x) than Gauss-Legendre quadrature.
- **Singularities**: Functions like `sqrt(1-x^2)` and `log(1-x)` are handled handled very efficiently, often matching or outperforming Gaussian quadrature due to the double exponential clustering of nodes.
- **Runge Function**: `1/(1+25x^2)` also shows competitive performance, especially with SIMD.

In summary, for smooth analytic functions, Gaussian quadrature (standard non-SIMD `integrate` vs GQ) is faster due to fewer nodes required for exactness using polynomials. However, `FastTanhSinhQuadrature`'s **SIMD implementation** often bridges or exceeds this gap, and it is the superior choice for **singular integrands**.
