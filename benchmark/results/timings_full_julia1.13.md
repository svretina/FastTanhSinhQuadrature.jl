# Benchmark Timings at Fixed Accuracy Target

- Tolerances: `rtol=1.0e-6`, `atol=1.0e-8`
- Max evaluations (external adaptive solvers): `200000`
- Timing method: `@belapsed` with interpolation (`samples=3`, `evals=1`).
- One warm call is executed before each timed benchmark (reduces first-call compilation effects).
- Adaptive cache construction is excluded from timed regions (caches are prebuilt).
- `FTS avx` is run at the `N` inferred from adaptive FTS convergence.
- `FastGauss` (1D) is independently calibrated to the tolerance target.

| Dim | Function | Method | Time (ns) | Abs Error | Rel Error | Status | Notes |
| :-- | :------- | :----- | --------: | --------: | --------: | :----- | :---- |
| 1D | 1/(1+25x^2) | CubatureJLh | 1.150e+04 | 9.469e-13 | 1.724e-12 | ok |  |
| 1D | 1/(1+25x^2) | CubatureJLp | 1.393e+04 | 1.110e-16 | 2.021e-16 | ok |  |
| 1D | 1/(1+25x^2) | FTS adaptive (typed) | 863.0000 | 1.665e-14 | 3.031e-14 | ok |  |
| 1D | 1/(1+25x^2) | FTS avx (precomputed) | 442.0000 | 1.688e-14 | 3.072e-14 | ok | N=256, adaptive_level=6 |
| 1D | 1/(1+25x^2) | FTS quad (convenience) | 545.0000 | 1.665e-14 | 3.031e-14 | ok |  |
| 1D | 1/(1+25x^2) | FastGauss (precomputed) | 185.0000 | 9.279e-12 | 1.689e-11 | ok | N=64, independently calibrated for tolerance |
| 1D | 1/(1+25x^2) | HCubature | 1.264e+04 | 1.894e-12 | 3.447e-12 | ok |  |
| 1D | 1/(1+25x^2) | QuadGK | 458.0000 | 1.894e-12 | 3.447e-12 | ok |  |
| 1D | 1/sqrt(1-x^2) | CubatureJLh | 2.770e+05 | 1.488e-07 | 4.737e-08 | ok |  |
| 1D | 1/sqrt(1-x^2) | CubatureJLp | 955.0000 | Inf | Inf | above_tol |  |
| 1D | 1/sqrt(1-x^2) | FTS adaptive (typed) | 713.0000 | 2.916e-09 | 9.282e-10 | ok |  |
| 1D | 1/sqrt(1-x^2) | FTS avx (precomputed) | 324.0000 | 2.916e-09 | 9.282e-10 | ok | N=32, adaptive_level=3 |
| 1D | 1/sqrt(1-x^2) | FTS quad (convenience) | 525.0000 | 2.916e-09 | 9.282e-10 | ok |  |
| 1D | 1/sqrt(1-x^2) | FastGauss (precomputed) | 9829.0000 | 0.0004 | 0.0001 | above_tol | N=4096, tol_not_met |
| 1D | 1/sqrt(1-x^2) | HCubature | 2.265e+05 | 1.971e-06 | 6.275e-07 | ok |  |
| 1D | 1/sqrt(1-x^2) | QuadGK | 9372.0000 | 1.971e-06 | 6.275e-07 | ok |  |
| 1D | log(1-x) | CubatureJLh | 6.426e+04 | 8.071e-10 | 1.315e-09 | ok |  |
| 1D | log(1-x) | CubatureJLp | 1157.0000 | Inf | Inf | above_tol |  |
| 1D | log(1-x) | FTS adaptive (typed) | 573.0000 | 2.220e-16 | 3.618e-16 | ok |  |
| 1D | log(1-x) | FTS avx (precomputed) | 353.0000 | 2.220e-16 | 3.618e-16 | ok | N=32, adaptive_level=3 |
| 1D | log(1-x) | FTS quad (convenience) | 507.0000 | 2.220e-16 | 3.618e-16 | ok |  |
| 1D | log(1-x) | FastGauss (precomputed) | 1.266e+04 | 3.009e-07 | 4.903e-07 | ok | N=2048, independently calibrated for tolerance |
| 1D | log(1-x) | HCubature | 5.281e+04 | 1.033e-07 | 1.683e-07 | ok |  |
| 1D | log(1-x) | QuadGK | 2996.0000 | 1.033e-07 | 1.683e-07 | ok |  |
| 1D | sin^2(1000x) | CubatureJLh | 9.820e+04 | 1.510e-14 | 4.806e-15 | ok |  |
| 1D | sin^2(1000x) | CubatureJLp | 1209.0000 | 3.1416 | 1.0000 | above_tol |  |
| 1D | sin^2(1000x) | FTS adaptive (typed) | 2.249e+05 | 1.377e-14 | 4.382e-15 | ok |  |
| 1D | sin^2(1000x) | FTS avx (precomputed) | 2.787e+04 | 1.110e-14 | 3.534e-15 | ok | N=16384, adaptive_level=12, max_levels_reached |
| 1D | sin^2(1000x) | FTS quad (convenience) | 1.281e+05 | 1.377e-14 | 4.382e-15 | ok |  |
| 1D | sin^2(1000x) | FastGauss (precomputed) | 3.266e+04 | 3.553e-14 | 1.131e-14 | ok | N=4096, independently calibrated for tolerance |
| 1D | sin^2(1000x) | HCubature | 9.984e+04 | 4.086e-14 | 1.300e-14 | ok |  |
| 1D | sin^2(1000x) | QuadGK | 1.555e+04 | 2.665e-15 | 8.481e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | CubatureJLh | 2103.0000 | 0.0000 | 0.0000 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | CubatureJLp | 2.013e+04 | 2.220e-16 | 1.727e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FTS adaptive (typed) | 1560.0000 | 4.441e-16 | 3.454e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FTS avx (precomputed) | 270.0000 | 2.220e-16 | 1.727e-16 | ok | N=64, adaptive_level=4 |
| 1D | x^6 - 2x^3 + 0.5 | FTS quad (convenience) | 652.0000 | 4.441e-16 | 3.454e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FastGauss (precomputed) | 94.0000 | 2.220e-16 | 1.727e-16 | ok | N=4, independently calibrated for tolerance |
| 1D | x^6 - 2x^3 + 0.5 | HCubature | 4396.0000 | 2.220e-16 | 1.727e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | QuadGK | 425.0000 | 2.220e-16 | 1.727e-16 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaCuhre | 6.673e+07 | 7.874e-05 | 7.978e-06 | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaDivonne | 1.198e+08 | Inf | Inf | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaVegas | 7.077e+07 | 0.0062 | 0.0006 | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubatureJLh | 2.312e+07 | 1.436e-06 | 1.455e-07 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubatureJLp | 1610.0000 | Inf | Inf | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS adaptive (typed) | 2572.0000 | 1.832e-08 | 1.856e-09 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS avx (precomputed) | 1695.0000 | 1.832e-08 | 1.856e-09 | ok | N=32, adaptive_level=3 |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS quad (convenience) | 2467.0000 | 1.832e-08 | 1.856e-09 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | HCubature | 4.709e+07 | 1.436e-06 | 1.455e-07 | ok |  |
| 2D | exp(x+y) | CubaCuhre | 7.614e+04 | 2.753e-14 | 4.984e-15 | ok |  |
| 2D | exp(x+y) | CubaDivonne | 7.041e+07 | 6.290e-07 | 1.139e-07 | ok |  |
| 2D | exp(x+y) | CubaVegas | 6.760e+07 | 8.491e-05 | 1.537e-05 | above_tol |  |
| 2D | exp(x+y) | CubatureJLh | 5.376e+04 | 2.387e-10 | 4.321e-11 | ok |  |
| 2D | exp(x+y) | CubatureJLp | 3.300e+04 | 4.441e-15 | 8.039e-16 | ok |  |
| 2D | exp(x+y) | FTS adaptive (typed) | 1.931e+04 | 3.553e-15 | 6.431e-16 | ok |  |
| 2D | exp(x+y) | FTS avx (precomputed) | 4365.0000 | 0.0000 | 0.0000 | ok | N=64, adaptive_level=4 |
| 2D | exp(x+y) | FTS quad (convenience) | 1.483e+04 | 3.553e-15 | 6.431e-16 | ok |  |
| 2D | exp(x+y) | HCubature | 6.944e+04 | 2.387e-10 | 4.321e-11 | ok |  |
| 2D | x^2 + y^2 | CubaCuhre | 7.621e+04 | 8.882e-16 | 3.331e-16 | ok |  |
| 2D | x^2 + y^2 | CubaDivonne | 7.821e+07 | 3.738e-08 | 1.402e-08 | ok |  |
| 2D | x^2 + y^2 | CubaVegas | 9.022e+07 | 1.013e-05 | 3.799e-06 | above_tol |  |
| 2D | x^2 + y^2 | CubatureJLh | 2720.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 2D | x^2 + y^2 | CubatureJLp | 5182.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 2D | x^2 + y^2 | FTS adaptive (typed) | 3285.0000 | 1.332e-15 | 4.996e-16 | ok |  |
| 2D | x^2 + y^2 | FTS avx (precomputed) | 659.0000 | 4.441e-16 | 1.665e-16 | ok | N=64, adaptive_level=4 |
| 2D | x^2 + y^2 | FTS quad (convenience) | 1576.0000 | 1.332e-15 | 4.996e-16 | ok |  |
| 2D | x^2 + y^2 | HCubature | 6847.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaCuhre | 9.282e+07 | 30.4601 | 0.9824 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaDivonne | 1.439e+08 | Inf | Inf | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaVegas | 1.314e+08 | 0.0283 | 0.0009 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubatureJLh | 2.623e+07 | 0.3997 | 0.0129 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubatureJLp | 4299.0000 | Inf | Inf | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS adaptive (typed) | 7.764e+04 | 8.634e-08 | 2.785e-09 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS avx (precomputed) | 5.692e+04 | 8.634e-08 | 2.785e-09 | ok | N=32, adaptive_level=3 |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS quad (convenience) | 8.144e+04 | 8.634e-08 | 2.785e-09 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | HCubature | 2.997e+07 | 0.3997 | 0.0129 | above_tol |  |
| 3D | exp(x+y+z) | CubaCuhre | 2.063e+05 | 8.225e-09 | 6.335e-10 | ok |  |
| 3D | exp(x+y+z) | CubaDivonne | 1.288e+08 | 0.0004 | 3.159e-05 | above_tol |  |
| 3D | exp(x+y+z) | CubaVegas | 9.754e+07 | 0.0004 | 3.300e-05 | above_tol |  |
| 3D | exp(x+y+z) | CubatureJLh | 2.461e+05 | 9.109e-08 | 7.015e-09 | ok |  |
| 3D | exp(x+y+z) | CubatureJLp | 6.173e+05 | 2.842e-14 | 2.189e-15 | ok |  |
| 3D | exp(x+y+z) | FTS adaptive (typed) | 1.353e+06 | 1.990e-13 | 1.532e-14 | ok |  |
| 3D | exp(x+y+z) | FTS avx (precomputed) | 2.902e+05 | 4.086e-14 | 3.147e-15 | ok | N=64, adaptive_level=4 |
| 3D | exp(x+y+z) | FTS quad (convenience) | 9.363e+05 | 1.990e-13 | 1.532e-14 | ok |  |
| 3D | exp(x+y+z) | HCubature | 3.181e+05 | 9.109e-08 | 7.015e-09 | ok |  |
| 3D | x^2*y^2*z^2 | CubaCuhre | 3.220e+07 | 5.551e-17 | 1.874e-16 | ok |  |
| 3D | x^2*y^2*z^2 | CubaDivonne | 1.038e+08 | 4.789e-06 | 1.616e-05 | above_tol |  |
| 3D | x^2*y^2*z^2 | CubaVegas | 9.335e+07 | 1.889e-05 | 6.375e-05 | above_tol |  |
| 3D | x^2*y^2*z^2 | CubatureJLh | 7.136e+06 | 7.772e-16 | 2.623e-15 | ok |  |
| 3D | x^2*y^2*z^2 | CubatureJLp | 1.827e+04 | 1.665e-16 | 5.621e-16 | ok |  |
| 3D | x^2*y^2*z^2 | FTS adaptive (typed) | 5.486e+04 | 1.266e-14 | 4.272e-14 | ok |  |
| 3D | x^2*y^2*z^2 | FTS avx (precomputed) | 1.086e+04 | 4.441e-16 | 1.499e-15 | ok | N=64, adaptive_level=4 |
| 3D | x^2*y^2*z^2 | FTS quad (convenience) | 6.344e+04 | 1.266e-14 | 4.272e-14 | ok |  |
| 3D | x^2*y^2*z^2 | HCubature | 9.470e+06 | 9.992e-16 | 3.372e-15 | ok |  |
