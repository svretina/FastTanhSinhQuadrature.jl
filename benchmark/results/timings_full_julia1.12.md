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
| 1D | 1/(1+25x^2) | CubatureJLh | 1.873e+04 | 9.469e-13 | 1.724e-12 | ok |  |
| 1D | 1/(1+25x^2) | CubatureJLp | 1.369e+04 | 1.110e-16 | 2.021e-16 | ok |  |
| 1D | 1/(1+25x^2) | FTS adaptive (typed) | 640.0000 | 1.665e-14 | 3.031e-14 | ok |  |
| 1D | 1/(1+25x^2) | FTS avx (precomputed) | 370.0000 | 1.688e-14 | 3.072e-14 | ok | N=256, adaptive_level=6 |
| 1D | 1/(1+25x^2) | FTS quad (convenience) | 477.0000 | 1.665e-14 | 3.031e-14 | ok |  |
| 1D | 1/(1+25x^2) | FastGauss (precomputed) | 114.0000 | 9.279e-12 | 1.689e-11 | ok | N=64, independently calibrated for tolerance |
| 1D | 1/(1+25x^2) | HCubature | 1.237e+04 | 1.894e-12 | 3.447e-12 | ok |  |
| 1D | 1/(1+25x^2) | QuadGK | 679.0000 | 1.894e-12 | 3.447e-12 | ok |  |
| 1D | 1/sqrt(1-x^2) | CubatureJLh | 1.759e+05 | 1.488e-07 | 4.737e-08 | ok |  |
| 1D | 1/sqrt(1-x^2) | CubatureJLp | 1697.0000 | Inf | Inf | above_tol |  |
| 1D | 1/sqrt(1-x^2) | FTS adaptive (typed) | 489.0000 | 2.916e-09 | 9.282e-10 | ok |  |
| 1D | 1/sqrt(1-x^2) | FTS avx (precomputed) | 366.0000 | 2.916e-09 | 9.282e-10 | ok | N=32, adaptive_level=3 |
| 1D | 1/sqrt(1-x^2) | FTS quad (convenience) | 333.0000 | 2.916e-09 | 9.282e-10 | ok |  |
| 1D | 1/sqrt(1-x^2) | FastGauss (precomputed) | 9849.0000 | 0.0004 | 0.0001 | above_tol | N=4096, tol_not_met |
| 1D | 1/sqrt(1-x^2) | HCubature | 1.716e+05 | 1.971e-06 | 6.275e-07 | ok |  |
| 1D | 1/sqrt(1-x^2) | QuadGK | 6846.0000 | 1.971e-06 | 6.275e-07 | ok |  |
| 1D | log(1-x) | CubatureJLh | 5.886e+04 | 8.071e-10 | 1.315e-09 | ok |  |
| 1D | log(1-x) | CubatureJLp | 1100.0000 | Inf | Inf | above_tol |  |
| 1D | log(1-x) | FTS adaptive (typed) | 570.0000 | 2.220e-16 | 3.618e-16 | ok |  |
| 1D | log(1-x) | FTS avx (precomputed) | 418.0000 | 2.220e-16 | 3.618e-16 | ok | N=32, adaptive_level=3 |
| 1D | log(1-x) | FTS quad (convenience) | 472.0000 | 2.220e-16 | 3.618e-16 | ok |  |
| 1D | log(1-x) | FastGauss (precomputed) | 1.491e+04 | 3.009e-07 | 4.903e-07 | ok | N=2048, independently calibrated for tolerance |
| 1D | log(1-x) | HCubature | 4.669e+04 | 1.033e-07 | 1.683e-07 | ok |  |
| 1D | log(1-x) | QuadGK | 4191.0000 | 1.033e-07 | 1.683e-07 | ok |  |
| 1D | sin^2(1000x) | CubatureJLh | 9.340e+04 | 1.510e-14 | 4.806e-15 | ok |  |
| 1D | sin^2(1000x) | CubatureJLp | 1156.0000 | 3.1416 | 1.0000 | above_tol |  |
| 1D | sin^2(1000x) | FTS adaptive (typed) | 2.557e+05 | 1.377e-14 | 4.382e-15 | ok |  |
| 1D | sin^2(1000x) | FTS avx (precomputed) | 3.082e+04 | 1.110e-14 | 3.534e-15 | ok | N=16384, adaptive_level=12, max_levels_reached |
| 1D | sin^2(1000x) | FTS quad (convenience) | 1.923e+05 | 1.377e-14 | 4.382e-15 | ok |  |
| 1D | sin^2(1000x) | FastGauss (precomputed) | 3.110e+04 | 3.553e-14 | 1.131e-14 | ok | N=4096, independently calibrated for tolerance |
| 1D | sin^2(1000x) | HCubature | 1.514e+05 | 4.086e-14 | 1.300e-14 | ok |  |
| 1D | sin^2(1000x) | QuadGK | 9978.0000 | 2.665e-15 | 8.481e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | CubatureJLh | 3046.0000 | 0.0000 | 0.0000 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | CubatureJLp | 3767.0000 | 2.220e-16 | 1.727e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FTS adaptive (typed) | 1231.0000 | 4.441e-16 | 3.454e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FTS avx (precomputed) | 307.0000 | 2.220e-16 | 1.727e-16 | ok | N=64, adaptive_level=4 |
| 1D | x^6 - 2x^3 + 0.5 | FTS quad (convenience) | 1163.0000 | 4.441e-16 | 3.454e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FastGauss (precomputed) | 92.0000 | 2.220e-16 | 1.727e-16 | ok | N=4, independently calibrated for tolerance |
| 1D | x^6 - 2x^3 + 0.5 | HCubature | 4747.0000 | 2.220e-16 | 1.727e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | QuadGK | 423.0000 | 2.220e-16 | 1.727e-16 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaCuhre | 8.516e+07 | 7.874e-05 | 7.978e-06 | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaDivonne | 1.291e+08 | Inf | Inf | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaVegas | 8.431e+07 | 0.0062 | 0.0006 | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubatureJLh | 2.682e+07 | 1.436e-06 | 1.455e-07 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubatureJLp | 1612.0000 | Inf | Inf | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS adaptive (typed) | 2607.0000 | 1.832e-08 | 1.856e-09 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS avx (precomputed) | 1435.0000 | 1.832e-08 | 1.856e-09 | ok | N=32, adaptive_level=3 |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS quad (convenience) | 2591.0000 | 1.832e-08 | 1.856e-09 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | HCubature | 5.005e+07 | 1.436e-06 | 1.455e-07 | ok |  |
| 2D | exp(x+y) | CubaCuhre | 5.521e+04 | 2.753e-14 | 4.984e-15 | ok |  |
| 2D | exp(x+y) | CubaDivonne | 7.766e+07 | 6.290e-07 | 1.139e-07 | ok |  |
| 2D | exp(x+y) | CubaVegas | 6.967e+07 | 8.491e-05 | 1.537e-05 | above_tol |  |
| 2D | exp(x+y) | CubatureJLh | 3.576e+04 | 2.387e-10 | 4.321e-11 | ok |  |
| 2D | exp(x+y) | CubatureJLp | 3.326e+04 | 4.441e-15 | 8.039e-16 | ok |  |
| 2D | exp(x+y) | FTS adaptive (typed) | 1.715e+04 | 3.553e-15 | 6.431e-16 | ok |  |
| 2D | exp(x+y) | FTS avx (precomputed) | 4151.0000 | 0.0000 | 0.0000 | ok | N=64, adaptive_level=4 |
| 2D | exp(x+y) | FTS quad (convenience) | 1.641e+04 | 3.553e-15 | 6.431e-16 | ok |  |
| 2D | exp(x+y) | HCubature | 1.035e+05 | 2.387e-10 | 4.321e-11 | ok |  |
| 2D | x^2 + y^2 | CubaCuhre | 7.291e+04 | 8.882e-16 | 3.331e-16 | ok |  |
| 2D | x^2 + y^2 | CubaDivonne | 7.028e+07 | 3.738e-08 | 1.402e-08 | ok |  |
| 2D | x^2 + y^2 | CubaVegas | 7.965e+07 | 1.013e-05 | 3.799e-06 | above_tol |  |
| 2D | x^2 + y^2 | CubatureJLh | 2163.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 2D | x^2 + y^2 | CubatureJLp | 4381.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 2D | x^2 + y^2 | FTS adaptive (typed) | 2413.0000 | 1.332e-15 | 4.996e-16 | ok |  |
| 2D | x^2 + y^2 | FTS avx (precomputed) | 602.0000 | 4.441e-16 | 1.665e-16 | ok | N=64, adaptive_level=4 |
| 2D | x^2 + y^2 | FTS quad (convenience) | 1720.0000 | 1.332e-15 | 4.996e-16 | ok |  |
| 2D | x^2 + y^2 | HCubature | 4714.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaCuhre | 1.064e+08 | 30.4601 | 0.9824 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaDivonne | 1.534e+08 | Inf | Inf | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaVegas | 9.790e+07 | 0.0283 | 0.0009 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubatureJLh | 2.635e+07 | 0.3997 | 0.0129 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubatureJLp | 3631.0000 | Inf | Inf | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS adaptive (typed) | 8.387e+04 | 8.634e-08 | 2.785e-09 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS avx (precomputed) | 2.138e+05 | 8.634e-08 | 2.785e-09 | ok | N=32, adaptive_level=3 |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS quad (convenience) | 9.529e+04 | 8.634e-08 | 2.785e-09 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | HCubature | 3.484e+07 | 0.3997 | 0.0129 | above_tol |  |
| 3D | exp(x+y+z) | CubaCuhre | 1.311e+05 | 8.225e-09 | 6.335e-10 | ok |  |
| 3D | exp(x+y+z) | CubaDivonne | 1.062e+08 | 0.0004 | 3.159e-05 | above_tol |  |
| 3D | exp(x+y+z) | CubaVegas | 1.260e+08 | 0.0004 | 3.300e-05 | above_tol |  |
| 3D | exp(x+y+z) | CubatureJLh | 2.447e+05 | 9.109e-08 | 7.015e-09 | ok |  |
| 3D | exp(x+y+z) | CubatureJLp | 6.615e+05 | 2.842e-14 | 2.189e-15 | ok |  |
| 3D | exp(x+y+z) | FTS adaptive (typed) | 8.721e+05 | 1.990e-13 | 1.532e-14 | ok |  |
| 3D | exp(x+y+z) | FTS avx (precomputed) | 2.773e+05 | 4.086e-14 | 3.147e-15 | ok | N=64, adaptive_level=4 |
| 3D | exp(x+y+z) | FTS quad (convenience) | 8.828e+05 | 1.990e-13 | 1.532e-14 | ok |  |
| 3D | exp(x+y+z) | HCubature | 5.108e+05 | 9.109e-08 | 7.015e-09 | ok |  |
| 3D | x^2*y^2*z^2 | CubaCuhre | 4.246e+07 | 5.551e-17 | 1.874e-16 | ok |  |
| 3D | x^2*y^2*z^2 | CubaDivonne | 1.433e+08 | 4.789e-06 | 1.616e-05 | above_tol |  |
| 3D | x^2*y^2*z^2 | CubaVegas | 9.618e+07 | 1.889e-05 | 6.375e-05 | above_tol |  |
| 3D | x^2*y^2*z^2 | CubatureJLh | 9.384e+06 | 7.772e-16 | 2.623e-15 | ok |  |
| 3D | x^2*y^2*z^2 | CubatureJLp | 1.831e+04 | 1.665e-16 | 5.621e-16 | ok |  |
| 3D | x^2*y^2*z^2 | FTS adaptive (typed) | 1.071e+05 | 1.266e-14 | 4.272e-14 | ok |  |
| 3D | x^2*y^2*z^2 | FTS avx (precomputed) | 1.999e+04 | 4.441e-16 | 1.499e-15 | ok | N=64, adaptive_level=4 |
| 3D | x^2*y^2*z^2 | FTS quad (convenience) | 7.264e+04 | 1.266e-14 | 4.272e-14 | ok |  |
| 3D | x^2*y^2*z^2 | HCubature | 1.033e+07 | 9.992e-16 | 3.372e-15 | ok |  |
