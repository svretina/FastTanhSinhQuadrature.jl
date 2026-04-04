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
| 1D | 1/(1+25x^2) | CubatureJLh | 1.644e+04 | 9.469e-13 | 1.724e-12 | ok |  |
| 1D | 1/(1+25x^2) | CubatureJLp | 2.007e+04 | 1.110e-16 | 2.021e-16 | ok |  |
| 1D | 1/(1+25x^2) | FTS adaptive (typed) | 735.0000 | 1.665e-14 | 3.031e-14 | ok |  |
| 1D | 1/(1+25x^2) | FTS avx (precomputed) | 474.0000 | 1.688e-14 | 3.072e-14 | ok | N=256, adaptive_level=6 |
| 1D | 1/(1+25x^2) | FTS quad (convenience) | 628.0000 | 1.665e-14 | 3.031e-14 | ok |  |
| 1D | 1/(1+25x^2) | FastGauss (precomputed) | 176.0000 | 9.279e-12 | 1.689e-11 | ok | N=64, independently calibrated for tolerance |
| 1D | 1/(1+25x^2) | HCubature | 1.620e+04 | 1.894e-12 | 3.447e-12 | ok |  |
| 1D | 1/(1+25x^2) | QuadGK | 1612.0000 | 1.894e-12 | 3.447e-12 | ok |  |
| 1D | 1/sqrt(1-x^2) | CubatureJLh | 2.389e+05 | 1.488e-07 | 4.737e-08 | ok |  |
| 1D | 1/sqrt(1-x^2) | CubatureJLp | 2834.0000 | Inf | Inf | above_tol |  |
| 1D | 1/sqrt(1-x^2) | FTS adaptive (typed) | 574.0000 | 2.916e-09 | 9.282e-10 | ok |  |
| 1D | 1/sqrt(1-x^2) | FTS avx (precomputed) | 523.0000 | 2.916e-09 | 9.282e-10 | ok | N=32, adaptive_level=3 |
| 1D | 1/sqrt(1-x^2) | FTS quad (convenience) | 528.0000 | 2.916e-09 | 9.282e-10 | ok |  |
| 1D | 1/sqrt(1-x^2) | FastGauss (precomputed) | 1.146e+04 | 0.0004 | 0.0001 | above_tol | N=4096, tol_not_met |
| 1D | 1/sqrt(1-x^2) | HCubature | 2.189e+05 | 1.971e-06 | 6.275e-07 | ok |  |
| 1D | 1/sqrt(1-x^2) | QuadGK | 1.154e+04 | 1.971e-06 | 6.275e-07 | ok |  |
| 1D | log(1-x) | CubatureJLh | 7.402e+04 | 8.071e-10 | 1.315e-09 | ok |  |
| 1D | log(1-x) | CubatureJLp | 1356.0000 | Inf | Inf | above_tol |  |
| 1D | log(1-x) | FTS adaptive (typed) | 1013.0000 | 2.220e-16 | 3.618e-16 | ok |  |
| 1D | log(1-x) | FTS avx (precomputed) | 442.0000 | 2.220e-16 | 3.618e-16 | ok | N=32, adaptive_level=3 |
| 1D | log(1-x) | FTS quad (convenience) | 592.0000 | 2.220e-16 | 3.618e-16 | ok |  |
| 1D | log(1-x) | FastGauss (precomputed) | 1.208e+04 | 3.009e-07 | 4.903e-07 | ok | N=2048, independently calibrated for tolerance |
| 1D | log(1-x) | HCubature | 6.031e+04 | 1.033e-07 | 1.683e-07 | ok |  |
| 1D | log(1-x) | QuadGK | 5237.0000 | 1.033e-07 | 1.683e-07 | ok |  |
| 1D | sin^2(1000x) | CubatureJLh | 1.041e+05 | 1.510e-14 | 4.806e-15 | ok |  |
| 1D | sin^2(1000x) | CubatureJLp | 1301.0000 | 3.1416 | 1.0000 | above_tol |  |
| 1D | sin^2(1000x) | FTS adaptive (typed) | 1.770e+05 | 1.377e-14 | 4.382e-15 | ok |  |
| 1D | sin^2(1000x) | FTS avx (precomputed) | 2.909e+04 | 1.110e-14 | 3.534e-15 | ok | N=16384, adaptive_level=12, max_levels_reached |
| 1D | sin^2(1000x) | FTS quad (convenience) | 2.542e+05 | 1.377e-14 | 4.382e-15 | ok |  |
| 1D | sin^2(1000x) | FastGauss (precomputed) | 3.748e+04 | 3.553e-14 | 1.131e-14 | ok | N=4096, independently calibrated for tolerance |
| 1D | sin^2(1000x) | HCubature | 1.313e+05 | 4.086e-14 | 1.300e-14 | ok |  |
| 1D | sin^2(1000x) | QuadGK | 1.142e+04 | 2.665e-15 | 8.481e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | CubatureJLh | 2249.0000 | 0.0000 | 0.0000 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | CubatureJLp | 8949.0000 | 2.220e-16 | 1.727e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FTS adaptive (typed) | 1306.0000 | 4.441e-16 | 3.454e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FTS avx (precomputed) | 335.0000 | 2.220e-16 | 1.727e-16 | ok | N=64, adaptive_level=4 |
| 1D | x^6 - 2x^3 + 0.5 | FTS quad (convenience) | 803.0000 | 4.441e-16 | 3.454e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FastGauss (precomputed) | 118.0000 | 2.220e-16 | 1.727e-16 | ok | N=4, independently calibrated for tolerance |
| 1D | x^6 - 2x^3 + 0.5 | HCubature | 4813.0000 | 2.220e-16 | 1.727e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | QuadGK | 1898.0000 | 2.220e-16 | 1.727e-16 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaCuhre | 9.744e+07 | 7.874e-05 | 7.978e-06 | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaDivonne | 1.140e+08 | Inf | Inf | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaVegas | 9.357e+07 | 0.0062 | 0.0006 | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubatureJLh | 2.832e+07 | 1.436e-06 | 1.455e-07 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubatureJLp | 2433.0000 | Inf | Inf | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS adaptive (typed) | 3663.0000 | 1.832e-08 | 1.856e-09 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS avx (precomputed) | 1865.0000 | 1.832e-08 | 1.856e-09 | ok | N=32, adaptive_level=3 |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS quad (convenience) | 3660.0000 | 1.832e-08 | 1.856e-09 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | HCubature | 5.218e+07 | 1.436e-06 | 1.455e-07 | ok |  |
| 2D | exp(x+y) | CubaCuhre | 6.697e+04 | 2.753e-14 | 4.984e-15 | ok |  |
| 2D | exp(x+y) | CubaDivonne | 8.362e+07 | 6.290e-07 | 1.139e-07 | ok |  |
| 2D | exp(x+y) | CubaVegas | 9.402e+07 | 8.491e-05 | 1.537e-05 | above_tol |  |
| 2D | exp(x+y) | CubatureJLh | 4.420e+04 | 2.387e-10 | 4.321e-11 | ok |  |
| 2D | exp(x+y) | CubatureJLp | 3.861e+04 | 4.441e-15 | 8.039e-16 | ok |  |
| 2D | exp(x+y) | FTS adaptive (typed) | 1.864e+04 | 3.553e-15 | 6.431e-16 | ok |  |
| 2D | exp(x+y) | FTS avx (precomputed) | 5720.0000 | 0.0000 | 0.0000 | ok | N=64, adaptive_level=4 |
| 2D | exp(x+y) | FTS quad (convenience) | 2.438e+04 | 3.553e-15 | 6.431e-16 | ok |  |
| 2D | exp(x+y) | HCubature | 8.769e+04 | 2.387e-10 | 4.321e-11 | ok |  |
| 2D | x^2 + y^2 | CubaCuhre | 6.858e+04 | 8.882e-16 | 3.331e-16 | ok |  |
| 2D | x^2 + y^2 | CubaDivonne | 7.713e+07 | 3.738e-08 | 1.402e-08 | ok |  |
| 2D | x^2 + y^2 | CubaVegas | 8.497e+07 | 1.013e-05 | 3.799e-06 | above_tol |  |
| 2D | x^2 + y^2 | CubatureJLh | 5980.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 2D | x^2 + y^2 | CubatureJLp | 9098.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 2D | x^2 + y^2 | FTS adaptive (typed) | 2678.0000 | 1.332e-15 | 4.996e-16 | ok |  |
| 2D | x^2 + y^2 | FTS avx (precomputed) | 1968.0000 | 4.441e-16 | 1.665e-16 | ok | N=64, adaptive_level=4 |
| 2D | x^2 + y^2 | FTS quad (convenience) | 2153.0000 | 1.332e-15 | 4.996e-16 | ok |  |
| 2D | x^2 + y^2 | HCubature | 9986.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaCuhre | 1.111e+08 | 30.4601 | 0.9824 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaDivonne | 1.391e+08 | Inf | Inf | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaVegas | 1.320e+08 | 0.0283 | 0.0009 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubatureJLh | 3.151e+07 | 0.3997 | 0.0129 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubatureJLp | 5400.0000 | Inf | Inf | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS adaptive (typed) | 1.168e+05 | 8.634e-08 | 2.785e-09 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS avx (precomputed) | 7.367e+04 | 8.634e-08 | 2.785e-09 | ok | N=32, adaptive_level=3 |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS quad (convenience) | 1.153e+05 | 8.634e-08 | 2.785e-09 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | HCubature | 3.339e+07 | 0.3997 | 0.0129 | above_tol |  |
| 3D | exp(x+y+z) | CubaCuhre | 1.795e+05 | 8.225e-09 | 6.335e-10 | ok |  |
| 3D | exp(x+y+z) | CubaDivonne | 1.283e+08 | 0.0004 | 3.159e-05 | above_tol |  |
| 3D | exp(x+y+z) | CubaVegas | 1.178e+08 | 0.0004 | 3.300e-05 | above_tol |  |
| 3D | exp(x+y+z) | CubatureJLh | 3.170e+05 | 9.109e-08 | 7.015e-09 | ok |  |
| 3D | exp(x+y+z) | CubatureJLp | 7.207e+05 | 2.842e-14 | 2.189e-15 | ok |  |
| 3D | exp(x+y+z) | FTS adaptive (typed) | 1.049e+06 | 1.990e-13 | 1.532e-14 | ok |  |
| 3D | exp(x+y+z) | FTS avx (precomputed) | 3.401e+05 | 4.086e-14 | 3.147e-15 | ok | N=64, adaptive_level=4 |
| 3D | exp(x+y+z) | FTS quad (convenience) | 1.069e+06 | 1.990e-13 | 1.532e-14 | ok |  |
| 3D | exp(x+y+z) | HCubature | 3.619e+05 | 9.109e-08 | 7.015e-09 | ok |  |
| 3D | x^2*y^2*z^2 | CubaCuhre | 4.454e+07 | 5.551e-17 | 1.874e-16 | ok |  |
| 3D | x^2*y^2*z^2 | CubaDivonne | 1.318e+08 | 4.789e-06 | 1.616e-05 | above_tol |  |
| 3D | x^2*y^2*z^2 | CubaVegas | 1.055e+08 | 1.889e-05 | 6.375e-05 | above_tol |  |
| 3D | x^2*y^2*z^2 | CubatureJLh | 8.682e+06 | 7.772e-16 | 2.623e-15 | ok |  |
| 3D | x^2*y^2*z^2 | CubatureJLp | 2.441e+04 | 1.665e-16 | 5.621e-16 | ok |  |
| 3D | x^2*y^2*z^2 | FTS adaptive (typed) | 7.994e+04 | 1.266e-14 | 4.272e-14 | ok |  |
| 3D | x^2*y^2*z^2 | FTS avx (precomputed) | 2.207e+04 | 4.441e-16 | 1.499e-15 | ok | N=64, adaptive_level=4 |
| 3D | x^2*y^2*z^2 | FTS quad (convenience) | 7.634e+04 | 1.266e-14 | 4.272e-14 | ok |  |
| 3D | x^2*y^2*z^2 | HCubature | 1.345e+07 | 9.992e-16 | 3.372e-15 | ok |  |
