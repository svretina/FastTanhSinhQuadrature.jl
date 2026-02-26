# Benchmark Timings

- Tolerances: `rtol=1.0e-6`, `atol=1.0e-8`
- Max evaluations (external adaptive solvers): `200000`
- Timing method: `@belapsed` with interpolation (`samples=3`, `evals=1`).

| Dim | Function | Method | Time (ns) | Abs Error | Rel Error | Status | Notes |
| :-- | :------- | :----- | --------: | --------: | --------: | :----- | :---- |
| 1D | 1/(1+25x^2) | CubatureJLh | 1.385e+04 | 9.469e-13 | 1.724e-12 | ok |  |
| 1D | 1/(1+25x^2) | CubatureJLp | 1.246e+04 | 1.110e-16 | 2.021e-16 | ok |  |
| 1D | 1/(1+25x^2) | FTS adaptive (typed) | 4058.0000 | 1.665e-14 | 3.031e-14 | ok |  |
| 1D | 1/(1+25x^2) | FTS avx (precomputed) | 217.0000 | 3.954e-07 | 7.197e-07 | ok | n=120 |
| 1D | 1/(1+25x^2) | FTS quad (convenience) | 4126.0000 | 1.665e-14 | 3.031e-14 | ok |  |
| 1D | 1/(1+25x^2) | FastGauss (precomputed) | 97.0000 | 9.279e-12 | 1.689e-11 | ok | n=64 |
| 1D | 1/(1+25x^2) | HCubature | 1.180e+04 | 1.894e-12 | 3.447e-12 | ok |  |
| 1D | 1/(1+25x^2) | QuadGK | 375.0000 | 1.894e-12 | 3.447e-12 | ok |  |
| 1D | 1/sqrt(1-x^2) | CubatureJLh | 1.634e+05 | 1.488e-07 | 4.737e-08 | ok |  |
| 1D | 1/sqrt(1-x^2) | CubatureJLp | 766.0000 | Inf | Inf | above_tol |  |
| 1D | 1/sqrt(1-x^2) | FTS adaptive (typed) | 1102.0000 | 2.916e-09 | 9.282e-10 | ok |  |
| 1D | 1/sqrt(1-x^2) | FTS avx (precomputed) | 339.0000 | 1.358e-07 | 4.321e-08 | ok | n=12 |
| 1D | 1/sqrt(1-x^2) | FTS quad (convenience) | 988.0000 | 2.916e-09 | 9.282e-10 | ok |  |
| 1D | 1/sqrt(1-x^2) | FastGauss (precomputed) | 9401.0000 | 0.0004 | 0.0001 | above_tol | n=4096, tol_not_met |
| 1D | 1/sqrt(1-x^2) | HCubature | 1.836e+05 | 1.971e-06 | 6.275e-07 | ok |  |
| 1D | 1/sqrt(1-x^2) | QuadGK | 6646.0000 | 1.971e-06 | 6.275e-07 | ok |  |
| 1D | log(1-x) | CubatureJLh | 5.497e+04 | 8.071e-10 | 1.315e-09 | ok |  |
| 1D | log(1-x) | CubatureJLp | 732.0000 | Inf | Inf | above_tol |  |
| 1D | log(1-x) | FTS adaptive (typed) | 1496.0000 | 2.220e-16 | 3.618e-16 | ok |  |
| 1D | log(1-x) | FTS avx (precomputed) | 247.0000 | 2.541e-07 | 4.141e-07 | ok | n=16 |
| 1D | log(1-x) | FTS quad (convenience) | 1513.0000 | 2.220e-16 | 3.618e-16 | ok |  |
| 1D | log(1-x) | FastGauss (precomputed) | 9460.0000 | 3.009e-07 | 4.903e-07 | ok | n=2048 |
| 1D | log(1-x) | HCubature | 4.422e+04 | 1.033e-07 | 1.683e-07 | ok |  |
| 1D | log(1-x) | QuadGK | 3399.0000 | 1.033e-07 | 1.683e-07 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | CubatureJLh | 1591.0000 | 0.0000 | 0.0000 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | CubatureJLp | 3107.0000 | 2.220e-16 | 1.727e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FTS adaptive (typed) | 2330.0000 | 4.441e-16 | 3.454e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FTS avx (precomputed) | 307.0000 | 4.158e-10 | 3.234e-10 | ok | n=32 |
| 1D | x^6 - 2x^3 + 0.5 | FTS quad (convenience) | 2221.0000 | 4.441e-16 | 3.454e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | FastGauss (precomputed) | 77.0000 | 2.220e-16 | 1.727e-16 | ok | n=4 |
| 1D | x^6 - 2x^3 + 0.5 | HCubature | 7359.0000 | 2.220e-16 | 1.727e-16 | ok |  |
| 1D | x^6 - 2x^3 + 0.5 | QuadGK | 322.0000 | 2.220e-16 | 1.727e-16 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaCuhre | 2.308e+07 | 7.874e-05 | 7.978e-06 | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaDivonne | 3.780e+07 | Inf | Inf | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubaVegas | 2.470e+07 | 0.0062 | 0.0006 | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubatureJLh | 2.752e+07 | 1.436e-06 | 1.455e-07 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | CubatureJLp | 1468.0000 | Inf | Inf | above_tol |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS adaptive (typed) | 3240.0000 | 1.832e-08 | 1.856e-09 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS avx (precomputed) | 578.0000 | 8.530e-07 | 8.643e-08 | ok | n=12 |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | FTS quad (convenience) | 3383.0000 | 1.832e-08 | 1.856e-09 | ok |  |
| 2D | 1/sqrt((1-x^2)(1-y^2)) | HCubature | 4.455e+07 | 1.436e-06 | 1.455e-07 | ok |  |
| 2D | exp(x+y) | CubaCuhre | 2.523e+04 | 2.753e-14 | 4.984e-15 | ok |  |
| 2D | exp(x+y) | CubaDivonne | 2.283e+07 | 6.290e-07 | 1.139e-07 | ok |  |
| 2D | exp(x+y) | CubaVegas | 2.408e+07 | 8.491e-05 | 1.537e-05 | above_tol |  |
| 2D | exp(x+y) | CubatureJLh | 3.199e+04 | 2.387e-10 | 4.321e-11 | ok |  |
| 2D | exp(x+y) | CubatureJLp | 3.125e+04 | 4.441e-15 | 8.039e-16 | ok |  |
| 2D | exp(x+y) | FTS adaptive (typed) | 1.602e+04 | 3.553e-15 | 6.431e-16 | ok |  |
| 2D | exp(x+y) | FTS avx (precomputed) | 933.0000 | 3.787e-09 | 6.855e-10 | ok | n=24 |
| 2D | exp(x+y) | FTS quad (convenience) | 1.712e+04 | 3.553e-15 | 6.431e-16 | ok |  |
| 2D | exp(x+y) | HCubature | 7.078e+04 | 2.387e-10 | 4.321e-11 | ok |  |
| 2D | x^2 + y^2 | CubaCuhre | 2.711e+04 | 8.882e-16 | 3.331e-16 | ok |  |
| 2D | x^2 + y^2 | CubaDivonne | 2.297e+07 | 3.738e-08 | 1.402e-08 | ok |  |
| 2D | x^2 + y^2 | CubaVegas | 2.351e+07 | 1.013e-05 | 3.799e-06 | above_tol |  |
| 2D | x^2 + y^2 | CubatureJLh | 1882.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 2D | x^2 + y^2 | CubatureJLp | 7069.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 2D | x^2 + y^2 | FTS adaptive (typed) | 3400.0000 | 1.332e-15 | 4.996e-16 | ok |  |
| 2D | x^2 + y^2 | FTS avx (precomputed) | 294.0000 | 3.212e-09 | 1.204e-09 | ok | n=24 |
| 2D | x^2 + y^2 | FTS quad (convenience) | 3621.0000 | 1.332e-15 | 4.996e-16 | ok |  |
| 2D | x^2 + y^2 | HCubature | 5857.0000 | 4.441e-16 | 1.665e-16 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaCuhre | 2.270e+07 | 30.4601 | 0.9824 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaDivonne | 3.176e+07 | Inf | Inf | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubaVegas | 2.452e+07 | 0.0283 | 0.0009 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubatureJLh | 2.071e+07 | 0.3997 | 0.0129 | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | CubatureJLp | 3250.0000 | Inf | Inf | above_tol |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS adaptive (typed) | 7.509e+04 | 8.634e-08 | 2.785e-09 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS avx (precomputed) | 748.0000 | 4.020e-06 | 1.296e-07 | ok | n=12 |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | FTS quad (convenience) | 7.803e+04 | 8.634e-08 | 2.785e-09 | ok |  |
| 3D | 1/sqrt((1-x^2)(1-y^2)(1-z^2)) | HCubature | 2.778e+07 | 0.3997 | 0.0129 | above_tol |  |
| 3D | exp(x+y+z) | CubaCuhre | 4.272e+04 | 8.225e-09 | 6.335e-10 | ok |  |
| 3D | exp(x+y+z) | CubaDivonne | 2.953e+07 | 0.0004 | 3.159e-05 | above_tol |  |
| 3D | exp(x+y+z) | CubaVegas | 2.682e+07 | 0.0004 | 3.300e-05 | above_tol |  |
| 3D | exp(x+y+z) | CubatureJLh | 2.271e+05 | 9.109e-08 | 7.015e-09 | ok |  |
| 3D | exp(x+y+z) | CubatureJLp | 6.014e+05 | 2.842e-14 | 2.189e-15 | ok |  |
| 3D | exp(x+y+z) | FTS adaptive (typed) | 8.888e+05 | 1.990e-13 | 1.532e-14 | ok |  |
| 3D | exp(x+y+z) | FTS avx (precomputed) | 1.106e+04 | 2.095e-06 | 1.613e-07 | ok | n=20 |
| 3D | exp(x+y+z) | FTS quad (convenience) | 8.946e+05 | 1.990e-13 | 1.532e-14 | ok |  |
| 3D | exp(x+y+z) | HCubature | 3.228e+05 | 9.109e-08 | 7.015e-09 | ok |  |
| 3D | x^2*y^2*z^2 | CubaCuhre | 9.269e+06 | 5.551e-17 | 1.874e-16 | ok |  |
| 3D | x^2*y^2*z^2 | CubaDivonne | 3.513e+07 | 4.789e-06 | 1.616e-05 | above_tol |  |
| 3D | x^2*y^2*z^2 | CubaVegas | 2.682e+07 | 1.889e-05 | 6.375e-05 | above_tol |  |
| 3D | x^2*y^2*z^2 | CubatureJLh | 7.470e+06 | 7.772e-16 | 2.623e-15 | ok |  |
| 3D | x^2*y^2*z^2 | CubatureJLp | 1.751e+04 | 1.665e-16 | 5.621e-16 | ok |  |
| 3D | x^2*y^2*z^2 | FTS adaptive (typed) | 5.562e+04 | 1.266e-14 | 4.272e-14 | ok |  |
| 3D | x^2*y^2*z^2 | FTS avx (precomputed) | 688.0000 | 1.812e-07 | 6.116e-07 | ok | n=20 |
| 3D | x^2*y^2*z^2 | FTS quad (convenience) | 5.563e+04 | 1.266e-14 | 4.272e-14 | ok |  |
| 3D | x^2*y^2*z^2 | HCubature | 9.873e+06 | 9.992e-16 | 3.372e-15 | ok |  |
