# Experimental 3D Fixed-Grid Results

Date: 2026-04-06

Source harness: [experimental_performance3D_stride.jl](/home/svretina/Codes/FastTanhSinhQuadrature.jl/test/experimental_performance3D_stride.jl)

## Goal

Compare fixed-grid 3D AVX kernels while keeping:

- zero hot-path allocations
- the analytic `exp(x+y+z)` integral correct
- the current package `integrate3D_avx` as the baseline

The analytic target was:

```math
\int_{[-1,1]^3} e^{x+y+z} \, dx\,dy\,dz = (e - e^{-1})^3
```

## Variants

- `package_avx`: current package `integrate3D_avx`
- `vector`: `Vector` nodes, weights, and scratch
- `ptr_vector`: `PtrArray(Vector)` nodes, weights, and scratch
- `stride`: `StrideArray` nodes, weights, and scratch
- `ptr_stride`: `PtrArray(StrideArray)` nodes, weights, and scratch
- `svector`: `SVector` nodes and weights with `Vector` scratch

## Benchmark Snapshot

All reported runs were `0` allocations and `0` bytes.

### N = 32

| Variant | Time | Vs `package_avx` | Error |
|---|---:|---:|---:|
| `package_avx` | 39.792 us | 1.000x | 2.1174173525650986e-12 |
| `vector` | 44.038 us | 0.904x | 2.1032064978498966e-12 |
| `ptr_vector` | 47.110 us | 0.845x | 2.1032064978498966e-12 |
| `stride` | 44.020 us | 0.904x | 2.1032064978498966e-12 |
| `ptr_stride` | 42.730 us | 0.931x | 2.1032064978498966e-12 |
| `svector` | 47.206 us | 0.843x | 2.1032064978498966e-12 |

Winner: `package_avx`

### N = 64

| Variant | Time | Vs `package_avx` | Error |
|---|---:|---:|---:|
| `package_avx` | 901.727 us | 1.000x | 3.552713678800501e-14 |
| `vector` | 438.720 us | 2.055x | 8.881784197001252e-15 |
| `ptr_vector` | 1181.302 us | 0.763x | 8.881784197001252e-15 |
| `stride` | 482.193 us | 1.870x | 8.881784197001252e-15 |
| `ptr_stride` | 369.937 us | 2.438x | 8.881784197001252e-15 |
| `svector` | 369.386 us | 2.441x | 8.881784197001252e-15 |

Winner: `svector`, with `ptr_stride` effectively tied

### N = 128

| Variant | Time | Vs `package_avx` | Error |
|---|---:|---:|---:|
| `package_avx` | 3298.539 us | 1.000x | 3.1796787425264483e-13 |
| `vector` | 3203.759 us | 1.030x | 3.907985046680551e-14 |
| `ptr_vector` | 3006.606 us | 1.097x | 3.907985046680551e-14 |
| `stride` | 2578.314 us | 1.279x | 3.907985046680551e-14 |
| `ptr_stride` | 3035.740 us | 1.087x | 3.907985046680551e-14 |
| `svector` | 2788.062 us | 1.183x | 3.907985046680551e-14 |

Winner: `stride`

### N = 256

| Variant | Time | Vs `package_avx` | Error |
|---|---:|---:|---:|
| `package_avx` | 29119.789 us | 1.000x | 3.3075764349632664e-12 |
| `vector` | 28525.967 us | 1.021x | 1.0462741784067475e-12 |
| `ptr_vector` | 27208.080 us | 1.070x | 1.0462741784067475e-12 |
| `stride` | 28627.978 us | 1.017x | 1.0462741784067475e-12 |
| `ptr_stride` | 29825.423 us | 0.976x | 1.0462741784067475e-12 |
| `svector` | skipped | skipped | skipped |

Winner: `ptr_vector`

## Takeaways

- Fixed-grid 3D still shows real wins from alternate storage/access patterns.
- `StrideArray` is strongest around medium node counts in this run.
- `PtrArray` is not universally good or bad. It helps in some cases, especially `PtrArray(StrideArray)` at `N=64`, but it also loses badly in others.
- `SVector` is excellent for smaller node counts, but that result does not transfer directly to adaptive refinement because adaptive levels have heterogeneous sizes.
- All tested variants preserved the analytic answer and stayed at `0` hot-path allocations.

## Interpretation

For fixed-grid kernels, there is still enough signal to justify experimenting with `StrideArray` or `PtrArray` integration into the package. The effect is size-dependent, so any adoption should be benchmark-driven rather than assumed.
