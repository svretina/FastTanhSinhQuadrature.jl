# Experimental 3D Adaptive Results

Date: 2026-04-06

Source harness: [experimental_adaptive_3d.jl](/home/svretina/Codes/FastTanhSinhQuadrature.jl/test/experimental_adaptive_3d.jl)

## Goal

Benchmark the adaptive 3D AVX path with a fairer baseline:

- keep the package-style AVX algorithm fixed
- keep hoisted coordinate scratch and `@turbo` loops
- vary only storage and access style

The analytic target was:

```math
\int_{[-1,1]^3} e^{x+y+z} \, dx\,dy\,dz = (e - e^{-1})^3
```

## Variants

- `package_scalar`: current scalar adaptive 3D package path
- `package_avx`: current package `adaptive_integrate_3D_avx`
- `vector`: package-style experimental kernel using `Vector` for nodes, weights, and scratch
- `ptr_vector`: package-style experimental kernel using `PtrArray(Vector)` for nodes, weights, and scratch
- `stride`: package-style experimental kernel using `StrideArray` for nodes, weights, and scratch
- `ptr_stride`: package-style experimental kernel using `PtrArray(StrideArray)` for nodes, weights, and scratch

## Important Notes

- This harness mirrors the package AVX algorithm and uses the same 3-set partitioning, hoisted coordinates, and `@turbo` structure.
- All reported runs were `0` allocations and `0` bytes.
- All variants matched the same analytic result for a given `max_levels`.
- This is a much fairer container comparison than the earlier unit-box-specialized adaptive experiment.

## Benchmark Snapshot

### max_levels = 3

| Variant | Time | Vs `package_avx` | Error |
|---|---:|---:|---:|
| `package_scalar` | 150.066 us | 0.592x | 2.1511681325137033e-12 |
| `package_avx` | 88.916 us | 1.000x | 2.1032064978498966e-12 |
| `vector` | 52.517 us | 1.693x | 2.1032064978498966e-12 |
| `ptr_vector` | 119.178 us | 0.746x | 2.1032064978498966e-12 |
| `stride` | 53.712 us | 1.655x | 2.1032064978498966e-12 |
| `ptr_stride` | 51.684 us | 1.720x | 2.1032064978498966e-12 |

Winner: `ptr_stride`

### max_levels = 4

| Variant | Time | Vs `package_avx` | Error |
|---|---:|---:|---:|
| `package_scalar` | 1099.890 us | 0.361x | 1.9895196601282805e-13 |
| `package_avx` | 396.717 us | 1.000x | 3.375077994860476e-14 |
| `vector` | 381.515 us | 1.040x | 3.375077994860476e-14 |
| `ptr_vector` | 911.289 us | 0.435x | 3.375077994860476e-14 |
| `stride` | 376.818 us | 1.053x | 3.375077994860476e-14 |
| `ptr_stride` | 459.519 us | 0.863x | 3.375077994860476e-14 |

Winner: `stride`

### max_levels = 5

| Variant | Time | Vs `package_avx` | Error |
|---|---:|---:|---:|
| `package_scalar` | 10482.404 us | 0.280x | 3.0091484859440243e-12 |
| `package_avx` | 2934.409 us | 1.000x | 1.7763568394002505e-14 |
| `vector` | 3718.083 us | 0.789x | 1.7763568394002505e-14 |
| `ptr_vector` | 9073.244 us | 0.323x | 1.7763568394002505e-14 |
| `stride` | 3361.064 us | 0.873x | 1.7763568394002505e-14 |
| `ptr_stride` | 2502.658 us | 1.173x | 1.7763568394002505e-14 |

Winner: `ptr_stride`

### max_levels = 6

| Variant | Time | Vs `package_avx` | Error |
|---|---:|---:|---:|
| `package_scalar` | 86167.589 us | 0.325x | 2.2319923687064147e-11 |
| `package_avx` | 28035.277 us | 1.000x | 5.435651928564766e-13 |
| `vector` | 34550.160 us | 0.811x | 5.435651928564766e-13 |
| `ptr_vector` | 28595.491 us | 0.980x | 5.435651928564766e-13 |
| `stride` | 27708.106 us | 1.012x | 5.435651928564766e-13 |
| `ptr_stride` | 28484.808 us | 0.984x | 5.435651928564766e-13 |

Winner: `stride`, but only by about 1.2%

## Takeaways

- `StrideArray` and `PtrArray` can help, but the gains are narrow and depth-dependent.
- `ptr_stride` is the best result at `max_levels=3` and `max_levels=5`.
- `stride` is best at `max_levels=4` and `max_levels=6`.
- `ptr_vector` is consistently weak in this package-style adaptive benchmark.
- At deeper refinement (`max_levels=6`), the best container variant is only about `1.012x` faster than the current package AVX baseline.

## Interpretation

The fairer adaptive benchmark says the current package AVX implementation is already close to the local optimum. `StrideArray` or `PtrArray(StrideArray)` might still be worth adopting for the adaptive path, but the upside looks modest compared with the fixed-grid case, and the best choice changes with refinement depth.

## Recursive Hybrid Note

I also prototyped a hybrid cache that uses:

- `SVector` for small read-only node and weight tables
- plain `Vector` for the mutable hoisted scratch buffers
- tuple-recursive level traversal to avoid the huge allocations seen in the naive mixed static/dynamic loop

Key finding:

- the tuple-recursive hybrid is now `0` allocation

Targeted results after the recursive rewrite:

| Case | Time | Allocs | Error |
|---|---:|---:|---:|
| `max_levels=3`, `hybrid64` | 69.492 us | 0 | 2.1032064978498966e-12 |
| `max_levels=3`, `hybrid128` | 45.590 us | 0 | 2.1032064978498966e-12 |
| `max_levels=4`, `hybrid64` | 366.938 us | 0 | 3.375077994860476e-14 |
| `max_levels=4`, `hybrid128` | 371.481 us | 0 | 3.375077994860476e-14 |
| `max_levels=5`, `hybrid64` | 2866.311 us | 0 | 1.7763568394002505e-14 |
| `max_levels=5`, `hybrid128` | 2892.829 us | 0 | 1.7763568394002505e-14 |

Interpretation:

- the recursive hybrid fixes the allocation problem completely
- `hybrid128` is very strong at shallow depth
- by `max_levels=5`, both hybrids beat the plain experimental `vector` run from the targeted check
- the remaining tradeoff is compile cost: deeper hybrid runs specialize much more heavily than the plain `Vector` path
