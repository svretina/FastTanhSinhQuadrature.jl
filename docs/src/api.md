# API Reference

```@index
Pages = ["api.md"]
```

## High-Level API

These functions provide the easiest interface for integration with automatic adaptation.

```@docs
quad
quad_split
quad_cmpl
```

## Tolerances and Fixed-Grid Interfaces

- `quad`, `quad_split`, and `quad_cmpl` are adaptive interfaces.
- Their stopping criterion is based on successive refinement levels:
  `err_est = abs(I_new - I_old)`, and refinement stops when
  `err_est <= max(atol, rtol * abs(I_new))`.
- If `rtol` is omitted and `atol == 0`, the default is `rtol = sqrt(eps(T))`,
  where `T` is the promoted floating-point endpoint type.
- For repeated calls, prebuild adaptive caches and pass `cache=...` to
  `quad` / `quad_split` / `quad_cmpl` or `adaptive_integrate_*` directly.
- `integrate1D`, `integrate1D_avx`, `integrate2D`, `integrate2D_avx`,
  `integrate3D`, and `integrate3D_avx` are fixed-grid interfaces.
- Fixed-grid interfaces do not estimate error internally; choose `N`
  externally, for example by doubling `N` until two successive results meet
  the same `max(atol, rtol * abs(I))` criterion.

## 1D Integration

Functions for integrating over 1D domains using pre-computed nodes and weights.

```@docs
integrate1D
integrate1D_avx
adaptive_integrate_1D
adaptive_integrate_1D_cmpl
```

## 2D Integration

Functions for integrating over 2D rectangular domains.

```@docs
integrate2D
integrate2D_avx
adaptive_integrate_2D
```

## 3D Integration

Functions for integrating over 3D box domains.

```@docs
integrate3D
integrate3D_avx
adaptive_integrate_3D
```

## Node Generation

```@docs
tanhsinh
```

## Adaptive Cache Utilities

```@docs
adaptive_cache_1D
adaptive_cache_2D
adaptive_cache_3D
```

## Bounds and Containers

- `quad` / `quad_split` / `quad_cmpl` accept scalar bounds as any `Real` type and convert internally to floating-point.
- For 2D/3D high-level calls, bounds can be `SVector` or regular vectors of reals (`length == 2` or `3`).
- Pre-computed multidimensional APIs (`integrate2D`, `integrate3D`, and `_avx` variants) support:
  - typed static bounds (`SVector{N,T}`),
  - mixed real static bounds (`SVector{N,<:Real}`),
  - vector bounds (`AbstractVector{<:Real}` with matching length).

## All Exported Symbols

```@autodocs
Modules = [FastTanhSinhQuadrature]
```
