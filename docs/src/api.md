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
