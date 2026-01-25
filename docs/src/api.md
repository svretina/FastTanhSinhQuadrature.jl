# API Reference

```@index
Pages = ["api.md"]
```

## High-Level API

These functions provide the easiest interface for integration with automatic adaptation.

```@docs
quad
quad_split
```

## 1D Integration

Functions for integrating over 1D domains using pre-computed nodes and weights.

```@docs
integrate1D
integrate1D_avx
adaptive_integrate_1D
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

## All Exported Symbols

```@autodocs
Modules = [FastTanhSinhQuadrature]
```
