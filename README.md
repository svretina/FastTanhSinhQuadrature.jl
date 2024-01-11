# FastTanhSinhQuadrature

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://svretina.github.io/FastTanhSinhQuadrature.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://svretina.github.io/FastTanhSinhQuadrature.jl/dev/)
[![Build Status](https://github.com/svretina/FastTanhSinhQuadrature.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/svretina/FastTanhSinhQuadrature.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/svretina/FastTanhSinhQuadrature.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/svretina/FastTanhSinhQuadrature.jl)

This library implements the [Tanh-sinh](https://en.wikipedia.org/wiki/Tanh-sinh_quadrature) quadrature rule in 1,2 and 3 dimensions. This method is capable of handling singular functions, but the singularity must be located at the endpoints of the integration domain.

!warning! this library is under heavy development. Use at your own risk

## Other TanhSinh quadrature libraries
    - [TanhSinhQuadrature.jl](https://github.com/eschnett/TanhSinhQuadrature.jl), D-dimensional quadrature but slower
    - [DoubleExponentialFormulas.jl](https://github.com/machakann/DoubleExponentialFormulas.jl) 1D quadrature
    
## TODOs
- add 2D and 3D tests
- add docs
- write benchmarks
