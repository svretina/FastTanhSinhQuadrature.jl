module FastTanhSinhQuadrature

using StaticArrays
using LambertW
using LoopVectorization

export tanhsinh, quad, quad_split, quad_cmpl,
    integrate1D, integrate2D, integrate3D,
    integrate1D_avx, integrate2D_avx, integrate3D_avx,
    adaptive_integrate_1D, adaptive_integrate_2D, adaptive_integrate_3D,
    adaptive_integrate_1D_cmpl

# Core transformation functions and node generation
include("core.jl")

# 1D Integration
include("integrate1D.jl")

# 2D Integration
include("integrate2D.jl")

# 3D Integration
include("integrate3D.jl")

# Adaptive integration
include("adaptive.jl")

# High-level quad interface
include("quad.jl")

end
