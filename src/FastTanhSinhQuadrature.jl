# FastTanhSinhQuadrature.jl: Fast and high-precision numerical integration using Tanh-Sinh (Double Exponential) quadrature.
# Copyright (C) 2024-2026 Stamatis Vretinaris
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
