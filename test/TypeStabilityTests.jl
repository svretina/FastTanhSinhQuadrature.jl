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

using FastTanhSinhQuadrature
using JET
using MultiFloats
using Quadmath
using Test
using StaticArrays

# Test Inputs (Global Scope)
const T_val = Float64
const x_vec, w_vec, h_val = tanhsinh(T_val, Val(10))
const T32_val = Float32
const x32_vec, w32_vec, h32_val = tanhsinh(T32_val, Val(10))

# 1D Bounds
const low1d_val, up1d_val = -1.0, 1.0

# 2D/3D Bounds (StaticArrays)
const low2d_val = SVector(-1.0, -1.0)
const up2d_val = SVector(1.0, 1.0)
const low3d_val = SVector(-1.0, -1.0, -1.0)
const up3d_val = SVector(1.0, 1.0, 1.0)
const low2d32_val = SVector(-1.0f0, -1.0f0)
const up2d32_val = SVector(1.0f0, 1.0f0)
const low3d32_val = SVector(-1.0f0, -1.0f0, -1.0f0)
const up3d32_val = SVector(1.0f0, 1.0f0, 1.0f0)

# Simple polynomial functions
func1(x) = x^2
func2(x, y) = x * y
func3(x, y, z) = x * y * z
func_cmpl(x, bmx, xma) = inv(sqrt(bmx * xma))

struct JetOffsetFunctor1D{C}
    c::C
end

(f::JetOffsetFunctor1D)(x) = x + f.c

struct JetAffineFunctor2D{C}
    c::C
end

(f::JetAffineFunctor2D)(x, y) = f.c * (x + y)

struct JetAffineFunctor3D{C}
    c::C
end

(f::JetAffineFunctor3D)(x, y, z) = f.c * (x + y + z)

struct JetEndpointFunctor1D{C}
    c::C
end

(f::JetEndpointFunctor1D)(x, bmx, xma) = f.c / sqrt(bmx * xma)

const jet_functor1_64 = JetOffsetFunctor1D(1.0)
const jet_functor2_32 = JetAffineFunctor2D(1.0f0)
const jet_functor3_mf = JetAffineFunctor3D(Float64x2(1))
const jet_functor_cmpl_32 = JetEndpointFunctor1D(1.0f0)

@testset "Type Stability" begin
    # JET 0.9 uses `target_defined_modules`, while JET 0.11 uses `target_modules`.
    # Keep compatibility across Julia 1.11/1.12 JET lines.
    config = if Base.pkgversion(JET) < v"0.11.0"
        (target_defined_modules=true,)
    else
        (target_modules=(FastTanhSinhQuadrature,),)
    end

    # Helper function to reduce boilerplate
    function check_opt(f, argtypes)
        JET.test_opt(f, argtypes; config...)
    end

    function check_call(f, argtypes)
        JET.test_call(f, argtypes; config...)
    end

    # --- Core ---
    check_opt(tanhsinh, (Type{Float64}, Int))
    check_call(tanhsinh, (Type{Float64}, Int))

    check_opt(tanhsinh, (Type{Float64}, Val{10}))
    check_call(tanhsinh, (Type{Float64}, Val{10}))

    check_opt(tanhsinh, (Type{Float32x2}, Int))
    check_call(tanhsinh, (Type{Float32x2}, Int))

    check_opt(tanhsinh, (Type{Float64x2}, Int))
    check_call(tanhsinh, (Type{Float64x2}, Int))

    check_opt(tanhsinh, (Type{BigFloat}, Int))
    check_call(tanhsinh, (Type{BigFloat}, Int))

    check_opt(tanhsinh, (Type{Double64}, Int))
    check_call(tanhsinh, (Type{Double64}, Int))

    check_opt(tanhsinh, (Type{Float128}, Int))
    check_call(tanhsinh, (Type{Float128}, Int))

    # --- Integrate 1D ---
    check_opt(integrate1D, (typeof(func1), Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate1D, (typeof(func1), Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate1D, (typeof(func1), Vector{Float32}, Vector{Float32}, Float32))
    check_call(integrate1D, (typeof(func1), Vector{Float32}, Vector{Float32}, Float32))

    check_opt(integrate1D, (typeof(func1), Float64, Float64, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate1D, (typeof(func1), Float64, Float64, Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate1D, (typeof(func1), Float32, Float32, Vector{Float32}, Vector{Float32}, Float32))
    check_call(integrate1D, (typeof(func1), Float32, Float32, Vector{Float32}, Vector{Float32}, Float32))

    check_opt(integrate1D, (typeof(func1), Float32x2, Float32x2, Vector{Float32x2}, Vector{Float32x2}, Float32x2))
    check_call(integrate1D, (typeof(func1), Float32x2, Float32x2, Vector{Float32x2}, Vector{Float32x2}, Float32x2))

    check_opt(integrate1D, (typeof(func1), Float64x2, Float64x2, Vector{Float64x2}, Vector{Float64x2}, Float64x2))
    check_call(integrate1D, (typeof(func1), Float64x2, Float64x2, Vector{Float64x2}, Vector{Float64x2}, Float64x2))

    check_opt(integrate1D, (typeof(func1), BigFloat, BigFloat, Vector{BigFloat}, Vector{BigFloat}, BigFloat))
    check_call(integrate1D, (typeof(func1), BigFloat, BigFloat, Vector{BigFloat}, Vector{BigFloat}, BigFloat))

    check_opt(integrate1D, (typeof(func1), Double64, Double64, Vector{Double64}, Vector{Double64}, Double64))
    check_call(integrate1D, (typeof(func1), Double64, Double64, Vector{Double64}, Vector{Double64}, Double64))

    check_opt(integrate1D, (typeof(func1), Float128, Float128, Vector{Float128}, Vector{Float128}, Float128))
    check_call(integrate1D, (typeof(func1), Float128, Float128, Vector{Float128}, Vector{Float128}, Float128))

    # --- Integrate 1D AVX ---
    check_opt(integrate1D_avx, (typeof(func1), Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate1D_avx, (typeof(func1), Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate1D_avx, (typeof(func1), Vector{Float32}, Vector{Float32}, Float32))
    check_call(integrate1D_avx, (typeof(func1), Vector{Float32}, Vector{Float32}, Float32))

    check_opt(integrate1D_avx, (typeof(func1), Float64, Float64, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate1D_avx, (typeof(func1), Float64, Float64, Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate1D_avx, (typeof(func1), Float32, Float32, Vector{Float32}, Vector{Float32}, Float32))
    check_call(integrate1D_avx, (typeof(func1), Float32, Float32, Vector{Float32}, Vector{Float32}, Float32))

    # --- Integrate 2D ---
    check_opt(integrate2D, (typeof(func2), Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate2D, (typeof(func2), Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate2D, (typeof(func2), Vector{Float32}, Vector{Float32}, Float32))
    check_call(integrate2D, (typeof(func2), Vector{Float32}, Vector{Float32}, Float32))

    check_opt(integrate2D, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate2D, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate2D, (typeof(func2), SVector{2,Float32}, SVector{2,Float32}, Vector{Float32}, Vector{Float32}, Float32))
    check_call(integrate2D, (typeof(func2), SVector{2,Float32}, SVector{2,Float32}, Vector{Float32}, Vector{Float32}, Float32))

    # --- Integrate 2D AVX ---
    check_opt(integrate2D_avx, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate2D_avx, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate2D_avx, (typeof(func2), SVector{2,Float32}, SVector{2,Float32}, Vector{Float32}, Vector{Float32}, Float32))
    check_call(integrate2D_avx, (typeof(func2), SVector{2,Float32}, SVector{2,Float32}, Vector{Float32}, Vector{Float32}, Float32))

    # --- Integrate 3D ---
    check_opt(integrate3D, (typeof(func3), Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate3D, (typeof(func3), Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate3D, (typeof(func3), Vector{Float32}, Vector{Float32}, Float32))
    check_call(integrate3D, (typeof(func3), Vector{Float32}, Vector{Float32}, Float32))

    check_opt(integrate3D, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate3D, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}, Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate3D, (typeof(func3), SVector{3,Float32}, SVector{3,Float32}, Vector{Float32}, Vector{Float32}, Float32))
    check_call(integrate3D, (typeof(func3), SVector{3,Float32}, SVector{3,Float32}, Vector{Float32}, Vector{Float32}, Float32))

    # --- Integrate 3D AVX ---
    check_opt(integrate3D_avx, (typeof(func3), Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate3D_avx, (typeof(func3), Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate3D_avx, (typeof(func3), Vector{Float32}, Vector{Float32}, Float32))
    check_call(integrate3D_avx, (typeof(func3), Vector{Float32}, Vector{Float32}, Float32))

    check_opt(integrate3D_avx, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate3D_avx, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}, Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate3D_avx, (typeof(func3), SVector{3,Float32}, SVector{3,Float32}, Vector{Float32}, Vector{Float32}, Float32))
    check_call(integrate3D_avx, (typeof(func3), SVector{3,Float32}, SVector{3,Float32}, Vector{Float32}, Vector{Float32}, Float32))

    # --- Adaptive Integration ---
    check_opt(adaptive_integrate_1D, (Type{Float64}, typeof(func1), Float64, Float64))
    check_call(adaptive_integrate_1D, (Type{Float64}, typeof(func1), Float64, Float64))

    check_opt(adaptive_integrate_1D, (Type{Float32}, typeof(func1), Float32, Float32))
    check_call(adaptive_integrate_1D, (Type{Float32}, typeof(func1), Float32, Float32))

    check_opt(adaptive_integrate_1D, (Type{Float32x2}, typeof(func1), Float32x2, Float32x2))
    check_call(adaptive_integrate_1D, (Type{Float32x2}, typeof(func1), Float32x2, Float32x2))

    check_opt(adaptive_integrate_1D, (Type{Float64x2}, typeof(func1), Float64x2, Float64x2))
    check_call(adaptive_integrate_1D, (Type{Float64x2}, typeof(func1), Float64x2, Float64x2))

    check_opt(adaptive_integrate_1D, (Type{BigFloat}, typeof(func1), BigFloat, BigFloat))
    check_call(adaptive_integrate_1D, (Type{BigFloat}, typeof(func1), BigFloat, BigFloat))

    check_opt(adaptive_integrate_1D, (Type{Double64}, typeof(func1), Double64, Double64))
    check_call(adaptive_integrate_1D, (Type{Double64}, typeof(func1), Double64, Double64))

    check_opt(adaptive_integrate_1D, (Type{Float128}, typeof(func1), Float128, Float128))
    check_call(adaptive_integrate_1D, (Type{Float128}, typeof(func1), Float128, Float128))

    check_opt(adaptive_integrate_1D, (Type{Float64}, typeof(jet_functor1_64), Float64, Float64))
    check_call(adaptive_integrate_1D, (Type{Float64}, typeof(jet_functor1_64), Float64, Float64))

    check_opt(adaptive_integrate_1D_avx, (Type{Float64}, typeof(func1), Float64, Float64))
    check_call(adaptive_integrate_1D_avx, (Type{Float64}, typeof(func1), Float64, Float64))

    check_opt(adaptive_integrate_1D_avx, (Type{Float32}, typeof(func1), Float32, Float32))
    check_call(adaptive_integrate_1D_avx, (Type{Float32}, typeof(func1), Float32, Float32))

    check_opt(adaptive_integrate_2D, (Type{Float64}, typeof(func2), SVector{2,Float64}, SVector{2,Float64}))
    check_call(adaptive_integrate_2D, (Type{Float64}, typeof(func2), SVector{2,Float64}, SVector{2,Float64}))

    check_opt(adaptive_integrate_2D_avx, (Type{Float64}, typeof(func2), SVector{2,Float64}, SVector{2,Float64}))
    check_call(adaptive_integrate_2D_avx, (Type{Float64}, typeof(func2), SVector{2,Float64}, SVector{2,Float64}))

    check_opt(adaptive_integrate_2D, (Type{Float32}, typeof(func2), SVector{2,Float32}, SVector{2,Float32}))
    check_call(adaptive_integrate_2D, (Type{Float32}, typeof(func2), SVector{2,Float32}, SVector{2,Float32}))

    check_opt(adaptive_integrate_2D_avx, (Type{Float32}, typeof(func2), SVector{2,Float32}, SVector{2,Float32}))
    check_call(adaptive_integrate_2D_avx, (Type{Float32}, typeof(func2), SVector{2,Float32}, SVector{2,Float32}))

    check_opt(adaptive_integrate_2D, (Type{Float32x2}, typeof(func2), SVector{2,Float32x2}, SVector{2,Float32x2}))
    check_call(adaptive_integrate_2D, (Type{Float32x2}, typeof(func2), SVector{2,Float32x2}, SVector{2,Float32x2}))

    check_opt(adaptive_integrate_2D, (Type{BigFloat}, typeof(func2), SVector{2,BigFloat}, SVector{2,BigFloat}))
    check_call(adaptive_integrate_2D, (Type{BigFloat}, typeof(func2), SVector{2,BigFloat}, SVector{2,BigFloat}))

    check_opt(adaptive_integrate_3D, (Type{Float64}, typeof(func3), SVector{3,Float64}, SVector{3,Float64}))
    check_call(adaptive_integrate_3D, (Type{Float64}, typeof(func3), SVector{3,Float64}, SVector{3,Float64}))

    check_opt(adaptive_integrate_3D_avx, (Type{Float64}, typeof(func3), SVector{3,Float64}, SVector{3,Float64}))
    check_call(adaptive_integrate_3D_avx, (Type{Float64}, typeof(func3), SVector{3,Float64}, SVector{3,Float64}))

    check_opt(adaptive_integrate_3D, (Type{Float32}, typeof(func3), SVector{3,Float32}, SVector{3,Float32}))
    check_call(adaptive_integrate_3D, (Type{Float32}, typeof(func3), SVector{3,Float32}, SVector{3,Float32}))

    check_opt(adaptive_integrate_3D_avx, (Type{Float32}, typeof(func3), SVector{3,Float32}, SVector{3,Float32}))
    check_call(adaptive_integrate_3D_avx, (Type{Float32}, typeof(func3), SVector{3,Float32}, SVector{3,Float32}))

    check_opt(adaptive_integrate_3D, (Type{Float64x2}, typeof(func3), SVector{3,Float64x2}, SVector{3,Float64x2}))
    check_call(adaptive_integrate_3D, (Type{Float64x2}, typeof(func3), SVector{3,Float64x2}, SVector{3,Float64x2}))

    check_opt(adaptive_integrate_3D, (Type{Double64}, typeof(func3), SVector{3,Double64}, SVector{3,Double64}))
    check_call(adaptive_integrate_3D, (Type{Double64}, typeof(func3), SVector{3,Double64}, SVector{3,Double64}))

    check_opt(adaptive_integrate_3D, (Type{Float128}, typeof(func3), SVector{3,Float128}, SVector{3,Float128}))
    check_call(adaptive_integrate_3D, (Type{Float128}, typeof(func3), SVector{3,Float128}, SVector{3,Float128}))

    check_opt(adaptive_integrate_1D_cmpl, (Type{Float32}, typeof(func_cmpl), Float32, Float32))
    check_call(adaptive_integrate_1D_cmpl, (Type{Float32}, typeof(func_cmpl), Float32, Float32))

    check_opt(adaptive_integrate_1D_cmpl_avx, (Type{Float64}, typeof(func_cmpl), Float64, Float64))
    check_call(adaptive_integrate_1D_cmpl_avx, (Type{Float64}, typeof(func_cmpl), Float64, Float64))

    check_opt(adaptive_integrate_1D_cmpl_avx, (Type{Float32}, typeof(func_cmpl), Float32, Float32))
    check_call(adaptive_integrate_1D_cmpl_avx, (Type{Float32}, typeof(func_cmpl), Float32, Float32))

    check_opt(adaptive_integrate_1D_cmpl, (Type{Float32}, typeof(jet_functor_cmpl_32), Float32, Float32))
    check_call(adaptive_integrate_1D_cmpl, (Type{Float32}, typeof(jet_functor_cmpl_32), Float32, Float32))

    # --- Quad Interface (1D) ---
    check_opt(quad, (typeof(func1), Float64, Float64))
    check_call(quad, (typeof(func1), Float64, Float64))

    check_opt(quad, (typeof(func1), Float32, Float32))
    check_call(quad, (typeof(func1), Float32, Float32))

    check_opt(quad, (typeof(func1), Float32x2, Float32x2))
    check_call(quad, (typeof(func1), Float32x2, Float32x2))

    check_opt(quad, (typeof(func1), Float64x2, Float64x2))
    check_call(quad, (typeof(func1), Float64x2, Float64x2))

    check_opt(quad, (typeof(func1), BigFloat, BigFloat))
    check_call(quad, (typeof(func1), BigFloat, BigFloat))

    check_opt(quad, (typeof(func1), Double64, Double64))
    check_call(quad, (typeof(func1), Double64, Double64))

    check_opt(quad, (typeof(func1), Float128, Float128))
    check_call(quad, (typeof(func1), Float128, Float128))

    check_opt(quad, (typeof(jet_functor1_64), Float64, Float64))
    check_call(quad, (typeof(jet_functor1_64), Float64, Float64))

    # --- Quad Interface (2D) ---
    check_opt(quad, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}))
    check_call(quad, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}))

    check_opt(quad, (typeof(func2), SVector{2,Float32}, SVector{2,Float32}))
    check_call(quad, (typeof(func2), SVector{2,Float32}, SVector{2,Float32}))

    check_opt(quad, (typeof(func2), SVector{2,Float32x2}, SVector{2,Float32x2}))
    check_call(quad, (typeof(func2), SVector{2,Float32x2}, SVector{2,Float32x2}))

    check_opt(quad, (typeof(func2), SVector{2,BigFloat}, SVector{2,BigFloat}))
    check_call(quad, (typeof(func2), SVector{2,BigFloat}, SVector{2,BigFloat}))

    check_opt(quad, (typeof(jet_functor2_32), SVector{2,Float32}, SVector{2,Float32}))
    check_call(quad, (typeof(jet_functor2_32), SVector{2,Float32}, SVector{2,Float32}))

    # Generic function for dynamic dispatch testing (accepts 2 or 3 args)
    # JET analyzes all branches of `quad(..., Vector, Vector)`, including 2D and 3D.
    # A 2-arg function like `func2` fails the 3D branch check.
    func_any(args...) = sum(args)

    # Quad with Vectors (error check only)
    # This allows checking that `quad` handles Vectors without immediate errors,
    # and that the branching logic is analyzed. We use `func_any` to satisfy both 2D and 3D branches.
    check_call(quad, (typeof(func_any), Vector{Float64}, Vector{Float64}))
    check_call(quad, (typeof(func_any), Vector{Float32}, Vector{Float32}))

    # --- Quad Interface (3D) ---
    check_opt(quad, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}))
    check_call(quad, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}))

    check_opt(quad, (typeof(func3), SVector{3,Float32}, SVector{3,Float32}))
    check_call(quad, (typeof(func3), SVector{3,Float32}, SVector{3,Float32}))

    check_opt(quad, (typeof(func3), SVector{3,Float64x2}, SVector{3,Float64x2}))
    check_call(quad, (typeof(func3), SVector{3,Float64x2}, SVector{3,Float64x2}))

    check_opt(quad, (typeof(func3), SVector{3,Double64}, SVector{3,Double64}))
    check_call(quad, (typeof(func3), SVector{3,Double64}, SVector{3,Double64}))

    check_opt(quad, (typeof(func3), SVector{3,Float128}, SVector{3,Float128}))
    check_call(quad, (typeof(func3), SVector{3,Float128}, SVector{3,Float128}))

    check_opt(quad, (typeof(jet_functor3_mf), SVector{3,Float64x2}, SVector{3,Float64x2}))
    check_call(quad, (typeof(jet_functor3_mf), SVector{3,Float64x2}, SVector{3,Float64x2}))

    # --- Quad Split ---
    check_opt(quad_split, (typeof(func1), Float64, Float64, Float64))
    check_call(quad_split, (typeof(func1), Float64, Float64, Float64))

    check_opt(quad_split, (typeof(func1), Float32, Float32, Float32))
    check_call(quad_split, (typeof(func1), Float32, Float32, Float32))

    check_opt(quad_split, (typeof(func1), Float32x2, Float32x2, Float32x2))
    check_call(quad_split, (typeof(func1), Float32x2, Float32x2, Float32x2))

    check_opt(quad_split, (typeof(func1), BigFloat, BigFloat, BigFloat))
    check_call(quad_split, (typeof(func1), BigFloat, BigFloat, BigFloat))

    check_opt(quad_split, (typeof(func1), Double64, Double64, Double64))
    check_call(quad_split, (typeof(func1), Double64, Double64, Double64))

    check_opt(quad_split, (typeof(func1), Float128, Float128, Float128))
    check_call(quad_split, (typeof(func1), Float128, Float128, Float128))

    check_opt(quad_split, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, SVector{2,Float64}))
    check_call(quad_split, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, SVector{2,Float64}))

    check_opt(quad_split, (typeof(func2), SVector{2,Float32}, SVector{2,Float32}, SVector{2,Float32}))
    check_call(quad_split, (typeof(func2), SVector{2,Float32}, SVector{2,Float32}, SVector{2,Float32}))

    check_opt(quad_split, (typeof(func3), SVector{3,Float32}, SVector{3,Float32}, SVector{3,Float32}))
    check_call(quad_split, (typeof(func3), SVector{3,Float32}, SVector{3,Float32}, SVector{3,Float32}))

    check_opt(quad_cmpl, (typeof(func_cmpl), Float32, Float32))
    check_call(quad_cmpl, (typeof(func_cmpl), Float32, Float32))

    check_opt(quad_cmpl, (typeof(func_cmpl), Float32x2, Float32x2))
    check_call(quad_cmpl, (typeof(func_cmpl), Float32x2, Float32x2))

    check_opt(quad_cmpl, (typeof(func_cmpl), Float64x2, Float64x2))
    check_call(quad_cmpl, (typeof(func_cmpl), Float64x2, Float64x2))

    check_opt(quad_cmpl, (typeof(func_cmpl), BigFloat, BigFloat))
    check_call(quad_cmpl, (typeof(func_cmpl), BigFloat, BigFloat))

    check_opt(quad_cmpl, (typeof(func_cmpl), Double64, Double64))
    check_call(quad_cmpl, (typeof(func_cmpl), Double64, Double64))

    check_opt(quad_cmpl, (typeof(func_cmpl), Float128, Float128))
    check_call(quad_cmpl, (typeof(func_cmpl), Float128, Float128))

    check_opt(quad_cmpl, (typeof(jet_functor_cmpl_32), Float32, Float32))
    check_call(quad_cmpl, (typeof(jet_functor_cmpl_32), Float32, Float32))
end
