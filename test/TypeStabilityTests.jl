using FastTanhSinhQuadrature
using JET
using Test
using StaticArrays

# Test Inputs (Global Scope)
const T_val = Float64
const x_vec, w_vec, h_val = tanhsinh(T_val, Val(10))

# 1D Bounds
const low1d_val, up1d_val = -1.0, 1.0

# 2D/3D Bounds (StaticArrays)
const low2d_val = SVector(-1.0, -1.0)
const up2d_val = SVector(1.0, 1.0)
const low3d_val = SVector(-1.0, -1.0, -1.0)
const up3d_val = SVector(1.0, 1.0, 1.0)

# Simple polynomial functions
func1(x) = x^2
func2(x, y) = x * y
func3(x, y, z) = x * y * z

@testset "Type Stability" begin
    # Configuration
    # target_defined_modules=true: focus analysis on this package
    # ignored_modules=[StaticArrays]: suppress false positives from StaticArrays internal dispatch
    config = (target_defined_modules=true,)

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

    # --- Integrate 1D ---
    check_opt(integrate1D, (typeof(func1), Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate1D, (typeof(func1), Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate1D, (typeof(func1), Float64, Float64, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate1D, (typeof(func1), Float64, Float64, Vector{Float64}, Vector{Float64}, Float64))

    # --- Integrate 1D AVX ---
    check_opt(integrate1D_avx, (typeof(func1), Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate1D_avx, (typeof(func1), Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate1D_avx, (typeof(func1), Float64, Float64, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate1D_avx, (typeof(func1), Float64, Float64, Vector{Float64}, Vector{Float64}, Float64))

    # --- Integrate 2D ---
    check_opt(integrate2D, (typeof(func2), Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate2D, (typeof(func2), Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate2D, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate2D, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, Vector{Float64}, Vector{Float64}, Float64))

    # --- Integrate 2D AVX ---
    check_opt(integrate2D_avx, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate2D_avx, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, Vector{Float64}, Vector{Float64}, Float64))

    # --- Integrate 3D ---
    check_opt(integrate3D, (typeof(func3), Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate3D, (typeof(func3), Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate3D, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate3D, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}, Vector{Float64}, Vector{Float64}, Float64))

    # --- Integrate 3D AVX ---
    check_opt(integrate3D_avx, (typeof(func3), Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate3D_avx, (typeof(func3), Vector{Float64}, Vector{Float64}, Float64))

    check_opt(integrate3D_avx, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}, Vector{Float64}, Vector{Float64}, Float64))
    check_call(integrate3D_avx, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}, Vector{Float64}, Vector{Float64}, Float64))

    # --- Adaptive Integration ---
    check_opt(adaptive_integrate_1D, (Type{Float64}, typeof(func1), Float64, Float64))
    check_call(adaptive_integrate_1D, (Type{Float64}, typeof(func1), Float64, Float64))

    check_opt(adaptive_integrate_2D, (Type{Float64}, typeof(func2), SVector{2,Float64}, SVector{2,Float64}))
    check_call(adaptive_integrate_2D, (Type{Float64}, typeof(func2), SVector{2,Float64}, SVector{2,Float64}))

    check_opt(adaptive_integrate_3D, (Type{Float64}, typeof(func3), SVector{3,Float64}, SVector{3,Float64}))
    check_call(adaptive_integrate_3D, (Type{Float64}, typeof(func3), SVector{3,Float64}, SVector{3,Float64}))

    # --- Quad Interface (1D) ---
    check_opt(quad, (typeof(func1), Float64, Float64))
    check_call(quad, (typeof(func1), Float64, Float64))

    # --- Quad Interface (2D) ---
    check_opt(quad, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}))
    check_call(quad, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}))

    # Generic function for dynamic dispatch testing (accepts 2 or 3 args)
    # JET analyzes all branches of `quad(..., Vector, Vector)`, including 2D and 3D.
    # A 2-arg function like `func2` fails the 3D branch check.
    func_any(args...) = sum(args)

    # Quad with Vectors (error check only)
    # This allows checking that `quad` handles Vectors without immediate errors,
    # and that the branching logic is analyzed. We use `func_any` to satisfy both 2D and 3D branches.
    check_call(quad, (typeof(func_any), Vector{Float64}, Vector{Float64}))

    # --- Quad Interface (3D) ---
    check_opt(quad, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}))
    check_call(quad, (typeof(func3), SVector{3,Float64}, SVector{3,Float64}))

    # --- Quad Split ---
    check_opt(quad_split, (typeof(func1), Float64, Float64, Float64))
    check_call(quad_split, (typeof(func1), Float64, Float64, Float64))

    check_opt(quad_split, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, SVector{2,Float64}))
    check_call(quad_split, (typeof(func2), SVector{2,Float64}, SVector{2,Float64}, SVector{2,Float64}))
end
