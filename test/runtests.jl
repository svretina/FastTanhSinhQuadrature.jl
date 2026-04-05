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
using StaticArrays
using DoubleFloats
using MultiFloats
using Quadmath
using Test
using Random
using LinearAlgebra

const Types = [Float32, Float64, Double64, BigFloat]
const rtol = Dict(Float32 => 10 * sqrt(eps(Float32)),
    Float64 => 100 * sqrt(eps(Float64)),
    Double64 => 10^7 * sqrt(eps(Double64)),
    BigFloat => 1.0e-18)
const ALLOW_AVX_FAILURES = VERSION >= v"1.13.0-0"

function run_avx_checks(f::Function)
    if ALLOW_AVX_FAILURES
        try
            f()
        catch err
            @info "AVX checks failed but are allowed on Julia $VERSION due to LoopVectorization compatibility: $(sprint(showerror, err))"
            @test_broken false
        end
    else
        f()
    end
end

inferred_return_type(f::Function) = only(code_typed(f, Tuple{}))[2]
typed_call_return_type(f, argtypes::Type{<:Tuple}) = only(code_typed(f, argtypes))[2]

@inline integral_monomial_1d(low::T, up::T, p::Int) where {T<:Real} =
    (up^(p + 1) - low^(p + 1)) / T(p + 1)

@inline integral_monomial_3d(low::SVector{3,T}, up::SVector{3,T},
    px::Int, py::Int, pz::Int) where {T<:Real} =
    integral_monomial_1d(low[1], up[1], px) *
    integral_monomial_1d(low[2], up[2], py) *
    integral_monomial_1d(low[3], up[3], pz)

f2_const(x, y) = one(x)
f2_xy(x, y) = x * y
f2_x2y2(x, y) = x^2 * y^2

f3_const(x, y, z) = one(x)
f3_xyz(x, y, z) = x * y * z
f3_x2y2z2(x, y, z) = x^2 * y^2 * z^2

struct OffsetFunctor1D{C}
    c::C
end

(f::OffsetFunctor1D)(x) = x + f.c

struct SingularFunctor1D{C}
    c::C
end

(f::SingularFunctor1D)(x) = f.c / sqrt(abs(x))

struct EndpointFunctor1D{C}
    c::C
end

(f::EndpointFunctor1D)(x, bmx, xma) = f.c / sqrt(bmx * xma)

struct AffineFunctor2D{C}
    c::C
end

(f::AffineFunctor2D)(x, y) = f.c * (x + y)

struct AffineFunctor3D{C}
    c::C
end

(f::AffineFunctor3D)(x, y, z) = f.c * (x + y + z)

@testset "FastTanhSinhQuadrature.jl" begin

    include("core.jl")
    include("caches.jl")

    @testset "Basic integration T=$T" for T in Types
        f0(x) = one(T)
        f1(x) = x + 1
        f2(x) = 3 * x^2
        x, w, h = tanhsinh(T, 80) # Increased N slightly for Float32 safety/generality
        # Using integrate1D(f, x, w, h)
        @test isapprox(integrate1D(f0, x, w, h), T(2), rtol=rtol[T])
        @test isapprox(integrate1D(f1, x, w, h), T(2), rtol=rtol[T])
        @test isapprox(integrate1D(f2, x, w, h), T(2), rtol=rtol[T])
    end

    Random.seed!(0)
    @testset "Integral bounds T=$T" for T in Types
        f0(x) = 1
        f1(x) = x + 1
        f2(x) = 3 * x^2

        xmin = -1 + T(rand(-5:5)) / 10
        xmax = +1 + T(rand(-5:5)) / 10

        x, w, h = tanhsinh(T, 80)
        # integrate1D(f, xmin, xmax, x, w, h)
        @test isapprox(integrate1D(f0, xmin, xmax, x, w, h), xmax - xmin, rtol=rtol[T])
        @test isapprox(integrate1D(f1, xmin, xmax, x, w, h), (xmax^2 - xmin^2) / 2 + xmax - xmin, rtol=rtol[T])
        @test isapprox(integrate1D(f2, xmin, xmax, x, w, h), (xmax^3 - xmin^3), rtol=rtol[T])
    end

    Random.seed!(0)
    @testset "Linearity T=$T" for T in Types
        a = 1 + T(rand(-5:5)) / 10
        b0 = 1 + T(rand(-5:5)) / 10
        b1 = 1 + T(rand(-5:5)) / 10
        b2 = 1 + T(rand(-5:5)) / 10
        c0 = 1 + T(rand(-5:5)) / 10
        c1 = 1 + T(rand(-5:5)) / 10
        c2 = 1 + T(rand(-5:5)) / 10

        x, w, h = tanhsinh(T, 80)
        f(x) = b0 + b1 * x + b2 * x^2
        g(x) = c0 + c1 * x + c2 * x^2

        F1 = integrate1D(f, x, w, h)
        G1 = integrate1D(g, x, w, h)

        # Test AVX variant only for float types where it's relevant, but API handles it generally? 
        # integrate1D_avx might error or fallback for BigFloat if not careful, 
        # but my implementation has @turbo which LoopVectorization handles or fastmath fallback.
        # Let's test it mostly for floats

        afg(x) = a * f(x) + g(x)
        @test integrate1D(afg, x, w, h) ≈ a * F1 + G1

        # Testing AVX parity if T is standard float
        if T <: Union{Float32,Float64}
            run_avx_checks() do
                F2 = integrate1D_avx(f, x, w, h)
                G2 = integrate1D_avx(g, x, w, h)
                @test integrate1D_avx(afg, x, w, h) ≈ a * F2 + G2
                @test integrate1D_avx(f, one(T), -one(T), x, w, h) ≈ -F2
            end
        end

        F01 = integrate1D(f, -one(T), a, x, w, h)
        F11 = integrate1D(f, a, one(T), x, w, h)
        @test isapprox(F01 + F11, F1, rtol=rtol[T])

        # Flipped bounds should match orientation with type-appropriate tolerance.
        @test isapprox(quad(f, T(1), T(-1)), -F1, rtol=rtol[T])
    end

    include("families_1d.jl")

    @testset "2D polynomials for [-1, 1], T=$T" for T in Types
        n = T == Float32 ? 40 : 80
        x, w, h = tanhsinh(T, n)
        ψ(x, y) = one(T)
        f(x, y) = x * y
        g(x, y) = x^2 * y^2

        low = SVector{2,T}(-1.0, -1.0)
        up = SVector{2,T}(1.0, 1.0)

        Ψ = integrate2D(ψ, low, up, x, w, h)
        # integrate2D assumes low < up usually, let's just test basic
        F = integrate2D(f, low, up, x, w, h)
        G = integrate2D(g, low, up, x, w, h)

        # Relax tolerance slightly for BigFloat due to accumulation in higher dims with fixed points
        # BigFloat convergence is sensitive to h choice. 1e-20 is still extremely high precision.
        tol_check = T == BigFloat ? 1e-20 : (Base.rtoldefault(T) * 10)

        @test isapprox(Ψ, 4one(T); atol=tol_check, rtol=sqrt(eps(T)))
        @test isapprox(F, zero(T); atol=max(1e-12, tol_check))
        @test isapprox(G, T(4) / T(9); atol=tol_check, rtol=sqrt(eps(T)))

        # Test quad helper 2D flip
        @test isapprox(quad(ψ, up, low), 4one(T); atol=tol_check, rtol=sqrt(eps(T)))
    end

    include("families_2d.jl")

    include("families_3d.jl")

    @testset "Adaptive integration 1D" begin
        f(x) = exp(x)
        val_true = exp(1.0) - exp(0.0)

        # Using quad which uses adaptive_integrate_1D
        val = quad(f, 0.0, 1.0; rtol=1e-8)
        @test isapprox(val, val_true, atol=1e-8)

        # Test with higher precision
        if Double64 in Types
            val_d = quad(f, Double64(0.0), Double64(1.0); rtol=1e-15)
            @test isapprox(val_d, val_true, atol=1e-15)
        end
    end

    @testset "Adaptive warning paths" begin
        @test_logs (:warn, r"adaptive_integrate_2D reached max_levels") adaptive_integrate_2D(
            Float64, (x, y) -> exp(x + y), SVector(-1.0, -1.0), SVector(1.0, 1.0);
            rtol=0.0, atol=0.0, max_levels=1, warn=true, cache=adaptive_cache_2D(Float64; max_levels=1)
        )
        @test_logs (:warn, r"adaptive_integrate_3D reached max_levels") adaptive_integrate_3D(
            Float64, (x, y, z) -> exp(x + y + z), SVector(-1.0, -1.0, -1.0), SVector(1.0, 1.0, 1.0);
            rtol=0.0, atol=0.0, max_levels=1, warn=true, cache=adaptive_cache_3D(Float64; max_levels=1)
        )
        @test_logs (:warn, r"adaptive_integrate_1D_cmpl reached max_levels") adaptive_integrate_1D_cmpl(
            Float64, (x, bmx, xma) -> exp(x) + bmx + xma, -1.0, 1.0;
            rtol=0.0, atol=0.0, max_levels=1, warn=true, cache=adaptive_cache_1D(Float64; max_levels=1, complement=true)
        )
    end

    @testset "Type preservation and inference" begin
        x32, w32, h32 = tanhsinh(Float32, 8)
        low2_32 = SVector(0.0f0, 0.0f0)
        up2_32 = SVector(1.0f0, 1.0f0)
        low3_32 = SVector(0.0f0, 0.0f0, 0.0f0)
        up3_32 = SVector(1.0f0, 1.0f0, 1.0f0)

        f1(x) = exp(x)
        f2(x, y) = x + y
        f3(x, y, z) = x + y + z
        f_cmpl(x, bmx, xma) = inv(sqrt(bmx * xma))
        f_split(x) = inv(sqrt(abs(x)))
        functor1 = OffsetFunctor1D(1.0f0)
        functor2 = AffineFunctor2D(1.0f0)
        functor3 = AffineFunctor3D(1.0f0)
        functor_cmpl = EndpointFunctor1D(1.0f0)
        functor_split = SingularFunctor1D(1.0f0)

        @test typed_call_return_type(integrate1D, Tuple{Type{Float32}, typeof(f1), Int}) === Float32
        @test typed_call_return_type(FastTanhSinhQuadrature.integrate1D_cmpl, Tuple{Type{Float32}, typeof(f_cmpl), Int}) === Float32
        @test inferred_return_type(() -> integrate1D(f1, x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate1D(f1, 0.0f0, 1.0f0, x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate1D_avx(f1, x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate1D_avx(f1, 0.0f0, 1.0f0, x32, w32, h32)) === Float32
        @test inferred_return_type(() -> adaptive_integrate_1D(Float32, f1, 0.0f0, 1.0f0; rtol=1f-5, max_levels=0)) === Float32
        @test inferred_return_type(() -> adaptive_integrate_1D_cmpl(Float32, f_cmpl, -1.0f0, 1.0f0; rtol=1f-5, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad(f1, 0.0f0, 1.0f0; max_levels=0)) === Float32
        @test inferred_return_type(() -> quad(f1, 0.0f0, 1.0f0; rtol=1f-5, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad_cmpl(f_cmpl, -1.0f0, 1.0f0; max_levels=0)) === Float32
        @test inferred_return_type(() -> quad_cmpl(f_cmpl, -1.0f0, 1.0f0; rtol=1f-5, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad_split(f_split, 0.0f0, -1.0f0, 1.0f0; max_levels=0)) === Float32
        @test inferred_return_type(() -> quad_split(f_split, 0.0f0, -1.0f0, 1.0f0; rtol=1f-4, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad_split(f_split, 0.0f0; rtol=1f-4, max_levels=0)) === Float32
        @test typed_call_return_type(quad, Tuple{typeof(functor1), Float32, Float32}) === Float32
        @test inferred_return_type(() -> quad(functor1, 0.0f0, 1.0f0; rtol=1f-5, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad_cmpl(functor_cmpl, -1.0f0, 1.0f0; rtol=1f-5, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad_split(functor_split, 0.0f0, -1.0f0, 1.0f0; rtol=1f-4, max_levels=0)) === Float32

        @test inferred_return_type(() -> integrate2D(f2, x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate2D(f2, low2_32, up2_32, x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate2D(f2, [0.0f0, 0.0f0], [1.0f0, 1.0f0], x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate2D_avx(f2, low2_32, up2_32, x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate2D_avx(f2, [0.0f0, 0.0f0], [1.0f0, 1.0f0], x32, w32, h32)) === Float32
        @test inferred_return_type(() -> adaptive_integrate_2D(Float32, f2, low2_32, up2_32; rtol=1f-4, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad(f2, low2_32, up2_32; rtol=1f-4, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad(f2, [0.0f0, 0.0f0], [1.0f0, 1.0f0]; rtol=1f-4, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad_split(f2, SVector(0.5f0, 0.5f0), low2_32, up2_32; rtol=1f-4, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad(functor2, low2_32, up2_32; rtol=1f-4, max_levels=0)) === Float32

        @test inferred_return_type(() -> integrate3D(f3, x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate3D(f3, low3_32, up3_32, x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate3D(f3, [0.0f0, 0.0f0, 0.0f0], [1.0f0, 1.0f0, 1.0f0], x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate3D_avx(f3, x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate3D_avx(f3, low3_32, up3_32, x32, w32, h32)) === Float32
        @test inferred_return_type(() -> integrate3D_avx(f3, [0.0f0, 0.0f0, 0.0f0], [1.0f0, 1.0f0, 1.0f0], x32, w32, h32)) === Float32
        @test inferred_return_type(() -> adaptive_integrate_3D(Float32, f3, low3_32, up3_32; rtol=1f-3, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad(f3, low3_32, up3_32; rtol=1f-3, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad(f3, [0.0f0, 0.0f0, 0.0f0], [1.0f0, 1.0f0, 1.0f0]; rtol=1f-3, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad_split(f3, SVector(0.5f0, 0.5f0, 0.5f0), low3_32, up3_32; rtol=1f-3, max_levels=0)) === Float32
        @test inferred_return_type(() -> quad(functor3, low3_32, up3_32; rtol=1f-3, max_levels=0)) === Float32
    end

    @testset "MultiFloat type preservation and inference" begin
        f1(x) = exp(x)
        f2(x, y) = x + y
        f3(x, y, z) = x + y + z
        f_cmpl(x, bmx, xma) = inv(sqrt(bmx * xma))
        f_split(x) = inv(sqrt(abs(x)))

        @testset "T=$T" for T in (Float32x2, Float64x2)
            x, w, h = tanhsinh(T, 8)
            z = zero(T)
            o = one(T)
            low2 = SVector(zero(T), zero(T))
            up2 = SVector(one(T), one(T))
            low3 = SVector(zero(T), zero(T), zero(T))
            up3 = SVector(one(T), one(T), one(T))
            rtol1 = T <: Float32x2 ? T(1e-4) : T(1e-6)
            rtol2 = T <: Float32x2 ? T(1e-3) : T(1e-5)

            @test eltype(x) === T
            @test eltype(w) === T
            @test h isa T

            @test typed_call_return_type(integrate1D, Tuple{Type{T}, typeof(f1), Int}) === T
            @test inferred_return_type(() -> integrate1D(f1, x, w, h)) === T
            @test inferred_return_type(() -> integrate1D(f1, z, o, x, w, h)) === T
            @test inferred_return_type(() -> adaptive_integrate_1D(typeof(z), f1, z, o; rtol=rtol1, max_levels=0)) === T
            @test inferred_return_type(() -> quad(f1, z, o; max_levels=0)) === T
            @test inferred_return_type(() -> quad(f1, z, o; rtol=rtol1, max_levels=0)) === T
            @test inferred_return_type(() -> quad_cmpl(f_cmpl, -o, o; rtol=rtol1, max_levels=0)) === T
            @test inferred_return_type(() -> quad_split(f_split, z, -o, o; rtol=rtol1, max_levels=0)) === T

            @test inferred_return_type(() -> integrate2D(f2, low2, up2, x, w, h)) === T
            @test inferred_return_type(() -> adaptive_integrate_2D(typeof(z), f2, low2, up2; rtol=rtol1, max_levels=0)) === T
            @test inferred_return_type(() -> quad(f2, low2, up2; rtol=rtol1, max_levels=0)) === T

            @test inferred_return_type(() -> integrate3D(f3, low3, up3, x, w, h)) === T
            @test inferred_return_type(() -> adaptive_integrate_3D(typeof(z), f3, low3, up3; rtol=rtol2, max_levels=0)) === T
            @test inferred_return_type(() -> quad(f3, low3, up3; rtol=rtol2, max_levels=0)) === T
        end
    end

    @testset "Extended precision type preservation and inference" begin
        f1(x) = exp(x)
        f2(x, y) = x + y
        f3(x, y, z) = x + y + z
        f_cmpl(x, bmx, xma) = inv(sqrt(bmx * xma))
        f_split(x) = inv(sqrt(abs(x)))

        @testset "T=$T" for T in (BigFloat, Double64, Float128)
            x, w, h = tanhsinh(T, 8)
            z = zero(T)
            o = one(T)
            low2 = SVector(zero(T), zero(T))
            up2 = SVector(one(T), one(T))
            low3 = SVector(zero(T), zero(T), zero(T))
            up3 = SVector(one(T), one(T), one(T))
            rtol1 = T === BigFloat ? BigFloat("1e-20") : T(1e-10)
            rtol2 = T === BigFloat ? BigFloat("1e-16") : T(1e-8)

            @test eltype(x) === T
            @test eltype(w) === T
            @test h isa T

            @test typed_call_return_type(integrate1D, Tuple{Type{T}, typeof(f1), Int}) === T
            @test inferred_return_type(() -> integrate1D(f1, x, w, h)) === T
            @test inferred_return_type(() -> integrate1D(f1, z, o, x, w, h)) === T
            @test inferred_return_type(() -> adaptive_integrate_1D(typeof(z), f1, z, o; rtol=rtol1, max_levels=0)) === T
            @test inferred_return_type(() -> quad(f1, z, o; max_levels=0)) === T
            @test inferred_return_type(() -> quad(f1, z, o; rtol=rtol1, max_levels=0)) === T
            @test inferred_return_type(() -> quad_cmpl(f_cmpl, -o, o; rtol=rtol1, max_levels=0)) === T
            @test inferred_return_type(() -> quad_split(f_split, z, -o, o; rtol=rtol1, max_levels=0)) === T

            @test inferred_return_type(() -> integrate2D(f2, low2, up2, x, w, h)) === T
            @test inferred_return_type(() -> adaptive_integrate_2D(typeof(z), f2, low2, up2; rtol=rtol1, max_levels=0)) === T
            @test inferred_return_type(() -> quad(f2, low2, up2; rtol=rtol1, max_levels=0)) === T

            @test inferred_return_type(() -> integrate3D(f3, low3, up3, x, w, h)) === T
            @test inferred_return_type(() -> adaptive_integrate_3D(typeof(z), f3, low3, up3; rtol=rtol2, max_levels=0)) === T
            @test inferred_return_type(() -> quad(f3, low3, up3; rtol=rtol2, max_levels=0)) === T
        end
    end

    @testset "Complement coordinates (quad_cmpl)" begin
        f_cmpl(x, bmx, xma) = 1 / sqrt(bmx * xma)

        # Default interval [-1, 1]
        @test isapprox(quad_cmpl(f_cmpl, -1.0, 1.0; rtol=1e-12), π, atol=1e-12)

        # General interval [a, b]
        a, b = -2.5, 3.25
        @test isapprox(quad_cmpl(f_cmpl, a, b; rtol=1e-12), π, atol=1e-12)
    end

    @testset "Logarithmic singularities T=$T" for T in Types
        exact = 2 * log(T(2)) - 2
        f1(x) = log(1 - x)
        f2(x) = log(1 + x)

        # High precision needs more points for singularities
        n = T == Float32 ? 40 : 120
        x, w, h = tanhsinh(T, n)

        # Pre-computed integration
        @test isapprox(integrate1D(f1, x, w, h), exact, rtol=rtol[T] * 10)
        @test isapprox(integrate1D(f2, x, w, h), exact, rtol=rtol[T] * 10)

        # Adaptive high-level quad
        # Relax tolerance for Float32 slightly or keep standard
        rtol_target = T == Float32 ? 1e-5 : 1e-12
        @test isapprox(quad(f1, -one(T), one(T); rtol=rtol_target), exact, atol=rtol_target * 10)
    end

    @testset "Internal Singularities (quad_split)" begin
        # 1/sqrt(|x|)
        f_abs(x) = 1 / sqrt(abs(x))
        # Exact is 4
        # Relax tolerance slightly for floating point noise near singularity
        @test isapprox(quad_split(f_abs, 0.0, -1.0, 1.0), 4.0, atol=1e-7)
        @test isapprox(quad_split(f_abs, 0.0), 4.0, atol=1e-7)
    end

    include("dispatch.jl")
end

# Aqua.jl quality assurance tests
using Aqua
@testset "Aqua.jl" begin
    Aqua.test_all(FastTanhSinhQuadrature)
end

if VERSION < v"1.11.0-0"
    @info "Skipping JET type stability tests on Julia $VERSION (enabled for Julia 1.11/1.12)."
elseif VERSION >= v"1.13.0-0"
    @info "Skipping JET type stability tests on Julia $VERSION (JET currently not compatible with Julia 1.13+)."
else
    using Pkg
    jet_series = VERSION < v"1.12.0-0" ? "0.9" : "0.11"
    @info "Installing JET $jet_series for Julia $VERSION and running type stability tests."
    Pkg.add(Pkg.PackageSpec(name="JET", version=jet_series))
    include("TypeStabilityTests.jl")
end
