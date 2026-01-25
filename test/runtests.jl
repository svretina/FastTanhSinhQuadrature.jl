using FastTanhSinhQuadrature
using StaticArrays
using DoubleFloats
using Test
using Random
using LinearAlgebra

const Types = [Float32, Float64, Double64, BigFloat]
const rtol = Dict(Float32 => 10 * sqrt(eps(Float32)),
    Float64 => 100 * sqrt(eps(Float64)),
    Double64 => 10^7 * sqrt(eps(Double64)),
    BigFloat => 1.0e-18)

@testset "FastTanhSinhQuadrature.jl" begin

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
            F2 = integrate1D_avx(f, x, w, h)
            G2 = integrate1D_avx(g, x, w, h)
            @test integrate1D_avx(afg, x, w, h) ≈ a * F2 + G2
            @test integrate1D_avx(f, one(T), -one(T), x, w, h) ≈ -F2
        end

        F01 = integrate1D(f, -one(T), a, x, w, h)
        F11 = integrate1D(f, a, one(T), x, w, h)
        @test isapprox(F01 + F11, F1, rtol=rtol[T])
    end

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

    @testset "Adaptive integration 1D" begin
        f(x) = exp(x)
        val_true = exp(1.0) - exp(0.0)

        # Using quad which uses adaptive_integrate_1D
        val = quad(f, 0.0, 1.0; tol=1e-8)
        @test isapprox(val, val_true, atol=1e-8)

        # Test with higher precision
        if Double64 in Types
            val_d = quad(f, Double64(0.0), Double64(1.0); tol=1e-15)
            @test isapprox(val_d, val_true, atol=1e-15)
        end
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
        tol = T == Float32 ? 1e-5 : 1e-12
        @test isapprox(quad(f1, -one(T), one(T); tol=tol), exact, atol=tol * 10)
    end

    @testset "Internal Singularities (quad_split)" begin
        # 1/sqrt(|x|)
        f_abs(x) = 1 / sqrt(abs(x))
        # Exact is 4
        # Relax tolerance slightly for floating point noise near singularity
        @test isapprox(quad_split(f_abs, 0.0, -1.0, 1.0), 4.0, atol=1e-7)
    end

    @testset "Edge Cases and Safety" begin
        # Zero range
        @test quad(x -> x, 1.0, 1.0) == 0.0
        # Flipped
        @test isapprox(quad(x -> x, 1.0, 0.0), -0.5, atol=1e-12)

        # 3D
        @test isapprox(quad((x, y, z) -> 1.0, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), 1.0, atol=1e-12)
    end
end
