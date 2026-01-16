using FastTanhSinhQuadrature
using Test
using DoubleFloats
using Random
using StaticArrays

const Types = [Float32, Float64, Double64, BigFloat]
const rtol = Dict(Float32 => 10 * sqrt(eps(Float32)),
    Float64 => 100 * sqrt(eps(Float64)),
    Double64 => 10^7 * sqrt(eps(Double64)),
    BigFloat => 1.0e-18)

@testset "Basic integration T=$T" for T in Types
    f0(x) = 1
    f1(x) = x + 1
    f2(x) = 3 * x^2
    x, w, h = tanhsinh(T, 40)
    @test isapprox(integrate(f0, x, w, h), T(2), rtol=rtol[T])
    @test isapprox(integrate(f1, x, w, h), T(2), rtol=rtol[T])
    @test isapprox(integrate(f2, x, w, h), T(2), rtol=rtol[T])
end

Random.seed!(0)
@testset "Integral bounds T=$T" for T in Types

    f0(x) = 1
    f1(x) = x + 1
    f2(x) = 3 * x^2

    xmin = -1 + T(rand(-5:5)) / 10
    xmax = +1 + T(rand(-5:5)) / 10

    x, w, h = tanhsinh(T, 40)
    @test isapprox(integrate(f0, xmin, xmax, x, w, h), xmax - xmin, rtol=rtol[T])
    @test isapprox(integrate(f1, xmin, xmax, x, w, h), (xmax^2 - xmin^2) / 2 + xmax - xmin, rtol=rtol[T])
    @test isapprox(integrate(f2, xmin, xmax, x, w, h), (xmax^3 - xmin^3), rtol=rtol[T])
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

    x, w, h = tanhsinh(T, 40)
    f(x) = b0 + b1 * x + b2 * x^2
    g(x) = c0 + c1 * x + c2 * x^2

    F1 = integrate(f, x, w, h)
    F2 = integrate_avx(f, x, w, h)
    G1 = integrate(g, x, w, h)
    G2 = integrate_avx(g, x, w, h)

    afg(x) = a * f(x) + g(x)
    @test integrate(afg, x, w, h) ≈ a * F1 + G1
    @test integrate_avx(afg, x, w, h) ≈ a * F2 + G2


    # @test integrate(f, one(T), -one(T), x, w, h) ≈ -F1
    @test integrate_avx(f, one(T), -one(T), x, w, h) ≈ -F2

    F01 = integrate(f, -one(T), a, x, w, h)
    F11 = integrate(f, a, one(T), x, w, h)
    F02 = integrate_avx(f, -one(T), a, x, w, h)
    F12 = integrate_avx(f, a, one(T), x, w, h)
    @test isapprox(F01 + F11, F1, rtol=rtol[T])
    @test isapprox(F02 + F12, F2, rtol=rtol[T])
end

@testset "2D polynomials for [-1, 1], T=$T" for T in Types
    # n=80 is too large for Float32
    n = T == Float32 ? 40 : 80
    x, w, h = tanhsinh(T, n)
    ψ(x, y) = one(T)
    f(x, y) = x * y
    g(x, y) = x^2 * y^2

    low = SVector{2,T}(-1.0, -1.0)
    up = SVector{2,T}(1.0, 1.0)

    Ψ = integrate(ψ, low, up, x, w, h)
    Ψ2 = integrate(ψ, up, low, x, w, h)
    F = integrate(f, low, up, x, w, h)
    G = integrate(g, low, up, x, w, h)
    @test Ψ ≈ 4one(T)
    @test Ψ2 ≈ 4one(T)
    @test F ≈ zero(T)
    @test G ≈ T(4) / T(9)
end

@testset "Adaptive integration" begin
    # Test function with known integral
    f(x) = exp(x)
    val_true = exp(1.0) - exp(0.0)
    
    val = adaptive_integrate(f, 0.0, 1.0, tol=1e-8)
    @test isapprox(val, val_true, atol=1e-8)
    
    # Test with higher precision if available
    if Double64 in Types
        val_d = adaptive_integrate(Double64, f, 0.0, 1.0, tol=1e-15)
        @test isapprox(val_d, val_true, atol=1e-15)
    end
end

@testset "Logarithmic singularities T=$T" for T in Types
    # Target value: 2*log(2) - 2 for integral from -1 to 1
    exact = 2 * log(T(2)) - 2
    
    f1(x) = log(1 - x)
    f2(x) = log(1 + x)
    
    # Use sufficient points for convergence on singularities
    # n=80 is robust for high precision, 40 for Float32
    n = T == Float32 ? 40 : 80
    x, w, h = tanhsinh(T, n)
    
    # Tanh-Sinh handles endpoint singularities naturally
    @test isapprox(integrate(f1, x, w, h), exact, rtol=rtol[T])
    @test isapprox(integrate(f2, x, w, h), exact, rtol=rtol[T])
end
