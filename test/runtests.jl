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
    x, w, h = tanhsinh(T, 80)
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
