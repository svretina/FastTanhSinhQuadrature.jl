using FastTanhSinhQuadrature
using Test
using DoubleFloats
using Random

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

    F = integrate(f, x, w, h)
    G = integrate(g, x, w, h)

    afg(x) = a * f(x) + g(x)
    @test integrate(afg, x, w, h) ≈ a * F + G

    d = T(rand(-9:9)) / 10

    #@test integrate(f, one(T), -one(T), x, w, h) ≈ -F

    F0 = integrate(f, -one(T), a, x, w, h)
    F1 = integrate(f, a, one(T), x, w, h)
    @test isapprox(F0 + F1, F, rtol=rtol[T])
end
