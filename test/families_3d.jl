@testset "3D fixed-grid and AVX on [-1,1], T=$T" for T in (Float32, Float64)
    n = T == Float32 ? 30 : 60
    x, w, h = tanhsinh(T, n)
    ψ = integrate3D(f3_const, x, w, h)
    odd = integrate3D(f3_xyz, x, w, h)
    even = integrate3D(f3_x2y2z2, x, w, h)

    @test isapprox(ψ, T(8), atol=T == Float32 ? T(2e-5) : T(1e-12))
    @test isapprox(odd, zero(T), atol=T == Float32 ? T(2e-6) : T(1e-14))
    @test isapprox(even, T(8) / T(27), atol=T == Float32 ? T(2e-6) : T(1e-13))

    run_avx_checks() do
        ψ_avx = integrate3D_avx(f3_const, x, w, h)
        odd_avx = integrate3D_avx(f3_xyz, x, w, h)
        even_avx = integrate3D_avx(f3_x2y2z2, x, w, h)
        @test isapprox(ψ_avx, ψ, atol=T == Float32 ? T(2e-5) : T(1e-12))
        @test isapprox(odd_avx, odd, atol=T == Float32 ? T(2e-6) : T(1e-14))
        @test isapprox(even_avx, even, atol=T == Float32 ? T(2e-6) : T(1e-13))
    end
end

@testset "3D anisotropic box regression (fixed-grid + AVX)" begin
    T = Float64
    low = SVector{3,T}(-2.0, -1.0, 0.0)
    up = SVector{3,T}(3.0, 2.0, 4.0)
    x, w, h = tanhsinh(T, 80)

    f_poly(x, y, z) = x^2 + 2y + 3z^2
    exact_poly =
        integral_monomial_3d(low, up, 2, 0, 0) +
        2 * integral_monomial_3d(low, up, 0, 1, 0) +
        3 * integral_monomial_3d(low, up, 0, 0, 2)

    val = integrate3D(f_poly, low, up, x, w, h)
    @test isapprox(val, exact_poly, atol=1e-10)

    run_avx_checks() do
        val_avx = integrate3D_avx(f_poly, low, up, x, w, h)
        @test isapprox(val_avx, exact_poly, atol=1e-10)
        @test isapprox(val_avx, val, atol=1e-10)
    end

    # Simple monomial catch for axis/plane bookkeeping: ∫ z² dV on anisotropic box.
    f_z2(x, y, z) = z^2
    exact_z2 = integral_monomial_3d(low, up, 0, 0, 2)
    @test isapprox(integrate3D(f_z2, low, up, x, w, h), exact_z2, atol=1e-10)
end

@testset "3D mathematical families (adaptive quad)" begin
    low = SVector(-1.0, -1.0, -1.0)
    up = SVector(1.0, 1.0, 1.0)

    exact_rat = (2 * atan(5.0) / 5) * (2 * atan(4.0) / 4) * (2 * atan(3.0) / 3)
    exact_exp = (exp(1.0) - exp(-1.0))^3
    exact_log_1d = 2 * log(2.0) - 2
    exact_log = exact_log_1d^3
    exact_sing = π^3

    rat3(x, y, z) = 1 / ((1 + 25x^2) * (1 + 16y^2) * (1 + 9z^2))
    exp3(x, y, z) = exp(x + y + z)
    log3(x, y, z) = log(1 - x) * log(1 - y) * log(1 - z)
    sing3(x, y, z) = 1 / sqrt((1 - x^2) * (1 - y^2) * (1 - z^2))

    @test isapprox(quad(rat3, low, up; tol=1e-9, max_levels=8), exact_rat, atol=1e-9)
    @test isapprox(quad(exp3, low, up; tol=1e-9, max_levels=8), exact_exp, atol=1e-9)
    @test isapprox(quad(log3, low, up; tol=1e-8, max_levels=8), exact_log, atol=5e-9)
    @test isapprox(quad(sing3, low, up; tol=1e-8, max_levels=8), exact_sing, atol=2e-7)
end

@testset "3D split-domain singular family" begin
    low = SVector(-1.0, -1.0, -1.0)
    up = SVector(1.0, 1.0, 1.0)
    c = SVector(0.0, 0.0, 0.0)
    f_abs3(x, y, z) = 1 / sqrt(abs(x * y * z))
    # Separable exact value: (∫_{-1}^1 |x|^{-1/2} dx)^3 = 4^3 = 64.
    @test isapprox(quad_split(f_abs3, c, low, up; tol=1e-6, max_levels=6), 64.0, atol=3e-6)
end
