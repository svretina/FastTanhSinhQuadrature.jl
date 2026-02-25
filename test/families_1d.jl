@testset "1D mathematical families" begin
    exact_poly = 9 / 7
    exact_rat = 2 * atan(5.0) / 5
    exact_exp = exp(1.0) - exp(-1.0)
    exact_log = 2 * log(2.0) - 2

    poly(x) = x^6 - 2x^3 + 0.5
    rat(x) = 1 / (1 + 25x^2)
    log_sing(x) = log(1 - x)
    cmpl_sing(x, bmx, xma) = 1 / sqrt(bmx * xma)

    x, w, h = tanhsinh(Float64, 120)
    @test isapprox(integrate1D(poly, x, w, h), exact_poly, atol=1e-12)
    @test isapprox(integrate1D(rat, x, w, h), exact_rat, atol=1e-6)
    @test isapprox(integrate1D(exp, x, w, h), exact_exp, atol=1e-12)
    @test isapprox(integrate1D(log_sing, x, w, h), exact_log, atol=1e-11)

    @test isapprox(integrate1D(Float64, poly, 120), exact_poly, atol=1e-12)
    @test isapprox(integrate1D(poly, 120), exact_poly, atol=1e-12)
    @test isapprox(FastTanhSinhQuadrature.integrate1D_cmpl(Float64, cmpl_sing, 120), π, atol=5e-8)

    @test isapprox(quad(poly, -1.0, 1.0; tol=1e-12), exact_poly, atol=1e-12)
    @test isapprox(quad(rat, -1.0, 1.0; tol=1e-12), exact_rat, atol=1e-12)
    @test isapprox(quad(exp, -1.0, 1.0; tol=1e-12), exact_exp, atol=1e-12)
    @test isapprox(quad(log_sing, -1.0, 1.0; tol=1e-12), exact_log, atol=1e-11)
    @test isapprox(quad_cmpl(cmpl_sing, -1.0, 1.0; tol=1e-12), π, atol=1e-12)
end
