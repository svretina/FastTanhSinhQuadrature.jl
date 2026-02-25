@testset "quad_split polynomial consistency in 2D" begin
    low = SVector(-1.0, -1.0)
    up = SVector(1.0, 1.0)
    c = SVector(0.0, 0.0)
    f2(x, y) = x^2 + y^2
    # ∫∫ (x²+y²) dx dy on [-1,1]² = 8/3.
    @test isapprox(quad_split(f2, c, low, up; tol=1e-9, max_levels=8), 8 / 3, atol=1e-9)
end

@testset "Adaptive base-level path (max_levels=0)" begin
    low2 = SVector(-1.0, -1.0)
    up2 = SVector(1.0, 1.0)
    low3 = SVector(-1.0, -1.0, -1.0)
    up3 = SVector(1.0, 1.0, 1.0)

    val1 = adaptive_integrate_1D(Float64, exp, 0.0, 1.0; max_levels=0)
    val2 = adaptive_integrate_2D(Float64, f2_const, low2, up2; max_levels=0)
    val3 = adaptive_integrate_3D(Float64, f3_const, low3, up3; max_levels=0)
    valc = adaptive_integrate_1D_cmpl(Float64, (x, bmx, xma) -> 1 / sqrt(bmx * xma), -1.0, 1.0; max_levels=0)

    @test isfinite(val1) && val1 > 0
    @test isfinite(val2) && val2 > 0
    @test isfinite(val3) && val3 > 0
    @test isfinite(valc) && valc > 0
end

@testset "Edge Cases and Safety" begin
    # Zero range
    @test quad(x -> x, 1.0, 1.0) == 0.0
    # Flipped
    @test isapprox(quad(x -> x, 1.0, 0.0), -0.5, atol=1e-12)

    # Pre-computed integration should accept irrational bounds (e.g. pi)
    x, w, h = tanhsinh(Float64, 80)
    @test isapprox(integrate1D(x -> sin(x)^2, 0.0, π, x, w, h), π / 2, atol=1e-12)
    @test isapprox(integrate1D_avx(x -> sin(x)^2, 0.0, π, x, w, h), π / 2, atol=1e-12)
    @test isapprox(integrate2D((x, y) -> 1.0, SVector(0.0, 0.0), SVector(π, π), x, w, h), π^2, atol=1e-12)
    @test isapprox(integrate2D_avx((x, y) -> 1.0, SVector(0.0, 0.0), SVector(π, π), x, w, h), π^2, atol=1e-12)
    @test isapprox(integrate2D((x, y) -> 1.0, [0.0, 0.0], [π, π], x, w, h), π^2, atol=1e-12)
    @test isapprox(integrate3D((x, y, z) -> 1.0, SVector(0.0, 0.0, 0.0), SVector(π, π, π), x, w, h), π^3, atol=1e-11)
    @test isapprox(integrate3D_avx((x, y, z) -> 1.0, SVector(0.0, 0.0, 0.0), SVector(π, π, π), x, w, h), π^3, atol=1e-11)
    @test isapprox(integrate3D((x, y, z) -> 1.0, [0.0, 0.0, 0.0], [π, π, π], x, w, h), π^3, atol=1e-11)

    # High-level interfaces should also accept mixed real bound types
    @test isapprox(quad(x -> x, 0.0, π), π^2 / 2, atol=1e-12)
    @test isapprox(quad((x, y) -> 1.0, [0.0, 0.0], [π, π]), π^2, atol=1e-11)
    @test isapprox(quad_split(x -> 1 / sqrt(abs(x)), 0, -1, 1.0), 4.0, atol=1e-7)

    # 2D/3D through AbstractVector dispatch
    @test isapprox(quad((x, y) -> x + y, [0.0, 0.0], [1.0, 2.0]), 3.0, atol=1e-10)
    @test isapprox(quad((x, y, z) -> 1.0, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]), 1.0, atol=1e-12)
    @test_throws ErrorException quad((args...) -> sum(args), zeros(4), ones(4))

    # Any degenerate dimension should return zero.
    @test quad((x, y) -> x + y, SVector(0.0, 1.0), SVector(0.0, 2.0)) == 0.0
    @test quad((x, y, z) -> x + y + z, SVector(0.0, 1.0, -2.0), SVector(0.0, 2.0, 4.0)) == 0.0

    # Odd number of flipped bounds => negative orientation.
    @test isapprox(quad((x, y) -> 1.0, SVector(1.0, -1.0), SVector(0.0, 2.0)), -3.0, atol=1e-12)
    @test isapprox(quad((x, y, z) -> 1.0, SVector(1.0, -1.0, -1.0), SVector(0.0, 2.0, 3.0)), -12.0, atol=1e-12)
    @test isapprox(quad((x, y, z) -> 1.0, SVector(0.0, 2.0, 3.0), SVector(1.0, -1.0, -1.0)), 12.0, atol=1e-12)

    f_cmpl(x, bmx, xma) = 1 / sqrt(bmx * xma)
    @test quad_cmpl(f_cmpl, 1.0, 1.0) == 0.0
    @test isapprox(quad_cmpl(f_cmpl, 1.0, -1.0; tol=1e-12), -π, atol=1e-12)
end
