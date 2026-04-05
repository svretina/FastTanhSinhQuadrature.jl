@testset "2D integration families and overloads" begin
    x, w, h = tanhsinh(Float64, 80)

    # Default-domain overload: integrate over [-1, 1]^2.
    f_default(x, y) = x^2 + y^2 + one(x)
    exact_default = 20 / 3
    val_default = integrate2D(f_default, x, w, h)
    @test isapprox(val_default, exact_default, atol=1e-12)
    run_avx_checks() do
        val_default_avx = integrate2D_avx(f_default, SVector(-1.0, -1.0), SVector(1.0, 1.0), x, w, h)
        @test isapprox(val_default_avx, exact_default, atol=1e-12)
        @test isapprox(val_default_avx, val_default, atol=1e-12)
    end

    # General bounded-domain SVector overload with an analytic reference.
    a, b = -2.0, 3.0
    c, d = -1.0, 2.0
    low = SVector(a, c)
    up = SVector(b, d)
    f_box(x, y) = x^2 + y
    exact_box = ((b^3 - a^3) / 3) * (d - c) + ((d^2 - c^2) / 2) * (b - a)
    val_box = integrate2D(f_box, low, up, x, w, h)
    @test isapprox(val_box, exact_box, atol=1e-11)

    run_avx_checks() do
        val_box_avx = integrate2D_avx(f_box, low, up, x, w, h)
        @test isapprox(val_box_avx, exact_box, atol=1e-11)
        @test isapprox(val_box_avx, val_box, atol=1e-11)
    end

    # AbstractVector overloads (including mixed real bounds with irrational pi).
    @test isapprox(integrate2D((x, y) -> one(x), [0, 0.0], [π, π], x, w, h), π^2, atol=1e-12)
    run_avx_checks() do
        @test isapprox(integrate2D_avx((x, y) -> one(x), [0, 0.0], [π, π], x, w, h), π^2, atol=1e-12)
    end

    # Orientation is encoded by bound ordering in fixed-grid APIs.
    low_flip = SVector(1.0, -1.0)
    up_flip = SVector(0.0, 2.0)
    exact_oriented = (up_flip[1] - low_flip[1]) * (up_flip[2] - low_flip[2]) # constant 1 integrand
    @test isapprox(integrate2D((x, y) -> one(x), low_flip, up_flip, x, w, h), exact_oriented, atol=1e-12)
    run_avx_checks() do
        @test isapprox(integrate2D_avx((x, y) -> one(x), low_flip, up_flip, x, w, h), exact_oriented, atol=1e-12)
    end

    # Length guards for vector-bound convenience overloads.
    @test_throws DimensionMismatch integrate2D((x, y) -> x + y, [0.0], [1.0, 2.0], x, w, h)
    @test_throws DimensionMismatch integrate2D((x, y) -> x + y, [0.0, 0.0], [1.0], x, w, h)
    @test_throws DimensionMismatch integrate2D_avx((x, y) -> x + y, [0.0], [1.0, 2.0], x, w, h)
    @test_throws DimensionMismatch integrate2D_avx((x, y) -> x + y, [0.0, 0.0], [1.0], x, w, h)

    @testset "2D split-domain singular family" begin
        low = SVector(-1.0, -1.0)
        up = SVector(1.0, 1.0)
        c = SVector(0.0, 0.0)
        f_abs2(x, y) = 1 / sqrt(abs(x * y))
        # Separable exact value: (∫_{-1}^1 |x|^{-1/2} dx)^2 = 4^2 = 16.
        @test isapprox(quad_split(f_abs2, c, low, up; rtol=1e-8, max_levels=8), 16.0, atol=5e-7)
    end
end
