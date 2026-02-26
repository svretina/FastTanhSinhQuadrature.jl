@testset "2D integration families and overloads" begin
    x, w, h = tanhsinh(Float64, 80)

    # Default-domain overload: integrate over [-1, 1]^2.
    f_default(x, y) = x^2 + y^2 + one(x)
    exact_default = 20 / 3
    @test isapprox(integrate2D(f_default, x, w, h), exact_default, atol=1e-12)

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
end
