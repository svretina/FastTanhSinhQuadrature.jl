@testset "Core helper identities and constructors" begin
    t = 1.75
    x_t = FastTanhSinhQuadrature.ordinate(t)
    x_mt = FastTanhSinhQuadrature.ordinate(-t)
    @test isapprox(x_mt, -x_t, atol=0.0, rtol=10 * eps(Float64))

    w_t = FastTanhSinhQuadrature.weight(t)
    w_mt = FastTanhSinhQuadrature.weight(-t)
    @test isapprox(w_mt, w_t, atol=0.0, rtol=10 * eps(Float64))

    c_t = FastTanhSinhQuadrature.ordinate_complement(t)
    @test isapprox(c_t, 1 - abs(x_t), atol=1e-16, rtol=1e-12)

    @test isapprox(FastTanhSinhQuadrature.inv_ordinate(x_t), t, atol=1e-13)

    # Branch guards near saturation/underflow
    @test FastTanhSinhQuadrature.ordinate(50.0) == prevfloat(1.0)
    @test FastTanhSinhQuadrature.ordinate(-50.0) == -prevfloat(1.0)
    @test FastTanhSinhQuadrature.weight(8.0) == 0.0

    tx = FastTanhSinhQuadrature.t_x_max(Float64)
    tw = FastTanhSinhQuadrature.t_w_max(Float64, 3)
    @test FastTanhSinhQuadrature.tmax(Float64, 3) == min(tx, tw)

    xodd, wodd, hodd = tanhsinh(Float64, 79)
    @test length(xodd) == 39
    @test length(wodd) == 39
    @test hodd > 0

    xval, wval, hval = tanhsinh(Float64, Val(11))
    @test xval isa SVector{5,Float64}
    @test wval isa SVector{5,Float64}
    @test hval > 0

    xval_even, wval_even, hval_even = tanhsinh(Float64, Val(10))
    @test xval_even isa SVector{5,Float64}
    @test wval_even isa SVector{5,Float64}
    @test hval_even > 0

    xdef, wdef, hdef = tanhsinh(80)
    x64, w64, h64 = tanhsinh(Float64, 80)
    @test xdef == x64
    @test wdef == w64
    @test hdef == h64

    xopt, wopt, hopt = FastTanhSinhQuadrature.tanhsinh_opt(Float64, 80)
    @test !isempty(xopt)
    @test length(xopt) == length(wopt)
    @test hopt > 0
    xopt_odd, wopt_odd, hopt_odd = FastTanhSinhQuadrature.tanhsinh_opt(Float64, 81)
    @test !isempty(xopt_odd)
    @test length(xopt_odd) == length(wopt_odd)
    @test hopt_odd > 0
    @test FastTanhSinhQuadrature.hopt(Float64, 80) > 0
    @test FastTanhSinhQuadrature.hmax(Float64, 80) > 0
    @test FastTanhSinhQuadrature.hmax(Float64, 81) > 0
end
