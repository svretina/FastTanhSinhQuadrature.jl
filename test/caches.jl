@testset "Adaptive cache APIs" begin
    cache1 = adaptive_cache_1D(Float64; max_levels=3)
    cache1_again = adaptive_cache_1D(Float64; max_levels=3)
    cache1_cmpl = adaptive_cache_1D(Float64; max_levels=3, complement=true)
    cache2 = adaptive_cache_2D(Float64; max_levels=3)
    cache3 = adaptive_cache_3D(Float64; max_levels=2)

    @test cache1 === cache1_again
    @test cache1 !== cache1_cmpl
    @test length(cache1.xs) == 3
    @test length(cache1.ws) == 3
    @test length(cache1.cs) == 3
    @test length(cache2.xs) == 3
    @test length(cache2.ws) == 3
    @test length(cache3.xs) == 2
    @test length(cache3.ws) == 2
    @test cache1.tm == FastTanhSinhQuadrature.tmax(Float64)
    @test cache1_cmpl.tm == FastTanhSinhQuadrature.t_w_max(Float64, 1)

    f1(x) = exp(x)
    f2(x, y) = exp(x + y)
    f3(x, y, z) = exp(x + y + z)

    @test_throws ArgumentError adaptive_integrate_1D(Float64, f1, -1.0, 1.0;
        cache=cache1, max_levels=4)
    @test_throws ArgumentError adaptive_integrate_2D(Float64, f2, SVector(-1.0, -1.0), SVector(1.0, 1.0);
        cache=cache2, max_levels=4)
    @test_throws ArgumentError adaptive_integrate_3D(Float64, f3, SVector(-1.0, -1.0, -1.0), SVector(1.0, 1.0, 1.0);
        cache=cache3, max_levels=3)
end
