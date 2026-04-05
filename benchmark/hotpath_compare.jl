using BenchmarkTools
using FastTanhSinhQuadrature
using Printf
using StaticArrays

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 0.05

const LOW2 = SVector(-1.0, -1.0)
const UP2 = SVector(1.0, 1.0)
const LOW3 = SVector(-1.0, -1.0, -1.0)
const UP3 = SVector(1.0, 1.0, 1.0)

const F1 = x -> exp(-x * x) + x^4 - 0.25 * x + sin(x)
const F2 = (x, y) -> exp(-(x * x + y * y)) + x * y + x^2 - y
const F3 = (x, y, z) -> exp(-(x * x + y * y + z * z)) + x * y * z + x - y + z^2

function bench_row(name::AbstractString, trial::BenchmarkTools.Trial)
    stats = median(trial)
    @printf("%-26s %12.2f ns  %4d allocs  %6d B\n",
        name, stats.time, stats.allocs, stats.memory)
end

function main()
    x128, w128, h128 = tanhsinh(Float64, 128)
    cache2 = adaptive_cache_2D(Float64; max_levels=9)
    cache3 = adaptive_cache_3D(Float64; max_levels=7)

    benches = [
        ("tanhsinh(Int)", @benchmark tanhsinh(Float64, 256)),
        ("tanhsinh(Val)", @benchmark tanhsinh(Float64, Val(256))),
        ("integrate1D typed", @benchmark integrate1D(Float64, F1, 256)),
        ("integrate1D_avx", @benchmark integrate1D_avx(F1, $x128, $w128, $h128)),
        ("integrate2D", @benchmark integrate2D(F2, LOW2, UP2, $x128, $w128, $h128)),
        ("integrate2D_avx", @benchmark integrate2D_avx(F2, LOW2, UP2, $x128, $w128, $h128)),
        ("integrate3D", @benchmark integrate3D(F3, LOW3, UP3, $x128, $w128, $h128)),
        ("integrate3D_avx", @benchmark integrate3D_avx(F3, LOW3, UP3, $x128, $w128, $h128)),
        ("adaptive2D", @benchmark adaptive_integrate_2D(Float64, F2, LOW2, UP2; rtol=1e-6, atol=1e-8, max_levels=9, warn=false, cache=$cache2)),
        ("adaptive3D", @benchmark adaptive_integrate_3D(Float64, F3, LOW3, UP3; rtol=1e-6, atol=1e-8, max_levels=7, warn=false, cache=$cache3)),
        ("quad2D", @benchmark quad(F2, LOW2, UP2; rtol=1e-6, atol=1e-8, max_levels=9, cache=$cache2)),
        ("quad3D", @benchmark quad(F3, LOW3, UP3; rtol=1e-6, atol=1e-8, max_levels=7, cache=$cache3)),
    ]

    println("Hot-path benchmarks")
    println("Julia: $(VERSION)")
    println("")
    for (name, trial) in benches
        bench_row(name, trial)
    end
end

main()
