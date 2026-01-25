using FastTanhSinhQuadrature
using FastGaussQuadrature
using BenchmarkTools
using Printf
using LinearAlgebra

println("| Function | Domain | Points | TS (ns) | TS SIMD (ns) | GQ (ns) | Ratio (TS/GQ) | Ratio (TS SIMD/GQ) |")
println("| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |")

const f1 = exp
const f2 = x -> sin(x)^2
const f3 = x -> 1.0 / (1.0 + 25.0 * x^2)
const f4 = x -> sqrt(1.0 - x^2)
const f5 = x -> x^2
const f6 = x -> log(1.0 - x)
const f7 = x -> x^3 # Odd
const f8 = x -> x^3 + x^2 + x + 1.0 # No symmetry

# All tests are on [-1, 1] for now
const domain = "[-1, 1]"

function run_benchmark()
    functions = [
        ("exp(x)", f1),
        ("sin(x)^2", f2),
        ("1/(1+25x^2)", f3),
        ("sqrt(1-x^2)", f4),
        ("x^2", f5),
        ("log(1-x)", f6),
        ("x^3", f7),
        ("x^3+x^2+x+1", f8)
    ]

    for (name, func) in functions
        for n_target in [10, 100, 1000]
            # We will use n_target directly
            level = n_target
            x_ts, w_ts, h_ts = tanhsinh(Float64, level)
            n_actual = length(x_ts)

            # FastGaussQuadrature - match the exact number of points
            x_gq, w_gq = gausslegendre(n_actual)

            # Benchmark TS
            b_ts = @belapsed integrate1D($func, $x_ts, $w_ts, $h_ts)

            # Benchmark TS SIMD
            # Explicitly use the exported integrate1D_avx
            b_ts_avx = @belapsed integrate1D_avx($func, $x_ts, $w_ts, $h_ts)

            # Benchmark GQ
            b_gq = @belapsed dot($w_gq, $func.($x_gq))

            ratio = b_ts / b_gq
            ratio_avx = b_ts_avx / b_gq
            @printf("| %s | %s | %d | %.2f | %.2f | %.2f | %.2f | %.2f |\n", name, domain, n_actual, b_ts * 1e9, b_ts_avx * 1e9, b_gq * 1e9, ratio, ratio_avx)
        end
    end
end