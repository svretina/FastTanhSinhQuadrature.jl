# FastTanhSinhQuadrature.jl: Fast and high-precision numerical integration using Tanh-Sinh (Double Exponential) quadrature.
# Copyright (C) 2024-2026 Stamatis Vretinaris
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

using BenchmarkTools
using FastTanhSinhQuadrature
using StaticArrays
using StrideArrays
using Printf

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 0.05

const low2 = SVector(-1.0, -1.0)
const up2 = SVector(1.0, 1.0)
const low3 = SVector(-1.0, -1.0, -1.0)
const up3 = SVector(1.0, 1.0, 1.0)

const f1 = x -> exp(-x * x) + x^4 - 0.25 * x + sin(x)
const f2 = (x, y) -> exp(-(x * x + y * y)) + x * y + x^2 - y
const f3 = (x, y, z) -> exp(-(x * x + y * y + z * z)) + x * y * z + x - y + z^2

function verify_same_results(x_static, w_static, h, x_stride, w_stride)
    checks = (
        ("integrate1D", integrate1D(f1, x_static, w_static, h), integrate1D(f1, x_stride, w_stride, h)),
        ("integrate1D_avx", integrate1D_avx(f1, x_static, w_static, h), integrate1D_avx(f1, x_stride, w_stride, h)),
        ("integrate2D", integrate2D(f2, low2, up2, x_static, w_static, h), integrate2D(f2, low2, up2, x_stride, w_stride, h)),
        ("integrate2D_avx", integrate2D_avx(f2, low2, up2, x_static, w_static, h), integrate2D_avx(f2, low2, up2, x_stride, w_stride, h)),
        ("integrate3D", integrate3D(f3, low3, up3, x_static, w_static, h), integrate3D(f3, low3, up3, x_stride, w_stride, h)),
        ("integrate3D_avx", integrate3D_avx(f3, low3, up3, x_static, w_static, h), integrate3D_avx(f3, low3, up3, x_stride, w_stride, h)),
    )

    for (name, a, b) in checks
        if !(a == b || isapprox(a, b; rtol=1e-14, atol=0.0))
            error("Result mismatch in $name: static=$a stride=$b")
        end
    end
end

function bench_pair(name::AbstractString, b_static::BenchmarkTools.Trial, b_stride::BenchmarkTools.Trial)
    s = median(b_static)
    t = median(b_stride)
    ratio = s.time / t.time

    @printf("%-18s | static %9.2f ns (%d allocs, %d B) | stride %9.2f ns (%d allocs, %d B) | static/stride %.3fx\n",
        name, s.time, s.allocs, s.memory, t.time, t.allocs, t.memory, ratio)
end

function run_case(N::Int)
    x_static, w_static, h = tanhsinh(Float64, Val(N))
    x_stride = StrideArray(collect(x_static))
    w_stride = StrideArray(collect(w_static))

    println("N = $N -> points = $(length(x_static))")
    println("")

    verify_same_results(x_static, w_static, h, x_stride, w_stride)
    println("Result check: PASS (same values within tolerance)")
    println("")

    b1s = @benchmark integrate1D($f1, $x_static, $w_static, $h)
    b1t = @benchmark integrate1D($f1, $x_stride, $w_stride, $h)
    b2s = @benchmark integrate1D_avx($f1, $x_static, $w_static, $h)
    b2t = @benchmark integrate1D_avx($f1, $x_stride, $w_stride, $h)

    b3s = @benchmark integrate2D($f2, $low2, $up2, $x_static, $w_static, $h)
    b3t = @benchmark integrate2D($f2, $low2, $up2, $x_stride, $w_stride, $h)
    b4s = @benchmark integrate2D_avx($f2, $low2, $up2, $x_static, $w_static, $h)
    b4t = @benchmark integrate2D_avx($f2, $low2, $up2, $x_stride, $w_stride, $h)

    b5s = @benchmark integrate3D($f3, $low3, $up3, $x_static, $w_static, $h)
    b5t = @benchmark integrate3D($f3, $low3, $up3, $x_stride, $w_stride, $h)
    b6s = @benchmark integrate3D_avx($f3, $low3, $up3, $x_static, $w_static, $h)
    b6t = @benchmark integrate3D_avx($f3, $low3, $up3, $x_stride, $w_stride, $h)

    println("Median runtimes from @benchmark:")
    bench_pair("integrate1D", b1s, b1t)
    bench_pair("integrate1D_avx", b2s, b2t)
    bench_pair("integrate2D", b3s, b3t)
    bench_pair("integrate2D_avx", b4s, b4t)
    bench_pair("integrate3D", b5s, b5t)
    bench_pair("integrate3D_avx", b6s, b6t)
end

function main()
    levels = isempty(ARGS) ? [64] : [parse(Int, a) for a in ARGS]

    println("StaticArrays vs StrideArrays benchmark")
    println("Julia: $(VERSION)")
    println("")

    for (idx, N) in pairs(levels)
        run_case(N)
        if idx != length(levels)
            println("")
            println("-"^110)
            println("")
        end
    end
end

main()
