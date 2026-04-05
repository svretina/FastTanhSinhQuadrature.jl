using BenchmarkTools
using FastTanhSinhQuadrature
using StaticArrays

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 0.03

const LOW2 = SVector(-1.0, -1.0)
const UP2 = SVector(1.0, 1.0)
const LOW3 = SVector(-1.0, -1.0, -1.0)
const UP3 = SVector(1.0, 1.0, 1.0)

const F2_REAL = (x, y) -> exp(-(x * x + y * y)) + x * y + x^2 - y
const F2_CHEAP = (x, y) -> inv(2 + x + y * y + 0.5 * x * y)
const F3_REAL = (x, y, z) -> exp(-(x * x + y * y + z * z)) + x * y * z + x - y + z^2
const F3_CHEAP = (x, y, z) -> inv(2 + x + y * y + z * z + 0.5 * x * y * z)

function summary(trial)
    s = BenchmarkTools.median(trial)
    return (time_ns = s.time, allocs = s.allocs, bytes = s.memory)
end

function adaptive_integrate_2D_splitmock(::Type{T}, f::S, low::SVector{2,T}, up::SVector{2,T};
    rtol=nothing, atol::Real=0, max_levels::Int=8,
    warn::Bool=true, cache=nothing) where {T<:Real,S}
    rtol_T, atol_T = FastTanhSinhQuadrature._resolve_tolerances(T; rtol=rtol, atol=atol)
    Δx, x₀ = FastTanhSinhQuadrature._midpoint_radius(low[1], up[1])
    Δy, y₀ = FastTanhSinhQuadrature._midpoint_radius(low[2], up[2])
    cache2d = cache === nothing ? adaptive_cache_2D(T; max_levels=max_levels) :
              FastTanhSinhQuadrature._require_cache_levels(cache, max_levels)
    half = FastTanhSinhQuadrature._half(T)
    h = cache2d.tm * half
    w0 = FastTanhSinhQuadrature._half_pi(T)

    @inline function eval_quadrants(xi, yi, wi, wj)
        dx, dy = Δx * xi, Δy * yi
        return wi * wj * (f(x₀ + dx, y₀ + dy) + f(x₀ - dx, y₀ + dy) +
                          f(x₀ + dx, y₀ - dy) + f(x₀ - dx, y₀ - dy))
    end

    @inline function eval_axes(val, wk)
        dx, dy = Δx * val, Δy * val
        return wk * w0 * (f(x₀ + dx, y₀) + f(x₀ - dx, y₀) +
                          f(x₀, y₀ + dy) + f(x₀, y₀ - dy))
    end

    s_total = (w0 * w0) * f(x₀, y₀)
    @inbounds for i in 1:2
        xi, wi = cache2d.initial_x[i], cache2d.initial_w[i]
        for j in 1:2
            xj, wj = cache2d.initial_x[j], cache2d.initial_w[j]
            s_total += eval_quadrants(xi, xj, wi, wj)
        end
        s_total += eval_axes(xi, wi)
    end

    old_res = Δx * Δy * (h * h) * s_total

    err_est = zero(T)
    for level in 1:max_levels
        h *= half
        s_new = zero(T)
        x_level = cache2d.xs[level]
        w_level = cache2d.ws[level]
        n = length(x_level)

        @inbounds for i in 1:n
            xi, wi = x_level[i], w_level[i]
            if isodd(i)
                for j in 1:n
                    s_new += eval_quadrants(xi, x_level[j], wi, w_level[j])
                end
                s_new += eval_axes(xi, wi)
            else
                for j in 1:2:n
                    s_new += eval_quadrants(xi, x_level[j], wi, w_level[j])
                end
            end
        end

        s_total += s_new
        new_res = Δx * Δy * (h * h) * s_total
        err_est = abs(new_res - old_res)

        if err_est <= FastTanhSinhQuadrature._error_target(new_res, rtol_T, atol_T)
            return new_res
        end
        old_res = new_res
    end
    if warn && max_levels > 0
        @warn "adaptive_integrate_2D_splitmock reached max_levels without meeting the requested tolerance." max_levels estimated_error=err_est
    end
    return old_res
end

function integrate2D_unit_direct_mock(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        zero_T = zero(T)
        w0 = FastTanhSinhQuadrature._half_pi(T)
        w0sq = w0 * w0
        s = w0sq * f(zero_T, zero_T)

        for i in eachindex(x)
            wi = w[i]
            xi = x[i]
            s += wi * w0 * (f(xi, zero_T) + f(-xi, zero_T) + f(zero_T, xi) + f(zero_T, -xi))
        end

        for i in eachindex(x)
            wi = w[i]
            xi = x[i]
            inner_s = zero(T)
            for j in eachindex(x)
                wj = w[j]
                xj = x[j]
                inner_s += wj * (f(-xi, -xj) + f(xi, -xj) + f(-xi, xj) + f(xi, xj))
            end
            s += wi * inner_s
        end
    end
    return (h * h) * s
end

function integrate3D_unit_direct_mock(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        zero_T = zero(T)
        w₀ = FastTanhSinhQuadrature._half_pi(T)
        w₀² = w₀ * w₀
        w₀³ = w₀² * w₀
        total_sum = w₀³ * f(zero_T, zero_T, zero_T)

        for i in eachindex(x)
            wi = w[i]
            xi = x[i]
            axis_sum = (f(xi, zero_T, zero_T) + f(-xi, zero_T, zero_T)) +
                       (f(zero_T, xi, zero_T) + f(zero_T, -xi, zero_T)) +
                       (f(zero_T, zero_T, xi) + f(zero_T, zero_T, -xi))
            total_sum += wi * w₀² * axis_sum
        end

        for i in eachindex(x)
            wi = w[i]
            xi = x[i]
            for j in eachindex(x)
                wj = w[j]
                wiwj = wi * wj
                xj = x[j]

                plane_sum = (f(xi, xj, zero_T) + f(-xi, xj, zero_T) + f(xi, -xj, zero_T) + f(-xi, -xj, zero_T)) +
                            (f(xi, zero_T, xj) + f(-xi, zero_T, xj) + f(xi, zero_T, -xj) + f(-xi, zero_T, -xj)) +
                            (f(zero_T, xi, xj) + f(zero_T, -xi, xj) + f(zero_T, xi, -xj) + f(zero_T, -xi, -xj))

                total_sum += wiwj * w₀ * plane_sum

                octant_sum = zero(T)
                for k in eachindex(x)
                    wk = w[k]
                    xk = x[k]
                    octant_sum += wk * (
                        (f(xi, xj, xk) + f(-xi, xj, xk) + f(xi, -xj, xk) + f(-xi, -xj, xk)) +
                        (f(xi, xj, -xk) + f(-xi, xj, -xk) + f(xi, -xj, -xk) + f(-xi, -xj, -xk))
                    )
                end
                total_sum += wiwj * octant_sum
            end
        end
    end
    return (h * h * h) * total_sum
end

function main()
    x128, w128, h128 = tanhsinh(Float64, 128)
    cache2 = adaptive_cache_2D(Float64; max_levels=9)

    correctness = (
        adaptive2d = adaptive_integrate_2D(Float64, F2_REAL, LOW2, UP2; rtol=1e-6, atol=1e-8, max_levels=9, warn=false, cache=cache2) ==
                     adaptive_integrate_2D_splitmock(Float64, F2_REAL, LOW2, UP2; rtol=1e-6, atol=1e-8, max_levels=9, warn=false, cache=cache2),
        integrate2d_unit = integrate2D(F2_CHEAP, x128, w128, h128) == integrate2D_unit_direct_mock(F2_CHEAP, x128, w128, h128),
        integrate3d_unit = integrate3D(F3_CHEAP, x128, w128, h128) == integrate3D_unit_direct_mock(F3_CHEAP, x128, w128, h128),
    )

    benches = Dict(
        "adaptive2D_prod" => summary(@benchmark adaptive_integrate_2D(Float64, $F2_REAL, $LOW2, $UP2; rtol=1e-6, atol=1e-8, max_levels=9, warn=false, cache=$cache2)),
        "adaptive2D_splitmock" => summary(@benchmark adaptive_integrate_2D_splitmock(Float64, $F2_REAL, $LOW2, $UP2; rtol=1e-6, atol=1e-8, max_levels=9, warn=false, cache=$cache2)),
        "integrate2D_unit_prod" => summary(@benchmark integrate2D($F2_CHEAP, $x128, $w128, $h128)),
        "integrate2D_unit_direct_mock" => summary(@benchmark integrate2D_unit_direct_mock($F2_CHEAP, $x128, $w128, $h128)),
        "integrate3D_unit_prod" => summary(@benchmark integrate3D($F3_CHEAP, $x128, $w128, $h128)),
        "integrate3D_unit_direct_mock" => summary(@benchmark integrate3D_unit_direct_mock($F3_CHEAP, $x128, $w128, $h128)),
    )

    println(correctness)
    println(benches)
end

main()
