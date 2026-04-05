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

# Adaptive integration functions

"""
    adaptive_integrate_1D(::Type{T}, f, a, b; rtol, atol, max_levels::Int=16, warn::Bool=true, cache=nothing)

Adaptive 1D Tanh-Sinh integration over `[a, b]`. Starts with a coarse grid (h ≈ tmax/2) and halves 
the step size at each level. Reuses function evaluations from previous levels by only computing
new (odd-indexed) nodes. Exploits symmetry around the center of the interval.
The error estimate is `abs(I_new - I_old)` between successive refinement levels,
and refinement stops when it is below `max(atol, rtol * abs(I_new))`.
"""
function adaptive_integrate_1D(::Type{T}, f::F, a_T::T, b_T::T,
    rtol_T::T, atol_T::T, max_levels::Int,
    warn::Bool, cache1d::_Adaptive1DCache{T}) where {T<:Real,F}
    Δx, x₀ = _midpoint_radius(a_T, b_T)
    half = _half(T)

    # Initial Grid (Level 0)
    h = cache1d.tm * half

    # Weight at t=0 is π/2
    w0 = T(π) * half
    s_total = w0 * f(x₀)

    # Initial Level evaluations
    @inbounds for i in 1:2
        xk = cache1d.initial_x[i]
        Δxxk = Δx * xk
        s_total += cache1d.initial_w[i] * (f(x₀ + Δxxk) + f(x₀ - Δxxk))
    end

    old_res = Δx * h * s_total

    err_est = zero(T)
    for level in 1:max_levels
        h *= half
        s_new = zero(T)
        x_level = cache1d.xs[level]
        w_level = cache1d.ws[level]
        @inbounds for i in eachindex(x_level)
            xk = x_level[i]
            xp = Δx * xk + x₀
            xm = -Δx * xk + x₀
            s_new += w_level[i] * (f(xp) + f(xm))
        end

        s_total += s_new
        new_res = Δx * h * s_total
        err_est = abs(new_res - old_res)

        if err_est <= _error_target(new_res, rtol_T, atol_T)
            return new_res
        end
        old_res = new_res
    end
    if warn && max_levels > 0
        @warn "adaptive_integrate_1D reached max_levels without meeting the requested tolerance." max_levels estimated_error = err_est target = _error_target(old_res, rtol_T, atol_T) value = old_res rtol = rtol_T atol = atol_T
    end
    return old_res
end

function adaptive_integrate_1D(::Type{T}, f::F, a, b;
    rtol=nothing, atol::Real=0, max_levels::Int=16,
    warn::Bool=true, cache=nothing) where {T<:Real,F}
    a_T, b_T = T(a), T(b)
    rtol_T, atol_T = _resolve_tolerances(T, rtol, atol)
    cache1d = cache === nothing ? adaptive_cache_1D(T; max_levels=max_levels) :
              _require_cache_levels(cache, max_levels)
    return adaptive_integrate_1D(T, f, a_T, b_T, rtol_T, atol_T, max_levels, warn, cache1d)
end

"""
    adaptive_integrate_2D(::Type{T}, f, low::SVector{2,T}, up::SVector{2,T}; rtol, atol, max_levels::Int=8, warn::Bool=true, cache=nothing)

Adaptive 2D Tanh-Sinh integration over a rectangle. Reuses indices by only evaluating new points 
where at least one coordinate corresponds to an odd multiple of the halved step size `h`. 
Exploits 4-way quadrant symmetry and 2-way axis symmetry.
The error estimate is `abs(I_new - I_old)` between successive refinement levels,
and refinement stops when it is below `max(atol, rtol * abs(I_new))`.
"""
function adaptive_integrate_2D(::Type{T}, f::S, low::SVector{2,T}, up::SVector{2,T},
    rtol_T::T, atol_T::T, max_levels::Int,
    warn::Bool, cache2d::_AdaptiveTensorCache{T}) where {T<:Real,S}
    Δx, x₀ = _midpoint_radius(low[1], up[1])
    Δy, y₀ = _midpoint_radius(low[2], up[2])
    half = _half(T)
    h = cache2d.tm * half
    w0 = _half_pi(T)
    w0² = w0 * w0

    # Helper to evaluate symmetric 4 quadrant points
    @inline function eval_quadrants(xi, yi, wi, wj)
        dx, dy = Δx * xi, Δy * yi
        xp, xm = x₀ + dx, x₀ - dx
        yp, ym = y₀ + dy, y₀ - dy
        return wi * wj * (f(xp, yp) + f(xm, yp) + f(xp, ym) + f(xm, ym))
    end

    # Helper to evaluate symmetric axis points
    @inline function eval_axes(val, wk)
        dx, dy = Δx * val, Δy * val
        xp, xm = x₀ + dx, x₀ - dx
        yp, ym = y₀ + dy, y₀ - dy
        return wk * w0 * (f(xp, y₀) + f(xm, y₀) + f(x₀, yp) + f(x₀, ym))
    end

    # Initial Level 0 (h, 2h)
    s_total = w0² * f(x₀, y₀)
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
            for j in 1:n
                # Point is new if at least one index is odd
                (iseven(i) && iseven(j)) && continue
                s_new += eval_quadrants(xi, x_level[j], wi, w_level[j])
            end
            if isodd(i)
                s_new += eval_axes(xi, wi)
            end
        end

        s_total += s_new
        new_res = Δx * Δy * (h * h) * s_total
        err_est = abs(new_res - old_res)

        if err_est <= _error_target(new_res, rtol_T, atol_T)
            return new_res
        end
        old_res = new_res
    end
    if warn && max_levels > 0
        @warn "adaptive_integrate_2D reached max_levels without meeting the requested tolerance." max_levels estimated_error = err_est target = _error_target(old_res, rtol_T, atol_T) value = old_res rtol = rtol_T atol = atol_T
    end
    return old_res
end

function adaptive_integrate_2D(::Type{T}, f::S, low::SVector{2,T}, up::SVector{2,T};
    rtol=nothing, atol::Real=0, max_levels::Int=8,
    warn::Bool=true, cache=nothing) where {T<:Real,S}
    rtol_T, atol_T = _resolve_tolerances(T, rtol, atol)
    cache2d = cache === nothing ? adaptive_cache_2D(T; max_levels=max_levels) :
              _require_cache_levels(cache, max_levels)
    return adaptive_integrate_2D(T, f, low, up, rtol_T, atol_T, max_levels, warn, cache2d)
end

function adaptive_integrate_2D_test(::Type{T}, f::S, low::SVector{2,T}, up::SVector{2,T};
    rtol=nothing, atol::Real=0, max_levels::Int=8,
    warn::Bool=true, cache=nothing) where {T<:Real,S}
    rtol_T, atol_T = _resolve_tolerances(T; rtol=rtol, atol=atol)
    Δx, x₀ = _midpoint_radius(low[1], up[1])
    Δy, y₀ = _midpoint_radius(low[2], up[2])
    cache2d = cache === nothing ? adaptive_cache_2D(T; max_levels=max_levels) :
              _require_cache_levels(cache, max_levels)
    half = _half(T)
    h = cache2d.tm * half
    w0 = _half_pi(T)
    w0² = w0 * w0

    # Helper to evaluate symmetric 4 quadrant points
    @inline function eval_quadrants(xi, yi, wi, wj)
        dx, dy = Δx * xi, Δy * yi
        xp, xm = x₀ + dx, x₀ - dx
        yp, ym = y₀ + dy, y₀ - dy
        return wi * wj * (f(xp, yp) + f(xm, yp) + f(xp, ym) + f(xm, ym))
    end

    # Helper to evaluate symmetric axis points
    @inline function eval_axes(val, wk)
        dx, dy = Δx * val, Δy * val
        xp, xm = x₀ + dx, x₀ - dx
        yp, ym = y₀ + dy, y₀ - dy
        return wk * w0 * (f(xp, y₀) + f(xm, y₀) + f(x₀, yp) + f(x₀, ym))
    end

    # Initial Level 0 (h, 2h)
    s_total = w0² * f(x₀, y₀)
    @inbounds for i in 1:2
        xi, wi = cache2d.initial_x[i], cache2d.initial_w[i]
        @simd for j in 1:2
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

        @inbounds begin
            # odd i: all j contribute, plus axis terms
            for i in 1:2:n
                xi = x_level[i]
                wi = w_level[i]

                @simd for j in 1:n
                    s_new += eval_quadrants(xi, x_level[j], wi, w_level[j])
                end

                s_new += eval_axes(xi, wi)
            end

            # even i: only odd j contribute
            for i in 2:2:n
                xi = x_level[i]
                wi = w_level[i]

                @simd for j in 1:2:n
                    s_new += eval_quadrants(xi, x_level[j], wi, w_level[j])
                end
            end
        end

        s_total += s_new
        new_res = Δx * Δy * (h * h) * s_total
        err_est = abs(new_res - old_res)

        if err_est <= _error_target(new_res, rtol_T, atol_T)
            return new_res
        end
        old_res = new_res
    end
    if warn && max_levels > 0
        @warn "adaptive_integrate_2D reached max_levels without meeting the requested tolerance." max_levels estimated_error = err_est target = _error_target(old_res, rtol_T, atol_T) value = old_res rtol = rtol_T atol = atol_T
    end
    return old_res
end


"""
    adaptive_integrate_3D(::Type{T}, f, low::SVector{3,T}, up::SVector{3,T}; rtol, atol, max_levels::Int=5, warn::Bool=true, cache=nothing)

Adaptive 3D Tanh-Sinh integration over a box. Reuses old points and exploits 8-way octant 
symmetry, 4-way plane symmetry, and 2-way axis symmetry to minimize function evaluations.
The error estimate is `abs(I_new - I_old)` between successive refinement levels,
and refinement stops when it is below `max(atol, rtol * abs(I_new))`.
"""
function adaptive_integrate_3D(::Type{T}, f::S, low::SVector{3,T}, up::SVector{3,T},
    rtol_T::T, atol_T::T, max_levels::Int,
    warn::Bool, cache3d::_AdaptiveTensorCache{T}) where {T<:Real,S}
    Δx, x₀ = _midpoint_radius(low[1], up[1])
    Δy, y₀ = _midpoint_radius(low[2], up[2])
    Δz, z₀ = _midpoint_radius(low[3], up[3])
    half = _half(T)
    h = cache3d.tm * half
    w₀ = _half_pi(T)
    w₀² = w₀ * w₀
    w₀³ = w₀² * w₀

    # Evaluate a single point in the octant (8 reflections)
    @inline function add_octant(vi, vj, vk, wi, wj, wk)
        dx, dy, dz = Δx * vi, Δy * vj, Δz * vk
        xp, xm = x₀ + dx, x₀ - dx
        yp, ym = y₀ + dy, y₀ - dy
        zp, zm = z₀ + dz, z₀ - dz
        w = wi * wj * wk
        return w * (
            (f(xp, yp, zp) + f(xm, yp, zp) +
             f(xp, ym, zp) + f(xm, ym, zp)) +
            (f(xp, yp, zm) + f(xm, yp, zm) +
             f(xp, ym, zm) + f(xm, ym, zm))
        )
    end

    # Evaluate points on the 3 planes (XY, XZ, YZ) - 4 reflections each
    @inline function add_planes(vi, vj, wi, wj)
        dx_i, dy_i, dz_i = Δx * vi, Δy * vi, Δz * vi
        dx_j, dy_j, dz_j = Δx * vj, Δy * vj, Δz * vj
        xp, xm = x₀ + dx_i, x₀ - dx_i
        yp_i, ym_i = y₀ + dy_i, y₀ - dy_i
        yp_j, ym_j = y₀ + dy_j, y₀ - dy_j
        zp_j, zm_j = z₀ + dz_j, z₀ - dz_j
        w = wi * wj * w₀
        return w * (
            (f(xp, yp_j, z₀) + f(xm, yp_j, z₀) + f(xp, ym_j, z₀) + f(xm, ym_j, z₀)) +
            (f(xp, y₀, zp_j) + f(xm, y₀, zp_j) + f(xp, y₀, zm_j) + f(xm, y₀, zm_j)) +
            (f(x₀, yp_i, zp_j) + f(x₀, ym_i, zp_j) + f(x₀, yp_i, zm_j) + f(x₀, ym_i, zm_j))
        )
    end

    # Evaluate points on the 3 axes (X, Y, Z) - 2 reflections each
    @inline function add_axes(vi, wi)
        dx, dy, dz = Δx * vi, Δy * vi, Δz * vi
        xp, xm = x₀ + dx, x₀ - dx
        yp, ym = y₀ + dy, y₀ - dy
        zp, zm = z₀ + dz, z₀ - dz
        w = wi * w₀²
        return w * (
            (f(xp, y₀, z₀) + f(xm, y₀, z₀)) +
            (f(x₀, yp, z₀) + f(x₀, ym, z₀)) +
            (f(x₀, y₀, zp) + f(x₀, y₀, zm))
        )
    end

    # Initial Level 0 (k=1, 2)
    s_total = w₀³ * f(x₀, y₀, z₀)
    @inbounds for i in 1:2
        vi, wi = cache3d.initial_x[i], cache3d.initial_w[i]
        for j in 1:2
            vj, wj = cache3d.initial_x[j], cache3d.initial_w[j]
            for k in 1:2
                vk, wk = cache3d.initial_x[k], cache3d.initial_w[k]
                s_total += add_octant(vi, vj, vk, wi, wj, wk)
            end
            s_total += add_planes(vi, vj, wi, wj)
        end
        s_total += add_axes(vi, wi)
    end

    old_res = Δx * Δy * Δz * (h * h * h) * s_total

    err_est = zero(T)
    for level in 1:max_levels
        h *= half
        s_new = zero(T)
        x_level = cache3d.xs[level]
        w_level = cache3d.ws[level]
        n = length(x_level)

        @inbounds for i in 1:n
            xi, wi = x_level[i], w_level[i]
            if isodd(i)
                for j in 1:n
                    xj, wj = x_level[j], w_level[j]
                    for k in 1:n
                        s_new += add_octant(xi, xj, x_level[k], wi, wj, w_level[k])
                    end
                    s_new += add_planes(xi, xj, wi, wj)
                end
                s_new += add_axes(xi, wi)
            else
                for j in 1:n
                    xj, wj = x_level[j], w_level[j]
                    if isodd(j)
                        for k in 1:n
                            s_new += add_octant(xi, xj, x_level[k], wi, wj, w_level[k])
                        end
                        s_new += add_planes(xi, xj, wi, wj)
                    else
                        for k in 1:2:n
                            s_new += add_octant(xi, xj, x_level[k], wi, wj, w_level[k])
                        end
                    end
                end
            end
        end

        s_total += s_new
        new_res = Δx * Δy * Δz * (h * h * h) * s_total
        err_est = abs(new_res - old_res)

        if err_est <= _error_target(new_res, rtol_T, atol_T)
            return new_res
        end
        old_res = new_res
    end
    if warn && max_levels > 0
        @warn "adaptive_integrate_3D reached max_levels without meeting the requested tolerance." max_levels estimated_error = err_est target = _error_target(old_res, rtol_T, atol_T) value = old_res rtol = rtol_T atol = atol_T
    end
    return old_res
end

function adaptive_integrate_3D(::Type{T}, f::S, low::SVector{3,T}, up::SVector{3,T};
    rtol=nothing, atol::Real=0, max_levels::Int=5,
    warn::Bool=true, cache=nothing) where {T<:Real,S}
    rtol_T, atol_T = _resolve_tolerances(T, rtol, atol)
    cache3d = cache === nothing ? adaptive_cache_3D(T; max_levels=max_levels) :
              _require_cache_levels(cache, max_levels)
    return adaptive_integrate_3D(T, f, low, up, rtol_T, atol_T, max_levels, warn, cache3d)
end

"""
    adaptive_integrate_1D_cmpl(::Type{T}, f, a, b; rtol, atol, max_levels::Int=16, warn::Bool=true, cache=nothing)

Adaptive 1D Tanh-Sinh integration for endpoint-distance-aware integrands.
For interval `[a, b]`, `f` should accept `f(x, b_minus_x, x_minus_a)`,
where `b_minus_x = b - x` and `x_minus_a = x - a`.
For the default interval `[-1, 1]`, this is `f(x, 1-x, 1+x)`.
The error estimate is `abs(I_new - I_old)` between successive refinement levels,
and refinement stops when it is below `max(atol, rtol * abs(I_new))`.
"""
function adaptive_integrate_1D_cmpl(::Type{T}, f::F, a_T::T, b_T::T,
    rtol_T::T, atol_T::T, max_levels::Int,
    warn::Bool, cache1d::_Adaptive1DCache{T}) where {T<:Real,F}
    Δx, x₀ = _midpoint_radius(a_T, b_T)
    half = _half(T)
    one_T = one(T)

    # Complement coordinates remain accurate well beyond t_x_max(T),
    # so use the weight-based window to avoid truncating endpoint tails.
    h = cache1d.tm * half
    w0 = _half_pi(T)
    s_total = w0 * f(x₀, Δx, Δx)

    @inbounds for i in 1:2
        xk = cache1d.initial_x[i]
        ck = cache1d.initial_c[i]
        wk = cache1d.initial_w[i]
        Δxxk = Δx * xk
        Δxck = Δx * ck
        Δx1pxk = Δx * (one_T + xk)
        # f(x, b-x, x-a)
        # At x = x₀ + Δx*xk:
        # b - (x₀ + Δx*xk) = (x₀ + Δx) - (x₀ + Δx*xk) = Δx * (1 - xk) = Δx * ck
        # (x₀ + Δx*xk) - a = (x₀ + Δx*xk) - (x₀ - Δx) = Δx * (1 + xk)
        s_total += wk * (f(x₀ + Δxxk, Δxck, Δx1pxk) +
                         f(x₀ - Δxxk, Δx1pxk, Δxck))
    end

    old_res = Δx * h * s_total

    err_est = zero(T)
    for level in 1:max_levels
        h *= half
        s_new = zero(T)
        x_level = cache1d.xs[level]
        w_level = cache1d.ws[level]
        c_level = cache1d.cs[level]
        @inbounds for i in eachindex(x_level)
            xk = x_level[i]
            ck = c_level[i]
            Δxxk = Δx * xk
            Δxck = Δx * ck
            Δx1pxk = Δx * (one_T + xk)
            s_new += w_level[i] * (f(x₀ + Δxxk, Δxck, Δx1pxk) +
                                   f(x₀ - Δxxk, Δx1pxk, Δxck))
        end

        s_total += s_new
        new_res = Δx * h * s_total
        err_est = abs(new_res - old_res)

        if err_est <= _error_target(new_res, rtol_T, atol_T)
            return new_res
        end
        old_res = new_res
    end
    if warn && max_levels > 0
        @warn "adaptive_integrate_1D_cmpl reached max_levels without meeting the requested tolerance." max_levels estimated_error = err_est target = _error_target(old_res, rtol_T, atol_T) value = old_res rtol = rtol_T atol = atol_T
    end
    return old_res
end

function adaptive_integrate_1D_cmpl(::Type{T}, f::F, a, b;
    rtol=nothing, atol::Real=0, max_levels::Int=16,
    warn::Bool=true, cache=nothing) where {T<:Real,F}
    a_T, b_T = T(a), T(b)
    rtol_T, atol_T = _resolve_tolerances(T, rtol, atol)
    cache1d = cache === nothing ? adaptive_cache_1D(T; max_levels=max_levels, complement=true) :
              _require_cache_levels(cache, max_levels)
    return adaptive_integrate_1D_cmpl(T, f, a_T, b_T, rtol_T, atol_T, max_levels, warn, cache1d)
end
