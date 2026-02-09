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
    adaptive_integrate_1D(::Type{T}, f::Function, a, b; tol::Real=1e-12, max_levels::Int=10)

Adaptive 1D Tanh-Sinh integration over `[a, b]`. Starts with a coarse grid (h ≈ tmax/2) and halves 
the step size at each level. Reuses function evaluations from previous levels by only computing
new (odd-indexed) nodes. Exploits symmetry around the center of the interval.
"""
function adaptive_integrate_1D(::Type{T}, f::Function, a, b;
    tol::Real=1e-12, max_levels::Int=10) where {T<:Real}
    a_T, b_T = T(a), T(b)
    Δx = 0.5 * (b_T - a_T)
    x₀ = 0.5 * (b_T + a_T)

    # Initial Grid (Level 0)
    tm = tmax(T)
    h = tm / 2

    # Weight at t=0 is π/2
    w0 = T(π) / 2
    s_total = w0 * f(x₀)

    # Initial Level evaluations
    for k in 1:2
        tk = k * h
        s_total += weight(tk) * (f(x₀ + Δx * ordinate(tk)) + f(x₀ - Δx * ordinate(tk)))
    end

    old_res = Δx * h * s_total

    # Pre-allocate caches for ordinates and weights to avoid $O(N)$ transcendental calls
    # Level 1 will add points at 1h, 3h... but h is halved.
    # We can pre-calculate the points needed for the next level.

    for level in 1:max_levels
        h /= 2
        s_new = zero(T)

        # New points are at ODD multiples of the new h: 1h, 3h, 5h...
        # We use a simple loop here. Performance for 1D is usually dominated by f.
        # But we can still optimize the weight/ordinate calls.
        k = 1
        while true
            tk = k * h
            tk > tm && break

            # Pre-calculation optimization: 
            # In a real high-performance scenario, we'd pre-generate these nodes.
            wk = weight(tk)
            xk = ordinate(tk)
            s_new += wk * (f(x₀ + Δx * xk) + f(x₀ - Δx * xk))
            k += 2
        end

        s_total += s_new
        new_res = Δx * h * s_total

        if abs(new_res - old_res) < tol * max(one(T), abs(new_res))
            return new_res
        end
        old_res = new_res
    end
    return old_res
end

"""
    adaptive_integrate_2D(::Type{T}, f::Function, low::SVector{2,T}, up::SVector{2,T}; tol::Real=1e-10, max_levels::Int=8)

Adaptive 2D Tanh-Sinh integration over a rectangle. Reuses indices by only evaluating new points 
where at least one coordinate corresponds to an odd multiple of the halved step size `h`. 
Exploits 4-way quadrant symmetry and 2-way axis symmetry.
"""
function adaptive_integrate_2D(::Type{T}, f::S, low::SVector{2,T}, up::SVector{2,T};
    tol::Real=1e-10, max_levels::Int=8) where {T<:Real,S}
    Δx = 0.5 * (up[1] - low[1])
    Δy = 0.5 * (up[2] - low[2])
    x₀ = 0.5 * (up[1] + low[1])
    y₀ = 0.5 * (up[2] + low[2])
    tm = tmax(T, 2)
    h = tm / 2
    w0 = T(π) / 2

    # Helper to evaluate symmetric 4 quadrant points
    @inline function eval_quadrants(xi, yi, wi, wj)
        dx, dy = Δx * xi, Δy * yi
        return wi * wj * (f(x₀ + dx, y₀ + dy) + f(x₀ - dx, y₀ + dy) +
                          f(x₀ + dx, y₀ - dy) + f(x₀ - dx, y₀ - dy))
    end

    # Helper to evaluate symmetric axis points
    @inline function eval_axes(val, wk)
        dx, dy = Δx * val, Δy * val
        return wk * w0 * (f(x₀ + dx, y₀) + f(x₀ - dx, y₀) +
                          f(x₀, y₀ + dy) + f(x₀, y₀ - dy))
    end

    # Initial Level 0 (h, 2h)
    s_total = (w0^2) * f(x₀, y₀)
    for i in 1:2
        ti = i * h
        xi, wi = ordinate(ti), weight(ti)
        for j in 1:2
            tj = j * h
            xj, wj = ordinate(tj), weight(tj)
            s_total += eval_quadrants(xi, xj, wi, wj)
        end
        s_total += eval_axes(xi, wi)
    end

    old_res = Δx * Δy * h^2 * s_total

    # Cache for nodes and weights to avoid repeated transcendental calls
    cache_x = T[]
    cache_w = T[]

    for level in 1:max_levels
        h /= 2
        s_new = zero(T)
        max_k = floor(Int, tm / h)

        # Populate cache for this level
        resize!(cache_x, max_k)
        resize!(cache_w, max_k)
        for i in 1:max_k
            ti = i * h
            cache_x[i] = ordinate(ti)
            cache_w[i] = weight(ti)
        end

        for i in 1:max_k
            xi, wi = cache_x[i], cache_w[i]
            for j in 1:max_k
                # Point is new if at least one index is odd
                (iseven(i) && iseven(j)) && continue
                s_new += eval_quadrants(xi, cache_x[j], wi, cache_w[j])
            end
            if isodd(i)
                s_new += eval_axes(xi, wi)
            end
        end

        s_total += s_new
        new_res = Δx * Δy * h^2 * s_total

        if abs(new_res - old_res) < tol * max(one(T), abs(new_res))
            return new_res
        end
        old_res = new_res
    end
    return old_res
end

"""
    adaptive_integrate_3D(::Type{T}, f::Function, low::SVector{3,T}, up::SVector{3,T}; tol::Real=1e-8, max_levels::Int=5)

Adaptive 3D Tanh-Sinh integration over a box. Reuses old points and exploits 8-way octant 
symmetry, 4-way plane symmetry, and 2-way axis symmetry to minimize function evaluations.
"""
function adaptive_integrate_3D(::Type{T}, f::S, low::SVector{3,T}, up::SVector{3,T};
    tol::Real=1e-8, max_levels::Int=5) where {T<:Real,S}
    Δx = 0.5 * (up[1] - low[1])
    Δy = 0.5 * (up[2] - low[2])
    Δz = 0.5 * (up[3] - low[3])
    x₀ = 0.5 * (up[1] + low[1])
    y₀ = 0.5 * (up[2] + low[2])
    z₀ = 0.5 * (up[3] + low[3])
    tm = tmax(T, 3)
    h = tm / 2
    w₀ = T(π) / 2

    # Evaluate a single point in the octant (8 reflections)
    @inline function add_octant(vi, vj, vk, wi, wj, wk)
        dx, dy, dz = Δx * vi, Δy * vj, Δz * vk
        w = wi * wj * wk
        return w * (
            (f(x₀ + dx, y₀ + dy, z₀ + dz) + f(x₀ - dx, y₀ + dy, z₀ + dz) +
             f(x₀ + dx, y₀ - dy, z₀ + dz) + f(x₀ - dx, y₀ - dy, z₀ + dz)) +
            (f(x₀ + dx, y₀ + dy, z₀ - dz) + f(x₀ - dx, y₀ + dy, z₀ - dz) +
             f(x₀ + dx, y₀ - dy, z₀ - dz) + f(x₀ - dx, y₀ - dy, z₀ - dz))
        )
    end

    # Evaluate points on the 3 planes (XY, XZ, YZ) - 4 reflections each
    @inline function add_planes(vi, vj, wi, wj)
        dx_i, dy_i, dz_i = Δx * vi, Δy * vi, Δz * vi
        dx_j, dy_j, dz_j = Δx * vj, Δy * vj, Δz * vj
        w = wi * wj * w₀
        return w * (
            (f(x₀ + dx_i, y₀ + dy_j, z₀) + f(x₀ - dx_i, y₀ + dy_j, z₀) + f(x₀ + dx_i, y₀ - dy_j, z₀) + f(x₀ - dx_i, y₀ - dy_j, z₀)) +
            (f(x₀ + dx_i, y₀, z₀ + dz_j) + f(x₀ - dx_i, y₀, z₀ + dz_j) + f(x₀ + dx_i, y₀, z₀ - dz_j) + f(x₀ - dx_i, y₀, z₀ - dz_j)) +
            (f(x₀, y₀ + dy_i, z₀ + dz_j) + f(x₀, y₀ - dy_i, z₀ + dz_j) + f(x₀, y₀ + dy_i, z₀ - dz_j) + f(x₀, y₀ - dy_i, z₀ - dz_j))
        )
    end

    # Evaluate points on the 3 axes (X, Y, Z) - 2 reflections each
    @inline function add_axes(vi, wi)
        dx, dy, dz = Δx * vi, Δy * vi, Δz * vi
        w = wi * w₀^2
        return w * (
            (f(x₀ + dx, y₀, z₀) + f(x₀ - dx, y₀, z₀)) +
            (f(x₀, y₀ + dy, z₀) + f(x₀, y₀ - dy, z₀)) +
            (f(x₀, y₀, z₀ + dz) + f(x₀, y₀, z₀ - dz))
        )
    end

    # Initial Level 0 (k=1, 2)
    s_total = (w₀^3) * f(x₀, y₀, z₀)
    for i in 1:2
        ti = i * h
        vi, wi = ordinate(ti), weight(ti)
        for j in 1:2
            tj = j * h
            vj, wj = ordinate(tj), weight(tj)
            for k in 1:2
                tk = k * h
                vk, wk = ordinate(tk), weight(tk)
                s_total += add_octant(vi, vj, vk, wi, wj, wk)
            end
            s_total += add_planes(vi, vj, wi, wj)
        end
        s_total += add_axes(vi, wi)
    end

    old_res = Δx * Δy * Δz * h^3 * s_total

    # Cache for nodes and weights
    cache_x = T[]
    cache_w = T[]

    for level in 1:max_levels
        h /= 2
        s_new = zero(T)
        max_k = floor(Int, tm / h)

        # Populate cache
        resize!(cache_x, max_k)
        resize!(cache_w, max_k)
        for i in 1:max_k
            ti = i * h
            cache_x[i] = ordinate(ti)
            cache_w[i] = weight(ti)
        end

        for i in 1:max_k
            xi, wi = cache_x[i], cache_w[i]
            for j in 1:max_k
                xj, wj = cache_x[j], cache_w[j]
                for k in 1:max_k
                    (iseven(i) && iseven(j) && iseven(k)) && continue
                    s_new += add_octant(xi, xj, cache_x[k], wi, wj, cache_w[k])
                end
                if !(iseven(i) && iseven(j))
                    s_new += add_planes(xi, xj, wi, wj)
                end
            end
            if isodd(i)
                s_new += add_axes(xi, wi)
            end
        end

        s_total += s_new
        new_res = Δx * Δy * Δz * h^3 * s_total

        if abs(new_res - old_res) < tol * max(one(T), abs(new_res))
            return new_res
        end
        old_res = new_res
    end
    return old_res
end

"""
    adaptive_integrate_1D_cmpl(::Type{T}, f::Function, a, b; tol::Real=1e-12, max_levels::Int=10)

Adaptive 1D Tanh-Sinh integration for functions sensitive near endpoints.
`f` should accept three arguments: `f(x, 1-x, 1+x)` for the interval `[-1, 1]`.
For general `[a, b]`, the arguments are `f(x, (b-x)/(b-a), (x-a)/(b-a))`.
"""
function adaptive_integrate_1D_cmpl(::Type{T}, f::Function, a, b;
    tol::Real=1e-12, max_levels::Int=10) where {T<:Real}
    a_T, b_T = T(a), T(b)
    Δx = 0.5 * (b_T - a_T)
    x₀ = 0.5 * (b_T + a_T)

    tm = tmax(T)
    h = tm / 2
    w0 = T(π) / 2
    s_total = w0 * f(x₀, T(0.5), T(0.5))

    for k in 1:2
        tk = k * h
        wk, xk, ck = weight(tk), ordinate(tk), ordinate_complement(tk)
        # f(x, (b-x)/(b-a), (x-a)/(b-a))
        # At x = x₀ + Δx*xk:
        # (b - (x₀ + Δx*xk)) / (b - a) = (Δx - Δx*xk) / (2Δx) = 0.5 * (1 - xk) = 0.5 * ck
        # (x₀ + Δx*xk - a) / (b - a) = (Δx + Δx*xk) / (2Δx) = 0.5 * (1 + xk)
        s_total += wk * (f(x₀ + Δx * xk, 0.5 * ck, 0.5 * (one(T) + xk)) +
                         f(x₀ - Δx * xk, 0.5 * (one(T) + xk), 0.5 * ck))
    end

    old_res = Δx * h * s_total

    for level in 1:max_levels
        h /= 2
        s_new = zero(T)
        k = 1
        while true
            tk = k * h
            tk > tm && break
            wk, xk, ck = weight(tk), ordinate(tk), ordinate_complement(tk)
            s_new += wk * (f(x₀ + Δx * xk, 0.5 * ck, 0.5 * (one(T) + xk)) +
                           f(x₀ - Δx * xk, 0.5 * (one(T) + xk), 0.5 * ck))
            k += 2
        end

        s_total += s_new
        new_res = Δx * h * s_total

        if abs(new_res - old_res) < tol * max(one(T), abs(new_res))
            return new_res
        end
        old_res = new_res
    end
    return old_res
end
