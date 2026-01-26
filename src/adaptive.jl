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
    s_origin = w0 * f(x₀)

    # Level 0 includes points at t=1h and t=2h
    s_weighted = zero(T)
    for k in 1:2
        tk = k * h
        s_weighted += weight(tk) * (f(x₀ + Δx * ordinate(tk)) + f(x₀ - Δx * ordinate(tk)))
    end

    total_weighted_sum = s_origin + s_weighted
    old_res = Δx * h * total_weighted_sum

    for level in 1:max_levels
        h /= 2
        s_new = zero(T)

        # New points are at ODD multiples of the new h: 1h, 3h, 5h...
        k = 1
        while true
            tk = k * h
            tk > tm && break
            s_new += weight(tk) * (f(x₀ + Δx * ordinate(tk)) + f(x₀ - Δx * ordinate(tk)))
            k += 2
        end

        total_weighted_sum += s_new
        new_res = Δx * h * total_weighted_sum

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
    function eval_quadrants(ti, tj, wi, wj)
        xi, yi = ordinate(ti), ordinate(tj)
        dx, dy = Δx * xi, Δy * yi
        return wi * wj * (f(x₀ + dx, y₀ + dy) + f(x₀ - dx, y₀ + dy) +
                          f(x₀ + dx, y₀ - dy) + f(x₀ - dx, y₀ - dy))
    end

    # Helper to evaluate symmetric axis points
    function eval_axes(tk, wk)
        val = ordinate(tk)
        dx, dy = Δx * val, Δy * val
        return wk * w0 * (f(x₀ + dx, y₀) + f(x₀ - dx, y₀) +
                          f(x₀, y₀ + dy) + f(x₀, y₀ - dy))
    end

    # Initial Level 0 (h, 2h)
    s_total = (w0^2) * f(x₀, y₀)
    for i in 1:2, j in 1:2
        s_total += eval_quadrants(i * h, j * h, weight(i * h), weight(j * h))
    end
    for k in 1:2
        s_total += eval_axes(k * h, weight(k * h))
    end

    old_res = Δx * Δy * h^2 * s_total

    for level in 1:max_levels
        h /= 2
        s_new = zero(T)
        max_k = floor(Int, tm / h)

        for i in 1:max_k
            wi, ti = weight(i * h), i * h
            for j in 1:max_k
                # Point is new if at least one index is odd
                (iseven(i) && iseven(j)) && continue
                s_new += eval_quadrants(ti, j * h, wi, weight(j * h))
            end
            if isodd(i)
                s_new += eval_axes(ti, wi)
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
    function add_octant(ti, tj, tk, wi, wj, wk)
        vi, vj, vk = ordinate(ti), ordinate(tj), ordinate(tk)
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
    function add_planes(ti, tj, wi, wj)
        vi, vj = ordinate(ti), ordinate(tj)
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
    function add_axes(ti, wi)
        vi = ordinate(ti)
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
    for i in 1:2, j in 1:2, k in 1:2
        s_total += add_octant(i * h, j * h, k * h, weight(i * h), weight(j * h), weight(k * h))
    end
    for i in 1:2, j in 1:2
        s_total += add_planes(i * h, j * h, weight(i * h), weight(j * h))
    end
    for i in 1:2
        s_total += add_axes(i * h, weight(i * h))
    end

    old_res = Δx * Δy * Δz * h^3 * s_total

    for level in 1:max_levels
        h /= 2
        s_new = zero(T)
        max_k = floor(Int, tm / h)

        for i in 1:max_k
            wi, ti = weight(i * h), i * h
            for j in 1:max_k
                wj, tj = weight(j * h), j * h
                for k in 1:max_k
                    (iseven(i) && iseven(j) && iseven(k)) && continue
                    s_new += add_octant(ti, tj, k * h, wi, wj, weight(k * h))
                end
                if !(iseven(i) && iseven(j))
                    s_new += add_planes(ti, tj, wi, wj)
                end
            end
            if isodd(i)
                s_new += add_axes(ti, wi)
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
