# 3D Integration functions

"""
    integrate3D(f, x, w, h)

Calculate the 3D integral of `f` over `[-1, 1]^3` using pre-computed nodes/weights.
"""
function integrate3D(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        # Tanh-Sinh weight at t=0 is π/2
        w0 = T(π) / 2
        w0sq = w0 * w0
        w0cub = w0sq * w0

        # 1. The Origin (0,0,0)
        total_sum = w0cub * f(zero(T), zero(T), zero(T))

        # 2. The Axes (One dimension non-zero)
        axis_acc = zero(T)
        for i in eachindex(x)
            val = x[i]
            axis_points = (f(val, zero(T), zero(T)) + f(-val, zero(T), zero(T))) +
                          (f(zero(T), val, zero(T)) + f(zero(T), -val, zero(T))) +
                          (f(zero(T), zero(T), val) + f(zero(T), zero(T), -val))
            axis_acc += w[i] * axis_points
        end
        total_sum += w0sq * axis_acc

        # 3. The Planes and Octants
        for i in eachindex(x)
            wi = w[i]
            xi = x[i]
            sub_sum_i = zero(T)

            for j in eachindex(x)
                wj = w[j]
                wiwj = wi * wj
                xj = x[j]

                planes = (f(xi, xj, zero(T)) + f(-xi, xj, zero(T)) + f(xi, -xj, zero(T)) + f(-xi, -xj, zero(T))) +
                         (f(xi, zero(T), xj) + f(-xi, zero(T), xj) + f(xi, zero(T), -xj) + f(-xi, zero(T), -xj)) +
                         (f(zero(T), xi, xj) + f(zero(T), -xi, xj) + f(zero(T), xi, -xj) + f(zero(T), -xi, -xj))

                total_sum += wiwj * w0 * planes

                octant_acc = zero(T)
                for k in eachindex(x)
                    xk = x[k]
                    octant_acc += w[k] * (
                        (f(xi, xj, xk) + f(-xi, xj, xk) + f(xi, -xj, xk) + f(-xi, -xj, xk)) +
                        (f(xi, xj, -xk) + f(-xi, xj, -xk) + f(xi, -xj, -xk) + f(-xi, -xj, -xk))
                    )
                end
                total_sum += wiwj * octant_acc
            end
        end
    end
    return (h^3) * total_sum
end

"""
    integrate3D_avx(f, x, w, h)

SIMD-accelerated 3D integral over `[-1, 1]^3`.
"""
function integrate3D_avx(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        w0 = T(π) / 2
        w0sq = w0 * w0
        w0cub = w0sq * w0

        total_sum = w0cub * f(zero(T), zero(T), zero(T))

        axis_acc = zero(T)
        @turbo for i in eachindex(x)
            val = x[i]
            axis_points = (f(val, zero(T), zero(T)) + f(-val, zero(T), zero(T))) +
                          (f(zero(T), val, zero(T)) + f(zero(T), -val, zero(T))) +
                          (f(zero(T), zero(T), val) + f(zero(T), zero(T), -val))
            axis_acc += w[i] * axis_points
        end
        total_sum += w0sq * axis_acc

        for i in eachindex(x)
            wi = w[i]
            xi = x[i]
            sub_sum_i = zero(T)

            for j in eachindex(x)
                wj = w[j]
                wiwj = wi * wj
                xj = x[j]

                planes = (f(xi, xj, zero(T)) + f(-xi, xj, zero(T)) + f(xi, -xj, zero(T)) + f(-xi, -xj, zero(T))) +
                         (f(xi, zero(T), xj) + f(-xi, zero(T), xj) + f(xi, zero(T), -xj) + f(-xi, zero(T), -xj)) +
                         (f(zero(T), xi, xj) + f(zero(T), -xi, xj) + f(zero(T), xi, -xj) + f(zero(T), -xi, -xj))

                total_sum += wiwj * w0 * planes

                octant_acc = zero(T)
                @turbo for k in eachindex(x)
                    xk = x[k]
                    octant_acc += w[k] * (
                        (f(xi, xj, xk) + f(-xi, xj, xk) + f(xi, -xj, xk) + f(-xi, -xj, xk)) +
                        (f(xi, xj, -xk) + f(-xi, xj, -xk) + f(xi, -xj, -xk) + f(-xi, -xj, -xk))
                    )
                end
                total_sum += wiwj * octant_acc
            end
        end
    end
    return (h^3) * total_sum
end

"""
    integrate3D(f, low, up, x, w, h)

Calculate the 3D integral of `f` over `[low, up]` (with `low, up::SVector{3}`).
"""
function integrate3D(f::S, low::SVector{3,T}, up::SVector{3,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        Δx = 0.5 * (up[1] - low[1])
        Δy = 0.5 * (up[2] - low[2])
        Δz = 0.5 * (up[3] - low[3])
        x₀ = 0.5 * (up[1] + low[1])
        y₀ = 0.5 * (up[2] + low[2])
        z₀ = 0.5 * (up[3] + low[3])
        w₀ = T(π) / 2
        w₀² = w₀^2
        w₀³ = w₀² * w₀
        total_sum = w₀³ * f(x₀, y₀, z₀)

        for i in eachindex(x)
            wi = w[i]
            dx, dy, dz = Δx * x[i], Δy * x[i], Δz * x[i]

            axis_sum = (f(x₀ + dx, y₀, z₀) + f(x₀ - dx, y₀, z₀)) +
                       (f(x₀, y₀ + dy, z₀) + f(x₀, y₀ - dy, z₀)) +
                       (f(x₀, y₀, z₀ + dz) + f(x₀, y₀, z₀ - dz))
            total_sum += wi * w₀² * axis_sum
        end

        for i in eachindex(x)
            wi = w[i]
            dx = Δx * x[i]
            xp, xm = x₀ + dx, x₀ - dx

            for j in eachindex(x)
                wj = w[j]
                wiwj = wi * wj
                dy = Δy * x[j]
                yp, ym = y₀ + dy, y₀ - dy
                dz_j = Δz * x[j]
                dx_j = Δx * x[j]

                plane_sum = (f(xp, yp, z₀) + f(xm, yp, z₀) + f(xp, ym, z₀) + f(xm, ym, z₀)) +
                            (f(xp, y₀, z₀ + dz_j) + f(xm, y₀, z₀ + dz_j) + f(xp, y₀, z₀ - dz_j) + f(xm, y₀, z₀ - dz_j)) +
                            (f(x₀, yp, z₀ + dx_j) + f(x₀, ym, z₀ + dx_j) + f(x₀, yp, z₀ - dx_j) + f(x₀, ym, z₀ - dx_j))

                total_sum += wiwj * w₀ * plane_sum

                octant_sum = zero(T)
                for k in eachindex(x)
                    wk = w[k]
                    dz = Δz * x[k]
                    zp, zm = z₀ + dz, z₀ - dz

                    octant_sum += wk * (
                        (f(xp, yp, zp) + f(xm, yp, zp) + f(xp, ym, zp) + f(xm, ym, zp)) +
                        (f(xp, yp, zm) + f(xm, yp, zm) + f(xp, ym, zm) + f(xm, ym, zm))
                    )
                end
                total_sum += wiwj * octant_sum
            end
        end
    end

    return (Δx * Δy * Δz * h^3) * total_sum
end

"""
    integrate3D_avx(f, low, up, x, w, h)

SIMD-accelerated 3D integral over `[low, up]`.
"""
function integrate3D_avx(f::S, low::SVector{3,T}, up::SVector{3,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        Δx = 0.5 * (up[1] - low[1])
        Δy = 0.5 * (up[2] - low[2])
        Δz = 0.5 * (up[3] - low[3])
        x₀ = 0.5 * (up[1] + low[1])
        y₀ = 0.5 * (up[2] + low[2])
        z₀ = 0.5 * (up[3] + low[3])
        w₀ = T(π) / 2
        w₀² = w₀^2
        w₀³ = w₀² * w₀
        total_sum = w₀³ * f(x₀, y₀, z₀)

        @turbo for i in eachindex(x)
            wi = w[i]
            dx, dy, dz = Δx * x[i], Δy * x[i], Δz * x[i]

            axis_sum = (f(x₀ + dx, y₀, z₀) + f(x₀ - dx, y₀, z₀)) +
                       (f(x₀, y₀ + dy, z₀) + f(x₀, y₀ - dy, z₀)) +
                       (f(x₀, y₀, z₀ + dz) + f(x₀, y₀, z₀ - dz))
            total_sum += wi * w₀² * axis_sum
        end

        for i in eachindex(x)
            wi = w[i]
            dx = Δx * x[i]
            xp, xm = x₀ + dx, x₀ - dx

            for j in eachindex(x)
                wj = w[j]
                wiwj = wi * wj
                dy = Δy * x[j]
                yp, ym = y₀ + dy, y₀ - dy
                dz_j = Δz * x[j]
                dx_j = Δx * x[j]

                plane_sum = (f(xp, yp, z₀) + f(xm, yp, z₀) + f(xp, ym, z₀) + f(xm, ym, z₀)) +
                            (f(xp, y₀, z₀ + dz_j) + f(xm, y₀, z₀ + dz_j) + f(xp, y₀, z₀ - dz_j) + f(xm, y₀, z₀ - dz_j)) +
                            (f(x₀, yp, z₀ + dx_j) + f(x₀, ym, z₀ + dx_j) + f(x₀, yp, z₀ - dx_j) + f(x₀, ym, z₀ - dx_j))

                total_sum += wiwj * w₀ * plane_sum

                octant_sum = zero(T)
                @turbo for k in eachindex(x)
                    wk = w[k]
                    dz = Δz * x[k]
                    zp, zm = z₀ + dz, z₀ - dz

                    octant_sum += wk * (
                        (f(xp, yp, zp) + f(xm, yp, zp) + f(xp, ym, zp) + f(xm, ym, zp)) +
                        (f(xp, yp, zm) + f(xm, yp, zm) + f(xp, ym, zm) + f(xm, ym, zm))
                    )
                end
                total_sum += wiwj * octant_sum
            end
        end
    end

    return (Δx * Δy * Δz * h^3) * total_sum
end
