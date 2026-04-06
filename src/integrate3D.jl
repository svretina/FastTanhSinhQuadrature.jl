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

# 3D Integration functions

@inline function _integrate3D_avx_unit(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        zero_T = zero(T)
        w₀ = _half_pi(T)
        w₀² = w₀ * w₀
        w₀³ = w₀² * w₀
        total_sum = w₀³ * f(zero_T, zero_T, zero_T)

        @turbo for i in eachindex(x)
            wi = w[i]
            xi = x[i]
            total_sum += wi * w₀² * (
                f(xi, zero_T, zero_T) + f(-xi, zero_T, zero_T) +
                f(zero_T, xi, zero_T) + f(zero_T, -xi, zero_T) +
                f(zero_T, zero_T, xi) + f(zero_T, zero_T, -xi)
            )
        end

        plane_sum = zero(T)
        @turbo for i in eachindex(x), j in eachindex(x)
            wi = w[i]
            wj = w[j]
            xi = x[i]
            xj = x[j]
            plane_sum += wi * wj * w₀ * (
                f(xi, xj, zero_T) + f(-xi, xj, zero_T) + f(xi, -xj, zero_T) + f(-xi, -xj, zero_T) +
                f(xi, zero_T, xj) + f(-xi, zero_T, xj) + f(xi, zero_T, -xj) + f(-xi, zero_T, -xj) +
                f(zero_T, xi, xj) + f(zero_T, -xi, xj) + f(zero_T, xi, -xj) + f(zero_T, -xi, -xj)
            )
        end

        octant_sum = zero(T)
        @turbo for i in eachindex(x), j in eachindex(x), k in eachindex(x)
            wi = w[i]
            wj = w[j]
            wk = w[k]
            xi = x[i]
            xj = x[j]
            xk = x[k]
            octant_sum += wi * wj * wk * (
                f(xi, xj, xk) + f(-xi, xj, xk) + f(xi, -xj, xk) + f(-xi, -xj, xk) +
                f(xi, xj, -xk) + f(-xi, xj, -xk) + f(xi, -xj, -xk) + f(-xi, -xj, -xk)
            )
        end

        return (h * h * h) * (total_sum + plane_sum + octant_sum)
    end
end

@inline function _integrate3D_avx_general(f::S, low::SVector{3,T}, up::SVector{3,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        Δx, x₀ = _midpoint_radius(low[1], up[1])
        Δy, y₀ = _midpoint_radius(low[2], up[2])
        Δz, z₀ = _midpoint_radius(low[3], up[3])
        w₀ = _half_pi(T)
        w₀² = w₀ * w₀
        w₀³ = w₀² * w₀
        total_sum = w₀³ * f(x₀, y₀, z₀)

        @turbo for i in eachindex(x)
            wi = w[i]
            xi = x[i]
            total_sum += wi * w₀² * (
                f(Δx * xi + x₀, y₀, z₀) + f(x₀ - Δx * xi, y₀, z₀) +
                f(x₀, Δy * xi + y₀, z₀) + f(x₀, y₀ - Δy * xi, z₀) +
                f(x₀, y₀, Δz * xi + z₀) + f(x₀, y₀, z₀ - Δz * xi)
            )
        end

        plane_sum = zero(T)
        @turbo for i in eachindex(x), j in eachindex(x)
            wi = w[i]
            wj = w[j]
            xi = x[i]
            xj = x[j]
            plane_sum += wi * wj * w₀ * (
                f(Δx * xi + x₀, Δy * xj + y₀, z₀) + f(x₀ - Δx * xi, Δy * xj + y₀, z₀) +
                f(Δx * xi + x₀, y₀ - Δy * xj, z₀) + f(x₀ - Δx * xi, y₀ - Δy * xj, z₀) +
                f(Δx * xi + x₀, y₀, Δz * xj + z₀) + f(x₀ - Δx * xi, y₀, Δz * xj + z₀) +
                f(Δx * xi + x₀, y₀, z₀ - Δz * xj) + f(x₀ - Δx * xi, y₀, z₀ - Δz * xj) +
                f(x₀, Δy * xi + y₀, Δz * xj + z₀) + f(x₀, y₀ - Δy * xi, Δz * xj + z₀) +
                f(x₀, Δy * xi + y₀, z₀ - Δz * xj) + f(x₀, y₀ - Δy * xi, z₀ - Δz * xj)
            )
        end

        octant_sum = zero(T)
        @turbo for i in eachindex(x), j in eachindex(x), k in eachindex(x)
            wi = w[i]
            wj = w[j]
            wk = w[k]
            xi = x[i]
            xj = x[j]
            xk = x[k]
            octant_sum += wi * wj * wk * (
                f(Δx * xi + x₀, Δy * xj + y₀, Δz * xk + z₀) + f(x₀ - Δx * xi, Δy * xj + y₀, Δz * xk + z₀) +
                f(Δx * xi + x₀, y₀ - Δy * xj, Δz * xk + z₀) + f(x₀ - Δx * xi, y₀ - Δy * xj, Δz * xk + z₀) +
                f(Δx * xi + x₀, Δy * xj + y₀, z₀ - Δz * xk) + f(x₀ - Δx * xi, Δy * xj + y₀, z₀ - Δz * xk) +
                f(Δx * xi + x₀, y₀ - Δy * xj, z₀ - Δz * xk) + f(x₀ - Δx * xi, y₀ - Δy * xj, z₀ - Δz * xk)
            )
        end

        return (Δx * Δy * Δz * (h * h * h)) * (total_sum + plane_sum + octant_sum)
    end
end

"""
    integrate3D(f, x, w, h)

Calculate the 3D integral of `f` over `[-1, 1]^3` using pre-computed nodes/weights.
"""
function integrate3D(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    low = SVector{3,T}(-one(T), -one(T), -one(T))
    up = SVector{3,T}(one(T), one(T), one(T))
    return integrate3D(f, low, up, x, w, h)
end

"""
    integrate3D_avx(f, x, w, h)

SIMD-accelerated 3D integral over `[-1, 1]^3`.
"""
function integrate3D_avx(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    return _integrate3D_avx_unit(f, x, w, h)
end

"""
    integrate3D(f, low, up, x, w, h)

Calculate the 3D integral of `f` over `[low, up]` (with `low, up::SVector{3}`).
"""
function integrate3D(f::S, low::SVector{3,T}, up::SVector{3,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        Δx, x₀ = _midpoint_radius(low[1], up[1])
        Δy, y₀ = _midpoint_radius(low[2], up[2])
        Δz, z₀ = _midpoint_radius(low[3], up[3])
        w₀ = _half_pi(T)
        w₀² = w₀ * w₀
        w₀³ = w₀² * w₀
        total_sum = w₀³ * f(x₀, y₀, z₀)

        for i in eachindex(x)
            wi = w[i]
            axis_sum = (f(Δx * x[i] + x₀, y₀, z₀) + f(x₀ - Δx * x[i], y₀, z₀)) +
                       (f(x₀, Δy * x[i] + y₀, z₀) + f(x₀, y₀ - Δy * x[i], z₀)) +
                       (f(x₀, y₀, Δz * x[i] + z₀) + f(x₀, y₀, z₀ - Δz * x[i]))
            total_sum += wi * w₀² * axis_sum
        end

        for i in eachindex(x)
            wi = w[i]
            xp = Δx * x[i] + x₀
            xm = x₀ - Δx * x[i]
            yp_i = Δy * x[i] + y₀
            ym_i = y₀ - Δy * x[i]

            for j in eachindex(x)
                wj = w[j]
                wiwj = wi * wj
                yp = Δy * x[j] + y₀
                ym = y₀ - Δy * x[j]
                zp_j = Δz * x[j] + z₀
                zm_j = z₀ - Δz * x[j]

                plane_sum = (f(xp, yp, z₀) + f(xm, yp, z₀) + f(xp, ym, z₀) + f(xm, ym, z₀)) +
                            (f(xp, y₀, zp_j) + f(xm, y₀, zp_j) + f(xp, y₀, zm_j) + f(xm, y₀, zm_j)) +
                            (f(x₀, yp_i, zp_j) + f(x₀, ym_i, zp_j) + f(x₀, yp_i, zm_j) + f(x₀, ym_i, zm_j))

                total_sum += wiwj * w₀ * plane_sum

                octant_sum = zero(T)
                for k in eachindex(x)
                    wk = w[k]
                    zp = Δz * x[k] + z₀
                    zm = z₀ - Δz * x[k]

                    octant_sum += wk * (
                        (f(xp, yp, zp) + f(xm, yp, zp) + f(xp, ym, zp) + f(xm, ym, zp)) +
                        (f(xp, yp, zm) + f(xm, yp, zm) + f(xp, ym, zm) + f(xm, ym, zm))
                    )
                end
                total_sum += wiwj * octant_sum
            end
        end
    end

    return (Δx * Δy * Δz * (h * h * h)) * total_sum
end

"""
    integrate3D(f, low::SVector{3,<:Real}, up::SVector{3,<:Real}, x, w, h)

Calculate the 3D integral of `f` over `[low, up]` using pre-computed nodes/weights.
Bounds are converted to the node type `T` to support mixed real types such as `π`.
"""
function integrate3D(f::S, low::SVector{3,L}, up::SVector{3,U},
    x::X, w::W, h::T) where {T<:Real,L<:Real,U<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    return integrate3D(
        f,
        SVector{3,T}(T(low[1]), T(low[2]), T(low[3])),
        SVector{3,T}(T(up[1]), T(up[2]), T(up[3])),
        x, w, h
    )
end

"""
    integrate3D(f, low::AbstractVector{<:Real}, up::AbstractVector{<:Real}, x, w, h)

Convenience overload for bounds given as vectors of length 3.
"""
function integrate3D(f::S, low::AbstractVector{<:Real}, up::AbstractVector{<:Real},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    length(low) == 3 || throw(DimensionMismatch("low must have length 3"))
    length(up) == 3 || throw(DimensionMismatch("up must have length 3"))
    return integrate3D(
        f,
        SVector{3,T}(T(low[1]), T(low[2]), T(low[3])),
        SVector{3,T}(T(up[1]), T(up[2]), T(up[3])),
        x, w, h
    )
end

"""
    integrate3D_avx(f, low, up, x, w, h)

SIMD-accelerated 3D integral over `[low, up]`.
"""
function integrate3D_avx(f::S, low::SVector{3,T}, up::SVector{3,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    if low[1] == -one(T) && low[2] == -one(T) && low[3] == -one(T) &&
       up[1] == one(T) && up[2] == one(T) && up[3] == one(T)
        return _integrate3D_avx_unit(f, x, w, h)
    end
    return _integrate3D_avx_general(f, low, up, x, w, h)
end

"""
    integrate3D_avx(f, low::SVector{3,<:Real}, up::SVector{3,<:Real}, x, w, h)

SIMD-accelerated 3D integration where bounds are converted to node type `T`.
"""
function integrate3D_avx(f::S, low::SVector{3,L}, up::SVector{3,U},
    x::X, w::W, h::T) where {T<:Real,L<:Real,U<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    return integrate3D_avx(
        f,
        SVector{3,T}(T(low[1]), T(low[2]), T(low[3])),
        SVector{3,T}(T(up[1]), T(up[2]), T(up[3])),
        x, w, h
    )
end

"""
    integrate3D_avx(f, low::AbstractVector{<:Real}, up::AbstractVector{<:Real}, x, w, h)

Convenience SIMD overload for bounds given as vectors of length 3.
"""
function integrate3D_avx(f::S, low::AbstractVector{<:Real}, up::AbstractVector{<:Real},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    length(low) == 3 || throw(DimensionMismatch("low must have length 3"))
    length(up) == 3 || throw(DimensionMismatch("up must have length 3"))
    return integrate3D_avx(
        f,
        SVector{3,T}(T(low[1]), T(low[2]), T(low[3])),
        SVector{3,T}(T(up[1]), T(up[2]), T(up[3])),
        x, w, h
    )
end
