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

# 2D Integration functions

@inline function _integrate2D_avx_unit(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        zero_T = zero(T)
        w₀ = _half_pi(T)
        w₀² = w₀ * w₀
        s = w₀² * f(zero_T, zero_T)

        @turbo for i in eachindex(x)
            wi = w[i]
            xi = x[i]
            s += wi * w₀ * (f(xi, zero_T) + f(-xi, zero_T) + f(zero_T, xi) + f(zero_T, -xi))
        end

        quadrants = zero(T)
        @turbo for i in eachindex(x), j in eachindex(x)
            wi = w[i]
            wj = w[j]
            xi = x[i]
            yj = x[j]
            quadrants += wi * wj * (f(-xi, -yj) + f(xi, -yj) + f(-xi, yj) + f(xi, yj))
        end

        return (h * h) * (s + quadrants)
    end
end

@inline function _integrate2D_avx_general(f::S, low::SVector{2,T}, up::SVector{2,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        Δx, x₀ = _midpoint_radius(low[1], up[1])
        Δy, y₀ = _midpoint_radius(low[2], up[2])

        w₀ = _half_pi(T)
        w₀² = w₀ * w₀
        s = w₀² * f(x₀, y₀)

        @turbo for i in eachindex(x)
            wi = w[i]
            xi = x[i]
            s += wi * w₀ * (f(Δx * xi + x₀, y₀) + f(x₀ - Δx * xi, y₀) + f(x₀, Δy * xi + y₀) + f(x₀, y₀ - Δy * xi))
        end

        quadrants = zero(T)
        @turbo for i in eachindex(x), j in eachindex(x)
            wi = w[i]
            wj = w[j]
            xi = x[i]
            yj = x[j]
            quadrants += wi * wj * (f(x₀ - Δx * xi, y₀ - Δy * yj) + f(Δx * xi + x₀, y₀ - Δy * yj) +
                                     f(x₀ - Δx * xi, Δy * yj + y₀) + f(Δx * xi + x₀, Δy * yj + y₀))
        end

        return Δx * Δy * (h * h) * (s + quadrants)
    end
end

"""
    integrate2D(f, x, w, h)

Calculate the 2D integral of `f` over `[-1, 1]^2` using pre-computed nodes/weights.
"""
function integrate2D(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    low = SVector{2,T}(-one(T), -one(T))
    up = SVector{2,T}(one(T), one(T))
    return integrate2D(f, low, up, x, w, h)
end

"""
    integrate2D(f, low, up, x, w, h)

Calculate the 2D integral of `f` over `[low, up]` (with `low, up::SVector{2}`) using pre-computed nodes/weights.
"""
function integrate2D(f::S, low::SVector{2,T}, up::SVector{2,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        Δx, x₀ = _midpoint_radius(low[1], up[1])
        Δy, y₀ = _midpoint_radius(low[2], up[2])

        w₀ = _half_pi(T)
        w₀² = w₀ * w₀
        s = w₀² * f(x₀, y₀)
        for k in eachindex(x)
            wk = w[k]
            xk_p = Δx * x[k] + x₀
            xk_m = x₀ - Δx * x[k]
            yk_p = Δy * x[k] + y₀
            yk_m = y₀ - Δy * x[k]
            s += wk * w₀ * (f(xk_p, y₀) + f(xk_m, y₀) + f(x₀, yk_p) + f(x₀, yk_m))
        end

        for i in eachindex(x)
            wi = w[i]
            xi_p = Δx * x[i] + x₀
            xi_m = x₀ - Δx * x[i]
            inner_s = zero(T)
            for j in eachindex(x)
                wj = w[j]
                yj_p = Δy * x[j] + y₀
                yj_m = y₀ - Δy * x[j]
                inner_s += wj * (f(xi_m, yj_m) + f(xi_p, yj_m) + f(xi_m, yj_p) + f(xi_p, yj_p))
            end
            s += wi * inner_s
        end
    end
    return Δx * Δy * (h * h) * s
end

"""
    integrate2D(f, low::SVector{2,<:Real}, up::SVector{2,<:Real}, x, w, h)

Calculate the 2D integral of `f` over `[low, up]` using pre-computed nodes/weights.
Bounds are converted to the node type `T` to support mixed real types such as `π`.
"""
function integrate2D(f::S, low::SVector{2,L}, up::SVector{2,U},
    x::X, w::W, h::T) where {T<:Real,L<:Real,U<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    return integrate2D(
        f,
        SVector{2,T}(T(low[1]), T(low[2])),
        SVector{2,T}(T(up[1]), T(up[2])),
        x, w, h
    )
end

"""
    integrate2D(f, low::AbstractVector{<:Real}, up::AbstractVector{<:Real}, x, w, h)

Convenience overload for bounds given as vectors of length 2.
"""
function integrate2D(f::S, low::AbstractVector{<:Real}, up::AbstractVector{<:Real},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    length(low) == 2 || throw(DimensionMismatch("low must have length 2"))
    length(up) == 2 || throw(DimensionMismatch("up must have length 2"))
    return integrate2D(
        f,
        SVector{2,T}(T(low[1]), T(low[2])),
        SVector{2,T}(T(up[1]), T(up[2])),
        x, w, h
    )
end

"""
    integrate2D_avx(f, low, up, x, w, h)

SIMD-accelerated 2D integration over `[low, up]` Using `LoopVectorization`.
"""
function integrate2D_avx(f::S, low::SVector{2,T}, up::SVector{2,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    if low[1] == -one(T) && low[2] == -one(T) && up[1] == one(T) && up[2] == one(T)
        return _integrate2D_avx_unit(f, x, w, h)
    end
    return _integrate2D_avx_general(f, low, up, x, w, h)
end

"""
    integrate2D_avx(f, low::SVector{2,<:Real}, up::SVector{2,<:Real}, x, w, h)

SIMD-accelerated 2D integration where bounds are converted to node type `T`.
"""
function integrate2D_avx(f::S, low::SVector{2,L}, up::SVector{2,U},
    x::X, w::W, h::T) where {T<:Real,L<:Real,U<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    return integrate2D_avx(
        f,
        SVector{2,T}(T(low[1]), T(low[2])),
        SVector{2,T}(T(up[1]), T(up[2])),
        x, w, h
    )
end

"""
    integrate2D_avx(f, low::AbstractVector{<:Real}, up::AbstractVector{<:Real}, x, w, h)

Convenience SIMD overload for bounds given as vectors of length 2.
"""
function integrate2D_avx(f::S, low::AbstractVector{<:Real}, up::AbstractVector{<:Real},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    length(low) == 2 || throw(DimensionMismatch("low must have length 2"))
    length(up) == 2 || throw(DimensionMismatch("up must have length 2"))
    return integrate2D_avx(
        f,
        SVector{2,T}(T(low[1]), T(low[2])),
        SVector{2,T}(T(up[1]), T(up[2])),
        x, w, h
    )
end
