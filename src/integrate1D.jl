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

# 1D Integration functions

"""
    integrate1D(::Type{T}, f::Function, N::Int) where {T<:Real}

Calculate the integral of `f` over `[-1, 1]` using `N` Tanh-Sinh quadrature points in precision `T`.
"""
function integrate1D(::Type{T}, f::Function, N::Int) where {T<:Real}
    x, w, h = tanhsinh(T, N)
    s = T(π) / 2 * f(zero(T))
    @inbounds for i in eachindex(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
    end
    return h * s
end

"""
    integrate1D(f::Function, N::Int)

Calculate the integral of `f` over `[-1, 1]` using `N` Tanh-Sinh quadrature points in `Float64` precision.
"""
# 24-25: integrate1D
function integrate1D(f::Function, N::Int)
    return integrate1D(Float64, f, N)
end

"""
    integrate1D_cmpl(::Type{T}, f::Function, N::Int) where {T<:Real}

Calculate the integral of `f(x, 1-x, 1+x)` over `[-1, 1]` using `N` points.
`f` should accept three arguments: `f(x, 1-x, 1+x)`.
"""
function integrate1D_cmpl(::Type{T}, f::Function, N::Int) where {T<:Real}
    x, w, h = tanhsinh(T, N)
    # Origin
    s = T(π) / 2 * f(zero(T), one(T), one(T))
    for i in eachindex(x)
        xi = x[i]
        # 1 - xi = ordinate_complement(ti)
        # However, for N-point fixed grid, we might just use the ordinate_complement logic
        # But wait, tanhsinh(T, N) only returns x. 
        # For simplicity in this fixed-N version, we use the standard complement.
        # Fixed nodes don't easily give us 't' unless we recompute.
        # Let's assume user wants accuracy.
        one_minus_x = one(T) - xi
        one_plus_x = one(T) + xi
        s += w[i] * (f(xi, one_minus_x, one_plus_x) + f(-xi, one_plus_x, one_minus_x))
    end
    return h * s
end

# [-1,1] by default 1D
"""
    integrate1D(f, x, w, h)

Calculate the integral of `f` over `[-1, 1]` using pre-computed nodes `x`, weights `w`, and step size `h`.
"""
function integrate1D(f::X, x::AbstractVector{T}, w::AbstractVector{T},
    h::T) where {T<:Real,X}
    s = T(π) / 2 * f(zero(T))
    @inbounds for i in 1:length(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
    end
    return h * s
end

# [a,b] 1D
"""
    integrate1D(f, low, up, x, w, h)

Calculate the integral of `f` over `[low, up]` using pre-computed nodes `x`, weights `w`, and step size `h`.
"""
function integrate1D(f::S, low::T, up::T, x::X,
    w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        Δx = 0.5(up - low)
        x₀ = 0.5(up + low)
        s = T(π) / 2 * f(x₀)
        for i in eachindex(x)
            xp = x₀ + Δx * x[i]
            xm = x₀ - Δx * x[i]
            s += w[i] * (f(xm) + f(xp))
        end
    end
    return Δx * h * s
end

"""
    integrate1D_avx(f, x, w, h)

SIMD-accelerated 1D integration over `[-1, 1]` Using `LoopVectorization`. 
Requires `f` to be compatible with `@turbo`. Only beneficial for `Float32/Float64`.
"""
function integrate1D_avx(f::S, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}
    s = weight(zero(T)) * f(zero(T))
    @turbo for i in eachindex(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
    end
    return h * s
end

"""
    integrate1D_avx(f, low, up, x, w, h)

SIMD-accelerated 1D integration over `[low, up]` Using `LoopVectorization`.
"""
function integrate1D_avx(f::S, low::T, up::T, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}
    Δx = 0.5(up - low)
    x₀ = 0.5(up + low)
    s = T(π) / 2 * f(x₀)
    @turbo for i in 1:length(x)
        Δxxi = Δx * x[i]
        xp = x₀ + Δxxi
        xm = x₀ - Δxxi
        s += w[i] * (f(xm) + f(xp))
    end
    return @fastmath Δx * h * s
end
