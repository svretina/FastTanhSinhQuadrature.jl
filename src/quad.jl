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

# High-level quad interface
"""
    quad(f::Function, [low, up]; tol=1e-12, max_levels=10)

High-level interface for Tanh-Sinh quadrature. Automatically detects dimensions (1D, 2D, or 3D) 
and chooses the most efficient implementation.

- If `low` and `up` are absent, integrates over `[-1, 1]` (1D).
- Uses SIMD-accelerated (`_avx`) implementation for `Float32/Float64`.
- Uses higher-precision implementation for other types (e.g., `BigFloat`).
- Automatically converts `AbstractVector` to `SVector` for multi-dimensional integration.
"""
function quad(f::Function, low::T, up::T; tol::Real=1e-12, max_levels::Int=10) where {T<:Real}
    if low == up
        return zero(T)
    end
    if low > up
        return -quad(f, up, low; tol=tol, max_levels=max_levels)
    end

    # 1D Case
    return adaptive_integrate_1D(T, f, low, up; tol=tol, max_levels=max_levels)
end

# Default 1D
"""
    quad(f::Function; tol=1e-12, max_levels=10)

Adaptive 1D integration of `f` over the default interval `[-1, 1]`.
"""
quad(f::Function; tol::Real=1e-12, max_levels::Int=10) = quad(f, -1.0, 1.0; tol=tol, max_levels=max_levels)

# Multi-dimensional cases
function quad(f::Function, low::AbstractVector{T}, up::AbstractVector{T}; tol::Real=1e-10, max_levels::Int=8) where {T<:Real}
    n = length(low)
    if n == 2
        return quad(f, SVector{2,T}(low[1], low[2]), SVector{2,T}(up[1], up[2]); tol=tol, max_levels=max_levels)
    elseif n == 3
        return quad(f, SVector{3,T}(low[1], low[2], low[3]), SVector{3,T}(up[1], up[2], up[3]); tol=tol, max_levels=max_levels)
    else
        error("Higher than 3D integration is not yet supported.")
    end
end

function quad(f::Function, low::SVector{2,T}, up::SVector{2,T}; tol::Real=1e-10, max_levels::Int=8) where {T<:Real}
    if any(low .== up)
        return zero(T)
    end
    # Check bounds and flip if necessary
    sign_flip = one(T)
    x1_low, x1_up = low[1], up[1]
    x2_low, x2_up = low[2], up[2]

    if x1_low > x1_up
        x1_low, x1_up = x1_up, x1_low
        sign_flip *= -one(T)
    end
    if x2_low > x2_up
        x2_low, x2_up = x2_up, x2_low
        sign_flip *= -one(T)
    end

    if sign_flip == -one(T)
        return -quad(f, SVector(x1_low, x2_low), SVector(x1_up, x2_up); tol=tol, max_levels=max_levels)
    end

    return adaptive_integrate_2D(T, f, SVector(x1_low, x2_low), SVector(x1_up, x2_up); tol=tol, max_levels=max_levels)
end

function quad(f::Function, low::SVector{3,T}, up::SVector{3,T}; tol::Real=1e-8, max_levels::Int=5) where {T<:Real}
    if any(low .== up)
        return zero(T)
    end

    sign_flip = one(T)
    x1_low, x1_up = low[1], up[1]
    x2_low, x2_up = low[2], up[2]
    x3_low, x3_up = low[3], up[3]

    if x1_low > x1_up
        x1_low, x1_up = x1_up, x1_low
        sign_flip *= -one(T)
    end
    if x2_low > x2_up
        x2_low, x2_up = x2_up, x2_low
        sign_flip *= -one(T)
    end
    if x3_low > x3_up
        x3_low, x3_up = x3_up, x3_low
        sign_flip *= -one(T)
    end

    if sign_flip == -one(T)
        return -quad(f, SVector(x1_low, x2_low, x3_low), SVector(x1_up, x2_up, x3_up); tol=tol, max_levels=max_levels)
    end

    return adaptive_integrate_3D(T, f, SVector(x1_low, x2_low, x3_low), SVector(x1_up, x2_up, x3_up); tol=tol, max_levels=max_levels)
end

"""
    quad_split(f, c, [low, up]; tol=1e-12, max_levels=10)

Split the integration domain at point `c` (singularity) and integrate sub-domains separately.
"""
function quad_split(f::Function, c::T, low::T, up::T; tol::Real=1e-12, max_levels::Int=10) where {T<:Real}
    # [low, c] + [c, up]
    return quad(f, low, c; tol=tol / 2, max_levels=max_levels) +
           quad(f, c, up; tol=tol / 2, max_levels=max_levels)
end

"""
    quad_split(f::Function, c; tol=1e-12, max_levels=10)

Split the default interval `[-1, 1]` at singularity `c` and integrate both sides adaptively.
"""
quad_split(f::Function, c::T; tol::Real=1e-12, max_levels::Int=10) where {T<:Real} = quad_split(f, c, -one(T), one(T); tol=tol, max_levels=max_levels)

function quad_split(f::Function, c::SVector{2,T}, low::SVector{2,T}, up::SVector{2,T}; tol::Real=1e-10, max_levels::Int=8) where {T<:Real}
    # Split into 4 quadrants
    cx, cy = c[1], c[2]
    x1, x2 = low[1], up[1]
    y1, y2 = low[2], up[2]

    q1 = quad(f, SVector(x1, y1), SVector(cx, cy); tol=tol / 4, max_levels=max_levels)
    q2 = quad(f, SVector(cx, y1), SVector(x2, cy); tol=tol / 4, max_levels=max_levels)
    q3 = quad(f, SVector(x1, cy), SVector(cx, y2); tol=tol / 4, max_levels=max_levels)
    q4 = quad(f, SVector(cx, cy), SVector(x2, y2); tol=tol / 4, max_levels=max_levels)

    return q1 + q2 + q3 + q4
end

function quad_split(f::Function, c::SVector{3,T}, low::SVector{3,T}, up::SVector{3,T}; tol::Real=1e-8, max_levels::Int=5) where {T<:Real}
    # Split into 8 octants
    cx, cy, cz = c[1], c[2], c[3]
    x1, x2 = low[1], up[1]
    y1, y2 = low[2], up[2]
    z1, z2 = low[3], up[3]

    xs = (x1, cx, x2)
    ys = (y1, cy, y2)
    zs = (z1, cz, z2)

    total = zero(T)
    for i in 1:2, j in 1:2, k in 1:2
        total += quad(f, SVector(xs[i], ys[j], zs[k]), SVector(xs[i+1], ys[j+1], zs[k+1]); tol=tol / 8, max_levels=max_levels)
    end
    return total
end

"""
    quad_cmpl(f, [low, up]; tol=1e-12, max_levels=10)

High-level interface for Tanh-Sinh quadrature with endpoint sensitivity.
`f` should accept three arguments: `f(x, 1-x, 1+x)`.
"""
function quad_cmpl(f::Function, low::T, up::T; tol::Real=1e-12, max_levels::Int=10) where {T<:Real}
    if low == up
        return zero(T)
    end
    if low > up
        return -quad_cmpl(f, up, low; tol=tol, max_levels=max_levels)
    end
    return adaptive_integrate_1D_cmpl(T, f, low, up; tol=tol, max_levels=max_levels)
end

quad_cmpl(f::Function; tol::Real=1e-12, max_levels::Int=10) = quad_cmpl(f, -1.0, 1.0; tol=tol, max_levels=max_levels)
