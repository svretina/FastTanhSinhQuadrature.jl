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
@inline _promote_real_bounds(low::Real, up::Real) = promote(float(low), float(up))

"""
    quad(f, [low, up]; tol=1e-12, max_levels=10)

High-level interface for Tanh-Sinh quadrature. Automatically detects dimensions (1D, 2D, or 3D) 
and chooses the most efficient implementation.

- If `low` and `up` are absent, integrates over `[-1, 1]` (1D).
- Uses SIMD-accelerated (`_avx`) implementation for `Float32/Float64`.
- Uses higher-precision implementation for other types (e.g., `BigFloat`).
- Automatically converts `AbstractVector` to `SVector` for multi-dimensional integration.
"""
function quad(f::F, low::T, up::T; tol::Real=1e-12, max_levels::Int=10) where {F,T<:Real}
    if !(T <: AbstractFloat)
        lowf, upf = _promote_real_bounds(low, up)
        return quad(f, lowf, upf; tol=tol, max_levels=max_levels)
    end
    if low == up
        return zero(T)
    end
    if low > up
        return -quad(f, up, low; tol=tol, max_levels=max_levels)
    end

    # 1D Case
    return adaptive_integrate_1D(T, f, low, up; tol=tol, max_levels=max_levels)
end

function quad(f::F, low::Real, up::Real; tol::Real=1e-12, max_levels::Int=10) where {F}
    lowf, upf = _promote_real_bounds(low, up)
    return quad(f, lowf, upf; tol=tol, max_levels=max_levels)
end

# Default 1D
"""
    quad(f; tol=1e-12, max_levels=10)

Adaptive 1D integration of `f` over the default interval `[-1, 1]`.
"""
quad(f::F; tol::Real=1e-12, max_levels::Int=10) where {F} = quad(f, -1.0, 1.0; tol=tol, max_levels=max_levels)

# Multi-dimensional cases
function quad(f::F, low::AbstractVector{<:Real}, up::AbstractVector{<:Real}; tol::Real=1e-10, max_levels::Int=8) where {F}
    n = length(low)
    n == length(up) || throw(DimensionMismatch("low and up must have the same length"))
    if n == 2
        t1 = typeof(float(low[1]))
        t2 = typeof(float(low[2]))
        t3 = typeof(float(up[1]))
        t4 = typeof(float(up[2]))
        T = promote_type(t1, t2, t3, t4)
        return quad(f, SVector{2,T}(T(low[1]), T(low[2])), SVector{2,T}(T(up[1]), T(up[2])); tol=tol, max_levels=max_levels)
    elseif n == 3
        t1 = typeof(float(low[1]))
        t2 = typeof(float(low[2]))
        t3 = typeof(float(low[3]))
        t4 = typeof(float(up[1]))
        t5 = typeof(float(up[2]))
        t6 = typeof(float(up[3]))
        T = promote_type(t1, t2, t3, t4, t5, t6)
        return quad(f, SVector{3,T}(T(low[1]), T(low[2]), T(low[3])), SVector{3,T}(T(up[1]), T(up[2]), T(up[3])); tol=tol, max_levels=max_levels)
    else
        error("Higher than 3D integration is not yet supported.")
    end
end

function quad(f::F, low::SVector{2,T}, up::SVector{2,T}; tol::Real=1e-10, max_levels::Int=8) where {F,T<:Real}
    if !(T <: AbstractFloat)
        return quad(f, SVector(float(low[1]), float(low[2])), SVector(float(up[1]), float(up[2])); tol=tol, max_levels=max_levels)
    end
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

function quad(f::F, low::SVector{3,T}, up::SVector{3,T}; tol::Real=1e-8, max_levels::Int=5) where {F,T<:Real}
    if !(T <: AbstractFloat)
        return quad(
            f,
            SVector(float(low[1]), float(low[2]), float(low[3])),
            SVector(float(up[1]), float(up[2]), float(up[3]));
            tol=tol, max_levels=max_levels
        )
    end
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
function quad_split(f::F, c::T, low::T, up::T; tol::Real=1e-12, max_levels::Int=10) where {F,T<:Real}
    if !(T <: AbstractFloat)
        cf, lowf, upf = promote(float(c), float(low), float(up))
        return quad_split(f, cf, lowf, upf; tol=tol, max_levels=max_levels)
    end
    # [low, c] + [c, up]
    return quad(f, low, c; tol=tol / 2, max_levels=max_levels) +
           quad(f, c, up; tol=tol / 2, max_levels=max_levels)
end

function quad_split(f::F, c::Real, low::Real, up::Real; tol::Real=1e-12, max_levels::Int=10) where {F}
    cf, lowf, upf = promote(float(c), float(low), float(up))
    return quad_split(f, cf, lowf, upf; tol=tol, max_levels=max_levels)
end

"""
    quad_split(f, c; tol=1e-12, max_levels=10)

Split the default interval `[-1, 1]` at singularity `c` and integrate both sides adaptively.
"""
function quad_split(f::F, c::Real; tol::Real=1e-12, max_levels::Int=10) where {F}
    cf = float(c)
    return quad_split(f, cf, -one(cf), one(cf); tol=tol, max_levels=max_levels)
end

function quad_split(f::F, c::SVector{2,T}, low::SVector{2,T}, up::SVector{2,T}; tol::Real=1e-10, max_levels::Int=8) where {F,T<:Real}
    if !(T <: AbstractFloat)
        return quad_split(
            f,
            SVector(float(c[1]), float(c[2])),
            SVector(float(low[1]), float(low[2])),
            SVector(float(up[1]), float(up[2]));
            tol=tol, max_levels=max_levels
        )
    end
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

function quad_split(f::F, c::SVector{2,<:Real}, low::SVector{2,<:Real}, up::SVector{2,<:Real}; tol::Real=1e-10, max_levels::Int=8) where {F}
    t1 = typeof(float(c[1]))
    t2 = typeof(float(c[2]))
    t3 = typeof(float(low[1]))
    t4 = typeof(float(low[2]))
    t5 = typeof(float(up[1]))
    t6 = typeof(float(up[2]))
    T = promote_type(t1, t2, t3, t4, t5, t6)
    return quad_split(
        f,
        SVector{2,T}(T(c[1]), T(c[2])),
        SVector{2,T}(T(low[1]), T(low[2])),
        SVector{2,T}(T(up[1]), T(up[2]));
        tol=tol, max_levels=max_levels
    )
end

function quad_split(f::F, c::SVector{3,T}, low::SVector{3,T}, up::SVector{3,T}; tol::Real=1e-8, max_levels::Int=5) where {F,T<:Real}
    if !(T <: AbstractFloat)
        return quad_split(
            f,
            SVector(float(c[1]), float(c[2]), float(c[3])),
            SVector(float(low[1]), float(low[2]), float(low[3])),
            SVector(float(up[1]), float(up[2]), float(up[3]));
            tol=tol, max_levels=max_levels
        )
    end
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

function quad_split(f::F, c::SVector{3,<:Real}, low::SVector{3,<:Real}, up::SVector{3,<:Real}; tol::Real=1e-8, max_levels::Int=5) where {F}
    t1 = typeof(float(c[1]))
    t2 = typeof(float(c[2]))
    t3 = typeof(float(c[3]))
    t4 = typeof(float(low[1]))
    t5 = typeof(float(low[2]))
    t6 = typeof(float(low[3]))
    t7 = typeof(float(up[1]))
    t8 = typeof(float(up[2]))
    t9 = typeof(float(up[3]))
    T = promote_type(t1, t2, t3, t4, t5, t6, t7, t8, t9)
    return quad_split(
        f,
        SVector{3,T}(T(c[1]), T(c[2]), T(c[3])),
        SVector{3,T}(T(low[1]), T(low[2]), T(low[3])),
        SVector{3,T}(T(up[1]), T(up[2]), T(up[3]));
        tol=tol, max_levels=max_levels
    )
end

"""
    quad_cmpl(f, [low, up]; tol=1e-12, max_levels=10)

High-level 1D interface for endpoint-distance-aware integrands.
Call as `quad_cmpl(f, a, b)`. At each quadrature node `x` in `[a, b]`,
the callback is evaluated as `f(x, b_minus_x, x_minus_a)`, where
`b_minus_x = b - x` and `x_minus_a = x - a` (for `[-1, 1]`: `f(x, 1-x, 1+x)`).

This is useful when expressions like `log(b_minus_x)` or
`1/sqrt(b_minus_x*x_minus_a)` are sensitive to cancellation near endpoints.
"""
function quad_cmpl(f::F, low::T, up::T; tol::Real=1e-12, max_levels::Int=10) where {F,T<:Real}
    if !(T <: AbstractFloat)
        lowf, upf = _promote_real_bounds(low, up)
        return quad_cmpl(f, lowf, upf; tol=tol, max_levels=max_levels)
    end
    if low == up
        return zero(T)
    end
    if low > up
        return -quad_cmpl(f, up, low; tol=tol, max_levels=max_levels)
    end
    return adaptive_integrate_1D_cmpl(T, f, low, up; tol=tol, max_levels=max_levels)
end

function quad_cmpl(f::F, low::Real, up::Real; tol::Real=1e-12, max_levels::Int=10) where {F}
    lowf, upf = _promote_real_bounds(low, up)
    return quad_cmpl(f, lowf, upf; tol=tol, max_levels=max_levels)
end

quad_cmpl(f::F; tol::Real=1e-12, max_levels::Int=10) where {F} = quad_cmpl(f, -1.0, 1.0; tol=tol, max_levels=max_levels)
