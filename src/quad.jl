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
    quad(f, [low, up]; rtol, atol, max_levels, cache)

High-level interface for Tanh-Sinh quadrature. Automatically detects dimensions (1D, 2D, or 3D) 
and chooses the most efficient implementation.

- If `low` and `up` are absent, integrates over `[-1, 1]` (1D).
- Uses SIMD-accelerated (`_avx`) implementation for `Float32/Float64`.
- Uses higher-precision implementation for other types (e.g., `BigFloat`).
- Automatically converts `AbstractVector` to `SVector` for multi-dimensional integration.
- The adaptive error estimate is `abs(I_new - I_old)` for successive refinement levels.
- Refinement stops when `abs(I_new - I_old) <= max(atol, rtol * abs(I_new))`.
- If `rtol` is omitted and `atol == 0`, the default is `sqrt(eps(T))` for the promoted endpoint type `T`.
- `cache` optionally supplies precomputed adaptive nodes/weights; when omitted,
  the adaptive backend builds a cache automatically.
"""
function quad(f::F, low::T, up::T; rtol=nothing, atol::Real=0, max_levels::Int=16, cache=nothing) where {F,T<:Real}
    if !(T <: AbstractFloat)
        lowf, upf = _promote_real_bounds(low, up)
        return quad(f, lowf, upf; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
    end
    if low == up
        return zero(T)
    end
    if low > up
        return -quad(f, up, low; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
    end

    # 1D Case
    return adaptive_integrate_1D(T, f, low, up; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
end

function quad(f::F, low::Real, up::Real; rtol=nothing, atol::Real=0, max_levels::Int=16, cache=nothing) where {F}
    lowf, upf = _promote_real_bounds(low, up)
    return quad(f, lowf, upf; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
end

# Default 1D
"""
    quad(f; rtol, atol, max_levels, cache)

Adaptive 1D integration of `f` over the default interval `[-1, 1]`.
Uses the same stopping rule as `quad(f, low, up; rtol, atol, max_levels)`.
"""
quad(f::F; rtol=nothing, atol::Real=0, max_levels::Int=16, cache=nothing) where {F} =
    quad(f, -1.0, 1.0; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)

# Multi-dimensional cases
function quad(f::F, low::AbstractVector{<:Real}, up::AbstractVector{<:Real}; rtol=nothing, atol::Real=0, max_levels::Int=8, cache=nothing) where {F}
    n = length(low)
    n == length(up) || throw(DimensionMismatch("low and up must have the same length"))
    if n == 2
        t1 = typeof(float(low[1]))
        t2 = typeof(float(low[2]))
        t3 = typeof(float(up[1]))
        t4 = typeof(float(up[2]))
        T = promote_type(t1, t2, t3, t4)
        return quad(f, SVector{2,T}(T(low[1]), T(low[2])), SVector{2,T}(T(up[1]), T(up[2])); rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
    elseif n == 3
        t1 = typeof(float(low[1]))
        t2 = typeof(float(low[2]))
        t3 = typeof(float(low[3]))
        t4 = typeof(float(up[1]))
        t5 = typeof(float(up[2]))
        t6 = typeof(float(up[3]))
        T = promote_type(t1, t2, t3, t4, t5, t6)
        return quad(f, SVector{3,T}(T(low[1]), T(low[2]), T(low[3])), SVector{3,T}(T(up[1]), T(up[2]), T(up[3])); rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
    else
        error("Higher than 3D integration is not yet supported.")
    end
end

function quad(f::F, low::SVector{2,T}, up::SVector{2,T}; rtol=nothing, atol::Real=0, max_levels::Int=8, cache=nothing) where {F,T<:Real}
    if !(T <: AbstractFloat)
        return quad(f, SVector(float(low[1]), float(low[2])), SVector(float(up[1]), float(up[2])); rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
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
        return -quad(f, SVector(x1_low, x2_low), SVector(x1_up, x2_up); rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
    end

    return adaptive_integrate_2D(T, f, SVector(x1_low, x2_low), SVector(x1_up, x2_up); rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
end

function quad(f::F, low::SVector{3,T}, up::SVector{3,T}; rtol=nothing, atol::Real=0, max_levels::Int=5, cache=nothing) where {F,T<:Real}
    if !(T <: AbstractFloat)
        return quad(
            f,
            SVector(float(low[1]), float(low[2]), float(low[3])),
            SVector(float(up[1]), float(up[2]), float(up[3]));
            rtol=rtol, atol=atol, max_levels=max_levels, cache=cache
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
        return -quad(f, SVector(x1_low, x2_low, x3_low), SVector(x1_up, x2_up, x3_up); rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
    end

    return adaptive_integrate_3D(T, f, SVector(x1_low, x2_low, x3_low), SVector(x1_up, x2_up, x3_up); rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
end

"""
    quad_split(f, c, [low, up]; rtol, atol, max_levels, cache)

Split the integration domain at point `c` (singularity) and integrate sub-domains separately.
The requested `rtol` and `atol` are split evenly across the subdomains.
"""
function quad_split(f::F, c::T, low::T, up::T; rtol=nothing, atol::Real=0, max_levels::Int=16, cache=nothing) where {F,T<:Real}
    if !(T <: AbstractFloat)
        cf, lowf, upf = promote(float(c), float(low), float(up))
        return quad_split(f, cf, lowf, upf; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
    end
    # [low, c] + [c, up]
    sub_rtol = rtol === nothing ? nothing : rtol / 2
    sub_atol = atol / 2
    return quad(f, low, c; rtol=sub_rtol, atol=sub_atol, max_levels=max_levels, cache=cache) +
           quad(f, c, up; rtol=sub_rtol, atol=sub_atol, max_levels=max_levels, cache=cache)
end

function quad_split(f::F, c::Real, low::Real, up::Real; rtol=nothing, atol::Real=0, max_levels::Int=16, cache=nothing) where {F}
    cf, lowf, upf = promote(float(c), float(low), float(up))
    return quad_split(f, cf, lowf, upf; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
end

"""
    quad_split(f, c; rtol, atol, max_levels, cache)

Split the default interval `[-1, 1]` at singularity `c` and integrate both sides adaptively.
The requested `rtol` and `atol` are split evenly across the subdomains.
"""
function quad_split(f::F, c::Real; rtol=nothing, atol::Real=0, max_levels::Int=16, cache=nothing) where {F}
    cf = float(c)
    return quad_split(f, cf, -one(cf), one(cf); rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
end

function quad_split(f::F, c::SVector{2,T}, low::SVector{2,T}, up::SVector{2,T}; rtol=nothing, atol::Real=0, max_levels::Int=8, cache=nothing) where {F,T<:Real}
    if !(T <: AbstractFloat)
        return quad_split(
            f,
            SVector(float(c[1]), float(c[2])),
            SVector(float(low[1]), float(low[2])),
            SVector(float(up[1]), float(up[2]));
            rtol=rtol, atol=atol, max_levels=max_levels, cache=cache
        )
    end
    # Split into 4 quadrants
    cx, cy = c[1], c[2]
    x1, x2 = low[1], up[1]
    y1, y2 = low[2], up[2]

    sub_rtol = rtol === nothing ? nothing : rtol / 4
    sub_atol = atol / 4
    q1 = quad(f, SVector(x1, y1), SVector(cx, cy); rtol=sub_rtol, atol=sub_atol, max_levels=max_levels, cache=cache)
    q2 = quad(f, SVector(cx, y1), SVector(x2, cy); rtol=sub_rtol, atol=sub_atol, max_levels=max_levels, cache=cache)
    q3 = quad(f, SVector(x1, cy), SVector(cx, y2); rtol=sub_rtol, atol=sub_atol, max_levels=max_levels, cache=cache)
    q4 = quad(f, SVector(cx, cy), SVector(x2, y2); rtol=sub_rtol, atol=sub_atol, max_levels=max_levels, cache=cache)

    return q1 + q2 + q3 + q4
end

function quad_split(f::F, c::SVector{2,<:Real}, low::SVector{2,<:Real}, up::SVector{2,<:Real}; rtol=nothing, atol::Real=0, max_levels::Int=8, cache=nothing) where {F}
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
        rtol=rtol, atol=atol, max_levels=max_levels, cache=cache
    )
end

function quad_split(f::F, c::SVector{3,T}, low::SVector{3,T}, up::SVector{3,T}; rtol=nothing, atol::Real=0, max_levels::Int=5, cache=nothing) where {F,T<:Real}
    if !(T <: AbstractFloat)
        return quad_split(
            f,
            SVector(float(c[1]), float(c[2]), float(c[3])),
            SVector(float(low[1]), float(low[2]), float(low[3])),
            SVector(float(up[1]), float(up[2]), float(up[3]));
            rtol=rtol, atol=atol, max_levels=max_levels, cache=cache
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

    sub_rtol = rtol === nothing ? nothing : rtol / 8
    sub_atol = atol / 8
    total = zero(T)
    for i in 1:2, j in 1:2, k in 1:2
        total += quad(f, SVector(xs[i], ys[j], zs[k]), SVector(xs[i+1], ys[j+1], zs[k+1]); rtol=sub_rtol, atol=sub_atol, max_levels=max_levels, cache=cache)
    end
    return total
end

function quad_split(f::F, c::SVector{3,<:Real}, low::SVector{3,<:Real}, up::SVector{3,<:Real}; rtol=nothing, atol::Real=0, max_levels::Int=5, cache=nothing) where {F}
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
        rtol=rtol, atol=atol, max_levels=max_levels, cache=cache
    )
end

"""
    quad_cmpl(f, [low, up]; rtol, atol, max_levels, cache)

High-level 1D interface for endpoint-distance-aware integrands.
Call as `quad_cmpl(f, a, b)`. At each quadrature node `x` in `[a, b]`,
the callback is evaluated as `f(x, b_minus_x, x_minus_a)`, where
`b_minus_x = b - x` and `x_minus_a = x - a` (for `[-1, 1]`: `f(x, 1-x, 1+x)`).

This is useful when expressions like `log(b_minus_x)` or
`1/sqrt(b_minus_x*x_minus_a)` are sensitive to cancellation near endpoints.

The adaptive error estimate is based on successive refinement levels and uses
`abs(I_new - I_old) <= max(atol, rtol * abs(I_new))`, with the same default
`rtol = sqrt(eps(T))` behavior as `quad` when `atol == 0`.
"""
function quad_cmpl(f::F, low::T, up::T; rtol=nothing, atol::Real=0, max_levels::Int=16, cache=nothing) where {F,T<:Real}
    if !(T <: AbstractFloat)
        lowf, upf = _promote_real_bounds(low, up)
        return quad_cmpl(f, lowf, upf; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
    end
    if low == up
        return zero(T)
    end
    if low > up
        return -quad_cmpl(f, up, low; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
    end
    return adaptive_integrate_1D_cmpl(T, f, low, up; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
end

function quad_cmpl(f::F, low::Real, up::Real; rtol=nothing, atol::Real=0, max_levels::Int=16, cache=nothing) where {F}
    lowf, upf = _promote_real_bounds(low, up)
    return quad_cmpl(f, lowf, upf; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
end

quad_cmpl(f::F; rtol=nothing, atol::Real=0, max_levels::Int=16, cache=nothing) where {F} =
    quad_cmpl(f, -1.0, 1.0; rtol=rtol, atol=atol, max_levels=max_levels, cache=cache)
