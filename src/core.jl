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

# Core Tanh-Sinh quadrature functions: transformation, weights, and node generation

@inline _half(::Type{T}) where {T<:Real} = inv(one(T) + one(T))

@inline function _midpoint_radius(low::T, up::T) where {T<:Real}
    half = _half(T)
    return half * (up - low), half * (up + low)
end

@inline function _needs_generic_hyperbolics(::Type{T}) where {T<:Real}
    return nameof(T) === :MultiFloat && nameof(parentmodule(T)) === :MultiFloats
end

@inline function _asinh_generic(x::T) where {T<:Real}
    y = sqrt(x * x + one(T))
    return x >= zero(T) ? log(x + y) : -log(y - x)
end

@inline function _sinhcosh_generic(x::T) where {T<:Real}
    ex = exp(x)
    inv_ex = inv(ex)
    half = _half(T)
    return half * (ex - inv_ex), half * (ex + inv_ex)
end

@inline function _tanh_generic(x::T) where {T<:Real}
    two = one(T) + one(T)
    if x >= zero(T)
        z = exp(-two * x)
        return (one(T) - z) / (one(T) + z)
    else
        z = exp(two * x)
        return (z - one(T)) / (z + one(T))
    end
end

@inline function _asinh_compat(x::T) where {T<:Real}
    if _needs_generic_hyperbolics(T)
        return _asinh_generic(x)
    end
    return asinh(x)
end

@inline function _sinhcosh_compat(x::T) where {T<:Real}
    if _needs_generic_hyperbolics(T)
        return _sinhcosh_generic(x)
    end
    return sinh(x), cosh(x)
end

@inline function _tanh_compat(x::T) where {T<:Real}
    if _needs_generic_hyperbolics(T)
        return _tanh_generic(x)
    end
    return tanh(x)
end

@inline function ordinate(t::T) where {T<:Real}
    # x = tanh(π/2 * sinh(t))
    # Stability: For large t, tanh(u) -> 1 - 2exp(-2u). 
    # Directly using tanh is fine for most cases, but we guard against rounding to 1.0.
    sinh_t, = _sinhcosh_compat(t)
    val = _tanh_compat(T(π) / 2 * sinh_t)

    if val >= one(T)
        return prevfloat(one(T))
    elseif val <= -one(T)
        return -prevfloat(one(T))
    end
    return val
end

"""
    ordinate_complement(t::T) where {T<:Real}

Return 1 - |ordinate(t)| accurately. Useful for f(1-x).
"""
@inline function ordinate_complement(t::T) where {T<:Real}
    sinh_abs_t, = _sinhcosh_compat(abs(t))
    u = (T(π) / 2) * sinh_abs_t
    # 1 - tanh(u) = 2*exp(-2u) / (1 + exp(-2u))
    # This avoids overflow in exp(2u) while preserving endpoint accuracy.
    z = exp(-2u)
    return 2 * z / (one(T) + z)
end

@inline function weight(t::T) where {T<:Real}
    sinh_t, cosh_t = _sinhcosh_compat(t)
    arg = T(π) / 2 * sinh_t
    # Stability: cosh can overflow for large t.
    # If the denominator cosh^2(...) would overflow, the weight is effectively 0.
    # For Float64, cosh(710) overflows.
    if abs(arg) > T(700.0)
        return zero(T)
    end
    _, tmp = _sinhcosh_compat(arg)
    # weight = (π/2 * cosh(t)) / cosh^2(π/2 * sinh(t))
    return ((T(π) / 2) * cosh_t) / (tmp * tmp)
end

@inline function inv_ordinate(t::T) where {T<:Real}
    return _asinh_compat(log((one(T) + t) / (one(T) - t)) / T(π))
end

"""
    t_x_max(::Type{T}) where {T<:Real}

Calculate the maximum `t` to avoid abscissae underflow (Eq. 13 in arXiv:2007.15057).
Ensures `1 - |t| >= floatmin(T)`.
"""
@inline function t_x_max(::Type{T}) where {T<:Real}
    return inv_ordinate(prevfloat(one(T)))
end

"""
    t_w_max(::Type{T}, D::Int=1) where {T<:Real}

Calculate the maximum `t` to avoid weight underflow (Eq. 15 in arXiv:2007.15057).
Ensures `(ψ'(t))^calD >= floatmin(T)` where `calD = max(1, D-1)`.
"""
@inline function t_w_max(::Type{T}, D::Int=1) where {T<:Real}
    calD = T(D)
    fmin = floatmin(T)
    # Target: π * sinh(t) - t <= ln(π) - (1/calD) * ln(fmin)
    rhs = log(T(π)) - log(fmin) / calD
    # Solve π * sinh(t) - t = rhs using Newton iteration
    t = _asinh_compat(rhs / T(π))
    for _ in 1:10
        sinh_t, cosh_t = _sinhcosh_compat(t)
        f = T(π) * sinh_t - t - rhs
        df = T(π) * cosh_t - 1
        t = t - f / df
    end
    return t
end

"""
    tmax(::Type{T}, D::Int=1) where {T<:Real}

Find the optimal window limit `tmax` considering abscissae and weight underflow (Eq. 14 in arXiv:2007.15057).
"""
@inline function tmax(::Type{T}, D::Int=1) where {T<:Real}
    return min(t_x_max(T), t_w_max(T, D))
end

"""
    tanhsinh(::Type{T}, N::Val{N}) where {T<:AbstractFloat,N}

Generate Tanh-Sinh quadrature nodes `x`, weights `w`, and step size `h` for a given floating point type `T` and number of points `Val{N}`.
Use this for approximately N<128
"""
function tanhsinh(::Type{T}, ::Val{N}; D::Int=1) where {T<:AbstractFloat,N}
    if iseven(N)
        n = N ÷ 2
    else
        n = (N - 1) ÷ 2
    end
    tm = tmax(T, D)
    h = tm / n
    t = range(h, tm, length=n)
    x = ordinate.(t)
    w = weight.(t)
    return SVector{n,T}(x), SVector{n,T}(w), h
end

"""
    tanhsinh(::Type{T}, N::Int; D::Int=1) where {T<:AbstractFloat}

Generate Tanh-Sinh quadrature nodes `x`, weights `w`, and step size `h` for a given floating point type `T` and number of points `N`.
"""
function tanhsinh(::Type{T}, N::Int; D::Int=1) where {T<:AbstractFloat}
    if iseven(N)
        n = N ÷ 2
    else
        n = (N - 1) ÷ 2
    end
    tm = tmax(T, D)
    h = tm / n
    t = range(h, tm, length=n)
    x = ordinate.(t)
    w = weight.(t)
    return x, w, h
end

"""
    tanhsinh(N::Int)

Generate Tanh-Sinh quadrature nodes `x`, weights `w`, and step size `h` for `Float64` precision and `N` points.
"""
tanhsinh(N::Int) = tanhsinh(Float64, N)

function tanhsinh_opt(::Type{T}, N::Int, d::Real=π / 2; D::Int=1) where {T<:Real}
    if iseven(N)
        n = N ÷ 2
    else
        n = (N - 1) ÷ 2
    end
    tm = tmax(T, D)
    h = hopt(T, N, d)
    t = h:h:tm
    x = ordinate.(t)
    w = weight.(t)
    return x, w, h
end

@inline function hopt(::Type{T}, N::Int, d::Real=π / 2) where {T}
    (T(2) / T(N)) * lambertw(T(2d) * T(N))
end

@inline function hmax(::Type{T}, N::Int, D::Int=1) where {T}
    if iseven(N)
        n = N ÷ 2
    else
        n = (N - 1) ÷ 2
    end
    tm = tmax(T, D)
    return tm / n
end
