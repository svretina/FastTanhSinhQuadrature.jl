# Core Tanh-Sinh quadrature functions: transformation, weights, and node generation

@inline function ordinate(t::T) where {T<:Real}
    val = tanh(T(π) / 2 * sinh(t))
    # Safety guard: Ensure ordinate never reaches exactly 1.0 due to rounding, 
    # which would cause Inf in singular functions.
    if val >= one(T)
        return prevfloat(one(T))
    elseif val <= -one(T)
        return -prevfloat(one(T))
    end
    return val
end

@inline function weight(t::T) where {T<:Real}
    arg = T(π) / 2 * sinh(t)
    # Stability: cosh can overflow for large t, which would lead to NaN.
    # If the denominator cosh^2(...) would overflow, the weight is effectively 0.
    if abs(arg) > 700.0
        return zero(T)
    end
    tmp = cosh(arg)
    return ((T(π) / 2) * cosh(t)) / (tmp * tmp)
end

@inline function inv_ordinate(t::T) where {T<:Real}
    return asinh(log((one(T) + t) / (one(T) - t)) / T(π))
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
    t = asinh(rhs / T(π))
    for _ in 1:10
        f = T(π) * sinh(t) - t - rhs
        df = T(π) * cosh(t) - 1
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
    tanhsinh(::Type{T}, N::Int) where {T<:AbstractFloat}

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
    if N <= 100
        return SVector{n,T}(x), SVector{n,T}(w), h
    else
        return x, w, h
    end
end

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
