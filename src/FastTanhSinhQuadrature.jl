module FastTanhSinhQuadrature

using StaticArrays
using LambertW

export tanhsinh, integrate, quad

@inline function ordinate(t::T)::T where {T<:Real}
    return tanh(T(π) / 2 * sinh(t))
end
@inline function weight(t::T)::T where {T<:Real}
    return ((T(π) / 2) * cosh(t)) / cosh(T(π) / 2 * sinh(t))^2
end

@inline function inv_ordinate(t::T)::T where {T<:Real}
    return asinh(log((one(T) + t) / (one(T) - t)) / T(π))
end

function tanhsinh(::Type{T}, n::Int)::Tuple{<:AbstractVector{T},<:AbstractVector{T},
    T} where {T<:AbstractFloat}
    tmax = inv_ordinate(prevfloat(one(T)))
    h = tmax / n
    t = h:h:tmax
    x = ordinate.(t)
    w = weight.(t)
    return x, w, h
end

tanhsinh(n::Int) = tanhsinh(Float64, n)

function tanhsinh_opt(::Type{T}, n::Int, d::Real=π / 2) where {T<:Real}
    tmax = inv_ordinate(prevfloat(one(T)))
    h = hopt(T, n, d)
    if n * h > tmax
        @show n
        throw(BoundsError)
    end
    t = h:h:tmax
    x = ordinate.(t)
    w = weight.(t)
    return x, w, h
end

function hopt(::Type{T}, n::Int, d::Real=π / 2) where {T}
    (T(2) / T(2n + 1)) * lambertw(T(2d) * T(2n + 1))
end

function tanhsinh!(::Type{T}, x::AbstractVector{T}, w::AbstractVector{T},
    h::MVector{1,T}, n::Int)::Nothing where {T<:Real}
    tmax = inv_ordinate(prevfloat(one(T)))
    h[1] = tmax / n
    for i in 1:n
        t = i * h[1]
        x[i] = ordinate(t)
        w[i] = weight(t)
    end
    return nothing
end


function integrate(::Type{T}, f::Function, n::Int)::T where {T<:Real}
    x, w, h = tanhsinh(T, n)
    s = weight(zero(T)) * f(zero(T))
    for i in 1:length(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
    end
    return h * s
end

function integrate(f::Function, n::Int)::Float64
    x, w, h = tanhsinh(n)
    s = weight(zero(Float64)) * f(zero(Float64))
    for i in 1:length(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
    end
    return h * s
end

# [-1,1] by default 1D
function integrate(f::Function, x::AbstractVector{T}, w::AbstractVector{T},
    h::T)::T where {T<:Real}
    s = weight(zero(T)) * f(zero(T))
    for i in 1:length(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
    end
    return h * s
end

# [a,b] 1D
function integrate(f::Function, xmin::T, xmax::T, x::AbstractVector{T},
    w::AbstractVector{T}, h::T)::T where {T}
    Δx = (xmax - xmin) / 2
    x₀ = (xmax + xmin) / 2
    s = weight(zero(T)) * f(x₀)
    for i in 1:length(x)
        xp = x₀ + Δx * x[i]
        xm = x₀ - Δx * x[i]
        (xm > xmin) && (s += w[i] * f(xm))
        (xp < xmax) && (s += w[i] * f(xp))
    end
    return Δx * h * s
end

## 2D
function integrate(f::Function, xmin::SVector{2,T}, xmax::SVector{2,T}, x::AbstractVector{T},
    w::AbstractVector{T}, h::T)::T where {T<:Real}
    f2(x1) = integrate(y -> f(x1, y), xmin[2], xmax[2], x, w, h)
    f3() = integrate(x1 -> f2(x1), xmin[1], xmax[1], x, w, h)
    return f3()
end

## 3D
function integrate(f::Function, xmin::SVector{3,T}, xmax::SVector{3,T},
    x::AbstractVector{T},
    w::AbstractVector{T}, h::T)::T where {T<:Real}
    f1(x1, y1) = integrate(z -> f(x1, y1, z), xmin[3], xmax[3], x, w, h)
    f2(x1) = integrate(y -> f1(x1, y), xmin[2], xmax[2], x, w, h)
    f3() = integrate(x1 -> f2(x1), xmin[1], xmax[1], x, w, h)
    return f3()
end

# convenience function to convert AbsrtactVectors to SVectors
function integrate(f::Function, xmin::AbstractVector{S}, xmax::AbstractVector{S}, x::AbstractVector{T},
    w::AbstractVector{T}, h::T)::T where {T<:Real,S<:Real}
    n = length(xmin)
    return integrate(f, SVector{n,T}(xmin), SVector{n,T}(xmax), x, w, h)
end

# helper function for generality
function quad(f::Function, xmin::T, xmax::T, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real}
    if xmin == xmax
        return zero(T)
    end

    if xmin > xmax
        return -integrate(f, xmax, xmin, x, w, h)
    end

    if xmin == -one(T) && xmax == one(T)
        return integrate(f, x, w, h)
    else
        return integrate(f, xmin, xmax, x, w, h)
    end
end

# 2D
function quad(f::Function, xmin::SVector{2,T}, xmax::SVector{2,T}, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real}
    if xmin == xmax
        return zero(T)
    end

    if xmin > xmax
        return -integrate(f, xmax, xmin, x, w, h)
    end

    if xmin == -one(T) && xmax == one(T)
        return integrate(f, x, w, h)
    else
        return integrate(f, xmin, xmax, x, w, h)
    end
end

# 3D
function quad(f::Function, xmin::SVector{3,T}, xmax::SVector{3,T}, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real}
    if xmin == xmax
        return zero(T)
    end

    if xmin > xmax
        return -integrate(f, xmax, xmin, x, w, h)
    end

    if xmin == -one(T) && xmax == one(T)
        return integrate(f, x, w, h)
    else
        return integrate(f, xmin, xmax, x, w, h)
    end
end

function quad(f::Function, xmin::AbstractVector{T}, xmax::AbstractVector{T}, x::AbstractVector{T},
    w::AbstractVector{T}, h::T)::T where {T<:Real}
    n = length(xmin)
    return quad(f, SVector{n,T}(xmin), SVector{n,T}(xmax), x, w, h)
end

end
