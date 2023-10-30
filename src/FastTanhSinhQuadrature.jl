module FastTanhSinhQuadrature

using StaticArrays
using LambertW

export tanhsinh, integrate, quad

@inline function ordinate(t::T) where {T<:Real}
    return tanh(T(π) / 2 * sinh(t))
end
@inline function weight(t::T) where {T<:Real}
    return ((T(π) / 2) * cosh(t)) / cosh(T(π) / 2 * sinh(t))^2
end

@inline function inv_ordinate(t::T) where {T<:Real}
    return asinh(log((one(T) + t) / (one(T) - t)) / T(π))
end

function tanhsinh(::Type{T}, n::Int) where {T<:AbstractFloat}
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
    h::MVector{1,T}, n::Int) where {T<:Real}
    tmax = inv_ordinate(prevfloat(one(T)))
    h[1] = tmax / n
    for i in 1:n
        t = i * h[1]
        x[i] = ordinate(t)
        w[i] = weight(t)
    end
    return nothing
end


function integrate(::Type{T}, f::Function, n::Int) where {T<:Real}
    x, w, h = tanhsinh(T, n)
    s = weight(zero(T)) * f(zero(T))
    for i in 1:length(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
    end
    return h * s
end

function integrate(f::Function, n::Int)
    x, w, h = tanhsinh(n)
    s = weight(zero(Float64)) * f(zero(Float64))
    for i in 1:length(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
    end
    return h * s
end

# [-1,1] by default 1D
function integrate(f::Function, x::AbstractVector{T}, w::AbstractVector{T},
    h::T) where {T<:Real}
    s = weight(zero(T)) * f(zero(T))
    # ncalls[1] += 1
    for i in 1:length(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
        # ncalls[1] += 2
    end
    return h * s
end

# [a,b] 1D
function integrate(f::Function, xmin::T, xmax::T, x::AbstractVector{T},
    w::AbstractVector{T}, h::T) where {T<:Real}
    Δx = (xmax - xmin) / 2
    x₀ = (xmax + xmin) / 2
    s = weight(zero(T)) * f(x₀)
    #ncalls[1] += 1
    for i in 1:length(x)
        xp = x₀ + Δx * x[i]
        xm = x₀ - Δx * x[i]
        if xm > xmin
            s += w[i] * f(xm)
            #ncalls[1] += 1
        end
        if xp < xmax
            s += w[i] * f(xp)
            # ncalls[1] += 1
        end
    end
    return Δx * h * s
end

## 2D
function integrate(f::Function, xmin::SVector{2,T}, xmax::SVector{2,T}, x::AbstractVector{T},
    w::AbstractVector{T}, h::T) where {T<:Real}
    function f1(x1::T) where {T<:Real}
        g1(y::T) where {T} = f(x1, y)
        integrate(g1, xmin[2], xmax[2], x, w, h)
        return res, nc
    end
    g2(x1::T) where {T} = f1(x1)
    res = integrate(g2, xmin[1], xmax[1], x, w, h)
    return res
end

## 3D
function integrate(f::Function, xmin::SVector{3,T}, xmax::SVector{3,T},
    x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real}
    f1(x1::T, y1::T) where {T<:Real} = integrate(z -> f(x1, y1, z), xmin[3], xmax[3], x, w, h)
    f2(x1::T) where {T<:Real} = integrate(y -> f1(x1, y), xmin[2], xmax[2], x, w, h)
    return integrate(x1 -> f2(x1), xmin[1], xmax[1], x, w, h)
end

# convenience function to convert AbsrtactVectors to SVectors
function integrate(f::Function, xmin::AbstractVector{S}, xmax::AbstractVector{S}, x::AbstractVector{T},
    w::AbstractVector{T}, h::T) where {T<:Real,S<:Real}
    n = length(xmin)
    return integrate(f, SVector{n,T}(xmin), SVector{n,T}(xmax), x, w, h)
end

function _integrate(f::Function, D::Int, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real}
    if D == 2
        f2(x1) = quad(y -> f(x1, y), x, w, h)
        return quad(x1 -> f2(x1), x, w, h)
    elseif D == 3
        g1(x1, y1) = quad(z -> f(x1, y1, z), x, w, h)
        g2(x1) = quad(y -> g1(x1, y), x, w, h)
        return quad(x1 -> g2(x1), x, w, h)
    end
    return zero(T)
end

# helper function for generality
function quad(f::Function, xmin::T, xmax::T, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real}
    if xmin == xmax
        return zero(T)
    end

    #ncalls = @MVector [0]
    if xmin > xmax
        return -integrate(f, xmax, xmin, x, w, h)#, ncalls[1]
    end

    if xmin == -one(T) && xmax == one(T)
        return integrate(f, x, w, h)#, ncalls[1]
    else
        return integrate(f, xmin, xmax, x, w, h)#, ncalls[1]
    end
end

# 2D
function quad(f::Function, xmin::SVector{2,T}, xmax::SVector{2,T},
    x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real}
    if any(xmin .== xmax)
        return zero(T)
    end
    #ncalls = @MVector [0]
    if (xmin[1] == -1) && (xmin[2] == -1) && (xmax[1] == 1) && (xmax[2] == -1)
        return _integrate(f, 2, x, w, h)#, ncalls[1]
    else
        return integrate(f, xmin, xmax, x, w, h)#, ncalls[1]
    end
end

# 3D
function quad(f::Function, xmin::SVector{3,T}, xmax::SVector{3,T},
    x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real}
    if any(xmin .== xmax)
        return zero(T)
    end
    #ncalls = @MVector [0]
    if (xmin[1] == -1) && (xmin[2] == -1) && (xmin[3] == -1) && (xmax[1] == 1) && (xmax[2] == -1) && (xmax[3] == -1)
        return _integrate(f, 3, x, w, h)#, ncalls[1]
    else
        return integrate(f, xmin, xmax, x, w, h)#, ncalls[1]
    end
end

function quad(f::Function, xmin::AbstractVector{S}, xmax::AbstractVector{S}, x::AbstractVector{T},
    w::AbstractVector{T}, h::T) where {T<:Real,S<:Real}
    n = length(xmin)
    return quad(f, SVector{n,T}(xmin), SVector{n,T}(xmax), x, w, h)
end

end
