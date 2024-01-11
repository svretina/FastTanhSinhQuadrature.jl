module FastTanhSinhQuadrature

using StaticArrays
using LambertW
using LoopVectorization

export tanhsinh, integrate, quad

@inline function ordinate(t::T) where {T<:Real}
    return @fastmath tanh(T(π) / 2 * sinh(t))
end
@inline function weight(t::T) where {T<:Real}
    tmp = cosh(T(π) / 2 * sinh(t))
    return @fastmath ((T(π) / 2) * cosh(t)) / (tmp * tmp)
end

@inline function inv_ordinate(t::T) where {T<:Real}
    return @fastmath asinh(log((one(T) + t) / (one(T) - t)) / T(π))
end

function tanhsinh(::Type{T}, n::Int) where {T<:AbstractFloat}
    tmax = inv_ordinate(prevfloat(one(T)))
    h = tmax / n
    t = h:h:tmax
    x = ordinate.(t)
    w = weight.(t)
    if n < 100
        return SVector{n,T}(x), SVector{n,T}(w), h
    else
        return x, w, h
    end
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

@inline function remove_endpoints!(xmin::T, xmax::T, x, w) where {T}
    Δx = (xmax - xmin) / 2
    x₀ = (xmax + xmin) / 2
    @fastmath @inbounds for i in 1:length(x)
        Δxxi = Δx * x[i]
        xp = x₀ + Δxxi
        xm = x₀ - Δxxi
        if xp ≥ xmax || xm ≤ xmin
            x[i] = x₀
            w[i] = zero(T)
        end
    end
end

@inline function remove_left_endpoint!(xmin::T, xmax::T, x, w) where {T}
    Δx = (xmax - xmin) / 2
    x₀ = (xmax + xmin) / 2
    @fastmath @inbounds for i in 1:length(x)
        Δxxi = Δx * x[i]
        xp = x₀ + Δxxi
        xm = x₀ - Δxxi
        if xm ≤ xmin
            x[i] = x₀
            w[i] = zero(T)
        end
    end
end

@inline function remove_right_endpoint!(xmin::T, xmax::T, x, w) where {T}
    Δx = (xmax - xmin) / 2
    x₀ = (xmax + xmin) / 2
    @fastmath @inbounds for i in 1:length(x)
        Δxxi = Δx * x[i]
        xp = x₀ + Δxxi
        xm = x₀ - Δxxi
        if xp ≥ xmax
            x[i] = x₀
            w[i] = zero(T)
        end
    end
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
@inline function integrate(f::S, xmin::T, xmax::T, x::AbstractVector{T},
    w::AbstractVector{T}, h::T) where {T<:Real,S}
    @fastmath @inbounds begin
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
    end
    return Δx * h * s
end


## carefull, this is unsafe for a function with a singularity at the endpoints
## if you want to use this with a singular function, then first run
## remove_endpoints! on your weights and points and then use this function
@inline function integrate_avx(f::S, xmin::T, xmax::T, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}
    @fastmath Δx = (xmax - xmin) / 2
    @fastmath x₀ = (xmax + xmin) / 2
    @fastmath s = weight(zero(T)) * f(x₀)
    @turbo for i in 1:length(x)
        Δxxi = Δx * x[i]
        xp = x₀ + Δxxi
        xm = x₀ - Δxxi
        s += w[i] * (f(xm) + f(xp))
    end
    return @fastmath Δx * h * s
end

## 2D
@inline function integrate(f::S, xmin::SVector{2,T}, xmax::SVector{2,T},
    x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}
    @inline function f2(x1::T) where {T<:Real}
        return @inbounds integrate_avx(y -> f(x1, y), xmin[2], xmax[2], x, w, h)
    end
    return @inbounds integrate(x1 -> f2(x1), xmin[1], xmax[1], x, w, h)
end

## 3D
@inline function integrate(f::S, xmin::SVector{3,T}, xmax::SVector{3,T},
    x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}
    f1(x1::T, y1::T) where {T<:Real} = @inbounds integrate_avx(z -> f(x1, y1, z), xmin[3], xmax[3], x, w, h)
    f2(x1::T) where {T<:Real} = @inbounds integrate(y -> f1(x1, y), xmin[2], xmax[2], x, w, h)
    return @inbounds integrate(x1 -> f2(x1), xmin[1], xmax[1], x, w, h)
end

# convenience function to convert AbsrtactVectors to SVectors
function integrate(f::X, xmin::AbstractVector{S}, xmax::AbstractVector{S}, x::AbstractVector{T},
    w::AbstractVector{T}, h::T) where {T<:Real,S<:Real,X}
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
@inline function quad(f::Function, xmin::T, xmax::T, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real}
    if xmin == xmax
        return zero(T)
    end

    #ncalls = @MVector [0]
    if xmin > xmax
        return -integrate_avx(f, xmin, xmax, x, w, h)#, ncalls[1]
    end

    if xmin == -one(T) && xmax == one(T)
        return integrate_avx(f, x, w, h)#, ncalls[1]
    else
        return integrate_avx(f, xmin, xmax, x, w, h)#, ncalls[1]
    end
end

# 2D
@inline function quad(f::Function, xmin::SVector{2,T}, xmax::SVector{2,T},
    x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real}
    if (xmin[1] == xmax[1]) || (xmin[2] == xmax[2])
        return zero(T)
    end
    #ncalls = @MVector [0]
    # if (xmin[1] == -1) && (xmin[2] == -1) && (xmax[1] == 1) && (xmax[2] == -1)
    #     return _integrate(f, 2, x, w, h)#, ncalls[1]
    # else
    #     return integrate(f, xmin, xmax, x, w, h)#, ncalls[1]
    # end
    return integrate(f, xmin, xmax, x, w, h)#, ncalls[1]
end

# 3D
@inline function quad(f::S, xmin::SVector{3,T}, xmax::SVector{3,T},
    x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}
    if (xmin[1] == xmax[1]) || (xmin[2] == xmax[2]) || (xmin[3] == xmax[3])
        return zero(T)
    end
    #ncalls = @MVector [0]
    # if (xmin[1] == -1) && (xmin[2] == -1) && (xmin[3] == -1) && (xmax[1] == 1) && (xmax[2] == -1) && (xmax[3] == -1)
    #     return _integrate(f, 3, x, w, h)#, ncalls[1]
    # else
    #     return integrate(f, xmin, xmax, x, w, h)#, ncalls[1]
    # end
    return integrate(f, xmin, xmax, x, w, h)#, ncalls[1]
end

function quad(f::X, xmin::AbstractVector{S}, xmax::AbstractVector{S}, x::AbstractVector{T},
    w::AbstractVector{T}, h::T) where {T<:Real,S<:Real,X}
    n = length(xmin)
    return quad(f, SVector{n,T}(xmin), SVector{n,T}(xmax), x, w, h)
end

end
