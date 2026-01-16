module FastTanhSinhQuadrature

using StaticArrays
using LambertW
using LoopVectorization

export tanhsinh, integrate, integrate_avx, quad

@inline function ordinate(t::T) where {T<:Real}
    return tanh(T(π) / 2 * sinh(t))
end
@inline function weight(t::T) where {T<:Real}
    tmp = cosh(T(π) / 2 * sinh(t))
    return ((T(π) / 2) * cosh(t)) / (tmp * tmp)
end

@inline function inv_ordinate(t::T) where {T<:Real}
    return asinh(log((one(T) + t) / (one(T) - t)) / T(π))
end

"""
    tanhsinh(::Type{T}, n::Int) where {T<:AbstractFloat}

Generate Tanh-Sinh quadrature nodes `x`, weights `w`, and step size `h` for a given floating point type `T` and level `n`.
The number of points generated is approximately `2^n`.
"""
function tanhsinh(::Type{T}, n::Int) where {T<:AbstractFloat}
    tmax = inv_ordinate(prevfloat(one(T)))
    h = tmax / n
    t = h:h:tmax
    x = ordinate.(t)
    w = weight.(t)
    N = length(x)
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
    @inbounds for i in 1:n
        t = i * h[1]
        x[i] = ordinate(t)
        w[i] = weight(t)
    end
    return nothing
end

@inline function remove_endpoints!(xmin::T, xmax::T, x, w) where {T}
    Δx = 0.5(xmax - xmin)
    x₀ = 0.5(xmax + xmin)
    @fastmath @inbounds for i in eachindex(x)
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
    Δx = 0.5(xmax - xmin)
    x₀ = 0.5(xmax + xmin)
    @fastmath @inbounds for i in eachindex(x)
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
    Δx = 0.5(xmax - xmin)
    x₀ = 0.5(xmax + xmin)
    @fastmath @inbounds for i in eachindex(x)
        Δxxi = Δx * x[i]
        xp = x₀ + Δxxi
        xm = x₀ - Δxxi
        if xp ≥ xmax
            x[i] = x₀
            w[i] = zero(T)
        end
    end
end

"""
    integrate(::Type{T}, f::Function, n::Int) where {T<:Real}

Integrate function `f` over `[-1, 1]` using Tanh-Sinh quadrature with level `n` and precision `T`.
"""
function integrate(::Type{T}, f::Function, n::Int) where {T<:Real}
    x, w, h = tanhsinh(T, n)
    s = weight(zero(T)) * f(zero(T))
    for i in eachindex(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
    end
    return h * s
end

function integrate(f::Function, n::Int)
    x, w, h = tanhsinh(n)
    s = weight(zero(Float64)) * f(zero(Float64))
    for i in eachindex(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
    end
    return h * s
end

# [-1,1] by default 1D
"""
    integrate(f::X, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,X}

Integrate function `f` over `[-1, 1]` using pre-computed nodes `x`, weights `w`, and step size `h`.
"""
function integrate(f::X, x::AbstractVector{T}, w::AbstractVector{T},
    h::T) where {T<:Real,X}
    s = weight(zero(T)) * f(zero(T))
    # ncalls[1] += 1
    for i in eachindex(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
        # ncalls[1] += 2
    end
    return h * s
end

# [a,b] 1D
"""
    integrate(f::S, xmin::T, xmax::T, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}

Integrate function `f` over `[xmin, xmax]` using pre-computed nodes `x`, weights `w`, and step size `h`.
"""
@inline function integrate(f::S, xmin::T, xmax::T, x::AbstractVector{T},
    w::AbstractVector{T}, h::T) where {T<:Real,S}
    @fastmath @inbounds begin
        Δx = 0.5(xmax - xmin)
        x₀ = 0.5(xmax + xmin)
        s = weight(zero(T)) * f(x₀)
        #ncalls[1] += 1
        if xmin < xmax    
            for i in eachindex(x)
                xp = x₀ + Δx * x[i]
                xm = x₀ - Δx * x[i]
                if xm > xmin
                    s += w[i] * f(xm)
                end
                if xp < xmax
                    s += w[i] * f(xp)
                end
            end
        else
            for i in eachindex(x)
                xp = x₀ + Δx * x[i]
                xm = x₀ - Δx * x[i]
                if xm < xmin
                    s += w[i] * f(xm)
                end
                if xp > xmax
                    s += w[i] * f(xp)
                end
            end
        end
    end
    return Δx * h * s
end

## carefull, this is unsafe for a function with a singularity at the endpoints
## if you want to use this with a singular function, then first run
## remove_endpoints! on your weights and points and then use this function
"""
    integrate_avx(f::S, xmin::T, xmax::T, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}

SIMD-accelerated integration using `LoopVectorization`. Ensure `f` is compatible with `@turbo`.
"""
@inline function integrate_avx(f::S, xmin::T, xmax::T, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}
    @fastmath Δx = 0.5(xmax - xmin)
    @fastmath x₀ = 0.5(xmax + xmin)
    @fastmath s = weight(zero(T)) * f(x₀)
    @turbo for i in eachindex(x)
        Δxxi = Δx * x[i]
        xp = x₀ + Δxxi
        xm = x₀ - Δxxi
        s += w[i] * (f(xm) + f(xp))
    end
    return @fastmath Δx * h * s
end

# [-1, 1] by default
@inline function integrate_avx(f::S, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}
    μηδεν = zero(T)
    @fastmath s = weight(μηδεν) * f(μηδεν)
    @turbo for i in eachindex(x)
        s += w[i] * (f(-x[i]) + f(x[i]))
    end
    return @fastmath h * s
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

# helper function for generality
"""
    quad(f::X, xmin::T, xmax::T, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,X}

Wrapper for `integrate` that handles domain checks such as `xmin > xmax` (flips sign) or `xmin == xmax` (returns 0).
"""
@inline function quad(f::X, xmin::T, xmax::T, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,X}
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
@inline function quad(f::X, xmin::SVector{2,T}, xmax::SVector{2,T},
    x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,X}
    if (xmin[1] == xmax[1]) || (xmin[2] == xmax[2])
        return zero(T)
    end
    if all(xmin .< xmax)
        return integrate(f, xmin, xmax, x, w, h)
    else
        sign = 1
        @inbounds for i in 1:2
            low[i] = xmin[i]
            up[i] = xmax[i]
            if xmin[i] > xmax[i]
                sign *= -1
                tmp = xmin[i]
                xmin[i] = xmax[i]
                xmax[i] = tmp
            end
        end
        return sign * integrate(f, xmin, xmax, x, w, h)
    end
end

# 3D
@inline function quad(f::S, xmin::SVector{3,T}, xmax::SVector{3,T},
    x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}
    if (xmin[1] == xmax[1]) || (xmin[2] == xmax[2]) || (xmin[3] == xmax[3])
        return zero(T)
    end
    return integrate(f, xmin, xmax, x, w, h)#, ncalls[1]
end

function quad(f::X, xmin::AbstractVector{S}, xmax::AbstractVector{S}, x::AbstractVector{T},
    w::AbstractVector{T}, h::T) where {T<:Real,S<:Real,X}
    n = length(xmin)
    return quad(f, SVector{n,T}(xmin), SVector{n,T}(xmax), x, w, h)
end

"""
    adaptive_integrate(::Type{T}, f::Function, a, b; tol::Real=1e-8, max_n::Int=10) where {T<:Real}

Adaptive integration of `f` over `[a, b]` using Tanh-Sinh quadrature. 
It doubles the number of quadrature points until the relative error is below `tol`.
"""
function adaptive_integrate(::Type{T}, f::Function, a, b; tol::Real=1e-8, max_n::Int=11) where {T<:Real}
    n = 2
    # First estimate
    x, w, h = tanhsinh(T, n)
    a_T = T(a)
    b_T = T(b)
    old_res = integrate(f, a_T, b_T, x, w, h)
    
    for i in 1:(max_n-2)
        n *= 2
        x, w, h = tanhsinh(T, n)
        new_res = integrate(f, a_T, b_T, x, w, h)
        
        err = abs(new_res - old_res)
        if err < tol || (abs(old_res) > zero(T) && err / abs(old_res) < tol)
            return new_res
        end
        old_res = new_res
    end
    # @warn "Adaptive integration reached max levels ($n) without converging to tol=$tol"
    return old_res
end

adaptive_integrate(f::Function, a, b; tol::Float64=1e-8, max_n::Int=10) = adaptive_integrate(Float64, f, a, b, tol=tol, max_n=max_n)

export adaptive_integrate

end

