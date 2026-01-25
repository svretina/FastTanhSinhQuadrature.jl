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

function integrate1D(f::Function, N::Int)
    return integrate1D(Float64, f, N)
end

# [-1,1] by default 1D
"""
    integrate1D(f, x, w, h)

Calculate the integral of `f` over `[-1, 1]` using pre-computed nodes `x`, weights `w`, and step size `h`.
"""
function integrate1D(f::X, x::AbstractVector{T}, w::AbstractVector{T},
    h::T) where {T<:Real,X}
    s = T(π) / 2 * f(zero(T))
    @inbounds for i in eachindex(x)
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
    @turbo for i in eachindex(x)
        Δxxi = Δx * x[i]
        xp = x₀ + Δxxi
        xm = x₀ - Δxxi
        s += w[i] * (f(xm) + f(xp))
    end
    return @fastmath Δx * h * s
end
