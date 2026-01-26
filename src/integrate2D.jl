# 2D Integration functions

"""
    integrate2D(f, x, w, h)

Calculate the 2D integral of `f` over `[-1, 1]^2` using pre-computed nodes/weights.
"""
function integrate2D(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    s = T(π)^2 / 4 * f(zero(T), zero(T))
    @inbounds for i in eachindex(x)
        for j in eachindex(x)
            xp = x[i]
            xm = -xp
            yp = x[j]
            ym = -yp
            s += w[i] * w[j] * (f(xm, ym) + f(xp, ym) + f(xm, yp) + f(xp, yp))
        end
    end
    return h^2 * s
end

"""
    integrate2D(f, low, up, x, w, h)

Calculate the 2D integral of `f` over `[low, up]` (with `low, up::SVector{2}`) using pre-computed nodes/weights.
"""
function integrate2D(f::S, low::SVector{2,T}, up::SVector{2,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        Δx = 0.5 * (up[1] - low[1])
        Δy = 0.5 * (up[2] - low[2])
        x₀ = 0.5 * (up[1] + low[1])
        y₀ = 0.5 * (up[2] + low[2])

        w₀ = T(π) / 2
        s = w₀^2 * f(x₀, y₀)
        for k in eachindex(x)
            wk = w[k]
            Δxk = Δx * x[k]
            Δyk = Δy * x[k]
            xk_p, xk_m = x₀ + Δxk, x₀ - Δxk
            yk_p, yk_m = y₀ + Δyk, y₀ - Δyk
            s += wk * w₀ * (f(xk_p, y₀) + f(xk_m, y₀) +
                            f(x₀, yk_p) + f(x₀, yk_m))
        end

        for i in eachindex(x)
            wi = w[i]
            Δxxi = Δx * x[i]
            xi_p = x₀ + Δxxi
            xi_m = x₀ - Δxxi
            inner_s = zero(T)
            for j in eachindex(x)
                wj = w[j]
                Δxyi = Δy * x[j]
                yj_p = y₀ + Δxyi
                yj_m = y₀ - Δxyi
                inner_s += wj * (f(xi_m, yj_m) + f(xi_p, yj_m) + f(xi_m, yj_p) + f(xi_p, yj_p))
            end
            s += wi * inner_s
        end
    end
    return Δx * Δy * h^2 * s
end

"""
    integrate2D_avx(f, low, up, x, w, h)

SIMD-accelerated 2D integration over `[low, up]` Using `LoopVectorization`.
"""
function integrate2D_avx(f::S, low::SVector{2,T}, up::SVector{2,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    Δx = 0.5 * (up[1] - low[1])
    Δy = 0.5 * (up[2] - low[2])
    x₀ = 0.5 * (up[1] + low[1])
    y₀ = 0.5 * (up[2] + low[2])

    w₀ = T(π) / 2
    s = w₀^2 * f(x₀, y₀)
    @turbo for k in eachindex(x)
        wk = w[k]
        Δxk = Δx * x[k]
        Δyk = Δy * x[k]
        s += wk * w₀ * (f(x₀ + Δxk, y₀) + f(x₀ - Δxk, y₀) +
                        f(x₀, y₀ + Δyk) + f(x₀, y₀ - Δyk))
    end

    for i in eachindex(x)
        wi = w[i]
        Δxxi = Δx * x[i]
        xi_p = x₀ + Δxxi
        xi_m = x₀ - Δxxi
        inner_s = zero(T)
        @turbo for j in eachindex(x)
            wj = w[j]
            Δxyi = Δy * x[j]
            yj_p = y₀ + Δxyi
            yj_m = y₀ - Δxyi
            inner_s += wj * (f(xi_m, yj_m) + f(xi_p, yj_m) + f(xi_m, yj_p) + f(xi_p, yj_p))
        end
        s += wi * inner_s
    end
    return Δx * Δy * h^2 * s
end
