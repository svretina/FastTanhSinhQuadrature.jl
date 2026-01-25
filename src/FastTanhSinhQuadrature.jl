module FastTanhSinhQuadrature

using StaticArrays
using LambertW
using LoopVectorization

export tanhsinh, quad, quad_split,
    integrate1D, integrate2D, integrate3D,
    integrate1D_avx, integrate2D_avx, integrate3D_avx,
    adaptive_integrate_1D, adaptive_integrate_2D, adaptive_integrate_3D

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
    integrate1D(f, xmin, xmax, x, w, h)

Calculate the integral of `f` over `[xmin, xmax]` using pre-computed nodes `x`, weights `w`, and step size `h`.
"""
function integrate1D(f::S, xmin::T, xmax::T, x::X,
    w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        Δx = 0.5(xmax - xmin)
        x₀ = 0.5(xmax + xmin)
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
    integrate1D_avx(f, xmin, xmax, x, w, h)

SIMD-accelerated 1D integration over `[xmin, xmax]` Using `LoopVectorization`.
"""
function integrate1D_avx(f::S, xmin::T, xmax::T, x::AbstractVector{T}, w::AbstractVector{T}, h::T) where {T<:Real,S}
    Δx = 0.5(xmax - xmin)
    x₀ = 0.5(xmax + xmin)
    s = T(π) / 2 * f(x₀)
    @turbo for i in eachindex(x)
        Δxxi = Δx * x[i]
        xp = x₀ + Δxxi
        xm = x₀ - Δxxi
        s += w[i] * (f(xm) + f(xp))
    end
    return @fastmath Δx * h * s
end

## 2D

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
    integrate2D(f, xmin, xmax, x, w, h)

Calculate the 2D integral of `f` over `[xmin, xmax]` (with `xmin, xmax::SVector{2}`) using pre-computed nodes/weights.
"""
function integrate2D(f::S, xmin::SVector{2,T}, xmax::SVector{2,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        Δx, Δy = 0.5 .* (xmax .- xmin)
        x₀, y₀ = 0.5 .* (xmax .+ xmin)

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
    integrate2D_avx(f, xmin, xmax, x, w, h)

SIMD-accelerated 2D integration over `[xmin, xmax]` Using `LoopVectorization`.
"""
function integrate2D_avx(f::S, xmin::SVector{2,T}, xmax::SVector{2,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    Δx, Δy = 0.5 .* (xmax .- xmin)
    x₀, y₀ = 0.5 .* (xmax .+ xmin)

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

## 3D
"""
    integrate3D(f, x, w, h)

Calculate the 3D integral of `f` over `[-1, 1]^3` using pre-computed nodes/weights.
"""
function integrate3D(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        # Tanh-Sinh weight at t=0 is π/2
        w0 = T(π) / 2
        w0sq = w0 * w0
        w0cub = w0sq * w0

        # 1. The Origin (0,0,0)
        total_sum = w0cub * f(zero(T), zero(T), zero(T))

        # 2. The Axes (One dimension non-zero)
        # Symmetry allows us to evaluate all 3 axes in one pass
        axis_acc = zero(T)
        for i in eachindex(x)
            val = x[i]
            # Points: (val, 0, 0), (-val, 0, 0), (0, val, 0)...
            axis_points = (f(val, zero(T), zero(T)) + f(-val, zero(T), zero(T))) +
                          (f(zero(T), val, zero(T)) + f(zero(T), -val, zero(T))) +
                          (f(zero(T), zero(T), val) + f(zero(T), zero(T), -val))
            axis_acc += w[i] * axis_points
        end
        total_sum += w0sq * axis_acc

        # 3. The Planes and Octants
        for i in eachindex(x)
            wi = w[i]
            xi = x[i]

            # Sub-sum for planes and octants involving this 'i'
            sub_sum_i = zero(T)

            for j in eachindex(x)
                wj = w[j]
                wiwj = wi * wj
                xj = x[j]

                # --- Planes (Two dimensions non-zero) ---
                # We handle XY, XZ, and YZ planes simultaneously
                planes = (f(xi, xj, zero(T)) + f(-xi, xj, zero(T)) + f(xi, -xj, zero(T)) + f(-xi, -xj, zero(T))) +
                         (f(xi, zero(T), xj) + f(-xi, zero(T), xj) + f(xi, zero(T), -xj) + f(-xi, zero(T), -xj)) +
                         (f(zero(T), xi, xj) + f(zero(T), -xi, xj) + f(zero(T), xi, -xj) + f(zero(T), -xi, -xj))

                total_sum += wiwj * w0 * planes

                # --- Octants (Three dimensions non-zero) ---
                # This is the O(N^3) bottleneck
                octant_acc = zero(T)
                for k in eachindex(x)
                    xk = x[k]
                    # Evaluate all 8 corners of the symmetric octants
                    octant_acc += w[k] * (
                        (f(xi, xj, xk) + f(-xi, xj, xk) + f(xi, -xj, xk) + f(-xi, -xj, xk)) +
                        (f(xi, xj, -xk) + f(-xi, xj, -xk) + f(xi, -xj, -xk) + f(-xi, -xj, -xk))
                    )
                end
                total_sum += wiwj * octant_acc
            end
        end
    end
    return (h^3) * total_sum
end

"""
    integrate3D_avx(f, x, w, h)

SIMD-accelerated 3D integral over `[-1, 1]^3`.
"""
function integrate3D_avx(f::S, x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        # Tanh-Sinh weight at t=0 is π/2
        w0 = T(π) / 2
        w0sq = w0 * w0
        w0cub = w0sq * w0

        # 1. The Origin (0,0,0)
        total_sum = w0cub * f(zero(T), zero(T), zero(T))

        # 2. The Axes (One dimension non-zero)
        # Symmetry allows us to evaluate all 3 axes in one pass
        axis_acc = zero(T)
        @turbo for i in eachindex(x)
            val = x[i]
            # Points: (val, 0, 0), (-val, 0, 0), (0, val, 0)...
            axis_points = (f(val, zero(T), zero(T)) + f(-val, zero(T), zero(T))) +
                          (f(zero(T), val, zero(T)) + f(zero(T), -val, zero(T))) +
                          (f(zero(T), zero(T), val) + f(zero(T), zero(T), -val))
            axis_acc += w[i] * axis_points
        end
        total_sum += w0sq * axis_acc

        # 3. The Planes and Octants
        for i in eachindex(x)
            wi = w[i]
            xi = x[i]

            # Sub-sum for planes and octants involving this 'i'
            sub_sum_i = zero(T)

            for j in eachindex(x)
                wj = w[j]
                wiwj = wi * wj
                xj = x[j]

                # --- Planes (Two dimensions non-zero) ---
                # We handle XY, XZ, and YZ planes simultaneously
                planes = (f(xi, xj, zero(T)) + f(-xi, xj, zero(T)) + f(xi, -xj, zero(T)) + f(-xi, -xj, zero(T))) +
                         (f(xi, zero(T), xj) + f(-xi, zero(T), xj) + f(xi, zero(T), -xj) + f(-xi, zero(T), -xj)) +
                         (f(zero(T), xi, xj) + f(zero(T), -xi, xj) + f(zero(T), xi, -xj) + f(zero(T), -xi, -xj))

                total_sum += wiwj * w0 * planes

                # --- Octants (Three dimensions non-zero) ---
                # This is the O(N^3) bottleneck
                octant_acc = zero(T)
                @turbo for k in eachindex(x)
                    xk = x[k]
                    # Evaluate all 8 corners of the symmetric octants
                    octant_acc += w[k] * (
                        (f(xi, xj, xk) + f(-xi, xj, xk) + f(xi, -xj, xk) + f(-xi, -xj, xk)) +
                        (f(xi, xj, -xk) + f(-xi, xj, -xk) + f(xi, -xj, -xk) + f(-xi, -xj, -xk))
                    )
                end
                total_sum += wiwj * octant_acc
            end
        end
    end
    return (h^3) * total_sum
end

# [a,b]x[c,d]x[e,f]
"""
    integrate3D(f, xmin, xmax, x, w, h)

Calculate the 3D integral of `f` over `[xmin, xmax]` (with `xmin, xmax::SVector{3}`).
"""
function integrate3D(f::S, xmin::SVector{3,T}, xmax::SVector{3,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        # 1. Pre-calculate constants and Jacobians
        Δx, Δy, Δz = 0.5 .* (xmax .- xmin)
        x₀, y₀, z₀ = 0.5 .* (xmax .+ xmin)
        w₀ = T(π) / 2
        # 1. THE ORIGIN (0, 0, 0)
        w₀² = w₀^2
        w₀³ = w₀² * w₀
        total_sum = w₀³ * f(x₀, y₀, z₀)

        # 2. THE AXES (One dimension is non-zero, two are zero)
        for i in eachindex(x)
            wi = w[i]
            dx, dy, dz = Δx * x[i], Δy * x[i], Δz * x[i]

            # X-axis, Y-axis, Z-axis
            axis_sum = (f(x₀ + dx, y₀, z₀) + f(x₀ - dx, y₀, z₀)) +
                       (f(x₀, y₀ + dy, z₀) + f(x₀, y₀ - dy, z₀)) +
                       (f(x₀, y₀, z₀ + dz) + f(x₀, y₀, z₀ - dz))
            total_sum += wi * w₀² * axis_sum
        end

        # 5. Nested Loops: Planes and Octants
        for i in eachindex(x)
            wi = w[i]
            dx = Δx * x[i]
            xp, xm = x₀ + dx, x₀ - dx

            # Plane sums that only depend on i
            for j in eachindex(x)
                wj = w[j]
                wiwj = wi * wj
                dy = Δy * x[j]
                yp, ym = y₀ + dy, y₀ - dy
                dz_j = Δz * x[j] # For the planes where j is the Z-index
                dx_j = Δx * x[j] # For the planes where j is the X-index

                # Plane contributions (XY, XZ, YZ)
                plane_sum = (f(xp, yp, z₀) + f(xm, yp, z₀) + f(xp, ym, z₀) + f(xm, ym, z₀)) +
                            (f(xp, y₀, z₀ + dz_j) + f(xm, y₀, z₀ + dz_j) + f(xp, y₀, z₀ - dz_j) + f(xm, y₀, z₀ - dz_j)) +
                            (f(x₀, yp, z₀ + dx_j) + f(x₀, ym, z₀ + dx_j) + f(x₀, yp, z₀ - dx_j) + f(x₀, ym, z₀ - dx_j))

                total_sum += wiwj * w₀ * plane_sum

                # Inner Loop: Octants
                octant_sum = zero(T)
                for k in eachindex(x)
                    wk = w[k]
                    dz = Δz * x[k]
                    zp, zm = z₀ + dz, z₀ - dz

                    # 8 points of the symmetry
                    octant_sum += wk * (
                        (f(xp, yp, zp) + f(xm, yp, zp) + f(xp, ym, zp) + f(xm, ym, zp)) +
                        (f(xp, yp, zm) + f(xm, yp, zm) + f(xp, ym, zm) + f(xm, ym, zm))
                    )
                end
                total_sum += wiwj * octant_sum
            end
        end
    end

    return (Δx * Δy * Δz * h^3) * total_sum
end

"""
    integrate3D_avx(f, xmin, xmax, x, w, h)

SIMD-accelerated 3D integral over `[xmin, xmax]`.
"""
function integrate3D_avx(f::S, xmin::SVector{3,T}, xmax::SVector{3,T},
    x::X, w::W, h::T) where {T<:Real,S,X<:AbstractVector{T},W<:AbstractVector{T}}
    @inbounds begin
        # 1. Pre-calculate constants and Jacobians
        Δx, Δy, Δz = 0.5 .* (xmax .- xmin)
        x₀, y₀, z₀ = 0.5 .* (xmax .+ xmin)
        w₀ = T(π) / 2
        # 1. THE ORIGIN (0, 0, 0)
        w₀² = w₀^2
        w₀³ = w₀² * w₀
        total_sum = w₀³ * f(x₀, y₀, z₀)

        # 2. THE AXES (One dimension is non-zero, two are zero)
        @turbo for i in eachindex(x)
            wi = w[i]
            dx, dy, dz = Δx * x[i], Δy * x[i], Δz * x[i]

            # X-axis, Y-axis, Z-axis
            axis_sum = (f(x₀ + dx, y₀, z₀) + f(x₀ - dx, y₀, z₀)) +
                       (f(x₀, y₀ + dy, z₀) + f(x₀, y₀ - dy, z₀)) +
                       (f(x₀, y₀, z₀ + dz) + f(x₀, y₀, z₀ - dz))
            total_sum += wi * w₀² * axis_sum
        end

        # 5. Nested Loops: Planes and Octants
        for i in eachindex(x)
            wi = w[i]
            dx = Δx * x[i]
            xp, xm = x₀ + dx, x₀ - dx

            # Plane sums that only depend on i
            for j in eachindex(x)
                wj = w[j]
                wiwj = wi * wj
                dy = Δy * x[j]
                yp, ym = y₀ + dy, y₀ - dy
                dz_j = Δz * x[j] # For the planes where j is the Z-index
                dx_j = Δx * x[j] # For the planes where j is the X-index

                # Plane contributions (XY, XZ, YZ)
                plane_sum = (f(xp, yp, z₀) + f(xm, yp, z₀) + f(xp, ym, z₀) + f(xm, ym, z₀)) +
                            (f(xp, y₀, z₀ + dz_j) + f(xm, y₀, z₀ + dz_j) + f(xp, y₀, z₀ - dz_j) + f(xm, y₀, z₀ - dz_j)) +
                            (f(x₀, yp, z₀ + dx_j) + f(x₀, ym, z₀ + dx_j) + f(x₀, yp, z₀ - dx_j) + f(x₀, ym, z₀ - dx_j))

                total_sum += wiwj * w₀ * plane_sum

                # Inner Loop: Octants
                octant_sum = zero(T)
                @turbo for k in eachindex(x)
                    wk = w[k]
                    dz = Δz * x[k]
                    zp, zm = z₀ + dz, z₀ - dz

                    # 8 points of the symmetry
                    octant_sum += wk * (
                        (f(xp, yp, zp) + f(xm, yp, zp) + f(xp, ym, zp) + f(xm, ym, zp)) +
                        (f(xp, yp, zm) + f(xm, yp, zm) + f(xp, ym, zm) + f(xm, ym, zm))
                    )
                end
                total_sum += wiwj * octant_sum
            end
        end
    end

    return (Δx * Δy * Δz * h^3) * total_sum
end

"""
    quad(f::Function, [xmin, xmax]; tol=1e-12, max_levels=10)

High-level interface for Tanh-Sinh quadrature. Automatically detects dimensions (1D, 2D, or 3D) 
and chooses the most efficient implementation.

- If `xmin` and `xmax` are absent, integrates over `[-1, 1]` (1D).
- Uses SIMD-accelerated (`_avx`) implementation for `Float32/Float64`.
- Uses higher-precision implementation for other types (e.g., `BigFloat`).
- Automatically converts `AbstractVector` to `SVector` for multi-dimensional integration.
"""
function quad(f::Function, xmin::T, xmax::T; tol::Real=1e-12, max_levels::Int=10) where {T<:Real}
    if xmin == xmax
        return zero(T)
    end
    if xmin > xmax
        return -quad(f, xmax, xmin; tol=tol, max_levels=max_levels)
    end

    # 1D Case
    return adaptive_integrate_1D(T, f, xmin, xmax; tol=tol, max_levels=max_levels)
end

# Default 1D
quad(f::Function; tol::Real=1e-12, max_levels::Int=10) = quad(f, -1.0, 1.0; tol=tol, max_levels=max_levels)

# Multi-dimensional cases
function quad(f::Function, xmin::AbstractVector{T}, xmax::AbstractVector{T}; tol::Real=1e-10, max_levels::Int=8) where {T<:Real}
    n = length(xmin)
    if n == 2
        return quad(f, SVector{2,T}(xmin), SVector{2,T}(xmax); tol=tol, max_levels=max_levels)
    elseif n == 3
        return quad(f, SVector{3,T}(xmin), SVector{3,T}(xmax); tol=tol, max_levels=max_levels)
    else
        error("Higher than 3D integration is not yet supported.")
    end
end

function quad(f::Function, xmin::SVector{2,T}, xmax::SVector{2,T}; tol::Real=1e-10, max_levels::Int=8) where {T<:Real}
    if any(xmin .== xmax)
        return zero(T)
    end
    # Check bounds and flip if necessary
    sign_flip = one(T)
    x1_min, x1_max = xmin[1], xmax[1]
    x2_min, x2_max = xmin[2], xmax[2]

    if x1_min > x1_max
        x1_min, x1_max = x1_max, x1_min
        sign_flip *= -one(T)
    end
    if x2_min > x2_max
        x2_min, x2_max = x2_max, x2_min
        sign_flip *= -one(T)
    end

    if sign_flip == -one(T)
        return -quad(f, SVector(x1_min, x2_min), SVector(x1_max, x2_max); tol=tol, max_levels=max_levels)
    end

    # If we made swaps but sign is positive, we just call adaptive with swapped bounds
    # but strictly speaking adaptive integrate expects ordered bounds? 
    # Yes, typically. So let's construct ordered SVectors.
    return adaptive_integrate_2D(T, f, SVector(x1_min, x2_min), SVector(x1_max, x2_max); tol=tol, max_levels=max_levels)
end

function quad(f::Function, xmin::SVector{3,T}, xmax::SVector{3,T}; tol::Real=1e-8, max_levels::Int=5) where {T<:Real}
    if any(xmin .== xmax)
        return zero(T)
    end

    sign_flip = one(T)
    x1_min, x1_max = xmin[1], xmax[1]
    x2_min, x2_max = xmin[2], xmax[2]
    x3_min, x3_max = xmin[3], xmax[3]

    if x1_min > x1_max
        x1_min, x1_max = x1_max, x1_min
        sign_flip *= -one(T)
    end
    if x2_min > x2_max
        x2_min, x2_max = x2_max, x2_min
        sign_flip *= -one(T)
    end
    if x3_min > x3_max
        x3_min, x3_max = x3_max, x3_min
        sign_flip *= -one(T)
    end

    if sign_flip == -one(T)
        return -quad(f, SVector(x1_min, x2_min, x3_min), SVector(x1_max, x2_max, x3_max); tol=tol, max_levels=max_levels)
    end

    return adaptive_integrate_3D(T, f, SVector(x1_min, x2_min, x3_min), SVector(x1_max, x2_max, x3_max); tol=tol, max_levels=max_levels)
end

"""
    quad_split(f, c, [xmin, xmax]; tol=1e-12, max_levels=10)

Split the integration domain at point `c` (singularity) and integrate sub-domains separately.
"""
function quad_split(f::Function, c::T, xmin::T, xmax::T; tol::Real=1e-12, max_levels::Int=10) where {T<:Real}
    # [xmin, c] + [c, xmax]
    return quad(f, xmin, c; tol=tol / 2, max_levels=max_levels) +
           quad(f, c, xmax; tol=tol / 2, max_levels=max_levels)
end

quad_split(f::Function, c::T; tol::Real=1e-12, max_levels::Int=10) where {T<:Real} = quad_split(f, c, -one(T), one(T); tol=tol, max_levels=max_levels)

function quad_split(f::Function, c::SVector{2,T}, xmin::SVector{2,T}, xmax::SVector{2,T}; tol::Real=1e-10, max_levels::Int=8) where {T<:Real}
    # Split into 4 quadrants
    # Q1: [xmin[1], c[1]] x [xmin[2], c[2]]
    # Q2: [c[1], xmax[1]] x [xmin[2], c[2]]
    # Q3: [xmin[1], c[1]] x [c[2], xmax[2]]
    # Q4: [c[1], xmax[1]] x [c[2], xmax[2]]
    cx, cy = c[1], c[2]
    x1, x2 = xmin[1], xmax[1]
    y1, y2 = xmin[2], xmax[2]

    q1 = quad(f, SVector(x1, y1), SVector(cx, cy); tol=tol / 4, max_levels=max_levels)
    q2 = quad(f, SVector(cx, y1), SVector(x2, cy); tol=tol / 4, max_levels=max_levels)
    q3 = quad(f, SVector(x1, cy), SVector(cx, y2); tol=tol / 4, max_levels=max_levels)
    q4 = quad(f, SVector(cx, cy), SVector(x2, y2); tol=tol / 4, max_levels=max_levels)

    return q1 + q2 + q3 + q4
end

function quad_split(f::Function, c::SVector{3,T}, xmin::SVector{3,T}, xmax::SVector{3,T}; tol::Real=1e-8, max_levels::Int=5) where {T<:Real}
    # Split into 8 octants
    cx, cy, cz = c[1], c[2], c[3]
    x1, x2 = xmin[1], xmax[1]
    y1, y2 = xmin[2], xmax[2]
    z1, z2 = xmin[3], xmax[3]

    xs = (x1, cx, x2)
    ys = (y1, cy, y2)
    zs = (z1, cz, z2)

    total = zero(T)
    for i in 1:2, j in 1:2, k in 1:2
        total += quad(f, SVector(xs[i], ys[j], zs[k]), SVector(xs[i+1], ys[j+1], zs[k+1]); tol=tol / 8, max_levels=max_levels)
    end
    return total
end

"""
    adaptive_integrate_1D(::Type{T}, f::Function, a, b; tol::Real=1e-12, max_levels::Int=10)

Adaptive 1D Tanh-Sinh integration over `[a, b]`. Starts with a coarse grid (h ≈ tmax/2) and halves 
the step size at each level. Reuses function evaluations from previous levels by only computing
new (odd-indexed) nodes. Exploits symmetry around the center of the interval.
"""
function adaptive_integrate_1D(::Type{T}, f::S, a, b;
    tol::Real=1e-12, max_levels::Int=10) where {T<:Real,S}
    a_T, b_T = T(a), T(b)
    Δx = 0.5 * (b_T - a_T)
    x₀ = 0.5 * (b_T + a_T)

    # Initial Grid (Level 0)
    tm = tmax(T)
    h = tm / 2

    # Weight at t=0 is π/2
    w0 = T(π) / 2
    s_origin = w0 * f(x₀)

    # Level 0 includes points at t=1h and t=2h
    s_weighted = zero(T)
    for k in 1:2
        tk = k * h
        s_weighted += weight(tk) * (f(x₀ + Δx * ordinate(tk)) + f(x₀ - Δx * ordinate(tk)))
    end

    total_weighted_sum = s_origin + s_weighted
    old_res = Δx * h * total_weighted_sum

    for level in 1:max_levels
        h /= 2
        s_new = zero(T)

        # New points are at ODD multiples of the new h: 1h, 3h, 5h...
        k = 1
        while true
            tk = k * h
            tk > tm && break
            s_new += weight(tk) * (f(x₀ + Δx * ordinate(tk)) + f(x₀ - Δx * ordinate(tk)))
            k += 2
        end

        total_weighted_sum += s_new
        new_res = Δx * h * total_weighted_sum

        if abs(new_res - old_res) < tol * max(one(T), abs(new_res))
            return new_res
        end
        old_res = new_res
    end
    return old_res
end

"""
    adaptive_integrate_2D(::Type{T}, f::Function, xmin::SVector{2,T}, xmax::SVector{2,T}; tol::Real=1e-10, max_levels::Int=8)

Adaptive 2D Tanh-Sinh integration over a rectangle. Reuses indices by only evaluating new points 
where at least one coordinate corresponds to an odd multiple of the halved step size `h`. 
Exploits 4-way quadrant symmetry and 2-way axis symmetry.
"""
function adaptive_integrate_2D(::Type{T}, f::S, xmin::SVector{2,T}, xmax::SVector{2,T};
    tol::Real=1e-10, max_levels::Int=8) where {T<:Real,S}
    Δx, Δy = 0.5 .* (xmax .- xmin)
    x₀, y₀ = 0.5 .* (xmax .+ xmin)
    tm = tmax(T, 2)
    h = tm / 2
    w0 = T(π) / 2

    # Helper to evaluate symmetric 4 quadrant points
    function eval_quadrants(ti, tj, wi, wj)
        xi, yi = ordinate(ti), ordinate(tj)
        dx, dy = Δx * xi, Δy * yi
        return wi * wj * (f(x₀ + dx, y₀ + dy) + f(x₀ - dx, y₀ + dy) +
                          f(x₀ + dx, y₀ - dy) + f(x₀ - dx, y₀ - dy))
    end

    # Helper to evaluate symmetric axis points
    function eval_axes(tk, wk)
        val = ordinate(tk)
        dx, dy = Δx * val, Δy * val
        return wk * w0 * (f(x₀ + dx, y₀) + f(x₀ - dx, y₀) +
                          f(x₀, y₀ + dy) + f(x₀, y₀ - dy))
    end

    # Initial Level 0 (h, 2h)
    s_total = (w0^2) * f(x₀, y₀)
    for i in 1:2, j in 1:2
        s_total += eval_quadrants(i * h, j * h, weight(i * h), weight(j * h))
    end
    for k in 1:2
        s_total += eval_axes(k * h, weight(k * h))
    end

    old_res = Δx * Δy * h^2 * s_total

    for level in 1:max_levels
        h /= 2
        s_new = zero(T)
        max_k = floor(Int, tm / h)

        for i in 1:max_k
            wi, ti = weight(i * h), i * h
            for j in 1:max_k
                # Point is new if at least one index is odd
                (iseven(i) && iseven(j)) && continue
                s_new += eval_quadrants(ti, j * h, wi, weight(j * h))
            end
            if isodd(i)
                s_new += eval_axes(ti, wi)
            end
        end

        s_total += s_new
        new_res = Δx * Δy * h^2 * s_total

        if abs(new_res - old_res) < tol * max(one(T), abs(new_res))
            return new_res
        end
        old_res = new_res
    end
    return old_res
end

"""
    adaptive_integrate_3D(::Type{T}, f::Function, xmin::SVector{3,T}, xmax::SVector{3,T}; tol::Real=1e-8, max_levels::Int=5)

Adaptive 3D Tanh-Sinh integration over a box. Reuses old points and exploits 8-way octant 
symmetry, 4-way plane symmetry, and 2-way axis symmetry to minimize function evaluations.
"""
function adaptive_integrate_3D(::Type{T}, f::S, xmin::SVector{3,T}, xmax::SVector{3,T};
    tol::Real=1e-8, max_levels::Int=5) where {T<:Real,S}
    Δx, Δy, Δz = 0.5 .* (xmax .- xmin)
    x₀, y₀, z₀ = 0.5 .* (xmax .+ xmin)
    tm = tmax(T, 3)
    h = tm / 2
    w₀ = T(π) / 2

    # Evaluate a single point in the octant (8 reflections)
    function add_octant(ti, tj, tk, wi, wj, wk)
        vi, vj, vk = ordinate(ti), ordinate(tj), ordinate(tk)
        dx, dy, dz = Δx * vi, Δy * vj, Δz * vk
        w = wi * wj * wk
        return w * (
            (f(x₀ + dx, y₀ + dy, z₀ + dz) + f(x₀ - dx, y₀ + dy, z₀ + dz) +
             f(x₀ + dx, y₀ - dy, z₀ + dz) + f(x₀ - dx, y₀ - dy, z₀ + dz)) +
            (f(x₀ + dx, y₀ + dy, z₀ - dz) + f(x₀ - dx, y₀ + dy, z₀ - dz) +
             f(x₀ + dx, y₀ - dy, z₀ - dz) + f(x₀ - dx, y₀ - dy, z₀ - dz))
        )
    end

    # Evaluate points on the 3 planes (XY, XZ, YZ) - 4 reflections each
    function add_planes(ti, tj, wi, wj)
        vi, vj = ordinate(ti), ordinate(tj)
        dx_i, dy_i, dz_i = Δx * vi, Δy * vi, Δz * vi
        dx_j, dy_j, dz_j = Δx * vj, Δy * vj, Δz * vj
        w = wi * wj * w₀
        return w * (
            (f(x₀ + dx_i, y₀ + dy_j, z₀) + f(x₀ - dx_i, y₀ + dy_j, z₀) + f(x₀ + dx_i, y₀ - dy_j, z₀) + f(x₀ - dx_i, y₀ - dy_j, z₀)) +
            (f(x₀ + dx_i, y₀, z₀ + dz_j) + f(x₀ - dx_i, y₀, z₀ + dz_j) + f(x₀ + dx_i, y₀, z₀ - dz_j) + f(x₀ - dx_i, y₀, z₀ - dz_j)) +
            (f(x₀, y₀ + dy_i, z₀ + dz_j) + f(x₀, y₀ - dy_i, z₀ + dz_j) + f(x₀, y₀ + dy_i, z₀ - dz_j) + f(x₀, y₀ - dy_i, z₀ - dz_j))
        )
    end

    # Evaluate points on the 3 axes (X, Y, Z) - 2 reflections each
    function add_axes(ti, wi)
        vi = ordinate(ti)
        dx, dy, dz = Δx * vi, Δy * vi, Δz * vi
        w = wi * w₀^2
        return w * (
            (f(x₀ + dx, y₀, z₀) + f(x₀ - dx, y₀, z₀)) +
            (f(x₀, y₀ + dy, z₀) + f(x₀, y₀ - dy, z₀)) +
            (f(x₀, y₀, z₀ + dz) + f(x₀, y₀, z₀ - dz))
        )
    end

    # Initial Level 0 (k=1, 2)
    s_total = (w₀^3) * f(x₀, y₀, z₀)
    for i in 1:2, j in 1:2, k in 1:2
        s_total += add_octant(i * h, j * h, k * h, weight(i * h), weight(j * h), weight(k * h))
    end
    for i in 1:2, j in 1:2
        s_total += add_planes(i * h, j * h, weight(i * h), weight(j * h))
    end
    for i in 1:2
        s_total += add_axes(i * h, weight(i * h))
    end

    old_res = Δx * Δy * Δz * h^3 * s_total

    for level in 1:max_levels
        h /= 2
        s_new = zero(T)
        max_k = floor(Int, tm / h)

        for i in 1:max_k
            wi, ti = weight(i * h), i * h
            for j in 1:max_k
                wj, tj = weight(j * h), j * h
                for k in 1:max_k
                    (iseven(i) && iseven(j) && iseven(k)) && continue
                    s_new += add_octant(ti, tj, k * h, wi, wj, weight(k * h))
                end
                if !(iseven(i) && iseven(j))
                    s_new += add_planes(ti, tj, wi, wj)
                end
            end
            if isodd(i)
                s_new += add_axes(ti, wi)
            end
        end

        s_total += s_new
        new_res = Δx * Δy * Δz * h^3 * s_total

        if abs(new_res - old_res) < tol * max(one(T), abs(new_res))
            return new_res
        end
        old_res = new_res
    end
    return old_res
end

end
