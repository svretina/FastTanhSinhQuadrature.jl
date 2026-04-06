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

# Adaptive cache construction and validation

"""
    _require_cache_levels(cache, max_levels::Int)

Validate that a provided adaptive cache contains at least `max_levels`
refinement tables, returning the cache unchanged on success.
"""
@inline function _require_cache_levels(cache, max_levels::Int)
    length(cache.xs) >= max_levels ||
        throw(ArgumentError("provided adaptive cache supports $(length(cache.xs)) levels, but `max_levels=$(max_levels)` was requested."))
    return cache
end

struct _Adaptive1DCache{T}
    tm::T
    initial_x::NTuple{2,T}
    initial_w::NTuple{2,T}
    initial_c::NTuple{2,T}
    xs::Vector{Vector{T}}
    ws::Vector{Vector{T}}
    cs::Vector{Vector{T}}
end

"""
    _adaptive_1d_window(::Type{T}, kind::Symbol) where {T<:Real}

Return the tanh-sinh truncation window used when building 1D adaptive caches.
Regular caches use `tmax(T)`, while complement-aware caches use `t_w_max(T, 1)`.
"""
@inline _adaptive_1d_window(::Type{T}, kind::Symbol) where {T<:Real} =
    kind === :complement ? t_w_max(T, 1) : tmax(T)

const _ADAPTIVE_1D_CACHES = Dict{Tuple{DataType, Int, Symbol}, Any}()
const _ADAPTIVE_1D_CACHES_LOCK = ReentrantLock()

"""
    _build_adaptive_1d_cache(::Type{T}, max_levels::Int, kind::Symbol) where {T<:Real}

Construct the uncached 1D adaptive tables for `max_levels` refinement steps.
`kind` selects either the regular or complement-aware node window.
"""
function _build_adaptive_1d_cache(::Type{T}, max_levels::Int, kind::Symbol) where {T<:Real}
    max_levels >= 0 || throw(ArgumentError("`max_levels` must be nonnegative."))
    tm = _adaptive_1d_window(T, kind)
    half = _half(T)
    h = tm * half

    h1 = h
    h2 = h + h
    initial_x = (ordinate(h1), ordinate(h2))
    initial_w = (weight(h1), weight(h2))
    initial_c = (ordinate_complement(h1), ordinate_complement(h2))

    xs = Vector{Vector{T}}(undef, max_levels)
    ws = Vector{Vector{T}}(undef, max_levels)
    cs = Vector{Vector{T}}(undef, max_levels)
    for level in 1:max_levels
        h *= half
        n = 1 << level
        x = Vector{T}(undef, n)
        w = Vector{T}(undef, n)
        c = Vector{T}(undef, n)
        @inbounds for i in 1:n
            tk = (2i - 1) * h
            x[i] = ordinate(tk)
            w[i] = weight(tk)
            c[i] = ordinate_complement(tk)
        end
        xs[level] = x
        ws[level] = w
        cs[level] = c
    end
    return _Adaptive1DCache{T}(tm, initial_x, initial_w, initial_c, xs, ws, cs)
end

"""
    _adaptive_1d_cache(::Type{T}, max_levels::Int, kind::Symbol=:regular) where {T<:Real}

Look up or build the shared 1D adaptive cache for the given element type,
level count, and cache `kind`.
"""
function _adaptive_1d_cache(::Type{T}, max_levels::Int, kind::Symbol=:regular) where {T<:Real}
    key = (T, max_levels, kind)
    Base.@lock _ADAPTIVE_1D_CACHES_LOCK begin
        cache = get!(_ADAPTIVE_1D_CACHES, key) do
            _build_adaptive_1d_cache(T, max_levels, kind)
        end
        return cache::_Adaptive1DCache{T}
    end
end

"""
    adaptive_cache_1D(::Type{T}; max_levels::Int=16, complement::Bool=false) where {T<:Real}

Return a reusable cache for `adaptive_integrate_1D` / `adaptive_integrate_1D_cmpl`.

- Use `complement=false` (default) for `adaptive_integrate_1D`.
- Use `complement=true` for `adaptive_integrate_1D_cmpl`.

Passing the returned cache through the `cache=` keyword avoids rebuilding cache
tables during repeated integrations.
"""
function adaptive_cache_1D(::Type{T}; max_levels::Int=16, complement::Bool=false) where {T<:Real}
    kind = complement ? :complement : :regular
    return _adaptive_1d_cache(T, max_levels, kind)
end

struct _AdaptiveTensorCache{T}
    tm::T
    initial_x::NTuple{2,T}
    initial_w::NTuple{2,T}
    xs::Vector{Vector{T}}
    ws::Vector{Vector{T}}
    xp::Vector{Vector{T}}
    xm::Vector{Vector{T}}
    yp::Vector{Vector{T}}
    ym::Vector{Vector{T}}
    zp::Vector{Vector{T}}
    zm::Vector{Vector{T}}
end

const _ADAPTIVE_TENSOR_CACHES = Dict{Tuple{DataType, Int, Int}, Any}()
const _ADAPTIVE_TENSOR_CACHES_LOCK = ReentrantLock()

"""
    _build_adaptive_tensor_cache(::Type{T}, D::Int, max_levels::Int) where {T<:Real}

Construct the uncached tensor-product adaptive tables used by the 2D and 3D
adaptive integrators.
"""
function _build_adaptive_tensor_cache(::Type{T}, D::Int, max_levels::Int) where {T<:Real}
    max_levels >= 0 || throw(ArgumentError("`max_levels` must be nonnegative."))
    tm = tmax(T, D)
    half = _half(T)
    h = tm * half

    h1 = h
    h2 = h + h
    initial_x = (ordinate(h1), ordinate(h2))
    initial_w = (weight(h1), weight(h2))

    xs = Vector{Vector{T}}(undef, max_levels)
    ws = Vector{Vector{T}}(undef, max_levels)
    xp = Vector{Vector{T}}(undef, max_levels)
    xm = Vector{Vector{T}}(undef, max_levels)
    yp = Vector{Vector{T}}(undef, max_levels)
    ym = Vector{Vector{T}}(undef, max_levels)
    zp = Vector{Vector{T}}(undef, max_levels)
    zm = Vector{Vector{T}}(undef, max_levels)
    max_k = 2
    for level in 1:max_levels
        h *= half
        max_k *= 2
        current_x = Vector{T}(undef, max_k)
        current_w = Vector{T}(undef, max_k)
        @inbounds for i in 1:max_k
            ti = T(i) * h
            current_x[i] = ordinate(ti)
            current_w[i] = weight(ti)
        end
        xs[level] = current_x
        ws[level] = current_w
        xp[level] = Vector{T}(undef, max_k)
        xm[level] = Vector{T}(undef, max_k)
        yp[level] = Vector{T}(undef, max_k)
        ym[level] = Vector{T}(undef, max_k)
        zp[level] = Vector{T}(undef, max_k)
        zm[level] = Vector{T}(undef, max_k)
    end
    return _AdaptiveTensorCache{T}(tm, initial_x, initial_w, xs, ws, xp, xm, yp, ym, zp, zm)
end

"""
    _adaptive_tensor_cache(::Type{T}, D::Int, max_levels::Int) where {T<:Real}

Look up or build the shared tensor-product adaptive cache for dimension `D`
and refinement depth `max_levels`.
"""
function _adaptive_tensor_cache(::Type{T}, D::Int, max_levels::Int) where {T<:Real}
    key = (T, D, max_levels)
    Base.@lock _ADAPTIVE_TENSOR_CACHES_LOCK begin
        cache = get!(_ADAPTIVE_TENSOR_CACHES, key) do
            _build_adaptive_tensor_cache(T, D, max_levels)
        end
        return cache::_AdaptiveTensorCache{T}
    end
end

"""
    adaptive_cache_2D(::Type{T}; max_levels::Int=8) where {T<:Real}

Return a reusable cache for `adaptive_integrate_2D`.
"""
adaptive_cache_2D(::Type{T}; max_levels::Int=8) where {T<:Real} =
    _adaptive_tensor_cache(T, 2, max_levels)

"""
    adaptive_cache_3D(::Type{T}; max_levels::Int=5) where {T<:Real}

Return a reusable cache for `adaptive_integrate_3D`.
"""
adaptive_cache_3D(::Type{T}; max_levels::Int=5) where {T<:Real} =
    _adaptive_tensor_cache(T, 3, max_levels)
