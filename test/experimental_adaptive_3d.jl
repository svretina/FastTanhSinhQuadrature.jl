using BenchmarkTools
using FastTanhSinhQuadrature
using LoopVectorization
using StaticArrays
using StrideArrays

struct PackageStyleLevelData{X,W,XP,XM,YP,YM,ZP,ZM,O}
    x::X
    w::W
    xp::XP
    xm::XM
    yp::YP
    ym::YM
    zp::ZP
    zm::ZM
    owner::O
end

struct PackageStyleVariant{T,L}
    name::String
    tm::T
    initial_x::NTuple{2,T}
    initial_w::NTuple{2,T}
    levels::L
end

@inline to_svector(::Val{N}, v::AbstractVector{T}) where {N,T} =
    SVector{N,T}(ntuple(i -> v[i], Val(N)))

@inline level_node_count(::Val{L}) where {L} = 1 << (L + 1)

function build_level_data(::Val{:vector}, ::Val{L}, x_src, w_src) where {L}
    T = eltype(x_src)
    n = length(x_src)
    x = collect(x_src)
    w = collect(w_src)
    xp = Vector{T}(undef, n)
    xm = Vector{T}(undef, n)
    yp = Vector{T}(undef, n)
    ym = Vector{T}(undef, n)
    zp = Vector{T}(undef, n)
    zm = Vector{T}(undef, n)
    return PackageStyleLevelData(x, w, xp, xm, yp, ym, zp, zm, nothing)
end

function build_level_data(::Val{:ptr_vector}, ::Val{L}, x_src, w_src) where {L}
    T = eltype(x_src)
    n = length(x_src)
    x = collect(x_src)
    w = collect(w_src)
    xp = Vector{T}(undef, n)
    xm = Vector{T}(undef, n)
    yp = Vector{T}(undef, n)
    ym = Vector{T}(undef, n)
    zp = Vector{T}(undef, n)
    zm = Vector{T}(undef, n)
    owner = (x, w, xp, xm, yp, ym, zp, zm)
    return PackageStyleLevelData(
        PtrArray(x),
        PtrArray(w),
        PtrArray(xp),
        PtrArray(xm),
        PtrArray(yp),
        PtrArray(ym),
        PtrArray(zp),
        PtrArray(zm),
        owner,
    )
end

function build_level_data(::Val{:stride}, ::Val{L}, x_src, w_src) where {L}
    T = eltype(x_src)
    n = length(x_src)
    x = StrideArray(collect(x_src))
    w = StrideArray(collect(w_src))
    xp = StrideArray{T}(undef, (n,))
    xm = StrideArray{T}(undef, (n,))
    yp = StrideArray{T}(undef, (n,))
    ym = StrideArray{T}(undef, (n,))
    zp = StrideArray{T}(undef, (n,))
    zm = StrideArray{T}(undef, (n,))
    return PackageStyleLevelData(x, w, xp, xm, yp, ym, zp, zm, nothing)
end

function build_level_data(::Val{:ptr_stride}, ::Val{L}, x_src, w_src) where {L}
    T = eltype(x_src)
    n = length(x_src)
    x = StrideArray(collect(x_src))
    w = StrideArray(collect(w_src))
    xp = StrideArray{T}(undef, (n,))
    xm = StrideArray{T}(undef, (n,))
    yp = StrideArray{T}(undef, (n,))
    ym = StrideArray{T}(undef, (n,))
    zp = StrideArray{T}(undef, (n,))
    zm = StrideArray{T}(undef, (n,))
    owner = (x, w, xp, xm, yp, ym, zp, zm)
    return PackageStyleLevelData(
        PtrArray(x),
        PtrArray(w),
        PtrArray(xp),
        PtrArray(xm),
        PtrArray(yp),
        PtrArray(ym),
        PtrArray(zp),
        PtrArray(zm),
        owner,
    )
end

function build_level_data(::Val{:hybrid_svector64}, ::Val{L}, x_src, w_src) where {L}
    T = eltype(x_src)
    N = level_node_count(Val(L))
    xp = Vector{T}(undef, N)
    xm = Vector{T}(undef, N)
    yp = Vector{T}(undef, N)
    ym = Vector{T}(undef, N)
    zp = Vector{T}(undef, N)
    zm = Vector{T}(undef, N)
    if N <= 64
        x = to_svector(Val(N), x_src)
        w = to_svector(Val(N), w_src)
        return PackageStyleLevelData(x, w, xp, xm, yp, ym, zp, zm, nothing)
    end
    x = collect(x_src)
    w = collect(w_src)
    return PackageStyleLevelData(x, w, xp, xm, yp, ym, zp, zm, nothing)
end

function build_level_data(::Val{:hybrid_svector128}, ::Val{L}, x_src, w_src) where {L}
    T = eltype(x_src)
    N = level_node_count(Val(L))
    xp = Vector{T}(undef, N)
    xm = Vector{T}(undef, N)
    yp = Vector{T}(undef, N)
    ym = Vector{T}(undef, N)
    zp = Vector{T}(undef, N)
    zm = Vector{T}(undef, N)
    if N <= 128
        x = to_svector(Val(N), x_src)
        w = to_svector(Val(N), w_src)
        return PackageStyleLevelData(x, w, xp, xm, yp, ym, zp, zm, nothing)
    end
    x = collect(x_src)
    w = collect(w_src)
    return PackageStyleLevelData(x, w, xp, xm, yp, ym, zp, zm, nothing)
end

function build_variant(::Type{T}, name::String, style::Symbol, max_levels::Int) where {T}
    cache = adaptive_cache_3D(T; max_levels=max_levels)
    levels = ntuple(level -> build_level_data(Val(style), Val(level), cache.xs[level], cache.ws[level]), max_levels)
    return PackageStyleVariant{T,typeof(levels)}(
        name,
        cache.tm,
        cache.initial_x,
        cache.initial_w,
        levels,
    )
end

@inline function package_style_level_sum(
    f::F, level, Δx::T, x0::T, Δy::T, y0::T, Δz::T, z0::T, w0::T, w0sq::T,
) where {F,T}
    s_new = zero(T)
    xs = level.x
    ws = level.w
    xp = level.xp
    xm = level.xm
    yp = level.yp
    ym = level.ym
    zp = level.zp
    zm = level.zm
    n = length(xs)
    odd_count = (n + 1) >>> 1

    @turbo warn_check_args=false for i in 1:n
        vx = Δx * xs[i]
        xp[i] = x0 + vx
        xm[i] = x0 - vx
        vy = Δy * xs[i]
        yp[i] = y0 + vy
        ym[i] = y0 - vy
        vz = Δz * xs[i]
        zp[i] = z0 + vz
        zm[i] = z0 - vz
    end

    @inbounds begin
        tile = 32

        for i in 1:2:n
            wi = ws[i]
            xpi, xmi = xp[i], xm[i]
            ypi, ymi = yp[i], ym[i]
            zpi, zmi = zp[i], zm[i]
            for j_tile in 1:tile:n
                j_end = min(j_tile + tile - 1, n)
                for k_tile in 1:tile:n
                    k_end = min(k_tile + tile - 1, n)
                    for j in j_tile:j_end
                        wj = ws[j]
                        ypj, ymj = yp[j], ym[j]
                        wiwj = wi * wj
                        inner_sum = zero(T)
                        @turbo warn_check_args=false for k in k_tile:k_end
                            zpk = zp[k]
                            zmk = zm[k]
                            inner_sum += ws[k] * (
                                (f(xpi, ypj, zpk) + f(xmi, ypj, zpk) + f(xpi, ymj, zpk) + f(xmi, ymj, zpk)) +
                                (f(xpi, ypj, zmk) + f(xmi, ypj, zmk) + f(xpi, ymj, zmk) + f(xmi, ymj, zmk))
                            )
                        end
                        s_new += wiwj * inner_sum
                    end
                end
                for j in j_tile:j_end
                    wj = ws[j]
                    ypj, ymj = yp[j], ym[j]
                    zpj, zmj = zp[j], zm[j]
                    s_new += wi * wj * w0 * (
                        (f(xpi, ypj, z0) + f(xmi, ypj, z0) + f(xpi, ymj, z0) + f(xmi, ymj, z0)) +
                        (f(xpi, y0, zpj) + f(xmi, y0, zpj) + f(xpi, y0, zmj) + f(xmi, y0, zmj)) +
                        (f(x0, ypi, zpj) + f(x0, ymi, zpj) + f(x0, ypi, zmj) + f(x0, ymi, zmj))
                    )
                end
            end
            s_new += wi * w0sq * (
                (f(xpi, y0, z0) + f(xmi, y0, z0)) +
                (f(x0, ypi, z0) + f(x0, ymi, z0)) +
                (f(x0, y0, zpi) + f(x0, y0, zmi))
            )
        end

        for i in 2:2:n
            wi = ws[i]
            xpi, xmi = xp[i], xm[i]
            ypi, ymi = yp[i], ym[i]
            for jj in 1:odd_count
                j = (jj << 1) - 1
                wj = ws[j]
                ypj, ymj = yp[j], ym[j]
                zpj, zmj = zp[j], zm[j]
                wiwj = wi * wj
                for k_tile in 1:tile:n
                    k_end = min(k_tile + tile - 1, n)
                    inner_sum = zero(T)
                    @turbo warn_check_args=false for k in k_tile:k_end
                        zpk = zp[k]
                        zmk = zm[k]
                        inner_sum += ws[k] * (
                            (f(xpi, ypj, zpk) + f(xmi, ypj, zpk) + f(xpi, ymj, zpk) + f(xmi, ymj, zpk)) +
                            (f(xpi, ypj, zmk) + f(xmi, ypj, zmk) + f(xpi, ymj, zmk) + f(xmi, ymj, zmk))
                        )
                    end
                    s_new += wiwj * inner_sum
                end
                s_new += wi * wj * w0 * (
                    (f(xpi, ypj, z0) + f(xmi, ypj, z0) + f(xpi, ymj, z0) + f(xmi, ymj, z0)) +
                    (f(xpi, y0, zpj) + f(xmi, y0, zpj) + f(xpi, y0, zmj) + f(xmi, y0, zmj)) +
                    (f(x0, ypi, zpj) + f(x0, ymi, zpj) + f(x0, ypi, zmj) + f(x0, ymi, zmj))
                )
            end
        end

        for i in 2:2:n
            wi = ws[i]
            xpi, xmi = xp[i], xm[i]
            for jj in 1:(n >>> 1)
                j = jj << 1
                wj = ws[j]
                ypj, ymj = yp[j], ym[j]
                wiwj = wi * wj
                for kk in 1:tile:odd_count
                    k_tile_end = min(kk + tile - 1, odd_count)
                    inner_sum = zero(T)
                    @turbo warn_check_args=false for k_idx in kk:k_tile_end
                        k = (k_idx << 1) - 1
                        zpk = zp[k]
                        zmk = zm[k]
                        inner_sum += ws[k] * (
                            (f(xpi, ypj, zpk) + f(xmi, ypj, zpk) + f(xpi, ymj, zpk) + f(xmi, ymj, zpk)) +
                            (f(xpi, ypj, zmk) + f(xmi, ypj, zmk) + f(xpi, ymj, zmk) + f(xmi, ymj, zmk))
                        )
                    end
                    s_new += wiwj * inner_sum
                end
            end
        end
    end

    return s_new
end

@inline function run_package_style_levels(
    f::F, levels::Tuple{}, h::T, s_total::T, old_res::T, half::T,
    rtol::T, atol::T, Δx::T, x0::T, Δy::T, y0::T, Δz::T, z0::T, w0::T, w0sq::T,
) where {F,T}
    return old_res
end

@inline function run_package_style_levels(
    f::F, levels::Tuple{L,Vararg{Any}}, h::T, s_total::T, old_res::T, half::T,
    rtol::T, atol::T, Δx::T, x0::T, Δy::T, y0::T, Δz::T, z0::T, w0::T, w0sq::T,
) where {F,T,L}
    level = first(levels)
    h_next = h * half
    s_new = package_style_level_sum(f, level, Δx, x0, Δy, y0, Δz, z0, w0, w0sq)
    s_total_next = s_total + s_new
    new_res = Δx * Δy * Δz * (h_next * h_next * h_next) * s_total_next
    if abs(new_res - old_res) <= max(atol, rtol * abs(new_res))
        return new_res
    end
    return run_package_style_levels(
        f, Base.tail(levels), h_next, s_total_next, new_res, half,
        rtol, atol, Δx, x0, Δy, y0, Δz, z0, w0, w0sq,
    )
end

@inline function initial_octant_sum(f, x0, y0, z0, Δx, Δy, Δz, vi, vj, vk, wi, wj, wk)
    xp, xm = Δx * vi + x0, x0 - Δx * vi
    yp, ym = Δy * vj + y0, y0 - Δy * vj
    zp, zm = Δz * vk + z0, z0 - Δz * vk
    w = wi * wj * wk
    return w * (
        (f(xp, yp, zp) + f(xm, yp, zp) + f(xp, ym, zp) + f(xm, ym, zp)) +
        (f(xp, yp, zm) + f(xm, yp, zm) + f(xp, ym, zm) + f(xm, ym, zm))
    )
end

@inline function initial_plane_sum(f, x0, y0, z0, Δx, Δy, Δz, vi, vj, wi, wj, w0)
    xp, xm = Δx * vi + x0, x0 - Δx * vi
    yp_i, ym_i = Δy * vi + y0, y0 - Δy * vi
    yp_j, ym_j = Δy * vj + y0, y0 - Δy * vj
    zp_j, zm_j = Δz * vj + z0, z0 - Δz * vj
    w = wi * wj * w0
    return w * (
        (f(xp, yp_j, z0) + f(xm, yp_j, z0) + f(xp, ym_j, z0) + f(xm, ym_j, z0)) +
        (f(xp, y0, zp_j) + f(xm, y0, zp_j) + f(xp, y0, zm_j) + f(xm, y0, zm_j)) +
        (f(x0, yp_i, zp_j) + f(x0, ym_i, zp_j) + f(x0, yp_i, zm_j) + f(x0, ym_i, zm_j))
    )
end

@inline function initial_axis_sum(f, x0, y0, z0, Δx, Δy, Δz, vi, wi, w0sq)
    xp, xm = Δx * vi + x0, x0 - Δx * vi
    yp, ym = Δy * vi + y0, y0 - Δy * vi
    zp, zm = Δz * vi + z0, z0 - Δz * vi
    w = wi * w0sq
    return w * (
        (f(xp, y0, z0) + f(xm, y0, z0)) +
        (f(x0, yp, z0) + f(x0, ym, z0)) +
        (f(x0, y0, zp) + f(x0, y0, zm))
    )
end

function run_package_style_variant(
    f::F,
    low::SVector{3,T},
    up::SVector{3,T},
    variant::PackageStyleVariant{T,L};
    rtol::T=zero(T),
    atol::T=zero(T),
) where {F,T,L}
    Δx, x0 = FastTanhSinhQuadrature._midpoint_radius(low[1], up[1])
    Δy, y0 = FastTanhSinhQuadrature._midpoint_radius(low[2], up[2])
    Δz, z0 = FastTanhSinhQuadrature._midpoint_radius(low[3], up[3])
    half = T(0.5)
    h = variant.tm * half
    w0 = T(π) * half
    w0sq = w0 * w0
    w0cu = w0sq * w0

    s_total = w0cu * f(x0, y0, z0)
    @inbounds for i in 1:2
        vi, wi = variant.initial_x[i], variant.initial_w[i]
        for j in 1:2
            vj, wj = variant.initial_x[j], variant.initial_w[j]
            for k in 1:2
                vk, wk = variant.initial_x[k], variant.initial_w[k]
                s_total += initial_octant_sum(f, x0, y0, z0, Δx, Δy, Δz, vi, vj, vk, wi, wj, wk)
            end
            s_total += initial_plane_sum(f, x0, y0, z0, Δx, Δy, Δz, vi, vj, wi, wj, w0)
        end
        s_total += initial_axis_sum(f, x0, y0, z0, Δx, Δy, Δz, vi, wi, w0sq)
    end

    old_res = Δx * Δy * Δz * (h * h * h) * s_total
    return run_package_style_levels(
        f, variant.levels, h, s_total, old_res, half,
        rtol, atol, Δx, x0, Δy, y0, Δz, z0, w0, w0sq,
    )
end

function measure_variant(f, analytic, variant; samples::Int=20)
    low = SVector(-1.0, -1.0, -1.0)
    up = SVector(1.0, 1.0, 1.0)
    value = run_package_style_variant(f, low, up, variant)
    trial = @benchmark run_package_style_variant($f, $low, $up, $variant) samples=samples evals=1
    best = minimum(trial)
    return (
        name = variant.name,
        value = value,
        error = abs(value - analytic),
        time_ns = best.time,
        allocs = best.allocs,
        alloc_bytes = best.memory,
    )
end

function measure_package_scalar(f, analytic, max_levels::Int; samples::Int=20)
    low = SVector(-1.0, -1.0, -1.0)
    up = SVector(1.0, 1.0, 1.0)
    cache = adaptive_cache_3D(Float64; max_levels=max_levels)
    value = adaptive_integrate_3D(Float64, f, low, up; rtol=0.0, atol=0.0, max_levels=max_levels, warn=false, cache=cache)
    trial = @benchmark adaptive_integrate_3D(Float64, $f, $low, $up; rtol=0.0, atol=0.0, max_levels=$max_levels, warn=false, cache=$cache) samples=samples evals=1
    best = minimum(trial)
    return (
        name = "package_scalar",
        value = value,
        error = abs(value - analytic),
        time_ns = best.time,
        allocs = best.allocs,
        alloc_bytes = best.memory,
    )
end

function measure_package_avx(f, analytic, max_levels::Int; samples::Int=20)
    low = SVector(-1.0, -1.0, -1.0)
    up = SVector(1.0, 1.0, 1.0)
    cache = adaptive_cache_3D(Float64; max_levels=max_levels)
    value = adaptive_integrate_3D_avx(Float64, f, low, up; rtol=0.0, atol=0.0, max_levels=max_levels, warn=false, cache=cache)
    trial = @benchmark adaptive_integrate_3D_avx(Float64, $f, $low, $up; rtol=0.0, atol=0.0, max_levels=$max_levels, warn=false, cache=$cache) samples=samples evals=1
    best = minimum(trial)
    return (
        name = "package_avx",
        value = value,
        error = abs(value - analytic),
        time_ns = best.time,
        allocs = best.allocs,
        alloc_bytes = best.memory,
    )
end

function print_result(result; baseline_time=nothing, baseline_label="baseline")
    speedup = baseline_time === nothing ? 1.0 : baseline_time / result.time_ns
    println(
        rpad(result.name, 16),
        " time=", lpad(round(result.time_ns / 1e3; digits=3), 10), " us",
        "  vs_", baseline_label, "=", round(speedup; digits=3),
        "x",
        "  allocs=", result.allocs,
        "  bytes=", result.alloc_bytes,
        "  err=", result.error,
    )
end

f_exp(x, y, z) = exp(x + y + z)

function run_adaptive_benchmark(; max_levels_set=(3, 4, 5, 6), samples::Int=20)
    println("Adaptive 3D package-style comparison: Vector vs Ptr(Vector) vs StrideArray vs Ptr(StrideArray) vs hybrid SVector/Vector")
    println("Integrand: f(x,y,z) = exp(x+y+z)")
    println("Analytic integral on [-1,1]^3: (exp(1) - exp(-1))^3")

    analytic = (exp(1) - exp(-1))^3

    for max_levels in max_levels_set
        println("\n--- max_levels = $max_levels ---")
        vector_variant = build_variant(Float64, "vector", :vector, max_levels)
        ptr_vector_variant = build_variant(Float64, "ptr_vector", :ptr_vector, max_levels)
        stride_variant = build_variant(Float64, "stride", :stride, max_levels)
        ptr_stride_variant = build_variant(Float64, "ptr_stride", :ptr_stride, max_levels)
        hybrid64_variant = build_variant(Float64, "hybrid64", :hybrid_svector64, max_levels)
        hybrid128_variant = build_variant(Float64, "hybrid128", :hybrid_svector128, max_levels)
        package_scalar = measure_package_scalar(f_exp, analytic, max_levels; samples=samples)
        package_avx = measure_package_avx(f_exp, analytic, max_levels; samples=samples)
        vector_result = measure_variant(f_exp, analytic, vector_variant; samples=samples)
        ptr_vector_result = measure_variant(f_exp, analytic, ptr_vector_variant; samples=samples)
        stride_result = measure_variant(f_exp, analytic, stride_variant; samples=samples)
        ptr_stride_result = measure_variant(f_exp, analytic, ptr_stride_variant; samples=samples)
        hybrid64_result = measure_variant(f_exp, analytic, hybrid64_variant; samples=samples)
        hybrid128_result = measure_variant(f_exp, analytic, hybrid128_variant; samples=samples)
        baseline_time = package_avx.time_ns

        print_result(package_scalar)
        print_result(package_avx)
        print_result(vector_result; baseline_time=baseline_time, baseline_label="package_avx")
        print_result(ptr_vector_result; baseline_time=baseline_time, baseline_label="package_avx")
        print_result(stride_result; baseline_time=baseline_time, baseline_label="package_avx")
        print_result(ptr_stride_result; baseline_time=baseline_time, baseline_label="package_avx")
        print_result(hybrid64_result; baseline_time=baseline_time, baseline_label="package_avx")
        print_result(hybrid128_result; baseline_time=baseline_time, baseline_label="package_avx")
    end
end

run_adaptive_benchmark()
