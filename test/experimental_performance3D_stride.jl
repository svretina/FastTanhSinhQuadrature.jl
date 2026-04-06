using FastTanhSinhQuadrature
using StaticArrays
using BenchmarkTools
using LoopVectorization
using StrideArrays

struct Integration3DVectorCache{T}
    xp::Vector{T}
    xm::Vector{T}
    yp::Vector{T}
    ym::Vector{T}
    zp::Vector{T}
    zm::Vector{T}
end

function Integration3DVectorCache(::Type{T}, n::Int) where {T}
    return Integration3DVectorCache{T}(
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
    )
end

struct Integration3DStrideCache{A}
    xp::A
    xm::A
    yp::A
    ym::A
    zp::A
    zm::A
end

function Integration3DStrideCache(::Type{T}, n::Int) where {T}
    A = typeof(StrideArray{T}(undef, (n,)))
    return Integration3DStrideCache{A}(
        StrideArray{T}(undef, (n,)),
        StrideArray{T}(undef, (n,)),
        StrideArray{T}(undef, (n,)),
        StrideArray{T}(undef, (n,)),
        StrideArray{T}(undef, (n,)),
        StrideArray{T}(undef, (n,)),
    )
end

struct Integration3DPtrCache{A}
    xp::A
    xm::A
    yp::A
    ym::A
    zp::A
    zm::A
end

Integration3DPtrCache(cache) =
    Integration3DPtrCache(
        PtrArray(cache.xp),
        PtrArray(cache.xm),
        PtrArray(cache.yp),
        PtrArray(cache.ym),
        PtrArray(cache.zp),
        PtrArray(cache.zm),
    )

struct KernelVariant{X,W,C}
    name::String
    x::X
    w::W
    cache::C
end

@inline to_svector(::Val{N}, v::AbstractVector{T}) where {N,T} =
    SVector{N,T}(ntuple(i -> v[i], Val(N)))

function _integrate3D_core(
    f::F, low::SVector{3,T}, up::SVector{3,T}, h::T,
    x::AbstractVector{T}, w::AbstractVector{T},
    xp::AbstractVector{T}, xm::AbstractVector{T},
    yp::AbstractVector{T}, ym::AbstractVector{T},
    zp::AbstractVector{T}, zm::AbstractVector{T},
) where {T<:Real,F}
    Δx, x₀ = (up[1] - low[1]) / 2, (up[1] + low[1]) / 2
    Δy, y₀ = (up[2] - low[2]) / 2, (up[2] + low[2]) / 2
    Δz, z₀ = (up[3] - low[3]) / 2, (up[3] + low[3]) / 2
    n = length(x)

    @turbo for i in 1:n
        xi = x[i]
        vx = Δx * xi
        vy = Δy * xi
        vz = Δz * xi
        xp[i] = x₀ + vx
        xm[i] = x₀ - vx
        yp[i] = y₀ + vy
        ym[i] = y₀ - vy
        zp[i] = z₀ + vz
        zm[i] = z₀ - vz
    end

    w₀ = T(π) / 2
    w₀² = w₀ * w₀
    w₀³ = w₀² * w₀

    total_sum = w₀³ * f(x₀, y₀, z₀)

    axis_sum = zero(T)
    @turbo for i in 1:n
        wi = w[i]
        axis_sum += wi * w₀² * (
            f(xp[i], y₀, z₀) + f(xm[i], y₀, z₀) +
            f(x₀, yp[i], z₀) + f(x₀, ym[i], z₀) +
            f(x₀, y₀, zp[i]) + f(x₀, y₀, zm[i])
        )
    end

    plane_sum = zero(T)
    @turbo for i in 1:n, j in 1:n
        wiwj = w[i] * w[j]
        plane_sum += wiwj * w₀ * (
            f(xp[i], yp[j], z₀) + f(xm[i], yp[j], z₀) +
            f(xp[i], ym[j], z₀) + f(xm[i], ym[j], z₀) +
            f(xp[i], y₀, zp[j]) + f(xm[i], y₀, zp[j]) +
            f(xp[i], y₀, zm[j]) + f(xm[i], y₀, zm[j]) +
            f(x₀, yp[i], zp[j]) + f(x₀, ym[i], zp[j]) +
            f(x₀, yp[i], zm[j]) + f(x₀, ym[i], zm[j])
        )
    end

    octant_sum = zero(T)
    tile_size = 32
    @inbounds for i_tile in 1:tile_size:n
        i_end = min(i_tile + tile_size - 1, n)
        for j_tile in 1:tile_size:n
            j_end = min(j_tile + tile_size - 1, n)
            for k_tile in 1:tile_size:n
                k_end = min(k_tile + tile_size - 1, n)
                for i in i_tile:i_end
                    wi = w[i]
                    xpi, xmi = xp[i], xm[i]
                    for j in j_tile:j_end
                        wj = w[j]
                        ypj, ymj = yp[j], ym[j]
                        wiwj = wi * wj
                        inner_sum = zero(T)
                        @turbo for k in k_tile:k_end
                            wk = w[k]
                            zpk, zmk = zp[k], zm[k]
                            inner_sum += wk * (
                                (f(xpi, ypj, zpk) + f(xmi, ypj, zpk) +
                                 f(xpi, ymj, zpk) + f(xmi, ymj, zpk)) +
                                (f(xpi, ypj, zmk) + f(xmi, ypj, zmk) +
                                 f(xpi, ymj, zmk) + f(xmi, ymj, zmk))
                            )
                        end
                        octant_sum += wiwj * inner_sum
                    end
                end
            end
        end
    end

    return (Δx * Δy * Δz) * (h * h * h) * (total_sum + axis_sum + plane_sum + octant_sum)
end

@inline function run_kernel(
    f, low, up, h, variant::KernelVariant,
)
    cache = variant.cache
    return _integrate3D_core(
        f, low, up, h, variant.x, variant.w,
        cache.xp, cache.xm, cache.yp, cache.ym, cache.zp, cache.zm,
    )
end

function measure_variant(f, low, up, h, analytic, variant::KernelVariant; samples::Int=20)
    value = run_kernel(f, low, up, h, variant)
    trial = @benchmark run_kernel($f, $low, $up, $h, $variant) samples=samples evals=1
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

function measure_package_baseline(f, low, up, x, w, h, analytic; samples::Int=20)
    value = integrate3D_avx(f, low, up, x, w, h)
    trial = @benchmark integrate3D_avx($f, $low, $up, $x, $w, $h) samples=samples evals=1
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

function print_result(result; baseline_time=nothing)
    speedup = baseline_time === nothing ? 1.0 : baseline_time / result.time_ns
    println(
        rpad(result.name, 14),
        " time=", lpad(round(result.time_ns / 1e3; digits=3), 9), " us",
        "  speedup=", round(speedup; digits=3),
        "x",
        "  allocs=", result.allocs,
        "  bytes=", result.alloc_bytes,
        "  err=", result.error,
    )
end

function run_experimental_benchmark(; Ns=(32, 64, 128, 256), samples::Int=20, static_max_n::Int=64)
    println("3D kernel comparison: Vector vs Ptr(Vector) vs StrideArray vs Ptr(StrideArray) vs SVector")
    println("Integrand: f(x,y,z) = exp(x+y+z)")
    println("Analytic integral on [-1,1]^3: (exp(1) - exp(-1))^3")

    f(x, y, z) = exp(x + y + z)
    analytic = (exp(1) - exp(-1))^3
    low = SVector(-1.0, -1.0, -1.0)
    up = SVector(1.0, 1.0, 1.0)

    for N in Ns
        println("\n--- N = $N ---")
        x_vec, w_vec, h = tanhsinh(Float64, N)
        node_count = length(x_vec)

        x_stride = StrideArray(x_vec)
        w_stride = StrideArray(w_vec)
        vector_cache = Integration3DVectorCache(Float64, node_count)
        stride_cache = Integration3DStrideCache(Float64, node_count)

        vector_variant = KernelVariant(
            "vector",
            x_vec,
            w_vec,
            vector_cache,
        )
        ptr_vector_variant = KernelVariant(
            "ptr_vector",
            PtrArray(x_vec),
            PtrArray(w_vec),
            Integration3DPtrCache(Integration3DVectorCache(Float64, node_count)),
        )
        stride_variant = KernelVariant(
            "stride",
            x_stride,
            w_stride,
            stride_cache,
        )
        ptr_variant = KernelVariant(
            "ptr_stride",
            PtrArray(x_stride),
            PtrArray(w_stride),
            Integration3DPtrCache(Integration3DStrideCache(Float64, node_count)),
        )

        package_result = measure_package_baseline(f, low, up, x_vec, w_vec, h, analytic; samples=samples)
        vector_result = measure_variant(f, low, up, h, analytic, vector_variant; samples=samples)
        ptr_vector_result = measure_variant(f, low, up, h, analytic, ptr_vector_variant; samples=samples)
        stride_result = measure_variant(f, low, up, h, analytic, stride_variant; samples=samples)
        ptr_result = measure_variant(f, low, up, h, analytic, ptr_variant; samples=samples)

        print_result(package_result)
        print_result(vector_result; baseline_time=vector_result.time_ns)
        print_result(ptr_vector_result; baseline_time=vector_result.time_ns)
        print_result(stride_result; baseline_time=vector_result.time_ns)
        print_result(ptr_result; baseline_time=vector_result.time_ns)

        if node_count <= static_max_n
            x_static = to_svector(Val(node_count), x_vec)
            w_static = to_svector(Val(node_count), w_vec)
            svector_variant = KernelVariant(
                "svector",
                x_static,
                w_static,
                Integration3DVectorCache(Float64, node_count),
            )
            svector_result = measure_variant(f, low, up, h, analytic, svector_variant; samples=samples)
            print_result(svector_result; baseline_time=vector_result.time_ns)
        else
            println(rpad("svector", 14), " skipped for node_count > ", static_max_n, " to avoid runaway compile cost")
        end

        # PtrArray(SVector) is intentionally omitted. In quick probes it did not expose
        # stable numeric data, which suggests it is not a safe ownership model here.
    end
end

run_experimental_benchmark()
