using FastTanhSinhQuadrature
using StaticArrays
using BenchmarkTools
using LoopVectorization

# Cache for Optimization 2: Pre-calculation of coordinate pairs
struct Integration3DCache{T}
    xp::Vector{T}
    xm::Vector{T}
    yp::Vector{T}
    ym::Vector{T}
    zp::Vector{T}
    zm::Vector{T}
end

function Integration3DCache(T::Type, n::Int)
    return Integration3DCache{T}(
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n),
        Vector{T}(undef, n)
    )
end

# Optimization 3: Handle 8 reflections in a single manual SIMD-like call 
# by passing SVectors to the function if the function supports it.
# This exploits SVML (SIMD Vector Math Library) if LoopVectorization can 
# effectively unroll the SVector calls.
@inline function eval_8_octants(f::F, xp, xm, yp, ym, zp, zm) where {F}
    # These are 8 independent calls that @turbo should ideally group 
    # to use the widest possible SIMD registers.
    return (
        (f(xp, yp, zp) + f(xm, yp, zp) + f(xp, ym, zp) + f(xm, ym, zp)) +
        (f(xp, yp, zm) + f(xm, yp, zm) + f(xp, ym, zm) + f(xm, ym, zm))
    )
end

function integrate3D_experimental(f::F, low, up, x, w, h, cache) where {F}
    T = eltype(x)
    Δx, x₀ = (T(up[1]) - T(low[1])) / 2, (T(up[1]) + T(low[1])) / 2
    Δy, y₀ = (T(up[2]) - T(low[2])) / 2, (T(up[2]) + T(low[2])) / 2
    Δz, z₀ = (T(up[3]) - T(low[3])) / 2, (T(up[3]) + T(low[3])) / 2
    
    n = length(x)
    
    # Optimization 2: Pre-calculation (Heavily used here)
    @turbo for i in 1:n
        val_x = Δx * x[i]
        cache.xp[i] = x₀ + val_x
        cache.xm[i] = x₀ - val_x
        val_y = Δy * x[i]
        cache.yp[i] = y₀ + val_y
        cache.ym[i] = y₀ - val_y
        val_z = Δz * x[i]
        cache.zp[i] = z₀ + val_z
        cache.zm[i] = z₀ - val_z
    end
    
    w₀ = T(π) / 2
    w₀² = w₀ * w₀
    w₀³ = w₀² * w₀
    
    total_sum = w₀³ * f(x₀, y₀, z₀)
    
    axis_sum = zero(T)
    @turbo for i in 1:n
       axis_sum += w[i] * w₀² * (
           f(cache.xp[i], y₀, z₀) + f(cache.xm[i], y₀, z₀) +
           f(x₀, cache.yp[i], z₀) + f(x₀, y₀ - (cache.yp[i]-y₀), z₀) + # Using ym alternative
           f(x₀, y₀, cache.zp[i]) + f(x₀, y₀, cache.zm[i])
       )
    end
    
    plane_sum = zero(T)
    @turbo for i in 1:n, j in 1:n
        plane_sum += w[i] * w[j] * w₀ * (
            f(cache.xp[i], cache.yp[j], z₀) + f(cache.xm[i], cache.yp[j], z₀) + 
            f(cache.xp[i], cache.ym[j], z₀) + f(cache.xm[i], cache.ym[j], z₀) +
            f(cache.xp[i], y₀, cache.zp[j]) + f(cache.xm[i], y₀, cache.zp[j]) + 
            f(cache.xp[i], y₀, cache.zm[j]) + f(cache.xm[i], y₀, cache.zm[j]) +
            f(x₀, cache.yp[i], cache.zp[j]) + f(x₀, cache.ym[i], cache.zp[j]) + 
            f(x₀, cache.yp[i], cache.zm[j]) + f(x₀, cache.ym[i], cache.zm[j])
        )
    end
    
    # Optimization 1 & 3: Tiling + Manual-style Vectorization setup
    octant_sum = zero(T)
    
    # Using a tile size of 16 for better register management in 3D
    tile_size = 16 
    for i_tile in 1:tile_size:n
        i_end = min(i_tile + tile_size - 1, n)
        for j_tile in 1:tile_size:n
            j_end = min(j_tile + tile_size - 1, n)
            for k_tile in 1:tile_size:n
                k_end = min(k_tile + tile_size - 1, n)
                
                for i in i_tile:i_end
                    wi = w[i]
                    xp_i, xm_i = cache.xp[i], cache.xm[i]
                    for j in j_tile:j_end
                        wj = w[j]
                        yp_j, ym_j = cache.yp[j], cache.ym[j]
                        wiwj = wi * wj
                        
                        # Optimization 3: This loop is the target for maximal SIMD throughput
                        inner_s = zero(T)
                        @turbo for k in k_tile:k_end
                            wk = w[k]
                            zp_k, zm_k = cache.zp[k], cache.zm[k]
                            
                            # Evaluating all 8 octants in one tight vectorizable block
                            inner_s += wk * (
                                (f(xp_i, yp_j, zp_k) + f(xm_i, yp_j, zp_k) + f(xp_i, ym_j, zp_k) + f(xm_i, ym_j, zp_k)) +
                                (f(xp_i, yp_j, zm_k) + f(xm_i, yp_j, zm_k) + f(xp_i, ym_j, zm_k) + f(xm_i, ym_j, zm_k))
                            )
                        end
                        octant_sum += wiwj * inner_s
                    end
                end
            end
        end
    end
    
    return (Δx * Δy * Δz) * (h * h * h) * (total_sum + axis_sum + plane_sum + octant_sum)
end

function run_experimental_benchmark()
    f(x, y, z) = exp(x + y + z)
    analytic = (exp(1) - exp(-1))^3
    low = SVector(-1.0, -1.0, -1.0)
    up = SVector(1.0, 1.0, 1.0)
    
    for N in [64, 128]
        println("\n--- Testing N = $N ---")
        x, w, h = tanhsinh(Float64, N)
        cache = Integration3DCache(Float64, N)
        
        val_std = integrate3D_avx(f, low, up, x, w, h)
        val_exp = integrate3D_experimental(f, low, up, x, w, h, cache)
        
        println("Analytic Error (Standard):     ", abs(val_std - analytic))
        println("Analytic Error (Experimental): ", abs(val_exp - analytic))
        
        println("Timing Standard:")
        @btime integrate3D_avx($f, $low, $up, $x, $w, $h)
        
        println("Timing Experimental:")
        @btime integrate3D_experimental($f, $low, $up, $x, $w, $h, $cache)
    end
end

run_experimental_benchmark()