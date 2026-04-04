using BenchmarkTools
using Cuba
using Cubature
using FastGaussQuadrature
using FastTanhSinhQuadrature
using HCubature
using Printf
using QuadGK
using StaticArrays

const RTOL = 1e-6
const ATOL = 1e-8
const MAXEVALS = 200_000
const CUBA_MINEVALS = 2_000
const BENCH_SAMPLES = 3
const BENCH_EVALS = 1

const ADAPTIVE_MAX_LEVELS_1D = 12
const ADAPTIVE_MAX_LEVELS_2D = 9
const ADAPTIVE_MAX_LEVELS_3D = 7
const CANDIDATES_GQ_1D = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]

struct Problem{F,B}
    dim::Int
    name::String
    f::F
    exact::Float64
    low::B
    up::B
end

const PROBLEMS = [
    Problem(1, "x^6 - 2x^3 + 0.5", x -> x^6 - 2x^3 + 0.5, 9 / 7, -1.0, 1.0),
    Problem(1, "1/(1+25x^2)", x -> 1 / (1 + 25x^2), 2 * atan(5.0) / 5, -1.0, 1.0),
    Problem(1, "log(1-x)", x -> log(1 - x), 2 * log(2.0) - 2, -1.0, 1.0),
    Problem(1, "1/sqrt(1-x^2)", x -> 1 / sqrt(1 - x^2), Float64(pi), -1.0, 1.0),
    Problem(1, "sin^2(1000x)", x -> sin(1000 * x)^2, Float64(pi), -Float64(pi), Float64(pi)),
    Problem(2, "exp(x+y)", (x, y) -> exp(x + y), (exp(1.0) - exp(-1.0))^2, SVector(-1.0, -1.0), SVector(1.0, 1.0)),
    Problem(2, "x^2 + y^2", (x, y) -> x^2 + y^2, 8 / 3, SVector(-1.0, -1.0), SVector(1.0, 1.0)),
    Problem(2, "1/sqrt((1-x^2)(1-y^2))",
        (x, y) -> 1 / sqrt((1 - x^2) * (1 - y^2)), Float64(pi^2), SVector(-1.0, -1.0), SVector(1.0, 1.0)),
    Problem(3, "exp(x+y+z)", (x, y, z) -> exp(x + y + z), (exp(1.0) - exp(-1.0))^3, SVector(-1.0, -1.0, -1.0), SVector(1.0, 1.0, 1.0)),
    Problem(3, "x^2*y^2*z^2", (x, y, z) -> x^2 * y^2 * z^2, 8 / 27, SVector(-1.0, -1.0, -1.0), SVector(1.0, 1.0, 1.0)),
    Problem(3, "1/sqrt((1-x^2)(1-y^2)(1-z^2))",
        (x, y, z) -> 1 / sqrt((1 - x^2) * (1 - y^2) * (1 - z^2)), Float64(pi^3), SVector(-1.0, -1.0, -1.0), SVector(1.0, 1.0, 1.0)),
]

@inline dim_label(dim::Int) = dim == 1 ? "1D" : (dim == 2 ? "2D" : "3D")
@inline function tol_target(exact::Float64)
    return max(ATOL, RTOL * abs(exact))
end
@inline function rel_error(val::Float64, exact::Float64)
    denom = abs(exact)
    return denom == 0 ? abs(val) : abs(val - exact) / denom
end

function integrate_gauss_precomputed(f::Function, x::AbstractVector{T},
    w::AbstractVector{T}) where {T<:Real}
    s = zero(T)
    @inbounds for i in eachindex(x)
        s += w[i] * f(x[i])
    end
    return s
end

function integrate_gauss_precomputed(f::Function, low::T, up::T,
    x::AbstractVector{T}, w::AbstractVector{T}) where {T<:Real}
    Δx, x₀ = (up - low) / 2, (up + low) / 2
    s = zero(T)
    @inbounds for i in eachindex(x)
        s += w[i] * f(x₀ + Δx * x[i])
    end
    return Δx * s
end

@inline bounds_vec(prob::Problem{F,SVector{N,T}}) where {F,N,T} =
    (collect(prob.low), collect(prob.up))
@inline bounds_vec(prob::Problem{F,Float64}) where {F} = ([prob.low], [prob.up])

@inline adaptive_max_levels(dim::Int) = dim == 1 ? ADAPTIVE_MAX_LEVELS_1D :
                                         dim == 2 ? ADAPTIVE_MAX_LEVELS_2D :
                                         ADAPTIVE_MAX_LEVELS_3D

function adaptive_cache_for_problem(prob::Problem)
    max_levels = adaptive_max_levels(prob.dim)
    if prob.dim == 1
        return adaptive_cache_1D(Float64; max_levels=max_levels)
    elseif prob.dim == 2
        return adaptive_cache_2D(Float64; max_levels=max_levels)
    else
        return adaptive_cache_3D(Float64; max_levels=max_levels)
    end
end

function adaptive_result_at_level(prob::Problem, level::Int, cache)
    level >= 0 || throw(ArgumentError("level must be nonnegative"))
    if prob.dim == 1
        return adaptive_integrate_1D(
            Float64, prob.f, prob.low, prob.up; rtol=0.0, atol=0.0, max_levels=level, warn=false, cache=cache
        )
    elseif prob.dim == 2
        return adaptive_integrate_2D(
            Float64, prob.f, prob.low, prob.up; rtol=0.0, atol=0.0, max_levels=level, warn=false, cache=cache
        )
    else
        return adaptive_integrate_3D(
            Float64, prob.f, prob.low, prob.up; rtol=0.0, atol=0.0, max_levels=level, warn=false, cache=cache
        )
    end
end

function adaptive_converged_grid(prob::Problem; cache)
    max_levels = adaptive_max_levels(prob.dim)
    prev = adaptive_result_at_level(prob, 0, cache)
    for level in 1:max_levels
        curr = adaptive_result_at_level(prob, level, cache)
        target = max(ATOL, RTOL * abs(curr))
        if abs(curr - prev) <= target
            N = 1 << (level + 2)  # N=2n with n=2^(level+1)
            return (level=level, N=N, value=Float64(curr), converged=true)
        end
        prev = curr
    end
    N = 1 << (max_levels + 2)
    return (level=max_levels, N=N, value=Float64(prev), converged=false)
end

function avx_config_at_N(prob::Problem, N::Int)
    x, w, h = tanhsinh(Float64, N)
    val = if prob.dim == 1
        integrate1D_avx(prob.f, prob.low, prob.up, x, w, h)
    elseif prob.dim == 2
        integrate2D_avx(prob.f, prob.low, prob.up, x, w, h)
    else
        integrate3D_avx(prob.f, prob.low, prob.up, x, w, h)
    end
    return (N=N, x=x, w=w, h=h, val=Float64(val))
end

function calibrate_gq_1d(prob::Problem)
    best = nothing
    for N in CANDIDATES_GQ_1D
        x, w = gausslegendre(N)
        val = integrate_gauss_precomputed(prob.f, prob.low, prob.up, x, w)
        best = (N=N, x=x, w=w, val=Float64(val))
        if isfinite(best.val) && abs(best.val - prob.exact) <= tol_target(prob.exact)
            return best, true
        end
    end
    return best, false
end

function evaluate_method(dim::Int, fname::String, method::String, exact::Float64,
    runner::Function; notes::AbstractString="")
    try
        val = Float64(runner())
        # Warm call done above; benchmark runtime only.
        time_s = @belapsed $runner() samples=BENCH_SAMPLES evals=BENCH_EVALS
        abs_err = abs(val - exact)
        rel_err = rel_error(val, exact)
        status = (isfinite(val) && abs_err <= tol_target(exact)) ? "ok" : "above_tol"
        return (dim=dim_label(dim), function_name=fname, method=method,
            time_ns=time_s * 1e9, value=val, exact=exact,
            abs_error=abs_err, rel_error=rel_err, status=status, notes=String(notes))
    catch err
        return (dim=dim_label(dim), function_name=fname, method=method,
            time_ns=NaN, value=NaN, exact=exact, abs_error=Inf, rel_error=Inf,
            status="error", notes="$(typeof(err))")
    end
end

function benchmark_problem(prob::Problem)
    rows = NamedTuple[]

    if prob.dim == 1
        f = prob.f
        low, up = prob.low, prob.up
        max_levels = adaptive_max_levels(1)
        adapt_cache = adaptive_cache_for_problem(prob)

        run_adapt = () -> adaptive_integrate_1D(
            Float64, f, low, up; rtol=RTOL, atol=ATOL, max_levels=max_levels, cache=adapt_cache
        )
        push!(rows, evaluate_method(1, prob.name, "FTS adaptive (typed)", prob.exact, run_adapt))

        run_quad = () -> quad(f, low, up; rtol=RTOL, atol=ATOL, max_levels=max_levels, cache=adapt_cache)
        push!(rows, evaluate_method(1, prob.name, "FTS quad (convenience)", prob.exact, run_quad))

        adapt_grid = adaptive_converged_grid(prob; cache=adapt_cache)
        avx_cfg = avx_config_at_N(prob, adapt_grid.N)
        run_avx = () -> integrate1D_avx(f, low, up, avx_cfg.x, avx_cfg.w, avx_cfg.h)
        note_avx = adapt_grid.converged ?
                   "N=$(avx_cfg.N), adaptive_level=$(adapt_grid.level)" :
                   "N=$(avx_cfg.N), adaptive_level=$(adapt_grid.level), max_levels_reached"
        push!(rows, evaluate_method(1, prob.name, "FTS avx (precomputed)", prob.exact, run_avx; notes=note_avx))

        run_qgk = () -> first(quadgk(f, low, up; rtol=RTOL, atol=ATOL, maxevals=MAXEVALS))
        push!(rows, evaluate_method(1, prob.name, "QuadGK", prob.exact, run_qgk))

        low_vec, up_vec = bounds_vec(prob)
        fv_hc = v -> f(v[1])
        run_hcub = () -> first(HCubature.hcubature(fv_hc, low_vec, up_vec;
            rtol=RTOL, atol=ATOL, maxevals=MAXEVALS))
        push!(rows, evaluate_method(1, prob.name, "HCubature", prob.exact, run_hcub))

        fv_hcub = v -> f(v[1])
        run_hcub_jl = () -> first(Cubature.hcubature(fv_hcub, low_vec, up_vec;
            reltol=RTOL, abstol=ATOL, maxevals=MAXEVALS))
        push!(rows, evaluate_method(1, prob.name, "CubatureJLh", prob.exact, run_hcub_jl))

        fv_pcub = v -> f(v[1])
        run_pcub_jl = () -> first(Cubature.pcubature(fv_pcub, low_vec, up_vec;
            reltol=RTOL, abstol=ATOL, maxevals=MAXEVALS))
        push!(rows, evaluate_method(1, prob.name, "CubatureJLp", prob.exact, run_pcub_jl))

        gq_cfg, gq_ok = calibrate_gq_1d(prob)
        run_gq = () -> integrate_gauss_precomputed(f, low, up, gq_cfg.x, gq_cfg.w)
        note_gq = gq_ok ? "N=$(gq_cfg.N), independently calibrated for tolerance" :
                  "N=$(gq_cfg.N), tol_not_met"
        push!(rows, evaluate_method(1, prob.name, "FastGauss (precomputed)", prob.exact, run_gq; notes=note_gq))

    elseif prob.dim == 2
        f = prob.f
        low, up = prob.low, prob.up
        max_levels = adaptive_max_levels(2)
        adapt_cache = adaptive_cache_for_problem(prob)

        run_adapt = () -> adaptive_integrate_2D(
            Float64, f, low, up; rtol=RTOL, atol=ATOL, max_levels=max_levels, cache=adapt_cache
        )
        push!(rows, evaluate_method(2, prob.name, "FTS adaptive (typed)", prob.exact, run_adapt))

        run_quad = () -> quad(f, low, up; rtol=RTOL, atol=ATOL, max_levels=max_levels, cache=adapt_cache)
        push!(rows, evaluate_method(2, prob.name, "FTS quad (convenience)", prob.exact, run_quad))

        adapt_grid = adaptive_converged_grid(prob; cache=adapt_cache)
        avx_cfg = avx_config_at_N(prob, adapt_grid.N)
        run_avx = () -> integrate2D_avx(f, low, up, avx_cfg.x, avx_cfg.w, avx_cfg.h)
        note_avx = adapt_grid.converged ?
                   "N=$(avx_cfg.N), adaptive_level=$(adapt_grid.level)" :
                   "N=$(avx_cfg.N), adaptive_level=$(adapt_grid.level), max_levels_reached"
        push!(rows, evaluate_method(2, prob.name, "FTS avx (precomputed)", prob.exact, run_avx; notes=note_avx))

        low_vec, up_vec = bounds_vec(prob)
        fv_hc = v -> f(v[1], v[2])
        run_hcub = () -> first(HCubature.hcubature(fv_hc, low_vec, up_vec;
            rtol=RTOL, atol=ATOL, maxevals=MAXEVALS))
        push!(rows, evaluate_method(2, prob.name, "HCubature", prob.exact, run_hcub))

        fv_hcub = v -> f(v[1], v[2])
        run_hcub_jl = () -> first(Cubature.hcubature(fv_hcub, low_vec, up_vec;
            reltol=RTOL, abstol=ATOL, maxevals=MAXEVALS))
        push!(rows, evaluate_method(2, prob.name, "CubatureJLh", prob.exact, run_hcub_jl))

        fv_pcub = v -> f(v[1], v[2])
        run_pcub_jl = () -> first(Cubature.pcubature(fv_pcub, low_vec, up_vec;
            reltol=RTOL, abstol=ATOL, maxevals=MAXEVALS))
        push!(rows, evaluate_method(2, prob.name, "CubatureJLp", prob.exact, run_pcub_jl))

        run_cuba_vegas = () -> begin
            jac = (up[1] - low[1]) * (up[2] - low[2])
            function cuba_f(u, out)
                x = low[1] + (up[1] - low[1]) * u[1]
                y = low[2] + (up[2] - low[2]) * u[2]
                out[1] = jac * f(x, y)
            end
            out = Cuba.vegas(cuba_f, 2, 1; rtol=RTOL, atol=ATOL,
                maxevals=MAXEVALS, minevals=CUBA_MINEVALS, flags=0)
            out.integral[1]
        end
        push!(rows, evaluate_method(2, prob.name, "CubaVegas", prob.exact, run_cuba_vegas))

        run_cuba_divonne = () -> begin
            jac = (up[1] - low[1]) * (up[2] - low[2])
            function cuba_f(u, out)
                x = low[1] + (up[1] - low[1]) * u[1]
                y = low[2] + (up[2] - low[2]) * u[2]
                out[1] = jac * f(x, y)
            end
            out = Cuba.divonne(cuba_f, 2, 1; rtol=RTOL, atol=ATOL,
                maxevals=MAXEVALS, minevals=CUBA_MINEVALS, flags=0)
            out.integral[1]
        end
        push!(rows, evaluate_method(2, prob.name, "CubaDivonne", prob.exact, run_cuba_divonne))

        run_cuba_cuhre = () -> begin
            jac = (up[1] - low[1]) * (up[2] - low[2])
            function cuba_f(u, out)
                x = low[1] + (up[1] - low[1]) * u[1]
                y = low[2] + (up[2] - low[2]) * u[2]
                out[1] = jac * f(x, y)
            end
            out = Cuba.cuhre(cuba_f, 2, 1; rtol=RTOL, atol=ATOL,
                maxevals=MAXEVALS, minevals=0, flags=0)
            out.integral[1]
        end
        push!(rows, evaluate_method(2, prob.name, "CubaCuhre", prob.exact, run_cuba_cuhre))

    else
        f = prob.f
        low, up = prob.low, prob.up
        max_levels = adaptive_max_levels(3)
        adapt_cache = adaptive_cache_for_problem(prob)

        run_adapt = () -> adaptive_integrate_3D(
            Float64, f, low, up; rtol=RTOL, atol=ATOL, max_levels=max_levels, cache=adapt_cache
        )
        push!(rows, evaluate_method(3, prob.name, "FTS adaptive (typed)", prob.exact, run_adapt))

        run_quad = () -> quad(f, low, up; rtol=RTOL, atol=ATOL, max_levels=max_levels, cache=adapt_cache)
        push!(rows, evaluate_method(3, prob.name, "FTS quad (convenience)", prob.exact, run_quad))

        adapt_grid = adaptive_converged_grid(prob; cache=adapt_cache)
        avx_cfg = avx_config_at_N(prob, adapt_grid.N)
        run_avx = () -> integrate3D_avx(f, low, up, avx_cfg.x, avx_cfg.w, avx_cfg.h)
        note_avx = adapt_grid.converged ?
                   "N=$(avx_cfg.N), adaptive_level=$(adapt_grid.level)" :
                   "N=$(avx_cfg.N), adaptive_level=$(adapt_grid.level), max_levels_reached"
        push!(rows, evaluate_method(3, prob.name, "FTS avx (precomputed)", prob.exact, run_avx; notes=note_avx))

        low_vec, up_vec = bounds_vec(prob)
        fv_hc = v -> f(v[1], v[2], v[3])
        run_hcub = () -> first(HCubature.hcubature(fv_hc, low_vec, up_vec;
            rtol=RTOL, atol=ATOL, maxevals=MAXEVALS))
        push!(rows, evaluate_method(3, prob.name, "HCubature", prob.exact, run_hcub))

        fv_hcub = v -> f(v[1], v[2], v[3])
        run_hcub_jl = () -> first(Cubature.hcubature(fv_hcub, low_vec, up_vec;
            reltol=RTOL, abstol=ATOL, maxevals=MAXEVALS))
        push!(rows, evaluate_method(3, prob.name, "CubatureJLh", prob.exact, run_hcub_jl))

        fv_pcub = v -> f(v[1], v[2], v[3])
        run_pcub_jl = () -> first(Cubature.pcubature(fv_pcub, low_vec, up_vec;
            reltol=RTOL, abstol=ATOL, maxevals=MAXEVALS))
        push!(rows, evaluate_method(3, prob.name, "CubatureJLp", prob.exact, run_pcub_jl))

        run_cuba_vegas = () -> begin
            jac = (up[1] - low[1]) * (up[2] - low[2]) * (up[3] - low[3])
            function cuba_f(u, out)
                x = low[1] + (up[1] - low[1]) * u[1]
                y = low[2] + (up[2] - low[2]) * u[2]
                z = low[3] + (up[3] - low[3]) * u[3]
                out[1] = jac * f(x, y, z)
            end
            out = Cuba.vegas(cuba_f, 3, 1; rtol=RTOL, atol=ATOL,
                maxevals=MAXEVALS, minevals=CUBA_MINEVALS, flags=0)
            out.integral[1]
        end
        push!(rows, evaluate_method(3, prob.name, "CubaVegas", prob.exact, run_cuba_vegas))

        run_cuba_divonne = () -> begin
            jac = (up[1] - low[1]) * (up[2] - low[2]) * (up[3] - low[3])
            function cuba_f(u, out)
                x = low[1] + (up[1] - low[1]) * u[1]
                y = low[2] + (up[2] - low[2]) * u[2]
                z = low[3] + (up[3] - low[3]) * u[3]
                out[1] = jac * f(x, y, z)
            end
            out = Cuba.divonne(cuba_f, 3, 1; rtol=RTOL, atol=ATOL,
                maxevals=MAXEVALS, minevals=CUBA_MINEVALS, flags=0)
            out.integral[1]
        end
        push!(rows, evaluate_method(3, prob.name, "CubaDivonne", prob.exact, run_cuba_divonne))

        run_cuba_cuhre = () -> begin
            jac = (up[1] - low[1]) * (up[2] - low[2]) * (up[3] - low[3])
            function cuba_f(u, out)
                x = low[1] + (up[1] - low[1]) * u[1]
                y = low[2] + (up[2] - low[2]) * u[2]
                z = low[3] + (up[3] - low[3]) * u[3]
                out[1] = jac * f(x, y, z)
            end
            out = Cuba.cuhre(cuba_f, 3, 1; rtol=RTOL, atol=ATOL,
                maxevals=MAXEVALS, minevals=0, flags=0)
            out.integral[1]
        end
        push!(rows, evaluate_method(3, prob.name, "CubaCuhre", prob.exact, run_cuba_cuhre))
    end

    return rows
end

function format_float(x::Float64)
    if isnan(x)
        return "n/a"
    elseif isinf(x)
        return "Inf"
    elseif abs(x) >= 1e4 || (abs(x) > 0 && abs(x) < 1e-4)
        return @sprintf("%.3e", x)
    else
        return @sprintf("%.4f", x)
    end
end

function write_csv(path::AbstractString, rows)
    open(path, "w") do io
        println(io, "dim,function,method,time_ns,value,exact,abs_error,rel_error,status,notes")
        for r in rows
            notes = replace(r.notes, "," => ";")
            println(io,
                "$(r.dim),$(r.function_name),$(r.method),$(r.time_ns),$(r.value),$(r.exact)," *
                "$(r.abs_error),$(r.rel_error),$(r.status),$(notes)")
        end
    end
end

function write_full_markdown(path::AbstractString, rows)
    open(path, "w") do io
        println(io, "# Benchmark Timings at Fixed Accuracy Target")
        println(io)
        println(io, "- Tolerances: `rtol=$(RTOL)`, `atol=$(ATOL)`")
        println(io, "- Max evaluations (external adaptive solvers): `$(MAXEVALS)`")
        println(io, "- Timing method: `@belapsed` with interpolation (`samples=$(BENCH_SAMPLES)`, `evals=$(BENCH_EVALS)`).")
        println(io, "- One warm call is executed before each timed benchmark (reduces first-call compilation effects).")
        println(io, "- Adaptive cache construction is excluded from timed regions (caches are prebuilt).")
        println(io, "- `FTS avx` is run at the `N` inferred from adaptive FTS convergence.")
        println(io, "- `FastGauss` (1D) is independently calibrated to the tolerance target.")
        println(io)
        println(io, "| Dim | Function | Method | Time (ns) | Abs Error | Rel Error | Status | Notes |")
        println(io, "| :-- | :------- | :----- | --------: | --------: | --------: | :----- | :---- |")
        for r in rows
            println(io, "| $(r.dim) | $(r.function_name) | $(r.method) | $(format_float(r.time_ns)) | " *
                "$(format_float(r.abs_error)) | $(format_float(r.rel_error)) | $(r.status) | $(r.notes) |")
        end
    end
end

function write_summary_markdown(path::AbstractString, rows)
    method_order = Dict(
        "FTS adaptive (typed)" => 1,
        "FTS quad (convenience)" => 2,
        "FTS avx (precomputed)" => 3,
        "QuadGK" => 4,
        "HCubature" => 5,
        "CubatureJLh" => 6,
        "CubatureJLp" => 7,
        "CubaVegas" => 8,
        "CubaDivonne" => 9,
        "CubaCuhre" => 10,
        "FastGauss (precomputed)" => 11,
    )
    grouped = Dict{Tuple{String,String},Vector{NamedTuple}}()
    for r in rows
        key = (r.dim, r.function_name)
        push!(get!(grouped, key, NamedTuple[]), r)
    end

    open(path, "w") do io
        println(io, "## Benchmark Timing Table at Fixed Accuracy Target")
        println(io)
        println(io, "| Dim | Function | FTS adaptive | FTS quad | FTS avx | QuadGK | HCubature | Cubature h | Cubature p | Cuba Vegas | Cuba Divonne | Cuba Cuhre | FastGauss |")
        println(io, "| :-- | :------- | ----------: | -------: | ------: | -----: | --------: | ---------: | ---------: | ---------: | -----------: | ---------: | --------: |")
        for key in sort(collect(keys(grouped)))
            dim, fname = key
            bucket = grouped[key]
            sort!(bucket, by=r -> get(method_order, r.method, 999))
            byname = Dict(r.method => r for r in bucket)
            function cell(m::String)
                if !haskey(byname, m)
                    return "n/a"
                end
                r = byname[m]
                marker = r.status == "ok" ? "" : " *"
                return "$(format_float(r.time_ns))$(marker)"
            end
            c_adapt = cell("FTS adaptive (typed)")
            c_quad = cell("FTS quad (convenience)")
            c_avx = cell("FTS avx (precomputed)")
            c_qgk = cell("QuadGK")
            c_hcub = cell("HCubature")
            c_cubh = cell("CubatureJLh")
            c_cubp = cell("CubatureJLp")
            c_vegas = cell("CubaVegas")
            c_divonne = cell("CubaDivonne")
            c_cuhre = cell("CubaCuhre")
            c_gq = cell("FastGauss (precomputed)")

            println(io, "| $(dim) | $(fname) | $(c_adapt) | $(c_quad) | $(c_avx) | $(c_qgk) | " *
                "$(c_hcub) | $(c_cubh) | $(c_cubp) | $(c_vegas) | $(c_divonne) | $(c_cuhre) | $(c_gq) |")
        end
        println(io)
        println(io, "`*` indicates result did not meet the requested tolerance.")
        println(io, "`FTS avx` uses the `N` inferred from adaptive FTS convergence on each case.")
        println(io, "`FastGauss` (1D) uses its own minimal `N` that meets tolerance (if found).")
    end
end

function run_benchmark(; print_table::Bool=true)
    rows = NamedTuple[]
    for prob in PROBLEMS
        append!(rows, benchmark_problem(prob))
    end

    sort!(rows, by=r -> (r.dim, r.function_name, r.method))

    result_dir = joinpath(@__DIR__, "results")
    mkpath(result_dir)
    csv_path = joinpath(result_dir, "timings.csv")
    full_md_path = joinpath(result_dir, "timings_full.md")
    summary_md_path = joinpath(result_dir, "timings_summary.md")

    write_csv(csv_path, rows)
    write_full_markdown(full_md_path, rows)
    write_summary_markdown(summary_md_path, rows)

    println("Wrote benchmark results:")
    println("  - $(csv_path)")
    println("  - $(full_md_path)")
    println("  - $(summary_md_path)")

    if print_table
        println(read(summary_md_path, String))
    end

    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmark()
end
