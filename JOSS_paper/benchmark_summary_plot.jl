using CairoMakie

const QUAD_COLOR = RGBAf(0.0, 0.447, 0.698, 1.0)
const AVX_COLOR = RGBAf(0.835, 0.369, 0.0, 1.0)

struct BenchmarkRow
    dim::String
    fname::String
    method::String
    time_ns::Float64
    status::String
end

function read_rows(path::AbstractString)
    lines = readlines(path)
    rows = BenchmarkRow[]
    for line in Iterators.drop(lines, 1)
        fields = split(chomp(line), ',')
        length(fields) >= 9 || continue
        push!(rows, BenchmarkRow(
            fields[1],
            fields[2],
            fields[3],
            parse(Float64, fields[4]),
            fields[9],
        ))
    end
    return rows
end

function benchmark_order(case)
    order = Dict(
        ("1D", "1/(1+25x^2)") => 1,
        ("1D", "x^6 - 2x^3 + 0.5") => 2,
        ("1D", "log(1-x)") => 3,
        ("1D", "1/sqrt(1-x^2)") => 4,
        ("2D", "x^2 + y^2") => 5,
        ("2D", "exp(x+y)") => 6,
        ("2D", "1/sqrt((1-x^2)(1-y^2))") => 7,
        ("3D", "x^2*y^2*z^2") => 8,
        ("3D", "exp(x+y+z)") => 9,
        ("3D", "1/sqrt((1-x^2)(1-y^2)(1-z^2))") => 10,
    )
    return order[case]
end

function label_for(case, baseline)
    dim, fname = case
    baseline_label = Dict(
        "FastGauss (precomputed)" => "FastGauss",
        "QuadGK" => "QuadGK",
        "CubatureJLh" => "Cubature h",
        "CubatureJLp" => "Cubature p",
        "HCubature" => "HCubature",
        "CubaCuhre" => "Cuba Cuhre",
    )[baseline]
    fname_label = Dict(
        "1/(1+25x^2)" => "1/(1+25x^2)",
        "x^6 - 2x^3 + 0.5" => "polynomial",
        "log(1-x)" => "log(1-x)",
        "1/sqrt(1-x^2)" => "(1-x^2)^(-1/2)",
        "x^2 + y^2" => "x^2+y^2",
        "exp(x+y)" => "exp(x+y)",
        "1/sqrt((1-x^2)(1-y^2))" => "endpoint-singular",
        "x^2*y^2*z^2" => "x^2 y^2 z^2",
        "exp(x+y+z)" => "exp(x+y+z)",
    )[fname]
    return "$(dim) $(fname_label)\n[$(baseline_label)]"
end

function summarize(rows)
    grouped = Dict{Tuple{String, String}, Vector{BenchmarkRow}}()
    for row in rows
        key = (row.dim, row.fname)
        push!(get!(grouped, key, BenchmarkRow[]), row)
    end

    cases = sort(collect(keys(grouped)); by=benchmark_order)
    labels = String[]
    quad_speedups = Float64[]
    avx_speedups = Float64[]

    for case in cases
        case_rows = grouped[case]
        others = filter(r -> !startswith(r.method, "FTS") && r.status == "ok", case_rows)
        isempty(others) && continue
        baseline = findmin(r -> r.time_ns, others)[2]
        baseline_row = others[baseline]

        quad_row = only(filter(r -> r.method == "FTS quad (convenience)", case_rows))
        avx_row = only(filter(r -> r.method == "FTS avx (precomputed)", case_rows))

        push!(labels, label_for(case, baseline_row.method))
        push!(quad_speedups, baseline_row.time_ns / quad_row.time_ns)
        push!(avx_speedups, baseline_row.time_ns / avx_row.time_ns)
    end

    return labels, quad_speedups, avx_speedups
end

function plot_summary(labels, quad_speedups, avx_speedups)
    CairoMakie.activate!()
    fig = Figure(size=(1450, 900), fontsize=24, figure_padding=(20, 20, 10, 10))
    ax = Axis(
        fig[1, 1],
        xlabel="Speedup vs fastest accurate competing method",
        ylabel="Benchmark case",
        xscale=log10,
        xgridvisible=true,
        ygridvisible=true,
        yticks=(1:length(labels), labels),
        xminorticksvisible=true,
        xticks=([0.02, 0.1, 1, 10, 100, 1000, 10000, 100000],
            ["0.02", "0.1", "1", "10", "10^2", "10^3", "10^4", "10^5"]),
        xlabelsize=28,
        ylabelsize=28,
        xticklabelsize=20,
        yticklabelsize=18,
    )

    y = collect(1:length(labels))
    y_quad = y .+ 0.14
    y_avx = y .- 0.14

    hlines!(ax, y, color=(:gray90, 0.9), linewidth=1)
    vlines!(ax, [1.0], color=:gray55, linestyle=:dash, linewidth=1.5)
    scatter!(ax, quad_speedups, y_quad, color=:white, strokecolor=QUAD_COLOR, strokewidth=3,
        marker=:circle, markersize=20,
        label="FastTanhSinh quad")
    scatter!(ax, avx_speedups, y_avx, color=AVX_COLOR, strokecolor=:black, strokewidth=1.5,
        marker=:utriangle, markersize=24,
        label="FastTanhSinh avx")

    text!(ax, 1.1, length(labels) + 0.55, text="speedup > 1 means faster", color=:gray35,
        align=(:left, :bottom), fontsize=18)
    xlims!(ax, 0.02, 1.0e5)
    ylims!(ax, 0.5, length(labels) + 0.8)
    ax.yreversed = true
    Legend(fig[2, 1], ax;
        orientation=:horizontal,
        nbanks=2,
        framevisible=false,
        labelsize=20,
        patchsize=(44, 20),
        tellwidth=false)
    rowgap!(fig.layout, 8)

    return fig
end

function main()
    rows = read_rows(joinpath(@__DIR__, "..", "benchmark", "results", "timings.csv"))
    labels, quad_speedups, avx_speedups = summarize(rows)
    fig = plot_summary(labels, quad_speedups, avx_speedups)
    save(joinpath(@__DIR__, "benchmark_summary.svg"), fig)
    save(joinpath(@__DIR__, "benchmark_summary.png"), fig)
end

main()
