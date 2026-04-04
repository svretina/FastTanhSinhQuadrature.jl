using CairoMakie
using LaTeXStrings

const QUAD_COLOR = RGBAf(0.0, 0.447, 0.698, 1.0)
const AVX_COLOR = RGBAf(0.835, 0.369, 0.0, 1.0)
const LINK_COLOR = RGBAf(0.55, 0.55, 0.55, 0.7)

struct BenchmarkRow
    dim::String
    fname::String
    method::String
    time_ns::Float64
    status::String
end

struct CaseSummary
    dim::String
    fname::String
    baseline::String
    quad_speedup::Float64
    avx_speedup::Float64
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
        ("1D", "sin^2(1000x)") => 3,
        ("1D", "log(1-x)") => 4,
        ("1D", "1/sqrt(1-x^2)") => 5,
        ("2D", "x^2 + y^2") => 6,
        ("2D", "exp(x+y)") => 7,
        ("2D", "1/sqrt((1-x^2)(1-y^2))") => 8,
        ("3D", "x^2*y^2*z^2") => 9,
        ("3D", "exp(x+y+z)") => 10,
        ("3D", "1/sqrt((1-x^2)(1-y^2)(1-z^2))") => 11,
    )
    return order[case]
end

function baseline_tex(method::String)
    return Dict(
        "FastGauss (precomputed)" => raw"\mathrm{FastGauss}",
        "QuadGK" => raw"\mathrm{QuadGK}",
        "CubatureJLh" => raw"\mathrm{Cubature\ h}",
        "CubatureJLp" => raw"\mathrm{Cubature\ p}",
        "HCubature" => raw"\mathrm{HCubature}",
        "CubaCuhre" => raw"\mathrm{Cuba\ Cuhre}",
        "CubaDivonne" => raw"\mathrm{Cuba\ Divonne}",
        "CubaVegas" => raw"\mathrm{Cuba\ Vegas}",
    )[method]
end

function case_expr_tex(fname::String)
    return Dict(
        "1/(1+25x^2)" => raw"\frac{1}{1+25x^2}",
        "x^6 - 2x^3 + 0.5" => raw"x^6-2x^3+\frac{1}{2}",
        "sin^2(1000x)" => raw"\sin^2(1000x)",
        "log(1-x)" => raw"\log(1-x)",
        "1/sqrt(1-x^2)" => raw"\frac{1}{\sqrt{1-x^2}}",
        "x^2 + y^2" => raw"x^2+y^2",
        "exp(x+y)" => raw"e^{x+y}",
        "1/sqrt((1-x^2)(1-y^2))" => raw"\frac{1}{\sqrt{(1-x^2)(1-y^2)}}",
        "x^2*y^2*z^2" => raw"x^2y^2z^2",
        "exp(x+y+z)" => raw"e^{x+y+z}",
        "1/sqrt((1-x^2)(1-y^2)(1-z^2))" => raw"\frac{1}{\sqrt{(1-x^2)(1-y^2)(1-z^2)}}",
    )[fname]
end

function case_label(summary::CaseSummary)
    return latexstring(case_expr_tex(summary.fname), raw"\;[", baseline_tex(summary.baseline), raw"]")
end

function summarize(rows)
    grouped = Dict{Tuple{String, String}, Vector{BenchmarkRow}}()
    for row in rows
        key = (row.dim, row.fname)
        push!(get!(grouped, key, BenchmarkRow[]), row)
    end

    summaries = CaseSummary[]
    for case in sort(collect(keys(grouped)); by=benchmark_order)
        case_rows = grouped[case]
        others = filter(r -> !startswith(r.method, "FTS") && r.status == "ok", case_rows)
        isempty(others) && continue

        baseline_ix = findmin(r -> r.time_ns, others)[2]
        baseline = others[baseline_ix]
        quad_row = only(filter(r -> r.method == "FTS quad (convenience)", case_rows))
        avx_row = only(filter(r -> r.method == "FTS avx (precomputed)", case_rows))

        push!(summaries, CaseSummary(
            case[1],
            case[2],
            baseline.method,
            baseline.time_ns / quad_row.time_ns,
            baseline.time_ns / avx_row.time_ns,
        ))
    end
    return summaries
end

function mytheme_aps()
    return Theme(
        Axis=Attributes(;
            spinewidth=1.1,
            xgridvisible=true,
            ygridvisible=true,
            xminorticks=IntervalsBetween(5, true),
            xminorticksvisible=true,
            xminortickalign=1,
            xminorticksize=3,
            xminortickwidth=0.75,
            yminorticks=IntervalsBetween(2, true),
            yminorticksvisible=false,
            xtickalign=1,
            ytickalign=1,
            xticksize=5,
            yticksize=5,
            xtickwidth=0.8,
            ytickwidth=0.8,
        ),
        Legend=Attributes(;
            framevisible=false,
            patchsize=(36, 18),
            rowgap=4,
            colgap=8,
        ),
        fonts=Attributes(;
            bold="NewComputerModern10 Bold",
            bold_italic="NewComputerModern10 Bold Italic",
            italic="NewComputerModern10 Italic",
            regular="NewComputerModern Math Regular",
        ),
    )
end

function dim_title(dim::String)
    return dim == "1D" ? L"\mathrm{1D\ benchmarks}" :
           dim == "2D" ? L"\mathrm{2D\ benchmarks}" :
           L"\mathrm{3D\ benchmarks}"
end

function plot_summary(summaries::Vector{CaseSummary})
    CairoMakie.activate!()
    CairoMakie.set_theme!(mytheme_aps())

    fig = Figure(size=(1760, 1320), figure_padding=(24, 26, 14, 12), fontsize=26)
    dims = ("1D", "2D", "3D")
    tick_vals = [0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0]
    tick_labs = [L"10^{-1}", L"10^{0}", L"10^{1}", L"10^{2}", L"10^{3}", L"10^{4}"]
    legend_handles = nothing

    for (rowid, dim) in enumerate(dims)
        sub = filter(s -> s.dim == dim, summaries)
        isempty(sub) && continue

        labels = case_label.(sub)
        y = collect(1:length(sub))
        yq = y .+ 0.12
        ya = y .- 0.12

        ax = Axis(
            fig[rowid, 1],
            title=dim_title(dim),
            xscale=log10,
            xticks=(tick_vals, tick_labs),
            yticks=(y, labels),
            xgridvisible=true,
            ygridvisible=true,
            xticklabelsize=22,
            yticklabelsize=21,
            titlegap=8,
            titlealign=:left,
            xlabel=rowid == length(dims) ?
                   L"\mathrm{Speedup}\;=\;t_{\mathrm{baseline}}/t_{\mathrm{method}}" : "",
        )

        for yi in y
            hlines!(ax, [yi], color=(:gray92, 0.9), linewidth=1)
        end

        vlines!(ax, [1.0], color=:gray45, linestyle=:dash, linewidth=1.4)
        for (k, case) in enumerate(sub)
            lines!(ax, [case.quad_speedup, case.avx_speedup], [yq[k], ya[k]];
                color=LINK_COLOR, linewidth=2.4)
        end

        quad_handle = scatter!(ax, getfield.(sub, :quad_speedup), yq;
            color=:white, strokecolor=QUAD_COLOR, strokewidth=2.7,
            marker=:circle, markersize=20, label=L"\mathrm{FastTanhSinh\;quad}")
        avx_handle = scatter!(ax, getfield.(sub, :avx_speedup), ya;
            color=AVX_COLOR, strokecolor=:black, strokewidth=1.2,
            marker=:utriangle, markersize=23, label=L"\mathrm{FastTanhSinh\;integrate^{\ast}_{avx}}")

        legend_handles === nothing && (legend_handles = (quad_handle, avx_handle))

        xlims!(ax, 0.07, 1.0e4)
        ylims!(ax, 0.5, length(sub) + 0.6)
        ax.yreversed = true

        if rowid == 1
            text!(ax, 1.12, 0.55;
                text=L"\mathrm{speedup}>1\;\Rightarrow\;\mathrm{faster}",
                color=:gray35, fontsize=19, align=(:left, :bottom))
        end
    end

    if legend_handles !== nothing
        Legend(fig[4, 1], [legend_handles[1], legend_handles[2]],
            [L"\mathrm{Adaptive\ convenience\ path}", L"\mathrm{Precomputed\ SIMD\ path}"];
            orientation=:horizontal, nbanks=2, tellwidth=false, labelsize=21)
    end
    rowgap!(fig.layout, 8)
    return fig
end

function main()
    repo_root = normpath(joinpath(@__DIR__, ".."))
    csv_path = joinpath(repo_root, "benchmark", "results", "timings.csv")
    rows = read_rows(csv_path)
    summaries = summarize(rows)
    fig = plot_summary(summaries)

    docs_assets = joinpath(repo_root, "docs", "src", "assets")
    bench_assets = joinpath(repo_root, "benchmark", "results")
    mkpath(docs_assets)
    mkpath(bench_assets)

    save(joinpath(@__DIR__, "benchmark_summary.svg"), fig)
    save(joinpath(@__DIR__, "benchmark_summary.png"), fig, px_per_unit=2)
    save(joinpath(docs_assets, "benchmark_summary.svg"), fig)
    save(joinpath(docs_assets, "benchmark_summary.png"), fig, px_per_unit=2)
    save(joinpath(bench_assets, "benchmark_summary.svg"), fig)
    save(joinpath(bench_assets, "benchmark_summary.png"), fig, px_per_unit=2)
    println("Saved benchmark summary plot to JOSS_paper, docs/src/assets, and benchmark/results.")
end

main()
