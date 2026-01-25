using FastTanhSinhQuadrature
using DoubleFloats
using LaTeXStrings
using GaussQuadrature
using CairoMakie

# Set BigFloat precision for high-accuracy benchmarks
setprecision(BigFloat, 256)
const T = BigFloat

# Define test functions and their exact integrals over [-1, 1] using BigFloat
test_functions = [
    # (name=L"\exp{x}", f=x -> exp(x), exact=T(exp(one(T))) - T(exp(-one(T)))),
    (name=L"\log(1+x)", f=x -> log(one(T) + x), exact=2 * log(T(2.0)) - T(2.0)),
    (name=L"1/(1+25x^2)", f=x -> 1 / (1 + 25 * x^2), exact=(T(2) / 5) * atan(T(5))),
    # (name=L"\sqrt{1-x^2}", f=x -> sqrt(1 - x^2), exact=T(π) / 2),
    # (name=L"x^{-1/4}\log(1/x)", f=x -> x^(-1 / 4) * log(1 / x), exact=T(16 / 9), range=T[0, 1]),
    (name=L"1/((x-2)(1-x)^{1/4}(1+x)^{3/4})", f=x -> 1 / ((x - 2) * (1 - x)^(1 / 4) * (1 + x)^(3 / 4)),
        exact=T("-1.949054259166747153657919113305184895821287200233066621785270125453326989447444885652648474542392802962572023494565808280455227350179904961894167915612983522793337378241965509394184885218503432473597542343269605628029234544188489125493174675613978756949408472350567259309791818066327536231903217201754779798024608668935810330010916825646585937869564656598924388817496515245290999475576102379748663946800010918013813610755103525742537281693887080238073231089698062623921438536568514121793542784984754848794609547174888340437525644514849527933316215226932197574446566338550419210926373034326391710478581965058200796649082090043398648116397620152842649085759285786383816326181758973349850630607502231380")),
    (name=L"\cos(\pi x)/\sqrt{1-x}", f=x -> cos(T(π) * x) / sqrt(one(T) - x), exact=T("-0.6904945887466050171527986111031877731184655785803339772826084001010560017822902042337565206614402623577274398229770150780304083926437098116947542807515746484283997389463793613930062215356510096141561417612700140904036249076398342266906001017584532421377530627373501812827614569857500092160444297408756579607546779687747497082859607927892078705630782776343639829983638709983890737070577232368080443846227859640840075063358751926065976450952302300905754593246223245550366901186644406720576334897902831323529556227663465197241346272735679485008981091265919131861875954158442735486104285131383407943329053800452709543824162308072844356094437494123306819451016838051727406791012618268808253826429458291598"))
]

function compute_convergence_errors(tf, levels)
    ns = Int[]
    errors = T[]

    for n_points in levels
        x, w, h = tanhsinh(T, n_points)
        if hasproperty(tf, :range)
            val = integrate1D_avx(tf.f, tf.range[1], tf.range[2], x, w, h)
        else
            val = integrate1D_avx(tf.f, x, w, h)
        end

        push!(ns, n_points)
        err = abs((val - tf.exact))# / exact)
        push!(errors, err)
    end
    return ns, errors
end

function compute_gauss_errors(tf, levels)
    ns = Int[]
    errors = T[]

    for n_points in levels
        x, w = legendre(T, n_points)
        val = sum(w .* tf.f.(x))

        push!(ns, n_points)
        err = abs(val - tf.exact)
        push!(errors, err)
    end
    return ns, errors
end

function generate_convergence_plots()
    CairoMakie.activate!()
    CairoMakie.set_theme!(mytheme_aps())
    fig = Figure(size=(2 * 253, 2 * 200), figure_padding=(1, 5, 5, 1), fontsize=20)
    ax = Axis(fig[1, 1],
        xlabel=L"\textrm{Number of Points (N)}",
        ylabel=L"\textrm{Error} \, |I - I_h|",
        yscale=log10,
        # xscale=log10,
        title="Convergence of Tanh-Sinh Quadrature",
        xgridvisible=true,
        ygridvisible=true)

    levels = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700]

    for (i, tf) in enumerate(test_functions)
        println("Computing errors for $(tf.name) with analytical exact value $(Float64(tf.exact))")

        # Tanh-Sinh
        ns, errors = compute_convergence_errors(tf, levels)
        sl_ts = scatterlines!(ax, ns, errors, label=tf.name, linewidth=2, markersize=8)
        # color = sl_ts.color[] # Extract color for Gauss comparison

        # Gauss-Legendre
        # ns_g, errors_g = compute_gauss_errors(tf.f, tf.exact, levels)
        # scatterlines!(ax, ns_g, errors_g, color=color, linestyle=:dash, linewidth=2, markersize=8)
    end

    axislegend(ax; position=:rb, labelsize=10, rowgap=0.2, margin=(0, 0, 25, 0))
    rowsize!(fig.layout, 1, Auto())
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)
    # Save the plot
    save("docs/src/assets/convergence.svg", fig)
    save("docs/src/assets/convergence.png", fig)
    println("Saved convergence plots to docs/src/assets/convergence.{svg,png}")
end

function mytheme_aps()
    return Theme(
        # Axis attributes
        ;
        Axis=Attributes(; spinewidth=1.1,
            xgridvisible=true,
            xlabelpadding=-0,
            xlabelsize=12,
            xminortickalign=1,
            xminorticks=IntervalsBetween(5, true),
            xminorticksize=3,
            xminorticksvisible=true,
            xminortickwidth=0.75,
            xtickalign=1,
            xticklabelsize=8,
            xticksize=5,
            xticksmirrored=true,
            xtickwidth=0.8,
            ygridvisible=true,
            ylabelpadding=2,
            ylabelsize=12,
            yminortickalign=1,
            yminorticks=IntervalsBetween(5, true),
            yminorticksize=3,
            yminorticksvisible=true,
            yminortickwidth=0.75,
            ytickalign=1,
            yticklabelsize=10,
            yticksize=5,
            yticksmirrored=true,
            ytickwidth=0.8,
            xticklabelfont="cmr10",  # Upright Computer Modern
            yticklabelfont="cmr10",  # Upright Computer Modern
            xticklabelstyle=Attributes(; italic=false),
            yticklabelstyle=Attributes(; italic=false)),
        # General figure settings
        colgap=8,
        figure_padding=0,
        rowgap=8,
        size=(243, 165),
        # Colorbar attributes
        Colorbar=Attributes(; labelpadding=2,
            labelsize=10,
            minortickalign=1,
            minorticksize=3,
            minorticksvisible=true,
            minortickwidth=0.75,
            size=8,
            spinewidth=1.1,
            tickalign=1,
            ticklabelpad=2,
            ticklabelsize=8,
            ticksize=5,
            tickwidth=0.8),
        fonts=Attributes(; bold="NewComputerModern10 Bold",
            bold_italic="NewComputerModern10 Bold Italic",
            italic="NewComputerModern10 Italic",
            regular="NewComputerModern Math Regular"),
        # fonts=Attributes(; bold="ComputerModern Bold",
        #                  bold_italic="ComputerModern Bold Italic",
        #                  italic="ComputerModern Italic",
        #                  regular="ComputerModern Math Regular"),
        # Legend attributes
        Legend=Attributes(; colgap=4,
            framecolor=(:grey, 0.5),
            framevisible=false,
            labelsize=7.5,
            margin=(0, 0, 0, 0),
            nbanks=1,
            padding=(2, 2, 2, 2),
            rowgap=-10,
            #labelfont="cmr10"
        ),
        # Lines attributes
        Lines=Attributes(;
            cycle=Cycle([[:color] => :color],
                true)),
        # Scatter attributes
        Scatter=Attributes(;
            cycle=Cycle([[:color] => :color, [:marker] => :marker],
                true),
            markersize=7,
            strokewidth=0),
        markersize=7,
        # Palette attributes
        palette=Attributes(;
            color=[RGBAf(0.298039, 0.447059, 0.690196, 1.0),
                RGBAf(0.866667, 0.517647, 0.321569, 1.0),
                RGBAf(0.333333, 0.658824, 0.407843, 1.0),
                RGBAf(0.768627, 0.305882, 0.321569, 1.0),
                RGBAf(0.505882, 0.447059, 0.701961, 1.0),
                RGBAf(0.576471, 0.470588, 0.376471, 1.0),
                RGBAf(0.854902, 0.545098, 0.764706, 1.0),
                RGBAf(0.54902, 0.54902, 0.54902, 1.0),
                RGBAf(0.8, 0.72549, 0.454902, 1.0),
                RGBAf(0.392157, 0.709804, 0.803922, 1.0)],
            linestyle=[nothing, :dash, :dot, :dashdot, :dashdotdot],
            marker=[:circle, :rect, :dtriangle, :utriangle, :cross,
                :diamond, :ltriangle, :rtriangle, :pentagon,
                :xcross, :hexagon],
            markercolor=[RGBAf(0.298039, 0.447059, 0.690196, 1.0),
                RGBAf(0.866667, 0.517647, 0.321569, 1.0),
                RGBAf(0.333333, 0.658824, 0.407843, 1.0),
                RGBAf(0.768627, 0.305882, 0.321569, 1.0),
                RGBAf(0.505882, 0.447059, 0.701961, 1.0),
                RGBAf(0.576471, 0.470588, 0.376471, 1.0),
                RGBAf(0.854902, 0.545098, 0.764706, 1.0),
                RGBAf(0.54902, 0.54902, 0.54902, 1.0),
                RGBAf(0.8, 0.72549, 0.454902, 1.0),
                RGBAf(0.392157, 0.709804, 0.803922, 1.0)],
            patchcolor=[RGBAf(0.298039, 0.447059, 0.690196, 1.0),
                RGBAf(0.866667, 0.517647, 0.321569, 1.0),
                RGBAf(0.333333, 0.658824, 0.407843, 1.0),
                RGBAf(0.768627, 0.305882, 0.321569, 1.0),
                RGBAf(0.505882, 0.447059, 0.701961, 1.0),
                RGBAf(0.576471, 0.470588, 0.376471, 1.0),
                RGBAf(0.854902, 0.545098, 0.764706, 1.0),
                RGBAf(0.54902, 0.54902, 0.54902, 1.0),
                RGBAf(0.8, 0.72549, 0.454902, 1.0),
                RGBAf(0.392157, 0.709804, 0.803922, 1.0)]),
        # PolarAxis attributes
        PolarAxis=Attributes(; spinewidth=1.1))
end


generate_convergence_plots()
