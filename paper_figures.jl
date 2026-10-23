using DataFrames
using CSV
using Plots
include("coherent_model.jl")
include("incoherent_model.jl")
#----------------------------------------------------------------------------------------
# Figures of the article and its supplement, from the optimization results in
# ./opt_coh, ./opt_inc and ./opt_fixed_alpha (run the three *_optimization.jl
# scripts first). Both forward models are integrated over the ARIEL antenna
# pattern; Fig. 3 also shows both at a single incidence angle of 40°.
#
# Output (./figures):
#   fig03_tb_vs_snow_depth.png          Fig. 3   TB vs snow depth, coherent and incoherent, 40° and footprint
#   fig04_geophysical_coh.png           Fig. 4   priors and optimal geophysical variables, coherent
#   fig04_geophysical_inc.png           Fig. 4   same, incoherent
#   fig05_tb_prior_optimal.png          Fig. 5   observed, prior-state and optimal TBs per measurement
#   fig06_costs.png                     Fig. 6   minimized cost per measurement and vs optimal snow depth
#   fig07_cost_histograms.png           Fig. 7   histograms of the minimized cost
#   fig08_tb_vs_snow_depth_sites.png    Fig. 8   TB vs snow depth per site, observed and optimal
#   pd_prior_optimal.png                         measured vs modeled PD at the prior and optimal states
#   figS1_fixed_alpha.png               Fig. S1  observed vs simulated TB with α fixed at 5
#----------------------------------------------------------------------------------------

const OUT = joinpath(@__DIR__, "figures")

# Site priors (identical to coh/inc-model_optimization.jl)
const SITES = [
    (name = "CLO", key = "clo", temp = -9.69, tsnow = [-14.0],
     rho   = [355.0, 300, 300, 355, 300, 267, 300, 300, 300, 350],
     dice  = [0.945, 0.95, 0.94, 0.94, 0.97, 0.985, 0.98, 0.895, 0.895, 0.99],
     tice  = fill(-13.0, 10), sice = fill(5.32, 10), ar = 5.0),
    (name = "MID-1", key = "mid1", temp = -7.37, tsnow = [-15.0],
     rho   = [350.0, 365, 367, 350, 350, 350],
     dice  = [0.89, 0.855, 0.855, 0.89, 0.855, 0.855],
     tice  = fill(-11.0, 6), sice = fill(4.6, 6), ar = 5.0),
    (name = "MID-2", key = "mid2", temp = -12.95, tsnow = [-25.0],
     rho   = [350.0, 350, 300, 347],
     dice  = [0.93, 0.93, 0.93, 0.93],
     tice  = fill(-20.0, 4), sice = fill(4.5, 4), ar = 5.0),
    (name = "FAR", key = "far", temp = -13.86, tsnow = [-25.0],
     rho   = [385.0, 385, 385, 385, 365, 400, 400, 400, 385, 385, 365, 400, 400, 400, 400],
     dice  = [0.86, 0.86, 0.86, 0.86, 0.92, 0.84, 0.84, 0.845, 0.855, 0.855, 0.855, 0.855, 0.855, 0.855, 0.855],
     tice  = fill(-17.0, 15), sice = fill(4.8, 15), ar = 5.0),
]
const SITE_COLORS     = [:darkgreen, :mediumblue, :darkorchid4, :brown4]
const SITE_BOUNDARIES = [10.5, 16.5, 20.5]
const SITE_LABEL_X    = [5, 13.5, 18.5, 28]
const SITE_LABELS     = ["CLO", "MID-1", "MID-2", "FAR"]

site_data(data, s) = filter(row -> row[:temp] == s.temp, data)
read_opt(dir, s) = CSV.read(joinpath(@__DIR__, dir, "opt_res_$(s.key)_jl.csv"), DataFrame)

function stats(sim, obs)
    d = sim .- obs; n = length(d)
    mx = sum(obs) / n; my = sum(sim) / n
    r = sum((obs .- mx) .* (sim .- my)) / sqrt(sum((obs .- mx) .^ 2) * sum((sim .- my) .^ 2))
    return sum(d) / n, sqrt(sum(d .^ 2) / n), r
end

# Forward models evaluated at the prior state x_a
function prior_tb(site, dsnow_prior, model)
    tbh = Float64[]; tbv = Float64[]
    for i in eachindex(dsnow_prior)
        h, v = model == :coh ?
            wilheit_model_integrated(40, [dsnow_prior[i]], site.tsnow, [site.rho[i]], [site.dice[i]],
                                     [site.tice[i]], [site.sice[i]], "mix", 0, site.ar) :
            incoherent_model_integrated(40, dsnow_prior[i], site.tsnow[1], site.rho[i], site.dice[i],
                                        site.tice[i], site.sice[i], "mix", 0, site.ar)
        push!(tbh, h); push!(tbv, v)
    end
    return tbh, tbv
end

#----------------------------------------------------------------------------------------
# Fig. 3
#----------------------------------------------------------------------------------------
function fig03()
    ts, rho, di, ti, si, ar, tw = -15.0, 300.0, 0.9, -10.0, 5.0, 5.0, -1.8
    dsnows = 0:0.001:0.40
    coh40 = [wilheit_model(40, [d], [ts], [rho], [di], [ti], [si], "mix", 0, ar) for d in dsnows]
    cohi  = [wilheit_model_integrated(40, [d], [ts], [rho], [di], [ti], [si], "mix", 0, ar) for d in dsnows]
    d_inc = dsnows[2:end]
    inc40 = [incoherent_model(40, d, ts, rho, di, ti, si, "mix", 0, ar) for d in d_inc]
    inci  = [incoherent_model_integrated(40, d, ts, rho, di, ti, si, "mix", 0, ar) for d in d_inc]

    # snow-free incoherent TB (air / ice / water): the incoherent model with a zero-thickness
    # snow layer keeps the air/snow and snow/ice interfaces, so d = 0 is drawn separately
    function nosnow(theta)
        epre_ice, epim_ice = e_eff_mix(ar, 0, 0, si, ti)
        epre_w, epim_w = epsilon_water(1.4e9, tw, 33)
        n = [sqrt(complex(1, 0)), sqrt(complex(epre_ice, epim_ice)), sqrt(complex(epre_w, epim_w))]
        return incoherent(di, n, deg2rad(theta), 1.4, [ti + 273.15, tw + 273.15])
    end
    w_h, w_v = antenna_pattern_integration(40, 15.29, 14.87)
    ns = [nosnow(t) for t in range(-89.99, stop = 89.99, length = 1000)]
    ns40 = nosnow(40)
    nsi = (sum(getindex.(ns, 1) .* w_h), sum(getindex.(ns, 2) .* w_v))

    p = plot(legend = :outertop, legend_columns = 2, dpi = 500, size = (800, 520), framestyle = :box)
    for (k, col, ch) in ((2, :red, "TBV"), (1, :blue, "TBH"))
        plot!(p, dsnows, getindex.(coh40, k), color = col, label = "$ch - coherent, fixed 40°")
        plot!(p, dsnows, getindex.(cohi, k), color = col, linestyle = :dash, label = "$ch - coherent, ARIEL footprint")
        plot!(p, d_inc, getindex.(inc40, k), color = col, linestyle = :dot, label = "$ch - incoherent, fixed 40°")
        plot!(p, d_inc, getindex.(inci, k), color = col, linestyle = :dashdot, label = "$ch - incoherent, ARIEL footprint")
        plot!(p, [0.0, 0.0], [ns40[k], inc40[1][k]], color = col, linestyle = :dot, label = "")
        plot!(p, [0.0, 0.0], [nsi[k], inci[1][k]], color = col, linestyle = :dashdot, label = "")
        scatter!(p, [0.0, 0.0], [ns40[k], nsi[k]], color = col, marker = :circle, markersize = 1.5,
                 markerstrokewidth = 0, label = "")
    end
    xlabel!(p, "Snow depth (m)"); ylabel!(p, "Brightness temperature (K)")
    savefig(p, joinpath(OUT, "fig03_tb_vs_snow_depth.png"))
end

#----------------------------------------------------------------------------------------
# Fig. 4
#----------------------------------------------------------------------------------------
function fig04(data, opt_dir, outname)
    opt = Dict(s.key => read_opt(opt_dir, s) for s in SITES)
    col(c) = vcat([opt[s.key][!, c] for s in SITES]...)
    dsnow_prior = vcat([site_data(data, s)[!, :dsnow] .* 0.01 for s in SITES]...)
    p = plot(layout = (6, 1), size = (600, 1200), legend = :topright, legend_columns = 2, dpi = 500,
             bottom_margin = 3Plots.mm, left_margin = 12Plots.mm, right_margin = 1.5Plots.mm, top_margin = 3Plots.mm)
    function row!(sub, opt_vals, prior_vals, ylab, lims)
        plot!(p, opt_vals, seriestype = :scatter, markersize = 6, marker = :star5, framestyle = :box, xticks = false,
              ylabel = ylab, guidefontsize = 9, legendfontsize = 6, ytickfontsize = 9,
              color = :black, label = "optimal", subplot = sub)
        plot!(p, prior_vals, seriestype = :scatter, color = :black, alpha = 0.2, label = "prior", markersize = 4, subplot = sub)
        vline!(p, SITE_BOUNDARIES, color = :grey, label = "", subplot = sub)
        ylims!(p, lims..., subplot = sub)
    end
    row!(1, col(:opt_dsnow), dsnow_prior, "Snow depth (m)", (0, 0.22))
    row!(2, col(:opt_rho_snow), vcat([s.rho for s in SITES]...), "Snow density (kg·m⁻³)", (100, 600))
    row!(3, log10.(col(:opt_axis_ratio)), log10.(vcat([fill(s.ar, length(s.rho)) for s in SITES]...)), "log₁₀(α)", (0, 1.2))
    row!(4, col(:opt_tice), vcat([s.tice for s in SITES]...), "Sea ice temperature (°C)", (-24, -8))
    row!(5, col(:opt_sice), vcat([s.sice for s in SITES]...), "Sea ice salinity", (3.9, 5.5))
    row!(6, col(:opt_dice), vcat([s.dice for s in SITES]...), "Sea ice thickness (m)", (0.75, 1.05))
    xlabel!(p, "Measurement", subplot = 6)
    savefig(p, joinpath(OUT, outname))
end

#----------------------------------------------------------------------------------------
# Fig. 5
#----------------------------------------------------------------------------------------
function tb_panel!(p, subplot, th, tv, ph, pv, oh, ov, title_str)
    n = length(th)
    for i in 1:n
        plot!(p, [i, i], [ph[i], th[i]], color = :grey60, lw = 0.4, label = "", subplot = subplot)
        plot!(p, [i, i], [pv[i], tv[i]], color = :grey60, lw = 0.4, label = "", subplot = subplot)
        plot!(p, [i, i], [th[i] - 5, th[i] + 5], color = :blue, alpha = 0.3, lw = 1.2, label = "", subplot = subplot)
        plot!(p, [i, i], [tv[i] - 5, tv[i] + 5], color = :red, alpha = 0.3, lw = 1.2, label = "", subplot = subplot)
    end
    scatter!(p, 1:n, th, marker = :circle, markercolor = :white, markerstrokecolor = :blue,
             markerstrokewidth = 1.6, markersize = 6.5, label = "TBH - observed (ARIEL)", subplot = subplot)
    scatter!(p, 1:n, tv, marker = :circle, markercolor = :white, markerstrokecolor = :red,
             markerstrokewidth = 1.6, markersize = 6.5, label = "TBV - observed (ARIEL)", subplot = subplot)
    scatter!(p, 1:n, ph, marker = :diamond, markerstrokewidth = 0.6, markerstrokecolor = :black,
             color = :blue, markersize = 5, label = "TBH - prior (xₐ)", subplot = subplot)
    scatter!(p, 1:n, pv, marker = :diamond, markerstrokewidth = 0.6, markerstrokecolor = :black,
             color = :red, markersize = 5, label = "TBV - prior (xₐ)", subplot = subplot)
    scatter!(p, 1:n, oh, marker = :star6, markerstrokewidth = 0.7, markerstrokecolor = :black,
             color = :blue, markersize = 7.5, label = "TBH - optimal", subplot = subplot)
    scatter!(p, 1:n, ov, marker = :star6, markerstrokewidth = 0.7, markerstrokecolor = :black,
             color = :red, markersize = 7.5, label = "TBV - optimal", subplot = subplot)
    vline!(p, SITE_BOUNDARIES, color = :grey, label = "", subplot = subplot)
    ylo, yhi = extrema(vcat(th, tv, ph, pv, oh, ov))
    pad = (yhi - ylo) * 0.12
    for (x, lbl) in zip(SITE_LABEL_X, SITE_LABELS)
        annotate!(p, x, yhi + pad * 1.3, (lbl, 15), subplot = subplot)
    end
    ylims!(p, ylo - pad, yhi + pad * 1.9, subplot = subplot)
    xlabel!(p, "Measurement", subplot = subplot)
    ylabel!(p, "Brightness temperature (K)", subplot = subplot)
    title!(p, title_str, titlefontsize = 18, subplot = subplot)
end

function fig05(data, opt_c, opt_i)
    th = Float64[]; tv = Float64[]; pr = Dict(m => (Float64[], Float64[]) for m in (:coh, :inc))
    for s in SITES
        sd = site_data(data, s)
        append!(th, sd[!, :tbh]); append!(tv, sd[!, :tbv])
        for m in (:coh, :inc)
            h, v = prior_tb(s, sd[!, :dsnow] .* 0.01, m)
            append!(pr[m][1], h); append!(pr[m][2], v)
        end
    end
    col(o, c) = vcat([df[!, c] for df in o]...)
    p = plot(layout = (1, 2), size = (2000, 800), dpi = 500, legend = :outertop, legend_columns = 2,
             legendfontsize = 13, guidefontsize = 17, tickfontsize = 14,
             left_margin = 14Plots.mm, bottom_margin = 12Plots.mm, top_margin = 3Plots.mm, framestyle = :box)
    tb_panel!(p, 1, th, tv, pr[:coh]..., col(opt_c, :opt_tbh), col(opt_c, :opt_tbv), "Coherent model")
    tb_panel!(p, 2, th, tv, pr[:inc]..., col(opt_i, :opt_tbh), col(opt_i, :opt_tbv), "Incoherent model")
    savefig(p, joinpath(OUT, "fig05_tb_prior_optimal.png"))
    return th, tv, pr
end

#----------------------------------------------------------------------------------------
# Figs. 6 and 7
#----------------------------------------------------------------------------------------
function fig06_07(opt_c, opt_i)
    costs     = vcat([o[!, :cost] for o in opt_c]...)
    costs_inc = vcat([o[!, :cost] for o in opt_i]...)
    ylab_pos = 1.08 * maximum(vcat(costs, costs_inc))
    plot(layout = (1, 2), size = (1200, 300), legend = :topright, legend_columns = 2, dpi = 1000,
         bottom_margin = 8Plots.mm, left_margin = 12Plots.mm, right_margin = 1.5Plots.mm, top_margin = 3Plots.mm)
    plot!(costs, markersize = 4, marker = :circle, framestyle = :box, xticks = false, xlabel = "Measurement",
          ylabel = "Cost", guidefontsize = 12, legendfontsize = 8, ytickfontsize = 9,
          color = :orange, label = "Coherent model", subplot = 1)
    plot!(costs_inc, markersize = 4, marker = :circle, framestyle = :box, xticks = false,
          ylabel = "Cost", guidefontsize = 12, legendfontsize = 8, ytickfontsize = 9,
          color = :green, label = "Incoherent model", subplot = 1)
    vline!(SITE_BOUNDARIES, color = :grey, label = "", subplot = 1)
    for (x, lbl) in zip(SITE_LABEL_X, SITE_LABELS)
        annotate!(x, ylab_pos, lbl, 12, subplot = 1)
    end
    scatter!(vcat([o[!, :opt_dsnow] for o in opt_c]...), costs, markersize = 4, marker = :circle, framestyle = :box,
             xlabel = "Snow depth (m)", guidefontsize = 12, legendfontsize = 8, ytickfontsize = 9, xtickfontsize = 9,
             color = :orange, label = "Coherent model", subplot = 2)
    scatter!(vcat([o[!, :opt_dsnow] for o in opt_i]...), costs_inc, markersize = 4, marker = :circle, framestyle = :box,
             guidefontsize = 12, legendfontsize = 8, ytickfontsize = 9, xtickfontsize = 9,
             color = :green, label = "Incoherent model", subplot = 2)
    savefig(joinpath(OUT, "fig06_costs.png"))

    plot(layout = (1, 2), size = (1200, 300), dpi = 500,
         bottom_margin = 8Plots.mm, left_margin = 12Plots.mm, right_margin = 1.5Plots.mm, top_margin = 3Plots.mm)
    x_min, x_max = extrema(vcat(costs, costs_inc))
    for (sub, c, colr, ttl) in ((1, costs, :orange, "Coherent model"), (2, costs_inc, :green, "Incoherent model"))
        histogram!(c, bins = 15, color = colr, alpha = 0.7, framestyle = :box, xlabel = "Cost", ylabel = "Frequency",
                   guidefontsize = 12, title = ttl, titlefontsize = 12, label = "", subplot = sub)
        xlims!(x_min, x_max, subplot = sub)
        ylims!(0, 18, subplot = sub)
    end
    savefig(joinpath(OUT, "fig07_cost_histograms.png"))
end

#----------------------------------------------------------------------------------------
# Fig. 8
#----------------------------------------------------------------------------------------
function fig08(data, opt_c, opt_i)
    dsnows = 0:0.001:0.2
    panels = []
    for (k, s) in enumerate(SITES)
        opt, opt_inc, obs = opt_c[k], opt_i[k], site_data(data, s)
        p = plot(xlabel = "Snow depth (m)", ylabel = "Brightness Temperature (K)", margin = 5Plots.mm,
                 framestyle = :box, guidefontsize = 12, title = s.name)
        for i in 1:nrow(opt)
            tb = [wilheit_model_integrated(40, [d], s.tsnow, [opt[i, :opt_rho_snow]], [opt[i, :opt_dice]],
                                           [opt[i, :opt_tice]], [opt[i, :opt_sice]], "mix", 0, opt[i, :opt_axis_ratio]) for d in dsnows]
            plot!(p, dsnows, getindex.(tb, 1), color = :blue, alpha = 0.5, lw = 0.4, label = "")
            plot!(p, dsnows, getindex.(tb, 2), color = :red, alpha = 0.5, lw = 0.4, label = "")
        end
        for (o, mk, ms, sw) in ((opt, :star5, 5, 0.5), (opt_inc, :xcross, 3, 1))
            scatter!(p, o[!, :opt_dsnow], o[!, :opt_tbh], marker = mk, markersize = ms, markerstrokewidth = sw, color = :blue, alpha = 0.5, label = "")
            scatter!(p, o[!, :opt_dsnow], o[!, :opt_tbv], marker = mk, markersize = ms, markerstrokewidth = sw, color = :red, alpha = 0.5, label = "")
            for i in 1:nrow(o)
                plot!(p, [o[i, :opt_dsnow]; obs[i, :dsnow] * 0.01], [o[i, :opt_tbh]; obs[i, :tbh]], color = :black, linestyle = :dot, lw = 0.5, label = "")
                plot!(p, [o[i, :opt_dsnow]; obs[i, :dsnow] * 0.01], [o[i, :opt_tbv]; obs[i, :tbv]], color = :black, linestyle = :dot, lw = 0.5, label = "")
            end
        end
        scatter!(p, obs[!, :dsnow] * 0.01, obs[!, :tbh], marker = :circle, markersize = 3, markerstrokewidth = 0.5, color = :blue, label = "")
        scatter!(p, obs[!, :dsnow] * 0.01, obs[!, :tbv], marker = :circle, markersize = 3, markerstrokewidth = 0.5, color = :red, label = "")
        push!(panels, p)
    end
    p_legend = plot(legend = :top, legend_columns = 4, legendfontsize = 10, framestyle = :none, grid = false,
                    xticks = false, yticks = false, bottom_margin = -10Plots.mm)
    plot!(p_legend, [], [], color = :blue, lw = 0.4, alpha = 0.5, label = "TBH coh model")
    plot!(p_legend, [], [], color = :red, lw = 0.4, alpha = 0.5, label = "TBV coh model")
    scatter!(p_legend, [], [], marker = :star5, color = :blue, markerstrokewidth = 0.5, alpha = 0.5, label = "opt - TBH coh model")
    scatter!(p_legend, [], [], marker = :star5, color = :red, markerstrokewidth = 0.5, alpha = 0.5, label = "opt - TBV coh model")
    scatter!(p_legend, [], [], marker = :circle, color = :blue, markerstrokewidth = 0.5, label = "obs - TBH measured")
    scatter!(p_legend, [], [], marker = :circle, color = :red, markerstrokewidth = 0.5, label = "obs - TBV measured")
    scatter!(p_legend, [], [], marker = :xcross, color = :blue, alpha = 0.5, markerstrokewidth = 1, label = "opt - TBH inc model")
    scatter!(p_legend, [], [], marker = :xcross, color = :red, alpha = 0.5, markerstrokewidth = 1, label = "opt - TBV inc model")
    l = @layout [a{0.1h}; grid(2, 2)]
    savefig(plot(p_legend, panels..., layout = l, size = (1400, 950), dpi = 500), joinpath(OUT, "fig08_tb_vs_snow_depth_sites.png"))
end

#----------------------------------------------------------------------------------------
# PD at the prior and optimal states
#----------------------------------------------------------------------------------------
function fig_pd(th, tv, pr, opt_c, opt_i)
    pd_meas = tv .- th
    site_of = vcat([fill(k, length(s.rho)) for (k, s) in enumerate(SITES)]...)
    col(o, c) = vcat([df[!, c] for df in o]...)
    pd = Dict(:coh => (pr[:coh][2] .- pr[:coh][1], col(opt_c, :opt_tbv) .- col(opt_c, :opt_tbh)),
              :inc => (pr[:inc][2] .- pr[:inc][1], col(opt_i, :opt_tbv) .- col(opt_i, :opt_tbh)))
    all_pd = vcat(pd_meas, pd[:coh]..., pd[:inc]...)
    lims = (floor(minimum(all_pd)) - 2, ceil(maximum(all_pd)) + 2); span = lims[2] - lims[1]
    p = plot(layout = (1, 2), size = (1100, 520), dpi = 500, framestyle = :box,
             legend = :topleft, legendfontsize = 8, guidefontsize = 11, tickfontsize = 10,
             bottom_margin = 8Plots.mm, left_margin = 10Plots.mm, right_margin = 4Plots.mm, top_margin = 4Plots.mm)
    for (sub, m, ttl) in ((1, :coh, "Coherent model"), (2, :inc, "Incoherent model"))
        prior, opt = pd[m]
        plot!(p, collect(lims), collect(lims); color = :black, linestyle = :dash, lw = 1, label = "", subplot = sub)
        for i in eachindex(pd_meas)
            plot!(p, [pd_meas[i], pd_meas[i]], [prior[i], opt[i]]; color = :grey60, lw = 0.5, label = "", subplot = sub)
        end
        for (k, s) in enumerate(SITES)
            sel = site_of .== k
            scatter!(p, pd_meas[sel], prior[sel]; marker = :diamond, markercolor = :white, markerstrokecolor = SITE_COLORS[k],
                     markerstrokewidth = 1.4, markersize = 6, label = "", subplot = sub)
            scatter!(p, pd_meas[sel], opt[sel]; marker = :circle, color = SITE_COLORS[k], markerstrokewidth = 0,
                     markersize = 6, label = s.name, subplot = sub)
        end
        scatter!(p, [NaN], [NaN]; marker = :diamond, markercolor = :white, markerstrokecolor = :black,
                 markerstrokewidth = 1.4, markersize = 6, label = "prior (xₐ)", subplot = sub)
        scatter!(p, [NaN], [NaN]; marker = :circle, color = :black, markerstrokewidth = 0, markersize = 6,
                 label = "optimal (x̂)", subplot = sub)
        for (k, (lab, v)) in enumerate((("optimal (x̂)", opt), ("prior (xₐ)", prior)))
            bias, rmse, r = stats(v, pd_meas)
            annotate!(p, lims[2] - 0.03span, lims[1] + (0.04 + 0.055 * (k - 1)) * span,
                      text("$lab: bias = $(round(bias, digits = 1)) K, RMSE = $(round(rmse, digits = 1)) K, r = $(round(r, digits = 2))",
                           8, :right, :bottom), subplot = sub)
        end
        xlims!(p, lims..., subplot = sub); ylims!(p, lims..., subplot = sub)
        xlabel!(p, "Measured PD (K)", subplot = sub); ylabel!(p, "Modelled PD (K)", subplot = sub)
        title!(p, ttl, titlefontsize = 12, subplot = sub)
    end
    savefig(p, joinpath(OUT, "pd_prior_optimal.png"))
end

#----------------------------------------------------------------------------------------
# Fig. S1
#----------------------------------------------------------------------------------------
function fig_s1()
    load(m) = [CSV.read(joinpath(@__DIR__, "opt_fixed_alpha", "opt_res_$(s.key)_$(m).csv"), DataFrame) for s in SITES]
    res = Dict("coh" => load("coh"), "inc" => load("inc"))
    all_tb = vcat([vcat(df.target_tbh, df.target_tbv, df.opt_tbh, df.opt_tbv) for m in ("coh", "inc") for df in res[m]]...)
    lims = (floor(minimum(all_tb)) - 3, ceil(maximum(all_tb)) + 3); span = lims[2] - lims[1]
    p = plot(layout = (2, 2), size = (950, 900), dpi = 500, framestyle = :box,
             legendfontsize = 8, guidefontsize = 10, tickfontsize = 9, titlefontsize = 10,
             bottom_margin = 6Plots.mm, left_margin = 10Plots.mm, right_margin = 4Plots.mm, top_margin = 3Plots.mm)
    for (row, (m, mlab)) in enumerate((("coh", "Coherent"), ("inc", "Incoherent")))
        total_cost = sum(sum(df.cost) for df in res[m])
        for (c, (pol, plab)) in enumerate(((:h, "TBH"), (:v, "TBV")))
            sub = 2 * (row - 1) + c
            obs_c, sim_c = Symbol("target_tb$(pol)"), Symbol("opt_tb$(pol)")
            plot!(p, collect(lims), collect(lims); color = :black, linestyle = :dash, lw = 1, label = "", subplot = sub)
            for (k, lab) in enumerate(SITE_LABELS)
                scatter!(p, res[m][k][!, obs_c], res[m][k][!, sim_c]; color = SITE_COLORS[k], markerstrokewidth = 0,
                         markersize = 6, label = lab, subplot = sub)
            end
            bias, rmse, r = stats(vcat([df[!, sim_c] for df in res[m]]...), vcat([df[!, obs_c] for df in res[m]]...))
            annotate!(p, lims[2] - 0.03span, lims[1] + 0.04span,
                      text("bias = $(round(bias, digits = 1)) K, RMSE = $(round(rmse, digits = 1)) K, r = $(round(r, digits = 2))", 8, :right, :bottom),
                      subplot = sub)
            annotate!(p, lims[2] - 0.03span, lims[1] + 0.10span,
                      text("total cost ($mlab, α fixed) = $(round(total_cost, digits = 1))", 8, :right, :bottom), subplot = sub)
            xlims!(p, lims..., subplot = sub); ylims!(p, lims..., subplot = sub)
            xlabel!(p, "Observed $plab (K)", subplot = sub); ylabel!(p, "Simulated $plab, α fixed (K)", subplot = sub)
            title!(p, "$mlab model, ARIEL footprint — $plab", subplot = sub)
            plot!(p, legend = :topleft, subplot = sub)
        end
    end
    savefig(p, joinpath(OUT, "figS1_fixed_alpha.png"))
end

#----------------------------------------------------------------------------------------
function main()
    mkpath(OUT)
    data = CSV.read(joinpath(@__DIR__, "chars_data_filtered.csv"), DataFrame)
    opt_c = [read_opt("opt_coh", s) for s in SITES]
    opt_i = [read_opt("opt_inc", s) for s in SITES]
    fig03();                                  println("Fig. 3 done")
    fig04(data, "opt_coh", "fig04_geophysical_coh.png")
    fig04(data, "opt_inc", "fig04_geophysical_inc.png"); println("Fig. 4 done")
    th, tv, pr = fig05(data, opt_c, opt_i);   println("Fig. 5 done")
    fig06_07(opt_c, opt_i);                   println("Figs. 6-7 done")
    fig08(data, opt_c, opt_i);                println("Fig. 8 done")
    fig_pd(th, tv, pr, opt_c, opt_i);         println("PD figure done")
    fig_s1();                                 println("Fig. S1 done")
end
#----------------------------------------------------------------------------------------

main()
