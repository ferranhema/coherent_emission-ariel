# coherence_damping.jl
#
# Quantifies how the coherent (interference) oscillations of the antenna-
# integrated TBs (wilheit_model_integrated) are damped as snow depth increases.
#
# For fixed geophysical conditions, TBH, TBV and PD = TBV - TBH are computed on
# a fine snow-depth grid with
#   coherent, ARIEL footprint   -> wilheit_model_integrated
#   coherent, single angle 40°  -> wilheit_model
#   incoherent, ARIEL footprint -> incoherent_model_integrated (same permittivities,
#                                  temperatures and antenna weights; no phase)
#
# Envelope of the footprint oscillations: the residual R(d) = coherent - incoherent
# (both antenna-integrated) is evaluated at its local extrema, A_k = |R(d_k)|.
# This gives the amplitude from the first extremum (~λ/4 in snow) onwards. It
# assumes the incoherent TB is the centre of the oscillation, which is checked
# in two ways: (i) mean(R) over 0.5-1.0 m, where the oscillations have died out,
# and (ii) a coherent-only estimate at the same extrema,
#     A_k = | y_k - (y_{k-1} + y_{k+1}) / 2 | / 2,
# which needs no reference but is only defined from the second extremum on.
#
# Decay laws fitted to ln A over the main decay (extrema with A >= 1 K):
#   exponential  A = A0 exp(-d / L)
#   Gaussian     A = A0 exp(-d^2 / (2 L^2))
# For the Gaussian case an analytic estimate is also given: averaging exp(iφ)
# over a Gaussian beam of width σθ, with φ = 4π d sqrt(ε' - sin²θ) / λ, damps the
# interference term by exp(-σφ²/2), σφ = |dφ/dθ| σθ, so
#     L = λ sqrt(ε' - sin²θ0) / (4π sinθ0 cosθ0 σθ)
# (linearised around θ0 = 40°, σθ = mean of the H and V beam widths).
#
# Configurations: the reference case of Fig. 3 and the site-mean priors of
# CLO, MID-1, MID-2, FAR (as in coh-model_optimization.jl).
#
# Products:
#   figures/coherence_damping.png             -- TBs and envelopes
#   figures/coherence_damping_envelope.csv    -- every envelope point (config, model, channel, d, A)
#   console                             -- amplitudes at the extrema, fitted damping lengths,
#                                          depth where A < 1 K, checks

using Plots
using DataFrames
using CSV

include("coherent_model.jl")
include("incoherent_model.jl")

const BASEDIR = @__DIR__

const DSNOW = collect(0.0:0.001:1.0)    # m
const A_THRESHOLD = 1.0                 # K
const D_LIST_MAX = 0.35                 # m, extrema listed in the console table up to this depth
const SIGMA_H, SIGMA_V = 15.29, 14.87   # deg, ARIEL beam (as in antenna_pattern_integration)

# ─────────────────────────────────────────────────────────────────────────────
# Configurations (scalars; site values are the means of the priors in
# coh-model_optimization.jl minimization())
# ─────────────────────────────────────────────────────────────────────────────

mean_(x) = sum(x) / length(x)

const CONFIGS = [
    (name = "reference", tsnow = -15.0, rho = 300.0, dice = 0.9, tice = -10.0, sice = 5.0),
    (name = "CLO",   tsnow = -14.0,
     rho  = mean_([355.0, 300, 300, 355, 300, 267, 300, 300, 300, 350]),
     dice = mean_([0.945, 0.95, 0.94, 0.94, 0.97, 0.985, 0.98, 0.895, 0.895, 0.99]),
     tice = -13.0, sice = 5.32),
    (name = "MID-1", tsnow = -15.0,
     rho  = mean_([350.0, 365, 367, 350, 350, 350]),
     dice = mean_([0.89, 0.855, 0.855, 0.89, 0.855, 0.855]),
     tice = -11.0, sice = 4.6),
    (name = "MID-2", tsnow = -25.0,
     rho  = mean_([350.0, 350, 300, 347]),
     dice = 0.93,
     tice = -20.0, sice = 4.5),
    (name = "FAR",   tsnow = -25.0,
     rho  = mean_([385.0, 385, 385, 385, 365, 400, 400, 400, 385, 385, 365, 400, 400, 400, 400]),
     dice = mean_([0.86, 0.86, 0.86, 0.86, 0.92, 0.84, 0.84, 0.845, 0.855, 0.855, 0.855, 0.855, 0.855, 0.855, 0.855]),
     tice = -17.0, sice = 4.8),
]
const PERM, MR, AR = "mix", 0, 5.0
const CONFIG_COLORS = [:black, :darkgreen, :mediumblue, :darkorchid4, :brown4]
const CHANNELS = (:H, :V, :PD)

# ─────────────────────────────────────────────────────────────────────────────
# Forward-model sweeps over snow depth
# ─────────────────────────────────────────────────────────────────────────────

function sweep(cfg, fwd)
    tbh = Float64[]; tbv = Float64[]
    for d in DSNOW
        h, v = fwd(d, cfg)
        push!(tbh, h); push!(tbv, v)
    end
    return (H = tbh, V = tbv, PD = tbv .- tbh)
end

coh_int(d, c) = wilheit_model_integrated(40, [d], [c.tsnow], [c.rho], [c.dice], [c.tice], [c.sice], PERM, MR, AR)
coh_40(d, c)  = wilheit_model(40, [d], [c.tsnow], [c.rho], [c.dice], [c.tice], [c.sice], PERM, MR, AR)
inc_int(d, c) = incoherent_model_integrated(40, d, c.tsnow, c.rho, c.dice, c.tice, c.sice, PERM, MR, AR)

# ─────────────────────────────────────────────────────────────────────────────
# Envelopes, fits, analytic estimate
# ─────────────────────────────────────────────────────────────────────────────

local_extrema(y) = [i for i in 2:length(y)-1 if (y[i] - y[i-1]) * (y[i+1] - y[i]) < 0]

# envelope from the residual against a reference (coherent - incoherent)
function envelope_residual(d, y, yref)
    r = y .- yref
    ex = local_extrema(r)
    return d[ex], abs.(r[ex])
end

# reference-free envelope (cancels a linear trend), from the 2nd extremum on
function envelope_3pt(d, y)
    ex = local_extrema(y)
    dk = Float64[]; Ak = Float64[]
    for k in 2:length(ex)-1
        push!(dk, d[ex[k]])
        push!(Ak, abs(y[ex[k]] - (y[ex[k-1]] + y[ex[k+1]]) / 2) / 2)
    end
    return dk, Ak
end

# least squares ln A = a + b * x, returns (a, b, R^2)
function linfit(x, lnA)
    n = length(x)
    mx = sum(x) / n; my = sum(lnA) / n
    b = sum((x .- mx) .* (lnA .- my)) / sum((x .- mx) .^ 2)
    a = my - b * mx
    ss_res = sum((lnA .- (a .+ b .* x)) .^ 2)
    ss_tot = sum((lnA .- my) .^ 2)
    return a, b, 1 - ss_res / ss_tot
end

# fit over the main decay: from the first extremum up to the last one before A first drops below the threshold
function fit_main_decay(dk, Ak)
    iend = findfirst(<(A_THRESHOLD), Ak)
    iend = iend === nothing ? length(Ak) : iend - 1
    iend < 3 && return nothing
    x = dk[1:iend]; lnA = log.(Ak[1:iend])
    a_e, b_e, r2_e = linfit(x, lnA)
    a_g, b_g, r2_g = linfit(x .^ 2, lnA)
    return (n = iend, dmax = x[end],
            A0_exp = exp(a_e), L_exp = -1 / b_e, R2_exp = r2_e,
            A0_gauss = exp(a_g), L_gauss = b_g < 0 ? sqrt(-1 / (2 * b_g)) : NaN, R2_gauss = r2_g)
end

function depth_below(dk, Ak, thr)
    for i in 2:length(Ak)
        if Ak[i-1] >= thr && Ak[i] < thr
            w = (log(Ak[i-1]) - log(thr)) / (log(Ak[i-1]) - log(Ak[i]))
            return dk[i-1] + w * (dk[i] - dk[i-1])
        end
    end
    return NaN
end

snow_n(c) = sqrt(epsilon_snow(1.4e9, c.rho, c.tsnow)[1])
wavelength() = 3e8 / 1.4e9

function theoretical_period(c)
    n = snow_n(c)
    return wavelength() / (2 * sqrt(n^2 - sin(deg2rad(40))^2))
end

function analytic_L_gauss(c)
    ε = snow_n(c)^2
    θ0 = deg2rad(40); σθ = deg2rad((SIGMA_H + SIGMA_V) / 2)
    return wavelength() * sqrt(ε - sin(θ0)^2) / (4π * sin(θ0) * cos(θ0) * σθ)
end

# ─────────────────────────────────────────────────────────────────────────────
# Run
# ─────────────────────────────────────────────────────────────────────────────

fmt(x; digits = 2) = (x === nothing || isnan(x)) ? "-" : string(round(x, digits = digits))

function analyse()
    env_rows = DataFrame(config = String[], model = String[], channel = String[], dsnow = Float64[], amplitude = Float64[])
    res = Dict{Tuple{String,Symbol}, Any}()
    curves = Dict{Tuple{String,String}, Any}()

    for cfg in CONFIGS
        println("Sweeping: $(cfg.name)")
        ci = sweep(cfg, coh_int); c40 = sweep(cfg, coh_40); ii = sweep(cfg, inc_int)
        curves[(cfg.name, "integrated")] = ci
        curves[(cfg.name, "40deg")]      = c40
        curves[(cfg.name, "incoherent")] = ii

        tail = DSNOW .>= 0.5
        for ch in CHANNELS
            dk, Ak     = envelope_residual(DSNOW, ci[ch], ii[ch])
            dk3, Ak3   = envelope_3pt(DSNOW, ci[ch])
            dk40, Ak40 = envelope_3pt(DSNOW, c40[ch])
            r = ci[ch] .- ii[ch]
            res[(cfg.name, ch)] = (dk = dk, Ak = Ak, dk3 = dk3, Ak3 = Ak3, dk40 = dk40, Ak40 = Ak40,
                                   fit = fit_main_decay(dk, Ak),
                                   tail_mean = mean_(r[tail]), tail_max = maximum(abs.(r[tail])))
            for (m, dd, AA) in (("integrated_residual", dk, Ak), ("integrated_3pt", dk3, Ak3), ("40deg_3pt", dk40, Ak40))
                for (d, A) in zip(dd, AA)
                    push!(env_rows, (cfg.name, m, String(ch), d, A))
                end
            end
        end
    end
    return env_rows, res, curves
end

function print_summary(res)
    println("\n── Oscillation period in snow depth at 40° (m) ────────────────────────")
    for cfg in CONFIGS
        r = res[(cfg.name, :H)]
        sim = length(r.dk40) >= 2 ? 2 * (r.dk40[end] - r.dk40[1]) / (length(r.dk40) - 1) : NaN
        println(rpad(cfg.name, 11), "theory λ/(2√(ε'-sin²θ)) = ", fmt(theoretical_period(cfg), digits = 3),
                "   from 40° TBH extrema = ", fmt(sim, digits = 3))
    end

    println("\n── Footprint oscillation amplitude |coh - inc| at the extrema, d ≤ $(D_LIST_MAX) m  (depth m : A K) ──")
    for cfg in CONFIGS, ch in CHANNELS
        r = res[(cfg.name, ch)]
        sel = r.dk .<= D_LIST_MAX
        println(rpad(cfg.name, 11), rpad(String(ch), 4),
                join(["$(fmt(d, digits = 3)):$(fmt(A, digits = 1))" for (d, A) in zip(r.dk[sel], r.Ak[sel])], "  "))
    end

    println("\n── Single angle 40°: amplitude over 0-1 m (no antenna averaging) ──────")
    for cfg in CONFIGS, ch in CHANNELS
        r = res[(cfg.name, ch)]
        println(rpad(cfg.name, 11), rpad(String(ch), 4), "A = ", fmt(minimum(r.Ak40)), " to ", fmt(maximum(r.Ak40)), " K")
    end

    println("\n── Main-decay fits (extrema with A ≥ $(A_THRESHOLD) K) ──────────────────────────")
    println(rpad("config", 11), rpad("ch", 4), rpad("n", 3), rpad("to d", 7), rpad("L_exp", 8), rpad("R²", 6),
            rpad("L_gauss", 9), rpad("R²", 6), rpad("L_gauss analytic", 18), "d(A<$(A_THRESHOLD) K)")
    for cfg in CONFIGS, ch in CHANNELS
        r = res[(cfg.name, ch)]; f = r.fit
        println(rpad(cfg.name, 11), rpad(String(ch), 4),
                f === nothing ? rpad("-", 39) :
                    rpad(string(f.n), 3) * rpad(fmt(f.dmax, digits = 3), 7) *
                    rpad(fmt(f.L_exp, digits = 3), 8) * rpad(fmt(f.R2_exp), 6) *
                    rpad(fmt(f.L_gauss, digits = 3), 9) * rpad(fmt(f.R2_gauss), 6),
                rpad(fmt(analytic_L_gauss(cfg), digits = 3), 18),
                fmt(depth_below(r.dk, r.Ak, A_THRESHOLD), digits = 3))
    end

    println("\n── Checks ──────────────────────────────────────────────────────────────")
    println(rpad("config", 11), rpad("ch", 4), rpad("mean(coh-inc) 0.5-1 m", 23), rpad("max|coh-inc| 0.5-1 m", 22),
            "max |A_residual - A_3pt| for d ≤ $(D_LIST_MAX) m")
    for cfg in CONFIGS, ch in CHANNELS
        r = res[(cfg.name, ch)]
        diffs = Float64[]
        for (d, A3) in zip(r.dk3, r.Ak3)
            d > D_LIST_MAX && continue
            i = findfirst(x -> abs(x - d) <= 0.003, r.dk)
            i === nothing || push!(diffs, abs(r.Ak[i] - A3))
        end
        println(rpad(cfg.name, 11), rpad(String(ch), 4), rpad(fmt(r.tail_mean), 23), rpad(fmt(r.tail_max), 22),
                isempty(diffs) ? "-" : fmt(maximum(diffs)))
    end
    println("(A in K, depths and lengths in m)")
end

# ─────────────────────────────────────────────────────────────────────────────
# Figure
# ─────────────────────────────────────────────────────────────────────────────

function figure_damping(res, curves, dmax_measured)
    p = plot(layout = (2, 2), size = (1200, 950), dpi = 500, framestyle = :box,
             legendfontsize = 8, guidefontsize = 11, tickfontsize = 9, titlefontsize = 11,
             bottom_margin = 6Plots.mm, left_margin = 10Plots.mm, right_margin = 4Plots.mm, top_margin = 3Plots.mm)

    # (1) reference TBs: 40°, footprint coherent, footprint incoherent
    ref = "reference"
    vspan!(p, [0, dmax_measured]; color = :grey85, alpha = 0.5, label = "ARIEL snow depths", subplot = 1)
    for (ch, col) in ((:V, :red), (:H, :blue))
        plot!(p, DSNOW, curves[(ref, "40deg")][ch]; color = col, lw = 0.8, alpha = 0.45,
              label = "TB$(ch) - coherent, 40°", subplot = 1)
        plot!(p, DSNOW, curves[(ref, "integrated")][ch]; color = col, lw = 1.8,
              label = "TB$(ch) - coherent, ARIEL footprint", subplot = 1)
        plot!(p, DSNOW, curves[(ref, "incoherent")][ch]; color = col, lw = 1.2, linestyle = :dot,
              label = "TB$(ch) - incoherent, ARIEL footprint", subplot = 1)
    end
    xlabel!(p, "Snow depth (m)", subplot = 1); ylabel!(p, "Brightness temperature (K)", subplot = 1)
    title!(p, "Reference case", subplot = 1)
    plot!(p, legend = :bottomright, subplot = 1)

    # (2) reference envelopes: footprint (filled) vs 40° (open), Gaussian fit and analytic estimate
    vspan!(p, [0, dmax_measured]; color = :grey85, alpha = 0.5, label = "", subplot = 2)
    cfg_ref = CONFIGS[1]
    Lan = analytic_L_gauss(cfg_ref)
    dd = range(0, 0.6, length = 300)
    for (ch, col) in ((:H, :blue), (:V, :red), (:PD, :black))
        r = res[(ref, ch)]
        scatter!(p, r.dk40, r.Ak40; marker = :circle, markercolor = :white, markerstrokecolor = col,
                 markerstrokewidth = 1, markersize = 3.5, label = "$ch - 40°", subplot = 2)
        scatter!(p, r.dk, r.Ak; marker = :circle, color = col, markerstrokewidth = 0, markersize = 4.5,
                 label = "$ch - ARIEL footprint", subplot = 2)
        f = r.fit
        if f !== nothing && !isnan(f.L_gauss)
            plot!(p, dd, f.A0_gauss .* exp.(-dd .^ 2 ./ (2 * f.L_gauss^2)); color = col, lw = 1, linestyle = :dash,
                  label = "", subplot = 2)
            plot!(p, dd, f.A0_gauss .* exp.(-dd .^ 2 ./ (2 * Lan^2)); color = col, lw = 1, linestyle = :dot,
                  label = "", subplot = 2)
        end
    end
    hline!(p, [A_THRESHOLD]; color = :grey40, lw = 0.8, linestyle = :dot, label = "", subplot = 2)
    plot!(p, yscale = :log10, ylims = (1e-2, 60), legend = :topright, legend_columns = 2, subplot = 2)
    xlabel!(p, "Snow depth (m)", subplot = 2); ylabel!(p, "Oscillation amplitude (K)", subplot = 2)
    title!(p, "Reference: dashed = Gaussian fit, dotted = analytic L", subplot = 2)

    # (3) and (4): footprint envelopes of TBH and PD for all configurations
    for (sub, ch, lab) in ((3, :H, "TBH"), (4, :PD, "PD"))
        vspan!(p, [0, dmax_measured]; color = :grey85, alpha = 0.5, label = "", subplot = sub)
        for (k, cfg) in enumerate(CONFIGS)
            r = res[(cfg.name, ch)]
            plot!(p, r.dk, r.Ak; color = CONFIG_COLORS[k], lw = 1.2, marker = :circle, markersize = 3,
                  markerstrokewidth = 0, label = cfg.name, subplot = sub)
        end
        hline!(p, [A_THRESHOLD]; color = :grey40, lw = 0.8, linestyle = :dot, label = "", subplot = sub)
        plot!(p, yscale = :log10, ylims = (1e-2, 60), legend = :topright, subplot = sub)
        xlabel!(p, "Snow depth (m)", subplot = sub); ylabel!(p, "Oscillation amplitude (K)", subplot = sub)
        title!(p, "$lab envelope |coh - inc|, ARIEL footprint", subplot = sub)
    end

    outpath = joinpath(BASEDIR, "figures", "coherence_damping.png")
    savefig(p, outpath)
    println("Saved -> $outpath")
    return p
end

# ─────────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────────

function main_damping()
    data = CSV.read(joinpath(BASEDIR, "chars_data_filtered.csv"), DataFrame)
    dmax_measured = maximum(data[!, :dsnow]) * 0.01

    env_rows, res, curves = analyse()

    mkpath(joinpath(BASEDIR, "figures"))
    csvpath = joinpath(BASEDIR, "figures", "coherence_damping_envelope.csv")
    CSV.write(csvpath, env_rows)
    println("Saved -> $csvpath")

    print_summary(res)
    figure_damping(res, curves, dmax_measured)
end

main_damping()
