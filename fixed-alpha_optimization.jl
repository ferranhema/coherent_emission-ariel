using DataFrames
using CSV
using Optim
include("coherent_model.jl")
include("incoherent_model.jl")
#----------------------------------------------------------------------------------------
# Sensitivity test: optimal estimation with the brine-inclusion axis ratio fixed at
# its prior value (α = 5) for both models, i.e. not retrieved. Same priors, σ,
# σ_TB = 5 K, first guess and Nelder-Mead settings as coh-model_optimization.jl and
# inc-model_optimization.jl; the cost has no α term. Both forward models are
# integrated over the ARIEL antenna pattern.
#
# State vector: [dsnow, rho_snow, dice, tice, sice]
# Output: ./opt_fixed_alpha/opt_res_<site>_{coh,inc}.csv
#----------------------------------------------------------------------------------------

const AR_FIXED = 5.0
const SIGMA_TB = 5.0    # TB measurement uncertainty (K)

# Site priors and uncertainties (identical to coh/inc-model_optimization.jl)
const SITES = [
    (key = "clo", temp = -9.69, tsnow = -14.0,
     rho   = [355.0, 300, 300, 355, 300, 267, 300, 300, 300, 350],
     dice  = [0.945, 0.95, 0.94, 0.94, 0.97, 0.985, 0.98, 0.895, 0.895, 0.99],
     tice  = -13.0, sice = 5.32,
     sig_ds  = [0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01, 0.02, 0.02, 0.01],
     sig_rho = [25.0, 50, 50, 25, 50, 25, 50, 50, 50, 35],
     sig_di  = [0.01, 0.01, 0.05, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05]),
    (key = "mid1", temp = -7.37, tsnow = -15.0,
     rho   = [350.0, 365, 367, 350, 350, 350],
     dice  = [0.89, 0.855, 0.855, 0.89, 0.855, 0.855],
     tice  = -11.0, sice = 4.6,
     sig_ds  = fill(0.01, 6),
     sig_rho = [100.0, 50, 50, 100, 100, 100],
     sig_di  = fill(0.01, 6)),
    (key = "mid2", temp = -12.95, tsnow = -25.0,
     rho   = [350.0, 350, 300, 347],
     dice  = [0.93, 0.93, 0.93, 0.93],
     tice  = -20.0, sice = 4.5,
     sig_ds  = fill(0.02, 4),
     sig_rho = [100.0, 100, 50, 50],
     sig_di  = fill(0.05, 4)),
    (key = "far", temp = -13.86, tsnow = -25.0,
     rho   = [385.0, 385, 385, 385, 365, 400, 400, 400, 385, 385, 365, 400, 400, 400, 400],
     dice  = [0.86, 0.86, 0.86, 0.86, 0.92, 0.84, 0.84, 0.845, 0.855, 0.855, 0.855, 0.855, 0.855, 0.855, 0.855],
     tice  = -17.0, sice = 4.8,
     sig_ds  = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.01, 0.005, 0.02, 0.005],
     sig_rho = [50.0, 50, 50, 50, 50, 100, 100, 100, 50, 50, 50, 100, 100, 100, 100],
     sig_di  = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]),
]
const SIG_TICE, SIG_SICE = 2.0, 1.0

fwd_coh(ds, rho, di, ti, si, s) = wilheit_model_integrated(40, [ds], [s.tsnow], [rho], [di], [ti], [si], "mix", 0, AR_FIXED)
fwd_inc(ds, rho, di, ti, si, s) = incoherent_model_integrated(40, ds, s.tsnow, rho, di, ti, si, "mix", 0, AR_FIXED)

#----------------------------------------------------------------------------------------
"""
    optimize_site(fwd, s, site_data)

Nelder-Mead optimization of [dsnow, rho_snow, dice, tice, sice] for every
measurement of one site, with the axis ratio fixed at `AR_FIXED`.
"""
function optimize_site(fwd, s, site_data)
    out = DataFrame(opt_dsnow = Float64[], opt_rho_snow = Float64[], opt_dice = Float64[], opt_tice = Float64[],
                    opt_sice = Float64[], opt_tbh = Float64[], opt_tbv = Float64[], target_tbh = Float64[],
                    target_tbv = Float64[], cost = Float64[], iterations = Int[], converged = Bool[])
    for i in 1:nrow(site_data)
        ds_a = site_data[i, :dsnow] * 0.01
        tbh_t, tbv_t = site_data[i, :tbh], site_data[i, :tbv]
        cost = x -> begin
            tbh, tbv = fwd(x[1], x[2], x[3], x[4], x[5], s)
            ((tbh - tbh_t)^2 + (tbv - tbv_t)^2) / SIGMA_TB^2 +
            (x[1] - ds_a)^2 / s.sig_ds[i]^2 + (x[2] - s.rho[i])^2 / s.sig_rho[i]^2 +
            (x[3] - s.dice[i])^2 / s.sig_di[i]^2 + (x[4] - s.tice)^2 / SIG_TICE^2 + (x[5] - s.sice)^2 / SIG_SICE^2
        end
        res = optimize(cost, [ds_a, s.rho[i], s.dice[i], s.tice, s.sice], NelderMead())
        m = res.minimizer
        tbh, tbv = fwd(m..., s)
        push!(out, (m..., tbh, tbv, tbh_t, tbv_t, res.minimum, Optim.iterations(res), Optim.converged(res)))
        println("Optimization finished, nº: $i")
    end
    return out
end

#----------------------------------------------------------------------------------------
function main()
    data = CSV.read(joinpath(@__DIR__, "chars_data_filtered.csv"), DataFrame)
    mkpath(joinpath(@__DIR__, "opt_fixed_alpha"))
    for s in SITES, (lab, fwd) in (("coh", fwd_coh), ("inc", fwd_inc))
        site_data = filter(row -> row[:temp] == s.temp, data)
        out = optimize_site(fwd, s, site_data)
        CSV.write(joinpath(@__DIR__, "opt_fixed_alpha", "opt_res_$(s.key)_$(lab).csv"), out)
        println("$(s.key) $(lab): total cost = $(round(sum(out.cost), digits = 2)), all converged = $(all(out.converged))")
    end
end
#----------------------------------------------------------------------------------------

main()
