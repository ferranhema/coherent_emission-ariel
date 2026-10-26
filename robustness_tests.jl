using DataFrames
using CSV
using Optim
using Random
using LinearAlgebra
using Printf
include("coherent_model.jl")
include("incoherent_model.jl")
#----------------------------------------------------------------------------------------
# Robustness tests of the optimization, as reported in the paper. Both forward models are
# integrated over the ARIEL antenna pattern; priors, σ, σ_TB and the cost function are
# those of coh-model_optimization.jl and inc-model_optimization.jl.
#
#   baseline   Nelder-Mead with this script's cost function; must reproduce ./opt_coh and
#              ./opt_inc exactly (checks that the cost below is the one used in the paper)
#   footprint  snow depth varying within the footprint: TBs averaged over a Gaussian
#              distribution of snow depth (σ_d = 0.01 and 0.02 m, 15-point Gauss-Hermite
#              quadrature, depths below 0 set to 0), (i) evaluated at the optimal states
#              and (ii) with the optimization repeated
#   priors     each prior shifted by -1σ and +1σ, one variable at a time (snow depth not
#              below 0; axis ratio divided/multiplied by 1.3), optimization repeated from
#              the shifted prior
#   annealing  optimization repeated with simulated annealing (Optim SAMIN, rt = 0.5,
#              20,000 cost evaluations, seed 1000 + measurement number) within x_a ± 5σ
#              (axis ratio within 5 ×/÷ 1.3⁵), clipped to snow depth ≥ 0, snow density
#              ≥ 100 kg m⁻³, ice salinity ≥ 0.5 and ice temperature ≤ -2 °C
#
# Run with several threads, e.g.  julia --project=. -t 9 robustness_tests.jl
# Optional arguments: any of baseline, footprint, priors, annealing (default: all).
# Output: ./robustness/*.csv and ./robustness/summary.txt
#----------------------------------------------------------------------------------------

const SIG_TB, SIG_TICE, SIG_SICE, SIG_AR = 5.0, 2.0, 1.0, 1.3
const OUT = joinpath(@__DIR__, "robustness")

# Site priors and uncertainties (identical to coh/inc-model_optimization.jl)
const SITES = [
    (key = "clo", temp = -9.69, tsnow = -14.0,
     rho   = [355.0, 300, 300, 355, 300, 267, 300, 300, 300, 350],
     dice  = [0.945, 0.95, 0.94, 0.94, 0.97, 0.985, 0.98, 0.895, 0.895, 0.99],
     tice  = -13.0, sice = 5.32, ar = 5.0,
     sig_ds  = [0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01, 0.02, 0.02, 0.01],
     sig_rho = [25.0, 50, 50, 25, 50, 25, 50, 50, 50, 35],
     sig_di  = [0.01, 0.01, 0.05, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05]),
    (key = "mid1", temp = -7.37, tsnow = -15.0,
     rho   = [350.0, 365, 367, 350, 350, 350],
     dice  = [0.89, 0.855, 0.855, 0.89, 0.855, 0.855],
     tice  = -11.0, sice = 4.6, ar = 5.0,
     sig_ds  = fill(0.01, 6),
     sig_rho = [100.0, 50, 50, 100, 100, 100],
     sig_di  = fill(0.01, 6)),
    (key = "mid2", temp = -12.95, tsnow = -25.0,
     rho   = [350.0, 350, 300, 347],
     dice  = [0.93, 0.93, 0.93, 0.93],
     tice  = -20.0, sice = 4.5, ar = 5.0,
     sig_ds  = fill(0.02, 4),
     sig_rho = [100.0, 100, 50, 50],
     sig_di  = fill(0.05, 4)),
    (key = "far", temp = -13.86, tsnow = -25.0,
     rho   = [385.0, 385, 385, 385, 365, 400, 400, 400, 385, 385, 365, 400, 400, 400, 400],
     dice  = [0.86, 0.86, 0.86, 0.86, 0.92, 0.84, 0.84, 0.845, 0.855, 0.855, 0.855, 0.855, 0.855, 0.855, 0.855],
     tice  = -17.0, sice = 4.8, ar = 5.0,
     sig_ds  = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.01, 0.005, 0.02, 0.005],
     sig_rho = [50.0, 50, 50, 50, 50, 100, 100, 100, 50, 50, 50, 100, 100, 100, 100],
     sig_di  = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]),
]

#----------------------------------------------------------------------------------------
# Measurements: observations, prior x_a = [dsnow, rho, dice, axis_ratio, tice, sice] and σ
#----------------------------------------------------------------------------------------
struct Meas
    no::Int; site::String; tsnow::Float64
    y::Vector{Float64}          # [TBH, TBV]
    xa::Vector{Float64}         # [dsnow, rho, dice, axis_ratio, tice, sice]
    sig::Vector{Float64}        # [σ_dsnow, σ_rho, σ_dice, σ_axis_ratio (factor), σ_tice, σ_sice]
end

function measurements()
    data = CSV.read(joinpath(@__DIR__, "chars_data_filtered.csv"), DataFrame)
    out = Meas[]
    for s in SITES
        sd = filter(row -> row[:temp] == s.temp, data)
        for i in 1:nrow(sd)
            push!(out, Meas(length(out) + 1, s.key, s.tsnow, [sd[i, :tbh], sd[i, :tbv]],
                            [sd[i, :dsnow] * 0.01, s.rho[i], s.dice[i], s.ar, s.tice, s.sice],
                            [s.sig_ds[i], s.sig_rho[i], s.sig_di[i], SIG_AR, SIG_TICE, SIG_SICE]))
        end
    end
    return out
end

#----------------------------------------------------------------------------------------
# Forward models, optionally averaged over a Gaussian distribution of snow depth
#----------------------------------------------------------------------------------------
function gauss_hermite(n)
    J = SymTridiagonal(zeros(n), [sqrt(k / 2) for k in 1:n-1])
    nodes, vecs = eigen(J)
    return nodes, sqrt(π) .* vecs[1, :] .^ 2
end
const GH_X, GH_W = gauss_hermite(15)

fwd_point(model, x, ts) = model == :coh ?
    wilheit_model_integrated(40, [x[1]], [ts], [x[2]], [x[3]], [x[5]], [x[6]], "mix", 0, x[4]) :
    incoherent_model_integrated(40, x[1], ts, x[2], x[3], x[5], x[6], "mix", 0, x[4])

function fwd(model, x, ts; sig_d = 0.0)
    sig_d == 0 && return fwd_point(model, x, ts)
    tbh = 0.0; tbv = 0.0
    for (xk, wk) in zip(GH_X, GH_W)
        xd = copy(x); xd[1] = max(0.0, x[1] + sqrt(2) * sig_d * xk)
        h, v = fwd_point(model, xd, ts)
        tbh += wk * h; tbv += wk * v
    end
    return tbh / sqrt(π), tbv / sqrt(π)
end

# Cost of Eq. (2), as in coh/inc-model_optimization.jl (axis ratio term in log10)
function cost(model, x, m::Meas, xa; sig_d = 0.0)
    tbh, tbv = fwd(model, x, m.tsnow; sig_d = sig_d)
    return ((tbh - m.y[1])^2 + (tbv - m.y[2])^2) / SIG_TB^2 +
           (x[1] - xa[1])^2 / m.sig[1]^2 + (x[2] - xa[2])^2 / m.sig[2]^2 + (x[3] - xa[3])^2 / m.sig[3]^2 +
           (log10(x[4]) - log10(xa[4]))^2 / log10(m.sig[4])^2 +
           (x[5] - xa[5])^2 / m.sig[5]^2 + (x[6] - xa[6])^2 / m.sig[6]^2
end

function optimize_nm(model, m, xa; sig_d = 0.0)
    res = optimize(x -> cost(model, x, m, xa; sig_d = sig_d), copy(xa), NelderMead())
    return Optim.minimizer(res), Optim.minimum(res)
end

# Departure between two states in units of the prior σ (axis ratio in log10)
function z_diff(x1, x2, m)
    d = [(x1[k] - x2[k]) / m.sig[k] for k in (1, 2, 3, 5, 6)]
    insert!(d, 4, (log10(x1[4]) - log10(x2[4])) / log10(m.sig[4]))
    return d
end

const VARS = ["dsnow", "rho_snow", "dice", "axis_ratio", "tice", "sice"]

function bic_delta(cc, ci)
    K = length(cc)
    return (sum(ci) + 5 * K * log(8)) - (sum(cc) + 6 * K * log(8))
end

function load_reference(meas)
    ref = Dict{Symbol, Vector{Tuple{Vector{Float64}, Float64}}}()
    for (model, dir) in ((:coh, "opt_coh"), (:inc, "opt_inc"))
        v = Tuple{Vector{Float64}, Float64}[]
        for s in SITES
            o = CSV.read(joinpath(@__DIR__, dir, "opt_res_$(s.key)_jl.csv"), DataFrame)
            for r in eachrow(o)
                push!(v, ([r.opt_dsnow, r.opt_rho_snow, r.opt_dice, r.opt_axis_ratio, r.opt_tice, r.opt_sice], r.cost))
            end
        end
        ref[model] = v
    end
    return ref
end

function threaded_map(f, n)
    out = Vector{Any}(undef, n)
    Threads.@threads for i in 1:n
        out[i] = f(i)
    end
    return out
end

#----------------------------------------------------------------------------------------
# Tests
#----------------------------------------------------------------------------------------
function test_baseline(meas, ref, io)
    for model in (:coh, :inc)
        res = threaded_map(i -> optimize_nm(model, meas[i], meas[i].xa), length(meas))
        dx = maximum(maximum(abs.(res[i][1] .- ref[model][i][1])) for i in eachindex(meas))
        dc = maximum(abs(res[i][2] - ref[model][i][2]) for i in eachindex(meas))
        println(io, @sprintf("baseline %s: max |Δx| = %.1e, max |Δcost| = %.1e (vs ./opt_%s)", model, dx, dc, model))
    end
end

function test_footprint(meas, ref, io)
    rows = DataFrame(no = Int[], site = String[], model = String[], sig_d = Float64[],
                     dtbh_at_opt = Float64[], dtbv_at_opt = Float64[], cost_reopt = Float64[], cost_ref = Float64[])
    totals = Dict{Tuple{Symbol, Float64}, Vector{Float64}}()
    for sig_d in (0.01, 0.02), model in (:coh, :inc)
        res = threaded_map(length(meas)) do i
            m = meas[i]; xopt = ref[model][i][1]
            h0, v0 = fwd(model, xopt, m.tsnow)
            h1, v1 = fwd(model, xopt, m.tsnow; sig_d = sig_d)
            _, c = optimize_nm(model, m, m.xa; sig_d = sig_d)
            (h1 - h0, v1 - v0, c)
        end
        for (i, r) in enumerate(res)
            push!(rows, (meas[i].no, meas[i].site, String(model), sig_d, r[1], r[2], r[3], ref[model][i][2]))
        end
        totals[(model, sig_d)] = [r[3] for r in res]
        dtb = vcat([abs(r[1]) for r in res], [abs(r[2]) for r in res])
        println(io, @sprintf("footprint σ_d = %.2f m, %s: at the optimal states max |ΔTB| = %.2f K, RMS ΔTB = %.2f K; re-optimized total cost = %.2f (reference %.2f)",
                              sig_d, model, maximum(dtb), sqrt(sum(dtb .^ 2) / length(dtb)), sum(r[3] for r in res), sum(c for (_, c) in ref[model])))
    end
    for sig_d in (0.01, 0.02)
        cc, ci = totals[(:coh, sig_d)], totals[(:inc, sig_d)]
        println(io, @sprintf("footprint σ_d = %.2f m: coherent cost lower for %d of %d; ΔBIC = %.1f",
                              sig_d, count(cc .< ci), length(cc), bic_delta(cc, ci)))
    end
    CSV.write(joinpath(OUT, "footprint_snow_depth.csv"), rows)
end

function shifted_prior(m, k, sgn)
    xa = copy(m.xa)
    if k == 4
        xa[4] = sgn > 0 ? xa[4] * m.sig[4] : xa[4] / m.sig[4]
    else
        xa[k] += sgn * m.sig[k]
    end
    k == 1 && (xa[1] = max(0.0, xa[1]))
    return xa
end

function test_priors(meas, ref, io)
    rows = DataFrame(no = Int[], site = String[], model = String[], variable = String[], shift = Int[],
                     cost = Float64[], cost_ref = Float64[],
                     dz_dsnow = Float64[], dz_rho_snow = Float64[], dz_dice = Float64[],
                     dz_axis_ratio = Float64[], dz_tice = Float64[], dz_sice = Float64[])
    for k in 1:6, sgn in (-1, 1)
        tot = Dict{Symbol, Vector{Float64}}()
        for model in (:coh, :inc)
            res = threaded_map(length(meas)) do i
                m = meas[i]
                x, c = optimize_nm(model, m, shifted_prior(m, k, sgn))
                (c, z_diff(x, ref[model][i][1], m))
            end
            tot[model] = [r[1] for r in res]
            for (i, r) in enumerate(res)
                push!(rows, (meas[i].no, meas[i].site, String(model), VARS[k], sgn, r[1], ref[model][i][2], r[2]...))
            end
            zself = [abs(r[2][k]) for r in res]
            zother = [maximum(abs.(r[2][setdiff(1:6, k)])) for r in res]
            println(io, @sprintf("priors %-10s %+dσ, %s: total cost = %.2f (reference %.2f); shift of optimal %s: median %.2fσ, max %.2fσ; largest shift of the other variables: median %.2fσ, max %.2fσ",
                                  VARS[k], sgn, model, sum(tot[model]), sum(c for (_, c) in ref[model]),
                                  VARS[k], sort(zself)[(end + 1) ÷ 2], maximum(zself), sort(zother)[(end + 1) ÷ 2], maximum(zother)))
        end
        println(io, @sprintf("priors %-10s %+dσ: coherent cost lower for %d of %d; ΔBIC = %.1f",
                              VARS[k], sgn, count(tot[:coh] .< tot[:inc]), length(meas), bic_delta(tot[:coh], tot[:inc])))
    end
    CSV.write(joinpath(OUT, "prior_shifts.csv"), rows)
end

# x_a ± 5σ, clipped to physical limits (snow density ≥ 100 kg m⁻³, ice salinity ≥ 0.5,
# ice temperature ≤ −2 °C); without the clipping, ρ_snow can go negative (MID-1/2, σ = 100)
function bounds(m)
    lo = [max(0.0, m.xa[1] - 5m.sig[1]), max(100.0, m.xa[2] - 5m.sig[2]), m.xa[3] - 5m.sig[3], m.xa[4] / m.sig[4]^5, m.xa[5] - 5m.sig[5], max(0.5, m.xa[6] - 5m.sig[6])]
    hi = [m.xa[1] + 5m.sig[1], m.xa[2] + 5m.sig[2], m.xa[3] + 5m.sig[3], m.xa[4] * m.sig[4]^5, min(-2.0, m.xa[5] + 5m.sig[5]), m.xa[6] + 5m.sig[6]]
    return lo, hi
end

function test_annealing(meas, ref, io)
    rows = DataFrame(no = Int[], site = String[], model = String[], cost_sa = Float64[], cost_nm = Float64[],
                     dz_dsnow = Float64[], dz_rho_snow = Float64[], dz_dice = Float64[],
                     dz_axis_ratio = Float64[], dz_tice = Float64[], dz_sice = Float64[])
    for model in (:coh, :inc)
        res = threaded_map(length(meas)) do i
            m = meas[i]
            lo, hi = bounds(m)
            Random.seed!(1000 + m.no)
            r = optimize(x -> cost(model, x, m, m.xa), lo, hi, clamp.(copy(m.xa), lo, hi), SAMIN(rt = 0.5),
                         Optim.Options(iterations = 20_000))
            (Optim.minimum(r), z_diff(Optim.minimizer(r), ref[model][i][1], m))
        end
        for (i, r) in enumerate(res)
            push!(rows, (meas[i].no, meas[i].site, String(model), r[1], ref[model][i][2], r[2]...))
        end
        dc = [r[1] - ref[model][i][2] for (i, r) in enumerate(res)]
        dz = [maximum(abs.(r[2])) for r in res]
        println(io, @sprintf("annealing %s: total cost = %.2f (Nelder-Mead %.2f); cost SA − NM: min %.3f, max %.3f; largest state difference: median %.3fσ, max %.3fσ",
                              model, sum(r[1] for r in res), sum(c for (_, c) in ref[model]), minimum(dc), maximum(dc),
                              sort(dz)[(end + 1) ÷ 2], maximum(dz)))
    end
    CSV.write(joinpath(OUT, "annealing.csv"), rows)
end

#----------------------------------------------------------------------------------------
function main()
    mkpath(OUT)
    which = isempty(ARGS) ? ["baseline", "footprint", "priors", "annealing"] : ARGS
    meas = measurements()
    ref = load_reference(meas)
    open(joinpath(OUT, "summary.txt"), "a") do f
        println(f, "\n# robustness_tests.jl, threads = $(Threads.nthreads()), tests = $(join(which, ", "))")
        for t in which
            t == "baseline"  && test_baseline(meas, ref, f)
            t == "footprint" && test_footprint(meas, ref, f)
            t == "priors"    && test_priors(meas, ref, f)
            t == "annealing" && test_annealing(meas, ref, f)
            flush(f)
            println("$t done")
        end
    end
    print(read(joinpath(OUT, "summary.txt"), String))
end
#----------------------------------------------------------------------------------------

main()
