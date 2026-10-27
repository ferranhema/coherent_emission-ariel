# Coherent L-Band Emission from Snow-Covered Arctic Sea Ice

This repository contains the code and data for reproducing the analysis presented in the article **"Observations of Coherent L-Band Emission from Snow-Covered Arctic Sea Ice"**.

## Overview

This research investigates microwave emission from snow-covered Arctic sea ice at L-band frequencies (1.4 GHz), comparing coherent and incoherent radiative transfer models to explain observed brightness temperatures. The code implements both modeling approaches and the optimization that, for each measurement, adjusts the geophysical state within its prior uncertainties to reproduce the observed brightness temperatures.

## Repository Structure

```
.
├── README.md                      
├── chars_data_filtered.csv        # Filtered in situ observations
├── permittivity_models.jl         # Dielectric permittivity models
├── ellipsoids_perm.jl            # Ellipsoidal inclusion mixing formulas
├── coherent_model.jl             # Coherent radiative transfer model
├── incoherent_model.jl           # Incoherent radiative transfer model
├── coh-model_optimization.jl     # Parameter optimization (coherent)
├── inc-model_optimization.jl     # Parameter optimization (incoherent)
├── fixed-alpha_optimization.jl   # Optimization with the axis ratio fixed (both models)
├── paper_figures.jl              # Figures of the article (Figs. 3-9 and B1)
├── coherence_damping.jl          # Damping of the coherent oscillations with snow depth
├── robustness_tests.jl           # Snow depth variability, prior shifts, simulated annealing
├── Project.toml, Manifest.toml   # Julia environment (pinned package versions)
├── LICENSE                       # MIT licence
├── opt_coh/                      # Optimization results, coherent model
├── opt_inc/                      # Optimization results, incoherent model
├── opt_fixed_alpha/              # Optimization results with the axis ratio fixed
└── robustness/                   # Results of robustness_tests.jl
```

## Models

### 1. Permittivity Models (`permittivity_models.jl`)

Implements dielectric permittivity models for snow, seawater and sea ice. The paper uses `epsilon_snow`, `epsilon_water`, and the sea ice option `"mix"` (`e_eff_mix` in `ellipsoids_perm.jl`, below). The other sea ice options are kept in the code but are **not used** in the paper.

Used in the paper:
- **`epsilon_snow(freq, rho_snow, T)`**: Snow permittivity; real part from Mätzler (1996), imaginary part from Tiuri et al. (1984)
- **`epsilon_water(freq, T, S)`**: Seawater permittivity (Klein & Swift, 1977)

Not used in the paper:
- **`epsilon_Vant_ice(S, T)`**: Sea ice permittivity following Vant (1978) (option `"vant"`)
- **`epsilon_rn_ice(freq, S, T)`**: Sea ice permittivity with random needle brine inclusions (option `"rn"`)
- **`epsilon_sp_ice(freq, S, T)`**: Sea ice permittivity with spherical brine inclusions (option `"sp"`)
- **`pvs_mixing(mixing_ratio, freq, S, T)`**: Mixing of random needles and spheres (option `"pvs"`)

**Key Parameters:**
- `freq`: Frequency (Hz)
- `T`: Temperature (°C)
- `S`: Salinity 
- `rho_snow`: Snow density (kg/m³)

### 2. Ellipsoidal Inclusion Model (`ellipsoids_perm.jl`)

Implements effective medium theory for sea ice with ellipsoidal brine inclusions:

- **`depolarization_factor(aspect_ratio, c)`**: Calculates depolarization factors for ellipsoids
- **`e_eff_mix(axis_ratio, c, v, S, T)`**: Computes effective permittivity using:
  - Maxwell-Garnett formulation (v=0) — **used in the paper**
  - Coherent potential formulation (v=1) — not used
  - Polder-van Santen formulation (v=2) — not used

In the configuration used in the paper (option `"mix"`, v = 0), brine inclusions are randomly oriented spheroids of axis ratio α in a pure ice host. The brine permittivity follows Stogryn & Desargant (1985) and the pure ice permittivity follows Mätzler (2006), with ε' = 3.1884 + 9.1×10⁻⁴ T.

**Key Features:**
- Iterative solution for effective medium permittivity
- Accounts for brine pocket geometry via axis ratios
- Temperature and salinity-dependent brine volume

### 3. Coherent Radiative Transfer Model (`coherent_model.jl`)

Implements the Wilheit (1978) coherent plane-stratified dielectric model:

- **`hpld(n, sn, ths, theta, xlam)`**: Horizontal polarization layer-decomposed absorption
- **`vpld(n, sn, ths, theta, xlam)`**: Vertical polarization layer-decomposed absorption
- **`wilheit_model(theta, ...)`**: Model at a single incidence angle (used only in Fig. 3)
- **`wilheit_model_integrated(theta, ...)`**: Main model with antenna pattern integration

**Model Features:**
- Multi-layer snow-ice-water system
- Accounts for coherent interference effects
- Antenna pattern integration over ±90° viewing angles
- Returns brightness temperatures for H and V polarizations

**Input Parameters:**
- `theta`: Incidence angle (degrees)
- `dsnow`: Snow layer thickness(es) (m)
- `tsnow`: Snow temperature(s) (°C)
- `rho_snow`: Snow density(ies) (kg/m³)
- `dice`: Ice layer thickness(es) (m)
- `tice`: Ice temperature(s) (°C)
- `sice`: Ice salinity(ies)
- `ice_permittivity`: Permittivity model ("vant", "rn", "sp", "pvs", "mix"; the paper uses "mix")
- `mixing_ratio`: Fraction of random needle inclusions (for "pvs")
- `axis_ratio`: Brine ellipsoid axis ratio (for "mix")

### 4. Incoherent Radiative Transfer Model (`incoherent_model.jl`)

Implements an incoherent emission model following MEMLS approach:

- **`fresnel(n1, n2, θ)`**: Fresnel reflection coefficients
- **`incoherent(da, N, θ, freq, T)`**: Incoherent emission with infinite reflections
- **`incoherent_model(theta, ...)`**: Incoherent model at a single incidence angle (used only in Fig. 3)
- **`incoherent_model_integrated(theta, ...)`**: Main incoherent model, integrated over the same antenna pattern as the coherent model

**Key Differences from Coherent Model:**
- Neglects phase information between layers
- Uses matrix formulation for multiple reflections
- Computationally faster but less accurate for thin layers

## Optimization Routines

Both `coh-model_optimization.jl` and `inc-model_optimization.jl` adjust, for each measurement, the geophysical state within its prior uncertainties so that the forward model reproduces the observed brightness temperatures.

### Priors

The priors and their uncertainties are hard-coded in the `minimization()` function of each optimization script (identical in both, and repeated in `fixed-alpha_optimization.jl` and `paper_figures.jl`); they are listed in Table A1 of the paper. The snow depth and ice thickness priors are the in situ values in `chars_data_filtered.csv`. That file has no snow density column, so the density priors exist only in the scripts.

### Optimization Settings

- Nelder–Mead (Optim.jl 1.13.3) with Optim's default options: affine initial simplex around the starting point, at most 1000 iterations, convergence tolerance `g_abstol` = 10⁻⁸
- Unconstrained (no bounds on any parameter); for example, a snow-free measurement can end at a snow depth of about −1×10⁻⁶ m
- The optimization starts at the prior
- Antenna-pattern integration: TBs are computed at 1000 incidence angles evenly spaced from −89.99° to 89.99° and weighted by a Gaussian in incidence angle centred at 40° (σ_H = 15.29°, σ_V = 14.87°) times the angle step. The weights are not renormalised: the part of the Gaussian beyond ±90° is dropped, so they sum to 0.9995 (H) and 0.9996 (V). Negative angles are included; both models are symmetric in incidence angle, so these act as contributions at |θ| and carry 0.4% of the weight

### Cost Function

The objective function minimizes:

```julia
Cost = Σ[(TB_modeled - TB_observed)²/σ_TB²] + Σ[(param - prior)²/σ_param²]
```

Where:
- Brightness temperature uncertainty: σ_TB = 5 K
- Prior terms for all parameters of the state (snow depth, density, ice thickness, etc.); the axis-ratio term is evaluated in log10
- Both forward models are integrated over the ARIEL antenna pattern (Gaussian in incidence angle centred at 40°, σ_H = 15.29°, σ_V = 14.87°), so the two optimizations differ only in the treatment of phase

### Optimization Variables

State adjusted in the optimization (6 parameters per measurement):
1. `dsnow`: Snow depth (m)
2. `rho_snow`: Snow density (kg/m³)
3. `dice`: Ice thickness (m)
4. `axis_ratio`: Brine pocket axis ratio
5. `tice`: Ice temperature (°C)
6. `sice`: Ice salinity (ppt)

### Sampling Sites

Data from four sites with different environmental conditions: CLO, MID-1, MID-2, and FAR.

### Axis ratio fixed (`fixed-alpha_optimization.jl`)

Sensitivity test (Fig. B1 of the paper) in which the axis ratio is fixed at its prior value (α = 5) for both models and removed from the state vector ([dsnow, rho_snow, dice, tice, sice]); the cost has no axis-ratio term and all other settings are as above. Results are written to `./opt_fixed_alpha/opt_res_<site>_{coh,inc}.csv`, which also record the number of iterations and whether the optimization converged.

### Robustness tests (`robustness_tests.jl`)

Tests of the optimization reported in the discussion of the paper. Both forward models are antenna-integrated, and the priors, σ values and cost function are those of the optimization scripts.

- **baseline**: Nelder–Mead with the script's own cost function; reproduces `opt_coh/` and `opt_inc/` exactly (largest cost difference 1.8 × 10⁻¹⁵), which checks that the tests use the cost function of the paper.
- **footprint**: snow depth varying within the footprint. TBs are averaged over a Gaussian distribution of snow depth (σ_d = 0.01 and 0.02 m, 15-point Gauss–Hermite quadrature, depths below 0 set to 0), evaluated at the optimal states and with the optimization repeated.
- **priors**: each prior shifted by −1σ and +1σ, one variable at a time (snow depth not below 0, axis ratio divided or multiplied by 1.3), with the optimization repeated from the shifted prior.
- **annealing**: optimization repeated with simulated annealing (Optim `SAMIN`, rt = 0.5, 20,000 cost evaluations, seed 1000 + measurement number) within x_a ± 5σ (axis ratio within 5 ×/÷ 1.3⁵), clipped to snow depth ≥ 0, snow density ≥ 100 kg m⁻³, ice salinity ≥ 0.5 and ice temperature ≤ −2 °C.

Results are written to `./robustness/`: `footprint_snow_depth.csv`, `prior_shifts.csv`, `annealing.csv` (per measurement and model) and `summary.txt` (totals, number of measurements with the lower coherent cost, and ΔBIC = Δcost − 35 ln 8). The four tests take about 2 hours on 9 threads; the simulated annealing depends on the random number generator, so its results may differ with other Julia versions.

### Output

Optimization results are saved as CSV files in `./opt_coh/` and `./opt_inc/` directories:
- `opt_res_clo_jl.csv`
- `opt_res_mid1_jl.csv`
- `opt_res_mid2_jl.csv`
- `opt_res_far_jl.csv`

Each file contains:
- Optimized parameters
- Modeled brightness temperatures (H and V polarizations)
- Target (observed) brightness temperatures
- Final cost function value

## Data File

### `chars_data_filtered.csv`

Contains the in situ observations from Arctic sea ice measurements.

**Columns:**
- `index`: Measurement index
- `tbh`: Horizontal polarization brightness temperature (K)
- `tbv`: Vertical polarization brightness temperature (K)
- `pd`: Polarization difference (K)
- `tsurf`: Surface temperature (K)
- `sal`: Ice salinity
- `temp`: Air/ice temperature identifier (°C)
- `dsnow`: Snow depth (cm)
- `dice`: Ice thickness (cm)

**Total Observations:** 35 measurements across 4 sites

## Reproducing the Analysis

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. coh-model_optimization.jl      # -> opt_coh/
julia --project=. inc-model_optimization.jl      # -> opt_inc/
julia --project=. fixed-alpha_optimization.jl    # -> opt_fixed_alpha/
julia --project=. paper_figures.jl               # -> figures/
julia --project=. coherence_damping.jl           # -> figures/ and console summary
julia --project=. -t 9 robustness_tests.jl       # -> robustness/
```

The optimizations are deterministic, so rerunning them reproduces the CSV files provided. `paper_figures.jl` writes Figs. 3 to 9 and Fig. B1 (appendix) of the paper.

## Changes Since the First Release

October 2026, after the review:

- **Snow permittivity**: the conversion of snow density to g cm⁻³ was commented out for densities above 400 kg m⁻³ (`epsilon_snow`), which made the imaginary part about 3 × 10⁵ times too large in that range. Fixed; all results were recomputed.
- **Incoherent model**: integrated over the ARIEL antenna pattern in the optimization (`incoherent_model_integrated`), as the coherent model already was; previously it was evaluated at 40° only. Seawater temperature set to -1.8 °C, as in the coherent model (previously -2 °C).
- **Optimization scripts**: the signature `minimization(:DataFrame)` is not valid Julia and prevented the scripts from running; corrected to `minimization(data::DataFrame)`.
- **New**: optimization with the axis ratio fixed, scripts for the figures and for the damping analysis, and a pinned Julia environment.
- **New** (after v1.0.0): robustness tests of the optimization (`robustness_tests.jl`, results in `robustness/`).

## Requirements

Julia 1.12.3. The package versions used for the results are pinned in `Project.toml` and `Manifest.toml` (BandedMatrices 1.11.0, CSV 0.10.16, DataFrames 1.8.1, Optim 1.13.3, Plots 1.41.6); install them with

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Citation

The code is archived on Zenodo. To cite the version used in the article:

Hernández-Macià, F. and Huntemann, M.: ferranhema/coherent_emission-ariel: v1.1.0, Zenodo [code], https://doi.org/10.5281/zenodo.23665110, 2026.

The concept DOI https://doi.org/10.5281/zenodo.23597502 always resolves to the latest version.

## License

The code is released under the MIT licence (see `LICENSE`).
