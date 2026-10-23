# Coherent L-Band Emission from Snow-Covered Arctic Sea Ice

This repository contains the code and data for reproducing the analysis presented in the article **"Observations of Coherent L-Band Emission from Snow-Covered Arctic Sea Ice"**.

## Overview

This research investigates microwave emission from snow-covered Arctic sea ice at L-band frequencies (1.4 GHz), comparing coherent and incoherent radiative transfer models to explain observed brightness temperatures. The code implements both modeling approaches and includes optimization routines to retrieve geophysical parameters from brightness temperature observations.

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
├── paper_figures.jl              # Figures of the article and supplement
├── coherence_damping.jl          # Damping of the coherent oscillations with snow depth
├── Project.toml, Manifest.toml   # Julia environment (pinned package versions)
├── opt_coh/                      # Optimization results, coherent model
├── opt_inc/                      # Optimization results, incoherent model
└── opt_fixed_alpha/              # Optimization results with the axis ratio fixed
```

## Models

### 1. Permittivity Models (`permittivity_models.jl`)

Implements dielectric permittivity models for sea ice constituents:

- **`epsilon_Vant_ice(S, T)`**: Pure ice permittivity following Vant (1978)
- **`epsilon_snow(freq, rho_snow, T)`**: Snow permittivity as function of density and temperature
- **`epsilon_water(freq, T, S)`**: Seawater permittivity (Stogryn & Desargant, 1985)
- **`epsilon_rn_ice(freq, S, T)`**: Sea ice permittivity with random needle brine inclusions (Shokr, 1998)
- **`epsilon_sp_ice(freq, S, T)`**: Sea ice permittivity with spherical brine inclusions
- **`pvs_mixing(mixing_ratio, freq, S, T)`**: Mixing of random needles and spheres

**Key Parameters:**
- `freq`: Frequency (Hz)
- `T`: Temperature (°C)
- `S`: Salinity 
- `rho_snow`: Snow density (kg/m³)

### 2. Ellipsoidal Inclusion Model (`ellipsoids_perm.jl`)

Implements effective medium theory for sea ice with ellipsoidal brine inclusions:

- **`depolarization_factor(aspect_ratio, c)`**: Calculates depolarization factors for ellipsoids
- **`e_eff_mix(axis_ratio, c, v, S, T)`**: Computes effective permittivity using:
  - Maxwell-Garnet formulation (v=0)
  - Coherent potential formulation (v=1)
  - Polder-van Santen formulation (v=2)

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
- `ice_permittivity`: Permittivity model ("vant", "rn", "sp", "pvs", "mix")
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

Both `coh-model_optimization.jl` and `inc-model_optimization.jl` implement Nelder-Mead optimization (Optim.jl default options, starting from the prior) to retrieve geophysical parameters from observed brightness temperatures.

### Cost Function

The objective function minimizes:

```julia
Cost = Σ[(TB_modeled - TB_observed)²/σ_TB²] + Σ[(param - prior)²/σ_param²]
```

Where:
- Brightness temperature uncertainty: σ_TB = 5 K
- Prior constraints on all retrieved parameters (snow depth, density, ice thickness, etc.); the axis-ratio term is evaluated in log10
- Both forward models are integrated over the ARIEL antenna pattern (Gaussian in incidence angle centred at 40°, σ_H = 15.29°, σ_V = 14.87°), so the two optimizations differ only in the treatment of phase

### Optimization Variables

Retrieved parameters (6 per measurement):
1. `dsnow`: Snow depth (m)
2. `rho_snow`: Snow density (kg/m³)
3. `dice`: Ice thickness (m)
4. `axis_ratio`: Brine pocket axis ratio
5. `tice`: Ice temperature (°C)
6. `sice`: Ice salinity (ppt)

### Sampling Sites

Data from four sites with different environmental conditions: CLO, MID-1, MID-2, and FAR.

### Axis ratio fixed (`fixed-alpha_optimization.jl`)

Sensitivity test in which the axis ratio is fixed at its prior value (α = 5) for both models and removed from the state vector ([dsnow, rho_snow, dice, tice, sice]); the cost has no axis-ratio term and all other settings are as above. Results are written to `./opt_fixed_alpha/opt_res_<site>_{coh,inc}.csv`, which also record the number of iterations and whether the optimization converged.

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
```

The optimizations are deterministic, so rerunning them reproduces the CSV files provided. `paper_figures.jl` writes Figs. 3 to 8, the PD comparison at the prior and optimal states, and Fig. S1 of the supplement.

## Changes Since the First Release

October 2026, after the review:

- **Snow permittivity**: the conversion of snow density to g cm⁻³ was commented out for densities above 400 kg m⁻³ (`epsilon_snow`), which made the imaginary part about 3 × 10⁵ times too large in that range. Fixed; all results were recomputed.
- **Incoherent model**: integrated over the ARIEL antenna pattern in the optimization (`incoherent_model_integrated`), as the coherent model already was; previously it was evaluated at 40° only. Seawater temperature set to -1.8 °C, as in the coherent model (previously -2 °C).
- **Optimization scripts**: the signature `minimization(:DataFrame)` is not valid Julia and prevented the scripts from running; corrected to `minimization(data::DataFrame)`.
- **New**: optimization with the axis ratio fixed, scripts for the figures and for the damping analysis, and a pinned Julia environment.

## Requirements

Julia 1.12.3. The package versions used for the results are pinned in `Project.toml` and `Manifest.toml` (BandedMatrices 1.11.0, CSV 0.10.16, DataFrames 1.8.1, Optim 1.13.3, Plots 1.41.6); install them with

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```