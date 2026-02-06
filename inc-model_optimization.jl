using DataFrames
using CSV
using Optim
include("incoherent_model.jl")
#----------------------------------------------------------------------------------------
"""
    objective_function(x, target_tbh, target_tbv, ...)

Cost function combining brightness temperature misfit and prior parameter constraints
for the incoherent emission model.

# Arguments
- `x`: Optimization variables [dsnow, rho_snow, dice, axis_ratio, tice, sice]
- `target_tbh/tbv`: Observed brightness temperatures (K)
- Prior means (mu_*) and uncertainties (sigma_*) for each parameter

# Returns
Total cost: TB misfit + parameter prior penalties
"""
function objective_function(x, target_tbh, target_tbv, theta, tair, tsnow, ice_permittivity, mixing_ratio, mu_dsnow, sigma_dsnow, mu_rhosnow, 
                            sigma_rhosnow, mu_dice, sigma_dice, mu_axis_ratio, sigma_axis_ratio, mu_tice, sigma_tice, mu_sice, sigma_sice)
    # Extract parameters from optimization vector
    dsnow = x[1:length(tsnow)]                                          # Snow depth (m)
    rho_snow = x[length(tsnow)+1:2*length(tsnow)]                      # Snow density (kg/m³)
    dice = x[2*length(tsnow)+1:2*length(tsnow)+length(tsnow)]         # Ice thickness (m)
    axis_ratio = x[2*length(tsnow)+length(tsnow)+1:2*length(tsnow)+2*length(tsnow)]  # Brine axis ratio
    tice = x[2*length(tsnow)+2*length(tsnow)+1:2*length(tsnow)+3*length(tsnow)]      # Ice temp (°C)
    sice = x[2*length(tsnow)+3*length(tsnow)+1:end]                   # Ice salinity (ppt)
    
    # Compute modeled brightness temperatures using incoherent model
    tbh, tbv = incoherent_model(theta, dsnow[1], tsnow[1], rho_snow[1], dice[1], tice[1], sice[1], ice_permittivity, mixing_ratio, axis_ratio[1])
    
    # Cost function: TB misfit + prior penalties
    sigma_tb = 5  # TB measurement uncertainty (K)
    error = (((tbh - target_tbh)^2) / (sigma_tb^2)) + (((tbv - target_tbv)^2) / (sigma_tb^2)) + 
            (((dsnow[1] - mu_dsnow[1])^2) / (sigma_dsnow^2)) + 
            (((rho_snow[1] - mu_rhosnow[1])^2) / (sigma_rhosnow^2)) + 
            (((dice[1] - mu_dice[1])^2) / (sigma_dice^2)) + 
            (((log10(axis_ratio[1]) - log10(mu_axis_ratio[1]))^2) / (log10(sigma_axis_ratio)^2) +
            (((tice[1] - mu_tice[1])^2) / (sigma_tice^2)) +
            (((sice[1] - mu_sice[1])^2) / (sigma_sice^2)))
    println(error)
    return error
end

#----------------------------------------------------------------------------------------
"""
    minimization(:DataFrame)

Main optimization routine for four sampling sites (clo, mid1, mid2, far).
For each site and measurement, optimizes snow/ice parameters using Nelder-Mead 
with the incoherent emission model.
Outputs results to CSV files in ./opt_inc/ directory.
"""
function minimization(:DataFrame)
    # Filter data by sampling site
    clo = filter(row -> row[:temp] == -9.69, data)
    mid1 = filter(row -> row[:temp] == -7.37, data)
    mid2 = filter(row -> row[:temp] == -12.95, data)
    far = filter(row -> row[:temp] == -13.86, data)
    
    # Site parameters: air temp, snow temp, and prior values for optimization
    # CLOSE SITE (10 measurements)
    tair_clo = -13.5
    tsnow_clo = [-14]
    rho_snow_clo = [355, 300, 300, 355, 300, 267, 300, 300, 300, 350]
    dice_clo = [0.945, 0.95, 0.94, 0.94, 0.97, 0.985, 0.98, 0.895, 0.895, 0.99]
    tice_clo = [-13, -13, -13, -13, -13, -13, -13, -13, -13, -13]
    sice_clo = [5.32, 5.32, 5.32, 5.32, 5.32, 5.32, 5.32, 5.32, 5.32, 5.32]
    ice_perm_clo = "mix"
    mr_clo = 0
    ar_clo = 5
    sigma_dsnow_clo = [0.01, 0.01, 0.02, 0.02, 0.01, 0.01, 0.01, 0.02, 0.02, 0.01] 
    sigma_rhosnow_clo = [25, 50, 50, 25, 50, 25, 50, 50, 50, 35]
    sigma_dice_clo = [0.01, 0.01, 0.05, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05]
    sigma_axis_ratio_clo = [1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3]
    sigma_tice_clo = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    sigma_sice_clo = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    # MID1 SITE (6 measurements)
    tair_mid1 = -15
    tsnow_mid1 = [-15]
    rho_snow_mid1 = [350, 365, 367, 350, 350, 350]
    dice_mid1 = [0.89, 0.855, 0.855, 0.89, 0.855, 0.855]
    tice_mid1 = [-11, -11, -11, -11, -11, -11]
    sice_mid1 = [4.6, 4.6, 4.6, 4.6, 4.6, 4.6]
    ice_perm_mid1 = "mix"
    mr_mid1 = 0
    ar_mid1 = 5
    sigma_dsnow_mid1 = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    sigma_rhosnow_mid1 = [100, 50, 50, 100, 100, 100]
    sigma_dice_mid1 = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
    sigma_axis_ratio_mid1 = [1.3, 1.3, 1.3, 1.3, 1.3, 1.3]
    sigma_tice_mid1 = [2, 2, 2, 2, 2, 2]
    sigma_sice_mid1 = [1, 1, 1, 1, 1, 1]

    # MID2 SITE (4 measurements)
    tair_mid2 = -30
    tsnow_mid2 = [-25]
    rho_snow_mid2 = [350, 350, 300, 347]
    dice_mid2 = [0.93, 0.93, 0.93, 0.93]
    tice_mid2 = [-20, -20, -20, -20]
    sice_mid2 = [4.5, 4.5, 4.5, 4.5]
    ice_perm_mid2 = "mix"
    mr_mid2 = 0
    ar_mid2 = 5
    sigma_dsnow_mid2 = [0.02, 0.02, 0.02, 0.02]
    sigma_rhosnow_mid2 = [100, 100, 50, 50]
    sigma_dice_mid2 = [0.05, 0.05, 0.05, 0.05]
    sigma_axis_ratio_mid2 = [1.3, 1.3, 1.3, 1.3]
    sigma_tice_mid2 = [2, 2, 2, 2]
    sigma_sice_mid2 = [1, 1, 1, 1]

    # FAR SITE (15 measurements)
    tair_far = -27
    tsnow_far = [-25]
    rho_snow_far = [385, 385, 385, 385, 365, 400, 400, 400, 385, 385, 365, 400, 400, 400, 400]
    dice_far = [0.86, 0.86, 0.86, 0.86, 0.92, 0.84, 0.84, 0.845, 0.855, 0.855, 0.855, 0.855, 0.855, 0.855, 0.855]
    tice_far = [-17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17, -17]    
    sice_far = [4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8, 4.8]
    ice_perm_far = "mix"
    mr_far = 0
    ar_far = 5
    sigma_dsnow_far = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.01, 0.005, 0.02, 0.005]
    sigma_rhosnow_far = [50, 50, 50, 50, 50, 100, 100, 100, 50, 50, 50, 100, 100, 100, 100]
    sigma_dice_far = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
    sigma_axis_ratio_far = [1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3]
    sigma_tice_far = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    sigma_sice_far = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    # Organize parameters for iteration
    tair = [tair_clo, tair_mid1, tair_mid2, tair_far]
    tsnow = [tsnow_clo, tsnow_mid1, tsnow_mid2, tsnow_far]
    rho_snow = [rho_snow_clo, rho_snow_mid1, rho_snow_mid2, rho_snow_far]
    dice = [dice_clo, dice_mid1, dice_mid2, dice_far]
    tice = [tice_clo, tice_mid1, tice_mid2, tice_far]
    sice = [sice_clo, sice_mid1, sice_mid2, sice_far]
    ice_perm = [ice_perm_clo, ice_perm_mid1, ice_perm_mid2, ice_perm_far]
    mr = [mr_clo, mr_mid1, mr_mid2, mr_far]
    ar = [ar_clo, ar_mid1, ar_mid2, ar_far]
    sigma_dsnow = [sigma_dsnow_clo, sigma_dsnow_mid1, sigma_dsnow_mid2, sigma_dsnow_far]
    sigma_rhosnow = [sigma_rhosnow_clo, sigma_rhosnow_mid1, sigma_rhosnow_mid2, sigma_rhosnow_far]
    sigma_dice = [sigma_dice_clo, sigma_dice_mid1, sigma_dice_mid2, sigma_dice_far]
    sigma_axis_ratio = [sigma_axis_ratio_clo, sigma_axis_ratio_mid1, sigma_axis_ratio_mid2, sigma_axis_ratio_far]
    sigma_tice = [sigma_tice_clo, sigma_tice_mid1, sigma_tice_mid2, sigma_tice_far]
    sigma_sice = [sigma_sice_clo, sigma_sice_mid1, sigma_sice_mid2, sigma_sice_far]

    theta = 40  # Incidence angle (degrees)
    
    # Loop over four sites
    for k in 1:4
        opt_tbh = []
        opt_tbv = []
        opt_dsnow = []
        opt_rho_snow = []
        opt_dice = []
        opt_axis_ratio = []
        opt_tice = []
        opt_sice = []
        costs = []
        
        # Select target data based on site
        if k == 1
            tbh_target = clo[!, :tbh]
            tbv_target = clo[!, :tbv]
            initial_dsnow = clo[!, :dsnow] * 0.01  # cm to m
        elseif k == 2
            tbh_target = mid1[!, :tbh]
            tbv_target = mid1[!, :tbv]
            initial_dsnow = mid1[!, :dsnow] * 0.01
        elseif k == 3
            tbh_target = mid2[!, :tbh]
            tbv_target = mid2[!, :tbv]
            initial_dsnow = mid2[!, :dsnow] * 0.01
        else
            tbh_target = far[!, :tbh]
            tbv_target = far[!, :tbv]
            initial_dsnow = far[!, :dsnow] * 0.01
        end
        
        # Optimize each measurement
        for i in eachindex(initial_dsnow)
            # Initial parameter vector
            x0 = vcat(initial_dsnow[i], rho_snow[k][i], dice[k][i], ar[k], tice[k][i], sice[k][i])
            
            # Run Nelder-Mead optimization
            result = optimize(
                x -> objective_function(x, tbh_target[i], tbv_target[i], theta, tair[k], tsnow[k], ice_perm[k], mr[k], 
                     initial_dsnow[i], sigma_dsnow[k][i], rho_snow[k][i], sigma_rhosnow[k][i], dice[k][i], 
                     sigma_dice[k][i], [ar[k]], sigma_axis_ratio[k][i], tice[k][i], sigma_tice[k][i], sice[k][i], sigma_sice[k][i]),
                x0,
                NelderMead())
            
            # Store results
            opt_dsnow = [opt_dsnow; result.minimizer[1]]
            opt_rho_snow = [opt_rho_snow; result.minimizer[2]]
            opt_dice = [opt_dice; result.minimizer[3]]
            opt_axis_ratio = [opt_axis_ratio; result.minimizer[4]]
            opt_tice = [opt_tice; result.minimizer[5]]
            opt_sice = [opt_sice; result.minimizer[6]]
            opt_tb = incoherent_model(theta, result.minimizer[1], tsnow[k][1], result.minimizer[2], 
                     result.minimizer[3], result.minimizer[5], result.minimizer[6], ice_perm[k], mr[k], result.minimizer[4])
            opt_tbh = [opt_tbh; opt_tb[1]]
            opt_tbv = [opt_tbv; opt_tb[2]]
            costs = [costs; result.minimum]
            println("Optimization finished, nº: $i")
        end

        # Save results
        opt_results = DataFrame(
            opt_dsnow = opt_dsnow, opt_rho_snow = opt_rho_snow, opt_dice = opt_dice,
            opt_axis_ratio = opt_axis_ratio, opt_tice = opt_tice, opt_sice = opt_sice,
            opt_tbh = opt_tbh, opt_tbv = opt_tbv,
            target_tbh = tbh_target, target_tbv = tbv_target, cost = costs)
        
        filenames = ["./opt_inc/opt_res_clo_jl.csv", "./opt_inc/opt_res_mid1_jl.csv", 
                     "./opt_inc/opt_res_mid2_jl.csv", "./opt_inc/opt_res_far_jl.csv"]
        CSV.write(filenames[k], opt_results)
    end
    return
end

#----------------------------------------------------------------------------------------
function main()
    data = CSV.read("./chars_data_filtered.csv", DataFrame)
    minimization(data)
end
#----------------------------------------------------------------------------------------

main()
