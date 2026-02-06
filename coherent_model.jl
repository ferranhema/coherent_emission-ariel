#---------------------------------------------------------------------------!
using LinearAlgebra
include("permittivity_models.jl")
include("ellipsoids_perm.jl")
#---------------------------------------------------------------------------!
function hpld(n, sn, ths, theta, xlam)
	"""
	!-- CALCULATE SURFACE REFLECTIVITY AND LAYER ABSORPTION FRACTIONS ----------!
	!----------- AT HORIZONTAL POLARIZATION ------------------------------------!
	!                                                                           !            
	! Input variables:                                                          !
	!   Name    Description                                                     !
	!                                                                           !    
	!   theta   incidence angle [rad]                                           !    
	!                                                                           !    
	!   n       number of layers (= number of soil layers plus 1, layer 1 = air)!
	!                                                                           !    
	!   sn      index of refraction for air (sn[1]) and each soil layer         !
	!                                                                           !    
	!   ths     layer thicknesses [m] (the first and last layer are assumed to  !
	!            be semi-infinite, their thickness value is ignored)            !
	!                                                                           !        
	!   xlam    radiation wave length [m]                                       !    
	!---------------------------------------------------------------------------!
	"""
	s = sin(theta)
	sp = zeros(Complex{Float64}, n)
	sp[1] = 1 + 0im
	nmax = 1
	for i in 2 : n - 1
		nmax = i + 1
		ss = sn[1] * (s / sn[i])
		sc = sqrt(1 + 0im - ss^2)
		arg = ths[i] * 2 * π / xlam
		sarg = 2 * arg * sn[i] * sc * im
		sp[i] = exp(sarg) * sp[i - 1]
		if abs(sp[i]) < 0.0001
			break
		end
	end
	sep = zeros(Complex{Float64}, n)
	sep[nmax] = 1 + 0im
	sem = zeros(Complex{Float64}, n)
	sem[nmax] = 0 + 0im
	for jj in 2 : nmax
		j = (nmax + 1) - jj
		ssj = sn[1] * s / sn[j]
		scj = sqrt(1 + 0im - ssj^2)
		ssjp1 = sn[1] * s / sn[j + 1]
		scjp1 = sqrt(1 + 0im - ssjp1^2)
		sa = 2 * sn[j] * scj / (sn[j] * scj + sn[j + 1] * scjp1)
		sb = (sn[j] * scj - sn[j + 1] * scjp1) / ((sn[j] * scj + sn[j + 1] * scjp1) * sp[j])
		sep[j] = sep[j + 1] / sa + sb * sem[j + 1] / sa
		sem[j] = sem[j + 1] + (sep[j + 1] - sep[j]) * sp[j]
	end
	sx = sep[1]
	for j in 1 : nmax
		sep[j] = sep[j] / sx
		sem[j] = sem[j] / sx
	end
	for j in nmax : n
		sp[j] = 0 + 1im / 1e50
	end
	for jj in 1 : nmax - 1
		j = (nmax + 1) - jj
		ss = sin(theta) / sn[j]
		sc = sqrt(1 + 0im - ss^2)
		r = abs(sp[j])
		s = abs(sp[j - 1])
		e2 = (s - r) * abs(sep[j])^2 + (1 / r - 1 / s) * abs(sem[j])^2
		dp = e2 * real(sn[j] * sc) / cos(theta)
		sxp = sep[j] * conj(sem[j])
		x = 2 * imag(sn[j] * sc / cos(theta)) * (imag(sxp * sp[j - 1] / abs(sp[j - 1])) - imag(sxp * sp[j] / abs(sp[j])))
		dp = dp - x
		sp[j] = complex(dp, 0)
	end
	r = abs(sem[1])^2 * real(sn[1])
	sp[1] = complex(r, 0)
	return real(sp)
end
#---------------------------------------------------------------------------!
function vpld(n, sn, ths, theta, xlam)
	"""
	!-- CALCULATE SURFACE REFLECTIVITY AND LAYER ABSORPTION FRACTIONS ----------!
	!----------- AT VERTICAL POLARIZATION ------------------------------------!
	!                                                                           !            
	! Input variables:                                                          !
	!   Name    Description                                                     !
	!                                                                           !    
	!   theta   incidence angle [rad]                                           !    
	!                                                                           !    
	!   n       number of layers (= number of soil layers plus 1, layer 1 = air)!
	!                                                                           !    
	!   sn      index of refraction for air (sn[1]) and each soil layer         !
	!                                                                           !    
	!   ths     layer thicknesses [m] (the first and last layer are assumed to  !
	!            be semi-infinite, their thickness value is ignored)            !
	!                                                                           !        
	!   xlam    radiation wave length [m]                                       !    
	!---------------------------------------------------------------------------!
	"""
	s = sin(theta)
	sp = zeros(Complex{Float64}, n)
	sp[1] = 1 + 0im
	nmax = 1
	for i in 2 : n - 1
		nmax = i + 1
		ss = sn[1] * (s/sn[i])
		sc = sqrt((1 + 0im) - ss^2)
		arg = ths[i] * 2 * π / xlam
		sarg = 2 * arg * sn[i] * sc * (0 + 1im)
		sp[i] = exp(sarg) * sp[i - 1]
		if abs(sp[i]) < 0.0001
			break
		end
	end
	sep = zeros(Complex{Float64}, n)
	sep[nmax] = 1 + 0im
	sem = zeros(Complex{Float64}, n)
	sem[nmax] = 0 + 0im
	for jj in 2 : nmax
		j = (nmax + 1) - jj 
		ssj = sn[1] * s / sn[j]
		scj = sqrt((1 + 0im) - ssj^2)
		ssjp1 = sn[1] * s / sn[j + 1]
		scjp1 = sqrt((1 + 0im) - ssjp1^2)
		sd = 2 * sn[j] * scj
		sa = sn[j] * scjp1 + sn[j + 1] * scj
		sb = sn[j] * scjp1 - sn[j + 1] * scj
		sep[j] = sa * sep[j + 1] / sd + sb * sem[j + 1] / (sd * sp[j])
		sr = sn[j + 1] / sn[j]
		sem[j] = sr * sem[j + 1] + (sep[j] - sep[j + 1] * sr) * sp[j]
	end
	sx = sep[1]
	for j in 1 : nmax
		sep[j] = sep[j] / sx
		sem[j] = sem[j] / sx
	end
	for j in nmax : n 
		sp[j] = (0 + 1im) / 1e50
	end
	for jj in 1 : nmax - 1
		j = (nmax + 1) - jj
		ss = sin(theta) / sn[j]
		sc = sqrt((1 + 0im) - ss^2)
		r = abs(sp[j])
		s = abs(sp[j - 1])
		e2 = (s - r) * abs(sep[j])^2 + (1 / r - 1 / s) * abs(sem[j])^2
		dp = e2 * real(sn[j] * sc) / cos(theta)
		sxp = sep[j] * conj(sem[j])
		x = 2 * imag(sn[j] * sc / cos(theta)) * (imag(sxp * sp[j - 1] / abs(sp[j - 1])) - imag(sxp * sp[j] / abs(sp[j])))
		dp = dp - x
		sp[j] = complex(dp,0)
	end
	r = abs(sem[1])^2 * real(sn[1])
	sp[1] = complex(r, 0)
	return real(sp)
end
#---------------------------------------------------------------------------!
function antenna_pattern_integration(mu, sigma_h, sigma_v)
	x = LinRange(-89.99, 89.99, 1000)
	dx = x[2] - x[1]
	f_h = (1 ./ sqrt.(2 * π * sigma_h^2)) .* exp.(- (x .- mu).^2 ./ (2 * sigma_h^2))
	f_v = (1 ./ sqrt.(2 * π * sigma_v^2)) .* exp.(- (x .- mu).^2 ./ (2 * sigma_v^2))
	return dx * f_h, dx * f_v
end
#---------------------------------------------------------------------------!
function wilheit_model_integrated(theta, dsnow, tsnow, rho_snow, dice, tice, sice, ice_permittivity, mixing_ratio, axis_ratio)
    # Radiative Transfer in a Plane Stratified Dielectric, Wilheit 1978
    sw = 33
    tw = -1.8
    freq = 1.4e9
    c = 3e8
    wl = c / freq
    w_h, w_v = antenna_pattern_integration(theta, 15.29, 14.87)
    thetav = range(-89.99, stop=89.99, length=1000)
    tbhs = []
    tbvs = []
    for theta0 in thetav
        theta0 = deg2rad(theta0)
        n_snow = length(dsnow)
        n_ice = length(dice)

        epre_snow = []
        epim_snow = []
        for i in 1:n_snow
            epre_snow0, epim_snow0 = epsilon_snow(freq, rho_snow[i], tsnow[i])
            push!(epre_snow, epre_snow0)
            push!(epim_snow, epim_snow0)
        end

        epre_ice = []
        epim_ice = []
        for i in 1:n_ice
            if ice_permittivity == "vant"
                epre_ice0, epim_ice0 = epsilon_Vant_ice(sice[i], tice[i])
                push!(epre_ice, epre_ice0)
                push!(epim_ice, epim_ice0)
            elseif ice_permittivity == "rn"
                epre_ice0, epim_ice0 = epsilon_rn_ice(freq, sice[i], tice[i])
                push!(epre_ice, epre_ice0)
                push!(epim_ice, epim_ice0)
            elseif ice_permittivity == "sp"
                epre_ice0, epim_ice0 = epsilon_sp_ice(freq, sice[i], tice[i])
                push!(epre_ice, epre_ice0)
                push!(epim_ice, epim_ice0)
            elseif ice_permittivity == "pvs"
                epre_ice0, epim_ice0 = pvs_mixing(mixing_ratio, freq, sice[i], tice[i])
                push!(epre_ice, epre_ice0)
                push!(epim_ice, epim_ice0)
            elseif ice_permittivity == "mix" # The 2nd 0 corresponds to the formulation: 0 = Maxwell-Garnet ; 1 = Coherent potential ; 2 = Polder van Santen
                epre_ice0, epim_ice0 = e_eff_mix(axis_ratio, 0, 0, sice[i], tice[i])
                push!(epre_ice, epre_ice0)
                push!(epim_ice, epim_ice0)
            end
        end
        epre_water, epim_water = epsilon_water(freq, tw, sw)

        sn = [sqrt(complex(1, 0))]
        for i in 1:n_snow
            push!(sn, sqrt(complex(epre_snow[i], epim_snow[i])))
        end
        for i in 1:n_ice
            push!(sn, sqrt(complex(epre_ice[i], epim_ice[i])))
        end
        push!(sn, sqrt(complex(epre_water, epim_water)))

        layers = [1e50]
        for i in 1:n_snow
            push!(layers, dsnow[i])
        end
        for i in 1:n_ice
            push!(layers, dice[i])
        end
        push!(layers, 1e50)

        sp_h = hpld(n_snow + n_ice + 2, sn, layers, theta0, wl)
        sp_v = vpld(n_snow + n_ice + 2, sn, layers, theta0, wl)

        tbh = sp_h[1]
        j = 2
        for i in 1:n_snow
            tbh += sp_h[j] * (tsnow[i] + 273.15)
            j += 1
        end
        for i in 1:n_ice
            tbh += sp_h[j] * (tice[i] + 273.15)
            j += 1
        end
        tbh += last(sp_h) * (tw + 273.15)
        
        tbv = sp_v[1]
        j = 2
        for i in 1:n_snow
            tbv += sp_v[j] * (tsnow[i] + 273.15)
            j += 1
        end
        for i in 1:n_ice
            tbv += sp_v[j] * (tice[i] + 273.15)
            j += 1
        end
        tbv += last(sp_v) * (tw + 273.15)

        push!(tbhs, tbh)
        push!(tbvs, tbv)
    end
    tbh = sum(tbhs .* w_h)
    tbv = sum(tbvs .* w_v)

    return tbh, tbv
end
#---------------------------------------------------------------------------!