#---------------------------------------------------------------------------!
function depolarization_factor(aspect_ratio, c)
    Na = 1 / (1 + 1.6 * aspect_ratio + 0.4 * aspect_ratio^2)
    Nb = 0.5 * (1 - Na)
    Nc = Nb
    return Na, Nb, Nc
end
#---------------------------------------------------------------------------!
function e_ice_brine_vi(freq, S, T)
    if T <= 0 && T >= -8.2
        a = 1.725
        b = -18.756
        c = -0.3964
        d = 0
    elseif T <= -8.2 && T >= -22.9
        a = 57.041
        b = -9.929
        c = -0.16204
        d = -0.002396
    elseif T <= -22.9 && T >= -36.8
        a = 242.94
        b = 1.5299
        c = 0.0429
        d = 0
    elseif T <= -36.8 && T >= -43.2
        a = 508.18
        b = 14.535
        c = 0.2018
        d = 0
    end
    Sb = a + b * T + c * T^2 + d * T^3
    Nb = 1.707 * 10^(-2) * Sb + 1.205 * 10^(-5) * Sb^2 + 4.058 * 10^(-9) * Sb^3
    Nb = Nb * 0.9141 # correction for sea water
    # Stogryn and Desargant [1985]
    epsiwoo = 4.9
    e0 = 8.854e-12
    f = freq
    eps = 88.045 - 0.4147 * T + 6.295e-4 * T^2 + 1.075e-5 * T^3
    a = 1 - 0.255 * Nb + 5.15e-2 * Nb^2 - 6.89e-3 * Nb^3
    eb0 = eps * a
    rel = (1.1109e-10) - (3.824e-12) * T + (6.938e-14) * T^2 - (5.096e-16) * T^3
    b = 1 + (0.146e-2) * T * Nb - (4.89e-2) * Nb - (2.97e-2) * Nb^2 + (5.64e-3) * Nb^3
    relax = rel * b
    D = 25 - T
    sig = Nb * (10.39 - 2.378 * Nb + 0.683 * Nb^2 - 0.135 * Nb^3 + (1.01e-2) * Nb^4)
    c = 1 - (1.96e-2) * D + (8.08e-5) * D^2 - Nb * D * ((3.02e-5) + (3.92e-5) * D + Nb * ((1.72e-5) - (6.58e-6) * D))
    conb = c * sig
    eb = epsiwoo + ((eb0 - epsiwoo) / (1 + (relax * f)^2))
    ebi = ((relax * f * (eb0 - epsiwoo)) / (1 + (relax * f)^2)) + (conb / (6.28 * f * e0))
    EpReb = eb
    EpImb = ebi
    rho_ice = 0.917 - 0.1404 * 10^(-3) * T # in Mg/m^3, density of pure ice from Pounder, 1965
    # coefficients from Cox and Weeks, 1983
    if T <= 0 && T > -2
        a1 = -0.041221
        b1 = -18.407
        c1 = 0.58402
        d1 = 0.21454
        a2 = 0.090312
        b2 = -0.016111
        c2 = 1.2291 * 10^(-4)
        d2 = 1.3603 * 10^(-4)
    elseif T <= -2 && T >= -22.9
        a1 = -4.732
        b1 = -22.45
        c1 = -0.6397
        d1 = -0.01074
        a2 = 0.08903
        b2 = -0.01763
        c2 = -5.330 * 10^(-4)
        d2 = -8.801 * 10^(-6)
    elseif T < -22.9
        a1 = 9899
        b1 = 1309
        c1 = 55.27
        d1 = 0.7160
        a2 = 8.547
        b2 = 1.089
        c2 = 0.04518
        d2 = 5.819 * 10^(-4)
    end
    F1 = a1 + b1 * T + c1 * T^2 + d1 * T^3
    F2 = a2 + b2 * T + c2 * T^2 + d2 * T^3
    vi = (rho_ice * S) / (F1 - rho_ice * S * F2)
    # Shokr 1998
    epi = 3.1884 + 9.1e-4 * T
    theta = 300 / (T + 273.15) - 1
    alpha = (0.00504 + 0.0062 * theta) * exp(-22.1 * theta)
    beta = ((0.502 - 0.131 * theta) / (1 + theta)) * 1e-4 + (0.542e-6 * ((1 + theta) / (theta + 0.0073))^2)
    epii = (alpha / (freq / 1e9)) + (beta * (freq / 1e9))
    eice = complex(epi, epii)
    ei = complex(EpReb, EpImb)
    return vi, eice, ei
end
#---------------------------------------------------------------------------!
function e_eff_i(e_eff, v, Ni, S, T)
    vi, e0, e1 = e_ice_brine_vi(1.4e9, S, T)
    e_eff_i1 = (vi * (e1 - e0) * (e0 + v * (e_eff - e0))) / (3 * (e0 + v * (e_eff - e0) + Ni * (e1 - e0)))
    e_eff_i2 = (vi * Ni * (e1 - e0)) / (3 * (e0 + v * (e_eff - e0) + Ni * (e1 - e0)))
    return e_eff_i1, e_eff_i2
end
 #---------------------------------------------------------------------------! 
function e_eff_ps(e_eff, Na, Nb, Nc, S, T)
    vi, e0, e1 = e_ice_brine_vi(1.4e9, S, T)
    va = 1 - Na
    vb = 1 - Nb
    vc = 1 - Nc
    e_eff_a1 = (vi * (e1 - e0) * (e0 + va * (e_eff - e0))) / (3 * (e0 + va * (e_eff - e0) + Na * (e1 - e0)))
    e_eff_a2 = (vi * Na * (e1 - e0)) / (3 * (e0 + va * (e_eff - e0) + Na * (e1 - e0)))
    e_eff_b1 = (vi * (e1 - e0) * (e0 + vb * (e_eff - e0))) / (3 * (e0 + vb * (e_eff - e0) + Nb * (e1 - e0)))
    e_eff_b2 = (vi * Nb * (e1 - e0)) / (3 * (e0 + vb * (e_eff - e0) + Nb * (e1 - e0)))
    e_eff_c1 = (vi * (e1 - e0) * (e0 + vc * (e_eff - e0))) / (3 * (e0 + vc * (e_eff - e0) + Nc * (e1 - e0)))
    e_eff_c2 = (vi * Nc * (e1 - e0)) / (3 * (e0 + vc * (e_eff - e0) + Nc * (e1 - e0)))
    return e_eff_a1, e_eff_a2, e_eff_b1, e_eff_b2, e_eff_c1, e_eff_c2
end
#---------------------------------------------------------------------------!
function e_eff_mix(axis_ratio, c, v, S, T)
    if v == 0 || v == 1 # Maxwell-Garnet / Coherent potential
        Na, Nb, Nc = depolarization_factor(axis_ratio, c)
        m = 0
        measureA = 1.2
        measureB = 1.2
        vi, e_eff0, e_brine = e_ice_brine_vi(1.4e9, S, T)
        e0 = e_eff0
        while measureA > 0.001 || measureB > 0.001
            e_eff_a1, e_eff_a2 = e_eff_i(e_eff0, v, Na, S, T)
            e_eff_b1, e_eff_b2 = e_eff_i(e_eff0, v, Nb, S, T)
            e_eff_c1, e_eff_c2 = e_eff_i(e_eff0, v, Nc, S, T)
            est = e0 + ((e_eff_a1 + e_eff_b1 + e_eff_c1) * (1 - (e_eff_a2 + e_eff_b2 + e_eff_c2))^-1)
            measureA = abs(real(e_eff0) - real(est))
            measureB = abs(imag(e_eff0) - imag(est))
            e_eff0 = est
            m += 1
            # if vi < 0.1
            #     break
            if m > 30
                break
            end
        end
        EpRem = real(e_eff0)
        EpImm = imag(e_eff0)
        return EpRem, EpImm
    else # Polder van Santen (v == 2)
        Na, Nb, Nc = depolarization_factor(axis_ratio, c)
        m = 0
        measureA = 1.2
        measureB = 1.2
        vi, e_eff0, e_brine = e_ice_brine_vi(1.4e9, S, T)
        e0 = e_eff0
        while measureA > 0.001 || measureB > 0.001
            e_eff_a1, e_eff_a2, e_eff_b1, e_eff_b2, e_eff_c1, e_eff_c2 = e_eff_ps(e_eff0, Na, Nb, Nc, S, T)
            est = e0 + ((e_eff_a1 + e_eff_b1 + e_eff_c1) * (1 - (e_eff_a2 + e_eff_b2 + e_eff_c2))^-1)
            measureA = abs(real(e_eff0) - real(est))
            measureB = abs(imag(e_eff0) - imag(est))
            e_eff0 = est
            m += 1
            if vi < 0.1
                break
            end
            if m > 30
                break
            end
        end
        EpRem = real(e_eff0)
        EpImm = imag(e_eff0)
        return EpRem, EpImm
    end
end
#---------------------------------------------------------------------------!