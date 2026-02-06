#---------------------------------------------------------------------------!
function epsilon_Vant_ice(S, T)
    rho_ice = 0.917 - 0.1404 * (10^-3) * T # in Mg/m^3, density of pure ice from Pounder, 1965
    # coefficients from Cox and Weeks, 1983
    if T > 0
        a1 = -0.041221
        b1 = -18.407
        c1 = 0.58402
        d1 = 0.21454
        a2 = 0.090312
        b2 = -0.016111
        c2 = 1.2291 * (10^-4)
        d2 = 1.3603 * (10^-4)
    elseif T <= 0 && T > -2
        a1 = -0.041221
        b1 = -18.407
        c1 = 0.58402
        d1 = 0.21454
        a2 = 0.090312
        b2 = -0.016111
        c2 = 1.2291 * (10^-4)
        d2 = 1.3603 * (10^-4)
    elseif T <= -2 && T >= -22.9
        a1 = -4.732
        b1 = -22.45
        c1 = -0.6397
        d1 = -0.01074
        a2 = 0.08903
        b2 = -0.01763
        c2 = -5.330 * (10^-4)
        d2 = -8.801 * (10^-6)
    end
    c1 = 55.27
    d1 = 0.7160
    a2 = 8.547
    b2 = 1.089
    c2 = 0.04518
    d2 = 5.819 * (10^-4)
    F1 = a1 + b1 * T + c1 * (T^2) + d1 * (T^3)
    F2 = a2 + b2 * T + c2 * (T^2) + d2 * (T^3)
    Vb = (rho_ice * S) / (F1 - rho_ice * S * F2)
    EpRe = 3.1 + 0.0084 * Vb * (10^3)
    EpIm = 0.037 + 0.00445 * Vb * (10^3)
    return EpRe, EpIm
end
#---------------------------------------------------------------------------!
function epsilon_snow(freq, rho_snow, T)
    # freq in GHz    
    # rho_snow in kg/m^3
    # T in K
    T = T + 273.15 # Convert to Kelvin
    freq = freq / 1e9 # Convert to GHz
    if rho_snow <= 400
        rho_snow = rho_snow / 1e3 # Convert to Mg/m^3
        EpRe = 1 + 1.599 * rho_snow + 1.861 * rho_snow^3
    else
        v = rho_snow / 917 
        e_h = 1
        e_s = 3.215
        EpRe = (1 - v) * e_h + v * e_s
        #rho_snow = rho_snow / 1e3 # Convert to Mg/m^3
    end
    B1 = 0.0207
    b = 335
    B2 = 1.16e-11
    beta = (B1 / T) * (exp(b / T) / (exp(b / T) - 1))^2 + B2 * freq^2 + exp(-9.963 + 0.0372 * (T - 273.16))
    theta = (300 / T) - 1
    alpha = (0.00504 + 0.0062) * exp(-22.1 * theta)
    e_i = (alpha / freq) + beta * freq    
    EpIm = e_i * (0.52 * rho_snow + 0.62 * rho_snow^2)
    return EpRe, EpIm
end
#---------------------------------------------------------------------------!
function epsilon_water(freq, T, S)
    """
    freq in Hz
    T in ºC
    S in PSU
    """
    E0 = 8.854 * 10^(-12.0)
    Eswinf = 4.9
    Esw00 = 87.174 - 1.949 * 10^(-1.0) * T - 1.279 * 10^(-2.0) * (T^2.0) + 2.491 * 10^(-4.0) * (T^3.0)
    a = 1.0 + 1.613 * 10^(-5.0) * T * S - 3.656 * 10^(-3.0) * S + 3.21 * 10^(-5.0) * (S^2.0) - 4.232 * 10^(-7.0) * (S^3.0)
    Esw0 = Esw00 * a
    Tausw0 = 1.1109 * 10^(-10.0) - 3.824 * 10^(-12.0) * T + 6.238 * 10^(-14.0) * (T^2.0) - 5.096 * 10^(-16.0) * (T^3.0)
    b = 1.0 + 2.282 * 10^(-5.0) * T * S - 7.638 * 10^(-4.0) * S - 7.760 * 10^(-6.0) * (S^2.0) + 1.105 * 10^(-8.0) * (S^3.0)
    Tausw = (Tausw0 / (2 * π)) * b
    EpRe = Eswinf + ((Esw0 - Eswinf) / (1.0 + ((2 * π * freq * Tausw)^2.0)))
    si1 = S * (0.18252 - 1.4619 * 10^(-3.0) * S + 2.093 * 10^(-5.0) * (S^2.0) - 1.282 * 10^(-7.0) * (S^3.0))
    inc = 25.0 - T
    phi = inc * (2.033 * 10^(-2.0) + 1.266 * 10^(-4.0) * inc + 2.464 * 10^(-6.0) * (inc^2.0) - S * (1.849 * 10^(-5.0) - 2.551 * 10^(-7.0) * inc + 2.551 * 10^(-8.0) * (inc^2.0)))
    si = si1 * exp(-phi)
    Epimn = (2 * π * freq * Tausw * (Esw0 - Eswinf))
    Epimd = 1.0 + ((2 * π * freq * Tausw)^2.0)
    C = (si / (2 * π * E0 * freq))
    EpIm = ((Epimn / Epimd) + C)
return EpRe, EpIm
end
#---------------------------------------------------------------------------!
function epsilon_rn_ice(freq, S, T)
    # Coefficients from Assur 1960
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

    Sb = a + b * T + c * (T^2) + d * (T^3)
    Nb = 1.707 * 10^(-2) * Sb + 1.205 * 10^(-5) * (Sb^2) + 4.058 * 10^(-9) * (Sb^3)
    Nb = Nb * 0.9141 # correction for sea water

    # Stogryn and Desargant [1985]
    epsiwoo = 4.9
    e0 = 8.854e-12
    f = freq
    eps = 88.045 - 0.4147 * T + 6.295e-4 * (T^2) + 1.075e-5 * (T^3)
    a = 1 - 0.255 * Nb + 5.15e-2 * (Nb^2) - 6.89e-3 * (Nb^3)
    eb0 = eps * a
    rel = (1.1109e-10) - (3.824e-12) * T + (6.938e-14) * (T^2) - (5.096e-16) * (T^3)
    b = 1 + (0.146e-2) * T * Nb - (4.89e-2) * Nb - (2.97e-2) * (Nb^2) + (5.64e-3) * (Nb^3)
    relax = rel * b
    D = 25 - T
    sig = Nb * (10.39 - 2.378 * Nb + 0.683 * (Nb^2) - 0.135 * (Nb^3) + (1.01e-2) * (Nb^4))
    c = 1 - (1.96e-2) * D + (8.08e-5) * (D^2) - Nb * D * ((3.02e-5) + (3.92e-5) * D + Nb * ((1.72e-5) - (6.58e-6) * D))
    conb = c * sig
    eb = epsiwoo + ((eb0 - epsiwoo) / (1 + ((relax * f)^2)))
    ebi = ((relax * f * (eb0 - epsiwoo)) / (1 + ((relax * f)^2))) + (conb / (6.28 * f * e0))
    EpReb = eb
    EpImb = ebi

    rho_ice = 0.917 - 0.1404 * 10^(-3) * T # in Mg/m^3
    # coefficients from Cox and Weeks, 1983
    if T > 0
        a1 = -0.041221
        b1 = -18.407
        c1 = 0.58402
        d1 = 0.21454
        a2 = 0.090312
        b2 = -0.016111
        c2 = 1.2291 * 10^(-4)
        d2 = 1.3603 * 10^(-4)
    elseif T <= 0 && T > -2
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
    elseif T < -22.9 && T >= -30
        a1 = 9899
        b1 = 1309
        c1 = 55.27
        d1 = 0.7160
        a2 = 8.547
        b2 = 1.089
        c2 = 0.04518
        d2 = 5.819 * 10^(-4)
    end

    F1 = a1 + b1 * T + c1 * (T^2) + d1 * (T^3)
    F2 = a2 + b2 * T + c2 * (T^2) + d2 * (T^3)
    vi = (rho_ice * S) / (F1 - rho_ice * S * F2)

    epi = 3.1884 + 9.1e-4 * T
    theta = 300 / (T + 273.15) - 1
    alpha = (0.00504 + 0.0062 * theta) * exp(-22.1 * theta)
    beta = ((0.502 - 0.131 * theta) / (1 + theta)) * 1e-4 + (0.542e-6 * ((1 + theta) / (theta + 0.0073))^2)
    epii = (alpha / (freq / 1e9)) + (beta * (freq / 1e9))
    eice = complex(epi, epii)
    ei = complex(EpReb, EpImb)
    eh = eice
    emi = eh

    # Shokr 1998
    m = 0
    measureA = 1.2
    measureB = 1.2
    while measureA > 0.001 || measureB > 0.001
        f1 = (ei - eh) / (ei + emi)
        f2 = 5 * emi + ei
        est = eh + ((vi / 3) * f1 * f2)
        measureA = abs(real(emi) - real(est))
        measureB = abs(imag(emi) - imag(est))
        emi = est
        m += 1
        if m > 30
            break
        end
    end

    EpRem = real(emi)
    EpImm = imag(emi)
    return EpRem, EpImm
end
#---------------------------------------------------------------------------!
function epsilon_sp_ice(freq, S, T)
    """
    Parameters
    ----------
    Ep_i : permittivity of the inclusion
    Ep_h : permittivity of the host
    Vi : volume fraction of the inclusions
    Returns
    -------
    Ep_m : dielectric constant of the mixture
    """
    # Coefficients from Assur 1960
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

    Sb = a + b * T + c * (T^2) + d * (T^3)
    Nb = (1.707e-2 * Sb) + (1.205e-5 * Sb^2) + (4.058e-9 * Sb^3)
    Nb = Nb * 0.9141 # correction for sea water

    # Stogryn and Desargant [1985]
    epsiwoo = 4.9
    e0 = 8.854e-12
    f = freq
    eps = 88.045 - 0.4147 * T + 6.295e-4 * (T^2) + 1.075e-5 * (T^3)
    a = 1 - 0.255 * Nb + 5.15e-2 * (Nb^2) - 6.89e-3 * (Nb^3)
    eb0 = eps * a
    rel = (1.1109e-10) - (3.824e-12) * T + (6.938e-14) * (T^2) - (5.096e-16) * (T^3)
    b = 1 + (0.146e-2) * T * Nb - (4.89e-2) * Nb - (2.97e-2) * (Nb^2) + (5.64e-3) * (Nb^3)
    relax = rel * b
    D = 25 - T
    sig = Nb * (10.39 - 2.378 * Nb + 0.683 * (Nb^2) - 0.135 * (Nb^3) + (1.01e-2) * (Nb^4))
    c = 1 - (1.96e-2) * D + (8.08e-5) * (D^2) - Nb * D * ((3.02e-5) + (3.92e-5) * D + Nb * ((1.72e-5) - (6.58e-6) * D))
    conb = c * sig
    eb = epsiwoo + ((eb0 - epsiwoo) / (1 + ((relax * f)^2)))
    ebi = ((relax * f * (eb0 - epsiwoo)) / (1 + ((relax * f)^2))) + (conb / (6.28 * f * e0))
    EpReb = eb
    EpImb = ebi

    rho_ice = 0.917 - 0.1404 * 10^(-3) * T # in Mg/m^3, density of pure ice from Pounder, 1965
    # coefficients from Cox and Weeks, 1983
    if T > 0
        a1 = -0.041221
        b1 = -18.407
        c1 = 0.58402
        d1 = 0.21454
        a2 = 0.090312
        b2 = -0.016111
        c2 = 1.2291 * 10^(-4)
        d2 = 1.3603 * 10^(-4)
    elseif T <= 0 && T > -2
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
    elseif T < -22.9 && T >= -30
        a1 = 9899
        b1 = 1309
        c1 = 55.27
        d1 = 0.7160
        a2 = 8.547
        b2 = 1.089
        c2 = 0.04518
        d2 = 5.819 * 10^(-4)
    end

    F1 = a1 + b1 * T + c1 * (T^2) + d1 * (T^3)
    F2 = a2 + b2 * T + c2 * (T^2) + d2 * (T^3)
    vi = (rho_ice * S) / (F1 - rho_ice * S * F2)

    # Shokr 1998
    epi = 3.1884 + 9.1e-4 * T
    theta = 300 / (T + 273.15) - 1
    alpha = (0.00504 + 0.0062 * theta) * exp(-22.1 * theta)
    beta = ((0.502 - 0.131 * theta) / (1 + theta)) * 1e-4 + (0.542e-6 * ((1 + theta) / (theta + 0.0073))^2)
    epii = (alpha / (freq / 1e9)) + (beta * (freq / 1e9))
    eice = complex(epi, epii)
    ei = complex(EpReb, EpImb)
    eh = eice
    emi = eh
    m = 0
    measureA = 1.2
    measureB = 1.2
    while measureA > 0.001 || measureB > 0.001
        f1 = (ei - eh)
        f2 = (ei + 2 * emi)
        est = eh + (3 * vi * emi * (f1 / f2))
        measureA = abs(real(emi) - real(est))
        measureB = abs(imag(emi) - imag(est))
        emi = est
        m += 1
        if m > 30
            break
        end
    end

    EpRem = real(emi)
    EpImm = imag(emi)
    return EpRem, EpImm
end
#---------------------------------------------------------------------------!
function pvs_mixing(mixing_ratio, freq, S, T)
    """
    The mixing_ratio corresponds to the fraction, from 0 to 1, of random needles inclusions, thus the spheres correspond to 1 - mixing_ratio
    """
    EpRe_rn, EpIm_rn = epsilon_rn_ice(freq, S, T)
    EpRe_sp, EpIm_sp = epsilon_sp_ice(freq, S, T)
    EpRem = EpRe_rn * mixing_ratio + EpRe_sp * (1 - mixing_ratio)
    EpImm = EpIm_rn * mixing_ratio + EpIm_sp * (1 - mixing_ratio)
    return EpRem, EpImm
end
#---------------------------------------------------------------------------!