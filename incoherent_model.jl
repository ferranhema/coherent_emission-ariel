using BandedMatrices
using LinearAlgebra
include("permittivity_models.jl")
include("ellipsoids_perm.jl")
#----------------------------------------------------------------------------------------
"""
    fresnel(n1, n2, θ)

the fresnel reflection coefficient, i.e., the reflected amplitude of
the electric field for propagation from a medium with refractive index n2
into medium with refractive index n1. θ is the propagation angle in medium 1.
"""
function fresnel(n1,n2,θ)
    ctt=cos(θ)
    ct2=sqrt(1-(real(n1)*sin(θ)/real(n2))^2)

    rh=(n1*ctt -n2 * ct2)/(n1*ctt +n2 * ct2)
    rv=(n2*ctt -n1 * ct2)/(n2*ctt +n1 * ct2)
    return rh,rv
end
#----------------------------------------------------------------------------------------
"""
    fresnel(N, θ)

The fresnel reflection coefficient for propagation from the bottom
of a stack of layers with refractive indices N.
θ is the propagation angle in vacuum.
"""
function fresnel(N::AbstractArray{T},θ) where {T<: Union{AbstractFloat,Complex}}
    rhl=Array{promote_type(T,typeof(θ))}(undef,length(N)-1)
    rvl=Array{promote_type(T,typeof(θ))}(undef,length(N)-1)
    for i in 1:length(N)-1
        rhl[i],rvl[i]=fresnel(N[i],N[i+1],theta_in_medium(N[i],θ))
    end
    return rhl,rvl
end
#----------------------------------------------------------------------------------------
"""
    theta_in_medium(n, θ)

the local incidence angle in the given layer with refractive index n.
The incidence angle θ is the one in vacuum.
"""
function theta_in_medium(n, θ)
    return asin((1/real(n))*sin(θ))
end
#----------------------------------------------------------------------------------------
"""
    incoherent(da, N, θ, freq, T)
simple incoherent emission model with infinite reflections

da is the layer thickness array (air not included)
N is the refractive index array (complex, air included)
θ incidence anle in air
freq is frequency (used for absorption, phase is not relevant in inchoerent modeling)
T temperature of each layer in K (air not included)
"""
function incoherent(da,N,θ,freq,T)
    Tall=[0.0;T]
    Nall=N
    d= length(da)+2==length(N) ? da : da[2:end-1]

    θall=theta_in_medium.(Nall[1:end],θ)
    κ=imag(4 .*π .* freq ./ 0.3 .* Nall)[2:end-1]
    rs=fresnel(Nall,θ) .|> x->abs.(x).^2
    trans=exp.(-κ.*d.*sec.(θall[2:end-1]))
    out=Array{Float64}(undef,2)
    for (i,r) in enumerate(rs)
        #interpolated from MEMLS documentation for without scattering
        #changed order of layers, first layer is on top, last is on the bottom
        M4= BandedMatrix( 1=>(1 .-r[2:end-1]).*trans[1:end-1])
        M3= Diagonal(r[2:end].*trans)
        M2= Diagonal(r[1:end-1].*trans)
        M1= BandedMatrix( 0=>ones(length(trans)),  -1=>-trans[2:end] .*(1 .-r[2:end-1]))
        
        F=(1 .-trans).*Tall[2:end-1]
        F[end]+=(1-r[end])*Tall[end]*trans[end]
        E=(1 .-trans).*Tall[2:end-1]
        E[1]+=trans[1]*(1 .-r[1]).*Tall[1]

        D=(I-M3*(M1\M2)-M4)\(M3*(M1\E)+F)

        out[i]=D[1]*(1-r[1])
    end
    out
end
#----------------------------------------------------------------------------------------
function incoherent_model(theta, dsnow, tsnow, rho_snow, dice, tice, sice, ice_permittivity, mixing_ratio, axis_ratio)
    freq = 1.4
    T = [tsnow + 273.15, tice + 273.15, -2 + 273.15]
    epre_snow, epim_snow = epsilon_snow(freq*1e9, rho_snow, T[1]-273.15)
    epre_ice, epim_ice = e_eff_mix(axis_ratio,0,0,sice,T[2]-273.15)
	epre_water, epim_water = epsilon_water(freq*1e9, T[3]-273.15, 33)

    n = [sqrt(complex(1,0)), sqrt(complex(epre_snow, epim_snow)), sqrt(complex(epre_ice, epim_ice)), sqrt(complex(epre_water, epim_water))]
    da = [dsnow, dice]
    tbh, tbv = incoherent(da, n, deg2rad(theta), freq, T)
    return tbh, tbv
end