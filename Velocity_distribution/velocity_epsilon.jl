using CSV, DataFrames, Interpolations, Plots, QuadGK
using Statistics
using LaTeXStrings
using Plots.PlotMeasures


function plotstyle()
    display(plot!(
            guidefont=font(24), left_margin=35px, bottom_margin=20px, 
            right_margin=20px,
            size=(1000,800),frame=:box, 
            xtickfontsize=18, ytickfontsize=18,
            legendfontsize=14
        ))
end

vdata = CSV.read("v_vs_theta.csv", DataFrame)
mdata = CSV.read("fM_vs_theta.csv", DataFrame)

θv = reverse(deg2rad.(90 .- vdata.theta))
v   = reverse(vdata.v)

θm = reverse(deg2rad.(90 .- mdata.theta))
fM  = reverse(mdata.fM)*length(θm)/(π/2) # Normalize to the number of points

v_interp  = LinearInterpolation(θv, v, extrapolation_bc=Line())
fM_interp = LinearInterpolation(θm, fM, extrapolation_bc=Line())

θ_dense = range(minimum(θv), maximum(θv), length=500)
v_vals = v_interp.(θ_dense)
fM_vals = fM_interp.(θ_dense)

plot(θ_dense, v_vals, label="v(θ)", xlabel="θ [rad]", ylabel="v", lw=2)
plot!(θ_dense, fM_vals*1000, label="f(θ)", lw=2)

scatter(v_vals, fM_vals ./ sin.(θ_dense), label="f(θ)/sin(θ) vs v", xlabel="v", ylabel="f(θ)/sin(θ)", lw=2, ylims=(0,0.3))

int1 = 2*quadgk(θ -> v_interp(θ) * sin(θ), minimum(θv), maximum(θv))[1]

int2 = 2*quadgk(θ -> v_interp(0)^2 *sin(θ) / v_interp(θ), minimum(θv), maximum(θv))[1]

int2/int1

length(fM)

int_f = quadgk(θ -> fM_interp(θ), 0.0 , π/2 )[1]*length(fM)/(π/2)

int3 = quadgk(θ -> v_interp(0)^2 * fM_interp(θ) / v_interp(θ)^2, 0,π/2)[1] *length(fM)/(π/2)

int_f / int3
int3 / int_f

yr = 3.2e7
v0 = 650e5 #cm/s
NA = 6.022e23 #mol^-1
X = 0.55
Msun = 1.9885e33 
10.0^23.0*4*π*(14*yr*v0)*(14*yr*v0)/NA/X/Msun *0.212

1.0*10^23 

θ_cut = θ_dense[θ_dense .> 0.1]
v_cut = v_vals[θ_dense .> 0.1]
fM_cut = fM_vals[θ_dense .> 0.1]



plot(θ_cut, v_cut .^2 , label="v(θ)", xlabel="θ [rad]", ylabel="v", lw=2)

val0 = mean(fM_cut ./ sin.(θ_cut) ./ v_cut.^2)
plot(θ_cut, fM_cut ./ sin.(θ_cut) ./ v_cut.^2 ./ val0, label="f(θ)", lw=2)

plot(θ_dense, v_vals.^2 .* fM_vals./sin.(θ_dense), label="f(θ)", lw=2)

val0 = mean(fM_cut ./ sin.(θ_cut) ./ v_cut)
plot(θ_cut, fM_cut ./ sin.(θ_cut) ./ v_cut ./ val0, label="f(θ)", lw=2)

10^14/(3.2e7*650e5)

mH = 1.67e-24 # g
#stephan bolzmann constant
σ = 5.67e-5 # erg cm^-2 K^-4 s^-1
#mass of the star
M = 150 * Msun # g
#Eddington luminosity of the star
# L_Edd = 4πGMc/σ
σT = 6.6524e-25 # cm^2
# Thomson opacity 
κ = σT / mH /(1+ X) # cm^2/g
G = 6.67430e-8 # cm^3 g^-1 s^-
L_Edd = 4π * G * M * 3e10 / κ # erg/s

# Radius for a given effective temperature
Teff = 5500 # K
R = (5*L_Edd / (4π * σ * Teff^4))^(1/2) # cm

# Time to reach this distance at the given velocity 650 km/s
t = R / (650e5) # s
t / (365.25 * 24 * 3600)*12 # months
zo

begin
val0 = (fM_cut ./ sin.(θ_cut) ./ v_cut.^2)[1]
plot(θ_cut, fM_cut ./ sin.(θ_cut) ./ v_cut.^2 / val0, label=L"f_M(\theta)/f_v^2(\theta)", xlabel=L"\theta", ylabel=L"f's", lw=2)
val0 = (fM_cut ./ sin.(θ_cut) )[1]
plot!(θ_cut, fM_cut ./ sin.(θ_cut) / val0, label=L"f_M(\theta)", xlabel=L"\theta", ylabel=L"f's", lw=2)
val0 = v_cut[1]
plot!(θ_cut, v_cut / val0, label=L"f_v(\theta)", xlabel=L"\theta", ylabel=L"f's", lw=2)

plotstyle()
end
savefig("fM_vs_v2.pdf")