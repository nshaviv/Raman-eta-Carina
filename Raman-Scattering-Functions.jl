# Functions to calculate the Raman scattering and fit the data
# N.J. Shaviv, July 2023 (with last modifications in July 2025)

# Import the different packages we need to use
# Note, for the first time you need to download/install the packages, which you do through
# using Pkg
# Pkg.add("DataFrames") etc.

using DataFrames, CSV     # Import from CSV files, defines DataFrame structure (like excel sheets)
using Plots               # Standard plotting packages
using BenchmarkTools      # More options to time things
using Optim, ForwardDiff  # Minimizer package and differentiator 
using Interpolations      # Interpolation of functions 
using ProgressMeter       # Progress bar for long calculationusing Roots               # Find roots of functions
using Base.Threads        # To use multithreading for long calculations
using LaTeXStrings        # To use LaTeX in plots
using ColorSchemes        # To use nice color schemes in plots
using Plots.PlotMeasures  # To add measures to plots, such as the time it took to run
using PrettyTables        # To print tables in a nice format
using Measurements        # To propagate errors in calculations
using Statistics          # To calculate mean, median, etc.
using CubicSplines        # To do cubic spline interpolation of functions
using Printf              # To print formatted strings

# - Some constants
NA = 6.022e23      # Avogadro's number
Msun = 1.98e33     # Solar mass in grams
σT = 6.6e-25       # Thomson cross-section in c.g.s.

# -------------- First lets define the functions we need

# This calculates the Hα flux passed through a slab. 
# The analytic calculation was carried out with Mathematica 
# pσ is the log10 of the total cross-section in c.g.s for the Lyman photon *Raman* scattering
# pΣ is the log10 of the total column density of the shell
# BR is the branching ratio into the Balmer line
# df is the fraction of dust relative to standard ISM value
function ramancalc(pσ, pΣ, BR; df=0.01, fbound=0.05, reflect=false)
    # Background scattering cross-sections
    # fbound is the fraction of electrons
    # Note that we assume H- is negligible, which it should for these densities.
    # We also assume that the Rayleigh scattering off of n=2 electrons is negligible
    σB_back_Bal = (0.003*fbound+(1-fbound)*1.0)*σT # Ly wing @ Balmer from bound e + scattering on free electrons
    σB_back_Ly  = ((1-fbound)*1.0)*σT # scattering on free electrons
    aB = 0.75  # scattering fraction on dust (albedo) around the Balmer lines
    aL = 0.3   # scattering fraction on dust (albedo) around the Lyman lines
    Σ = 10.0^pΣ
    κB         =  df*(1-aB)*3.7e-22          # Balmer frequency absorption cross-section 
    σB         =  df*  aB  *3.7e-22+σB_back_Bal  # Balmer frequnecy scattering cross-section
    σ_raman_tot = 10.0^pσ                 # The Raman scattering cross-section for Lyman photons
    κL_dust    =  df*(1-aL)*2.2e-21          # Lyman frequency absorption cross-section
    σL_dust_th =  df*  aL  *2.2e-21+σB_back_Ly  # Lyman frequency scattering cross-section
    σ_rayleigh = (1-BR)*σ_raman_tot       # cross-section for Rayleigh scattering   
    σ_raman    =   (BR)*σ_raman_tot       # For the Lyman band, it all appears as absorption 
    κL         = σ_raman    + κL_dust        # Raman scattering is like absorption for Lyman
    σL         = σL_dust_th + σ_rayleigh
    χL = sqrt(κL*(κL+σL))
    χB = sqrt(κB*(κB+σB))
    cτL = cosh(Σ*χL) 
    sτL = sinh(Σ*χL) 
    cτB = cosh(Σ*χB) 
    sτB = sinh(Σ*χB) 
    sqκB = sqrt(κB)
    sqχB = sqrt(κB + σB)
    denom = (sqκB*κL*(sτB*(2*κB + σB) + 2*cτB*χB)*(χB - χL)*(χB + χL)*(cτL*(κL + σL) + sτL*χL))
    enum  = reflect == true ? 
       (sqχB*σ_raman*(κL + σL)*(κB*κL*(-κB + κL - σB + σL) + sτB*χB*(-(cτL*κL*(-κB + κL + σL)) +
            sτL*(κB - κL)*χL) + cτB*κB*(cτL*κL*(κB - κL + σB - σL) + sτL*(κB - κL + σB)*χL))) :
       (sqκB*sqχB*σ_raman*(κL + σL)*(sqχB*sτB*κL*(κB + κL + σL) + cτB*sqκB*κL*(κB + κL + σB + σL) 
           - cτL*sqκB*κL*(κB + κL + σB + σL) - sqκB*sτL*(κB + κL + σB)*χL))
    enum/denom
end

# Remove NaNs from results:
fixNaN(num) = isnan(num) ? 0 : num

# The expression for the pass through flux has very large exponents. The simplest solution is to 
# do the calculation with higher precision in those regions of phase space that require it.
# We use Julia's "native" arbitrary float precision option (There are packages, such as QuadMath which doesn't 
# run on ARM (new macs) or DoubleFloats, which doesn't have enough exponent bits!). 
setprecision(128)
function passthroughcombined(pσ, pΣin, BR; df=0.01, reflect=true)
    if reflect == true
        pΣ = min(pΣin,24-log10(df))
 #        if pσ + pΣ > 1.8 && pσ >-23 return ramancalc(-18.8, 20.0, BR; df=df) end 
        if pσ + pΣ < 0 return ramancalc(pσ, pΣ, BR; df=df)  end
        Float64(ramancalc(BigFloat(pσ), pΣ, BR; df=df)) |> fixNaN 
    else
        pΣ = min(pΣin,23-log10(df))
        if pσ + pΣ < 2.8 return ramancalc(pσ, pΣ, BR; df=df) |> fixNaN end
        if pσ + pΣ < 4.4 return Float64(ramancalc(BigFloat(pσ), pΣ, BR; df=df)) |> fixNaN end
        return(Float64(ramancalc(BigFloat(4.4-pΣ), pΣ, BR; df=df)) |> fixNaN)
    end
end


# simplest fit - one component

# One model parameters: Aback, Sback, Aly, pΣ, pdf
# Two model parameters: Aback, Sback, Aly, pΣ1, pΣ2, f, pdf

OneCompFit(v,Aback,Sback,Aly,pΣ,pdf) =
    Aback + (Sback*v/10000) +
    Aly*passthroughcombined(LyβCS(v),pΣ,LyβBr(v);df=10.0^pdf)

OneCompFitHβ(v,Aback,Sback,Aly,pΣ,pdf) =
    Aback + (Sback*v/10000) +
    Aly*passthroughcombined(LyγCS(v),pΣ,LyγBr(v);df=10.0^pdf)

function χ2_oneComp(x,datavel,dataflux)
    pred = map(v->OneCompFit(v,x...),datavel)
    sum(map( x->max(min(x,0.8),-0.8),(pred-dataflux).^2))+(max(25.0,x[4])-25.0)
#    sum(map( x->x,(pred-dataflux).^2))+(max(25.0,x[4])-25.0)
end

frac(f) = (1+tanh(f))/2.0

#TwoCompFit(v,Aback,Sback,Aly,pΣ1,pΣ2,f,pdf) =
#    Aback + (Sback*v/10000) +
#    frac(f)*Aly*passthroughcombined(LyβCS(v),pΣ1,LyβBr(v);df=10.0^pdf) +
#    (1-frac(f))*Aly*passthroughcombined(LyβCS(v),pΣ2,LyβBr(v);df=10.0^pdf)

TwoCompFit(v,Aback,Sback,Aly,pΣ1,pΣ2,f,pdf,Sly) =
    Aback + (Sback*v/10000) +
    frac(f)*(Aly+ (Sly*v/10000)) *passthroughcombined(LyβCS(v),pΣ1,LyβBr(v);df=10.0^(-7.0)) +
    (1-frac(f))*(Aly+ (Sly*v/10000))*passthroughcombined(LyβCS(v),pΣ2,LyβBr(v);df=10.0^(pdf))


function χ2_twoComp(x,datavel,dataflux)
    pred = map(v->TwoCompFit(v,x...),datavel)
    sum(map( x->max(min(x,0.8),-0.8),(pred-dataflux).^2))+
        (max(x[4],x[5],24.5)-24.5) + (max(abs(x[6]-0.5),0.4)-0.4) +
        (max(abs(x[4]-x[5]),1.0)-1.0) +
        (max(x[4]-x[5],0.0)-0.0) + 
        (1.0 + tanh((x[4]-x[5])*10))
end

function degradedData(dataΔv, dataflx)
    selection = rand(length(dataΔv)) .< 1-1/exp(1)
    (dataΔv[selection], dataflx[selection])
end

#-------

sReLU(x) = log(1 + exp(3*x))

function χ2_twoComp(x,datavel,dataflux)
    pred = map(v->TwoCompFit(v,x...),datavel)
    sum(map( x->x,(pred-dataflux).^2)) +
    1.0*(sReLU(x[4]-26.0)+sReLU(x[5]-26.0)) + 
    1.0*(sReLU(18-x[4])+sReLU(18-x[5])) +
 #  (x[4]-22.8)^2 + (x[5]-21.5)^2 + 
    0.5*sReLU(5.0*(x[5]-x[4])) #+
 #  0.25*sReLU(2.0*(x[4]-x[5])-4.0)
end

function degradedData(dataΔv, dataflx)
    (dataΔv , dataflx .* (1.0 .+ 0.02.*randn(length(dataΔv))))
end

# To estimate the mass we assume 
# a shell that expanded at 650km/s for the 
# following number of years:
t_eff = 17.0 
# Hydrogen mass fraction
X = 0.55
# Geometrical relation between mass and column density
ϵ = 4.5

# One model parameters: Aback, Sback, Aly, pΣ, pdf
# Two model parameters: Aback, Sback, Aly, pΣ1, pΣ2, f, pdf

# Derive the total mass or the mass of each component from the one/two compoennt solutions:
MassFromOneComp(x) =  (10.0^x[4])*4*π*(t_eff*3.2e7*650.0*1e5)^2 / (ϵ*X * Msun * NA)
MassFromTwoComp(x) =  (10.0^x[4]*frac(x[6])+10.0^x[5]*(1-frac(x[6])))*4*π*(t_eff*3.2e7*650.0*1e5)^2 / (ϵ*X * Msun * NA)
Mass1FromTwoComp(x) = (10.0^x[4]*frac(x[6]))*4*π*(t_eff*3.2e7*650.0*1e5)^2 / (ϵ*X * Msun * NA)
Mass2FromTwoComp(x) = (10.0^x[5]*(1-frac(x[6])))*4*π*(t_eff*3.2e7*650.0*1e5)^2 / (ϵ*X * Msun * NA)

function OneCompFitDriver(dataΔv, dataflx;SA=true,initguess=[1.0,-0.04,1.5,23.0,-1.5],fixeddust=false,degraded=false)

  if degraded==true
     dΔv, dflx = degradedData(dataΔv, dataflx)
    else
    dΔv, dflx = dataΔv, dataflx
  end

  if SA==true
      if fixeddust == true
          resSA = optimize(x->χ2_oneComp([x[1:4]...,initguess[5]],dΔv,dflx),
                           initguess[1:4],
                           SimulatedAnnealing(),
                           Optim.Options(iterations = 2000))
          newinit = [Optim.minimizer(resSA)...,initguess[5]]
        else
          resSA = optimize(x->χ2_oneComp(x[1:5],dΔv,dflx),
                           initguess,
                           SimulatedAnnealing(),
                           Optim.Options(iterations = 2000))
          newinit = Optim.minimizer(resSA)
      end
    else
      newinit = initguess
  end

  if fixeddust == true
          res = optimize(x->χ2_oneComp([x[1:4]...,newinit[5]],dΔv,dflx),
                           newinit[1:4],
                           Optim.Options(iterations = 2000))
          newinit = [Optim.minimizer(res)...,initguess[5]]
      else
          res = optimize(x->χ2_oneComp(x[1:5],dΔv,dflx),
                           newinit,
                           Optim.Options(iterations = 2000))
          newinit = Optim.minimizer(res)
   end

   (Optim.minimum(res)/length(dΔv), newinit, MassFromOneComp(newinit))
end

function TwoCompFitDriver(dataΔv, dataflx;SA=true,initguess=[1.0,-0.04,1.5,23.0,23.1,0.5,-1.5,0.0],fixeddust=false,degraded=false)

  if degraded==true
     dΔv, dflx = degradedData(dataΔv, dataflx)
    else
    dΔv, dflx = dataΔv, dataflx
  end

  if SA==true
      if fixeddust == true
          resSA = optimize(x->χ2_twoComp([x[1:6]...,initguess[7],x[7]],dΔv,dflx),
                           [initguess[1:6]...,initguess[8]],
                           SimulatedAnnealing(),
                           Optim.Options(iterations = 10000))
          newinit = [Optim.minimizer(resSA)[1:6]...,initguess[7],Optim.minimizer(resSA)[7]]
        else
          resSA = optimize(x->χ2_twoComp(x[1:8],dΔv,dflx),
                           initguess,
                           SimulatedAnnealing(),
                           Optim.Options(iterations = 10000))
          newinit = Optim.minimizer(resSA)
      end
    else
      newinit = initguess
  end

  if fixeddust == true
          res = optimize(x->χ2_twoComp([x[1:6]...,newinit[7],x[7]],dΔv,dflx),
                           [newinit[1:6]...,newinit[8]],
                           Optim.Options(iterations = 30000))
          newinit = [Optim.minimizer(res)[1:6]...,initguess[7],Optim.minimizer(res)[7]]
      else
          res = optimize(x->χ2_twoComp(x[1:8],dΔv,dflx),
                           newinit,
                           Optim.Options(iterations = 30000))
          newinit = Optim.minimizer(res)
   end

   (Optim.minimum(res)/length(dΔv), newinit, MassFromTwoComp(newinit))
end

function plotstyle()
    display(plot!(
            guidefont=font(24), left_margin=35px, bottom_margin=20px, 
            right_margin=20px,
            size=(1000,800),frame=:box, 
            xtickfontsize=18, ytickfontsize=18,
            legendfontsize=14
        ))
end

function plot_passthrough(dust_fraction)
    contourf(-25:0.125:-18,20:0.125:26,
        (x,y)->passthroughcombined((x),y,0.2;df=dust_fraction),
        levels=[0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1],
        legend=:none,
        xlabel = L"\sigma_{RR} (cm^2)",ylabel=L"\Sigma ~(H/cm^2)",
        seriescolor=cgrad(:heat,[0.0,1.1,0.925]),
        guidefont=font(24), left_margin=35px, bottom_margin=20px, 
        size=(1000,800),frame=:box, xtickfontsize=18, ytickfontsize=18)
    contour!(-25:0.125:-18,20:0.125:26,
        (x,y)->passthroughcombined((x),y,0.2;df=dust_fraction),
        levels=[0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1],
        legend=:none,
        contour_labels = true,
        seriescolor=:black)
    dfrac = dust_fraction > 0.000001 ? dust_fraction : 0
    annotate!(-24.8, 20.25, text("Dust/ISM = $dfrac", :blue, :left, 18))
    plotstyle()
    plot!([-18,-18,-25,-25,-18],[20,26,26,20,20],c=:black,xlims=(-25,-18),ylims=(20,26))
end 

function plot_model_results(pdf;)
    Δv = collect(-5e4:100:6.3e4)
    plot()
    for pΣ = 20:0.5:25.5
        plot!(Δv,(x->100*6.6^-2*(1+x/6.6/300000)^2*OneCompFit(x,0.0,0.0,1,pΣ,pdf)).(Δv), 
              label="$pΣ", legend=:topright)
    end
    p = plot!(frame=:box,
         xlabel=L"v (km/s)",ylabel=L"100 ~\times ~ dn/dλ~/~\left(dn/dλ\right)_\alpha")
    rangey = Plots.ylims(p)
    dfrac = pdf > -6.5 ? round(100000*10^pdf)/100000 : 0
    loc_dust_y = rangey[1]+0.95*(rangey[2]-rangey[1])
    plotstyle()
    annotate!(-4.5e4, loc_dust_y, text("Dust/ISM = $dfrac", :blue, :left, 18))
end

function plot_Hα_Hβ_compare(pdf, pΣ, aα, aβ ;)
    Δv = collect(-1e4:100:1.0e4)
    plot()
    plot!(Δv,1 .+ (x->100*6.6^-2*(1+x/6.6/300000)^2*OneCompFit(x,0.0,0.0,aα,pΣ,pdf)).(Δv), 
              label="Hα $pΣ", legend=:topright)
    plot!(Δv,1 .+ (x->100*6.6^-2*(1+x/6.6/300000)^2*OneCompFitHβ(x,0.0,0.0,aβ,pΣ,pdf)).(Δv), 
              label="Hβ $pΣ", legend=:topright)
    p = plot!(frame=:box,
         xlabel=L"v (km/s)",ylabel=L"(F_{cont} + F_{Raman})/F_{cont}", ylims = (0,1.34))
    rangey = Plots.ylims(p)
    dfrac = pdf > -6.5 ? round(100000*10^pdf)/100000 : 0
    loc_dust_y = rangey[1]+0.95*(rangey[2]-rangey[1])
    plotstyle()
    annotate!(-4.5e4, loc_dust_y, text("Dust/ISM = $dfrac", :blue, :left, 18))
end

function planck_radiation_wavelength(T::Float64, λ::Float64)
    h = 6.62607015e-34  # Planck's constant in J*s
    c = 3.0e8           # Speed of light in m/s
    kB = 1.380649e-23   # Boltzmann constant in J/K

    numerator = 2.0 * h * c^2 / λ^5
    denominator = exp(h * c / (λ * kB * T)) - 1.0

    return numerator / denominator
end


function plot_data_cuts(;dist=300)
    diffs=[0,[ dataΔv_plot[i+1]-dataΔv_plot[i]  for i=1:length(dataΔv_plot)-1]...] 
    vecs_id = cumsum(map( x-> x > dist ? 1 : 0 , diffs)) .+1
    for i=1:vecs_id[end]
        plot!( dataΔv_plot[vecs_id .== i], dataflx_plot[vecs_id .== i], c=RGB(0.1,0.1,0.7), label=:none) 
    end
    plot!(legend=:none)
end 


function plot_fit(;with_error=false, two_comp=false, plot_data=true, label=:none, color=:auto, linestyle=:auto, lw=1.0)
    if plot_data==true
        plot()
        plot_data_cuts()
    end
    range=-25000:250:24000
    if with_error == true
        for i ∈ 1:length(restable)
            if two_comp == true 
                plot!(range,map(v->TwoCompFit(v,restable[i][2]...),range),linewidth=1, label=:none, color=color, linestyle=linestyle)
            else 
                plot!(range,map(v->OneCompFit(v,restable[i][2]...),range),linewidth=1, label=:none, color=color, linestyle=linestyle)       
            end
        end
    end
    if two_comp == true 
        p = plot!(range,map(v->TwoCompFit(v,resbest[2]...),range),linewidth=lw, label=label, color=color, linestyle=linestyle)
    else
        p = plot!(range,map(v->OneCompFit(v,resbest[2]...),range),linewidth=lw, label=label, color=color, linestyle=linestyle)
    end
    plot!(frame=:box, xlims=(range[1],range[end]), legend=:topright)
    rangey = Plots.ylims(p)
    loc_data = rangey[1]+0.92*(rangey[2]-rangey[1])
    annotate!(-22e3 , loc_data, text(data_name, c=RGB(0.1,0.1,0.7), :left, 18))
    plotstyle()
end 

