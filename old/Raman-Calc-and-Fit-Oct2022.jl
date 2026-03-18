# Import the different packages we need to use
# Note, for the first time you need to download/install the packages, which you do through
# using Pkg
# Pkg.add("DataFrames") etc.

using DataFrames, CSV     # Import from CSV files, defines DataFrame structure (like excel sheets)
using Plots               # Standard plotting packages
using BenchmarkTools      # More options to time things
using Optim, ForwardDiff  # Minimizer package and differentiator 
using Interpolations      # Interpolation of functions 
using ProgressMeter       # Progress bar for long calculations
using Roots               # Find roots of functions
using Base.Threads        # To use multithreading for long calculations
using LaTeXStrings
using ColorSchemes
using Plots.PlotMeasures


# - Some constants

NA = 6.022e23      # Avogadro's number
Msun = 1.98e33     # Solar mass in grams

# -------------- First lets define the functions we need

# This calculates the Hα flux passed through a slab. 
# The analytic calculation was carried out with Mathematica 
# pσ is the log10 of the total cross-section in c.g.s for the Lyman photon *Raman* scattering
# pΣ is the log10 of the total column density of the shell
# BR is the branching ratio into the Balmer line
# df is the fraction of dust relative to standard ISM value
function passthrough(pσ, pΣ, BR; df=0.01)
    aH = 0.75  # scattering fraction (albedo) for the Balmer lines
    aL = 0.3   # scattering fraction for the Lyman lines
    Σ = 10.0^pΣ
    σ_raman_tot = 10.0^pσ                 # The Raman scattering cross-section for Lyman photons
    κH_dust =  df*(1-aH)*3.7e-22          # Balmer frequency absorption cross-section 
    σH_dust =  df*  aH  *3.7e-22+6.6e-25  # Balmer frequnecy scattering cross-section
    κL_dust =  df*(1-aL)*2.2e-21          # Lyman frequency absorption cross-section
    σL_dust =  df*  aL  *2.2e-21+6.6e-25  # Lyman frequency scattering cross-section
    σ_rayleigh = (1-BR)*σ_raman_tot       # cross-section for Rayleigh scattering   
    σ_raman    =   (BR)*σ_raman_tot       # For the Lyman band, it all appears as absorption 
    κL = σ_raman + κL_dust                # Raman scattering is like absorption for Lyman
    σL = σL_dust + σ_rayleigh

    χLeff = sqrt(κL*(κL+σL))
    χHeff = sqrt(κH_dust*(κH_dust+σH_dust))
    χLtot = κL + σL
    χHtot = κH_dust + σH_dust

    enum=2*σ_raman*χHtot*χLeff*(κH_dust*(χHtot + χLtot)*cosh(Σ*χHeff) -
       κH_dust*(χHtot + χLtot)*cosh(Σ*χLeff) + χHeff*(κH_dust + χLtot)*sinh(Σ*χHeff) -
       (κH_dust*(κL + χHtot)*sqrt(χLtot)*sinh(Σ*χLeff))/sqrt(κL))

    denom=sqrt(κH_dust)*sqrt(χHtot)*(κH_dust*χHtot - κL*χLtot)*(2*χHeff*cosh(Σ*χHeff) +
       (2*κH_dust + σH_dust)*sinh(Σ*χHeff))*(2*χLeff*cosh(Σ*χLeff) + (2*κL + σL)*sinh(Σ*χLeff))

    enum/denom
end


# Remove NaNs from results:
fixNaN(num) = isnan(num) ? 0 : num

# The expression for the pass through flux has very large exponents. The simplest solution is to 
# do the calculation with higher precision in those regions of phase space that require it.
# We use Julia's "native" arbitrary float precision option (There are packages, such as QuadMath which doesn't 
# run on ARM (new macs) or DoubleFloats, which doesn't have enough exponent bits). 
setprecision(128)
function passthroughcombined(pσ, pΣ, BR; df=0.01)
    if pσ + pΣ < 2.8 return passthrough(pσ, pΣ, BR; df=df) |> fixNaN end
    if pσ + pΣ < 4.4 return Float64(passthrough(BigFloat(pσ), pΣ, BR; df=df)) |> fixNaN end
    Float64(passthrough(BigFloat(4.4-pΣ), pΣ, BR; df=df)) |> fixNaN
end

# simplest fit - one component

OneCompFit(v,Aback,Sback,Aly,pΣ,pdf) =
    Aback + (Sback*v/10000) +
    Aly*passthroughcombined(LyβCS(v),pΣ,LyβBr(v);df=10.0^pdf)

function χ2_oneComp(x,datavel,dataflux)
    pred = map(v->OneCompFit(v,x...),datavel)
    sum(map( x->max(min(x,0.8),-0.8),(pred-dataflux).^2))
end

TwoCompFit(v,Aback,Sback,Aly,pΣ1,pΣ2,f,pdf) =
    Aback + (Sback*v/10000) +
    f*Aly*passthroughcombined(LyβCS(v),pΣ1,LyβBr(v);df=10.0^pdf) +
    (1-f)*Aly*passthroughcombined(LyβCS(v),pΣ2,LyβBr(v);df=10.0^pdf)


function χ2_twoComp(x,datavel,dataflux)
    pred = map(v->TwoCompFit(v,x...),datavel)
    sum(map( x->max(min(x,0.8),-0.8),(pred-dataflux).^2))+
        (max(x[4],x[5],24.5)-24.5) + (max(abs(x[6]-0.5),0.4)-0.4) +
        (max(abs(x[4]-x[5]),1.0)-1.0) +
        (max(x[4]-x[5],0.0)-0.0)
end


function degradedData(dataΔv, dataflx)
    selection = rand(length(dataΔv)) .< 1-1/exp(1)
    (dataΔv[selection], dataflx[selection])
end


# To estimate the mass we assume that the effecive are is of a shell that expanded at 500km/s for the 
# following number of years:
t_eff = 10.0
# Hydrogen mass fraction
X = 0.5

# Derive the total mass or the mass of each component from the one/two compoennt solutions:
MassFromOneComp(x) =  (x[3]*10.0^x[4])*4*π*(t_eff*3.2e7*500.0*1e5)^2 / (X * Msun * NA)
MassFromTwoComp(x) =  x[3]*(10.0^x[4]*x[6]+10.0^x[5]*(1-x[6]))*4*π*(t_eff*3.2e7*500.0*1e5)^2 / (X * Msun * NA)
Mass1FromTwoComp(x) = x[3]*(10.0^x[4]*x[6])*4*π*(t_eff*3.2e7*500.0*1e5)^2 / (X * Msun * NA)
Mass2FromTwoComp(x) = x[3]*(10.0^x[5]*(1-x[6]))*4*π*(t_eff*3.2e7*500.0*1e5)^2 / (X * Msun * NA)

function OneCompFitDriver(dataΔv, dataflx;SA=true,initguess=[1.0,-0.04,1.5,22.8,-1.5],fixeddust=false,degraded=false)

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


function TwoCompFitDriver(dataΔv, dataflx;SA=true,initguess=[1.0,-0.04,1.5,22.8,23.1,0.5,-1.5],fixeddust=false,degraded=false)

  if degraded==true
     dΔv, dflx = degradedData(dataΔv, dataflx)
    else
    dΔv, dflx = dataΔv, dataflx
  end

  if SA==true
      if fixeddust == true
          resSA = optimize(x->χ2_twoComp([x[1:6]...,initguess[7]],dΔv,dflx),
                           initguess[1:6],
                           SimulatedAnnealing(),
                           Optim.Options(iterations = 10000))
          newinit = [Optim.minimizer(resSA)...,initguess[7]]
        else
          resSA = optimize(x->χ2_twoComp(x[1:7],dΔv,dflx),
                           initguess,
                           SimulatedAnnealing(),
                           Optim.Options(iterations = 10000))
          newinit = Optim.minimizer(resSA)
      end
    else
      newinit = initguess
  end

  if fixeddust == true
          res = optimize(x->χ2_twoComp([x[1:6]...,newinit[7]],dΔv,dflx),
                           newinit[1:6],
                           Optim.Options(iterations = 30000))
          newinit = [Optim.minimizer(res)...,initguess[7]]
      else
          res = optimize(x->χ2_twoComp(x[1:7],dΔv,dflx),
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
        levels=[0.0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.45,0.5,1.0],
        legend=:none,
        xlabel = L"\sigma (cm^2/g)",ylabel=L"\Sigma ~(g/cm^2)",
        seriescolor=cgrad(:heat,[0.0,1.1,0.95]),
        guidefont=font(24), left_margin=35px, bottom_margin=20px, 
        size=(1000,800),frame=:box, xtickfontsize=18, ytickfontsize=18)
    contour!(-25:0.125:-18,20:0.125:26,
        (x,y)->passthroughcombined((x),y,0.2;df=dust_fraction),
        levels=[0.0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.45,0.5,1.0],
        legend=:none,
        contour_labels = true,
        seriescolor=:black)
    annotate!(-24.5, 25.5, text("Dust/ISM = $dust_fraction", :blue, :left, 18))
    plotstyle()
end 

function plot_model_results(pdf)
    Δv = collect(-5e4:100:5e4)
    plot()
    for pΣ = 20:0.5:25.5
        plot!(Δv,(x->100*6.6^-2*(1+x/6.6/300000)^2*OneCompFit(x,0.0,0.0,1,pΣ,pdf)).(Δv), label="$pΣ")
    end
    plot!(frame=:box,
         xlabel=L"v (km/s)",ylabel=L"100 ~\times ~ dn/dλ~/~\left(dn/dλ\right)_\alpha")
    plotstyle()
end


#=
Note that Julia doesn't require the functions to be defined separately
like it is here, however, it cleaner this way. Now we can start doing 
the analysis.
=#

#-------------------------------------------------------------------------- 
## First, plot the pass through function (that we derived from Mathematica)
#-------------------------------------------------------------------------- 


plot_passthrough(1)
savefig("figures/Pass-through-ISM-dust.pdf")
plot_passthrough(0.1)
savefig("figures/Pass-through-0.1-ISM-dust.pdf")
plot_passthrough(0.01)
savefig("figures/Pass-through-0.01-ISM-dust.pdf")
plot_passthrough(0.001)
savefig("figures/Pass-through-0.001-ISM-dust.pdf")
plot_passthrough(0.0000001)
savefig("figures/Pass-through-No-ISM-dust.pdf")

#-------------------------------------------------------------------------- 
## Next, plot solutions
#-------------------------------------------------------------------------- 

plot_model_results(-7.0)
savefig("figures/Model_Results_no-dust.pdf")
plot_model_results(-3.0)
savefig("figures/Model_Results_dust-pdf-m3.0.pdf")
plot_model_results(-2.0)
savefig("figures/Model_Results_dust-pdf-m2.0.pdf")
plot_model_results(-1.0)
savefig("figures/Model_Results_dust-pdf-m1.0.pdf")
plot_model_results(-0.0)
savefig("figures/Model_Results_dust-pdf-m0.0.pdf")

# Start with one component model

#OneCompFit(v,Aback,Sback,Aly,pΣ,pdf) 

savefig("figures/Pass-through-No-ISM-dust.pdf")


OneCompFit(+1000,1.0,0.0,0.1,23,-6.5)

for i ∈ 1:length(restable)
    plot!(dataΔv,map(v->OneCompFit(v,restable[i][2]...),dataΔv),linewidth=1,label=:none)
end


## read the data, uncomment one line

data_to_read = "Gemini-1"
#data_to_read = "Magelan-1"

if data_to_read == "Gemini-1"
    data = CSV.read("data/Gemini-etaCarinae-2014-11-fit.csv",
           DataFrame, header = ["Δv","flux"]);offset=-1
    data_name = "Gemini 11 2014"
end
if data_to_read == "Magelan-1"
    data = CSV.read("data/Magelan-etaCarinae-2015-01-fit.csv",
           DataFrame, header = ["Δv","flux"]);offset=0
    data_name = "Magellan 01 2015"
end

## We cut out the data near zero velocity.
dataΔv = 0 .+ 1000 .* data.Δv[abs.(data.Δv) .>1]
dataflx = 1.0 .+ data.flux[abs.(data.Δv) .>1] .+ offset

scatter(dataΔv,dataflx)
errornorm = mean((dataflx[2:end]-dataflx[1:end-1]).^2)/2

# We shall work in velocity space relative to the Balmer alpha frequency (Hα)
# The files have the BR and Cross-section in units of the Lyman β doppler.
# The ratio between the two is 32/5

LyβBRRaw = CSV.read("data/Lyβ-BranchingRatio-better.csv",DataFrame,  header = ["Δv0","br"])
LyβBr = LinearInterpolation(LyβBRRaw.Δv0 .* 32/5,LyβBRRaw.br)

LyβCSRaw = CSV.read("data/Lyβ-TotalCrossSection.csv",DataFrame,  header = ["Δv0","cs"])
LyβCS = LinearInterpolation(LyβCSRaw.Δv0 .* 32/5,LyβCSRaw.cs)

#------------------------------------
# Best fit with dust
#------------------------------------

resbest = OneCompFitDriver(dataΔv,dataflx; SA=false, degraded = false,
                           fixeddust = false, initguess=[1.0,-0.04,1.5,22.8,-1.5])
restable =  @showprogress [
          OneCompFitDriver(dataΔv,dataflx; SA=false, degraded = true,
                           fixeddust = false, initguess=resbest[2])
                          for i ∈ 1:200]

χ20   = resbest[1]
Mass0 = resbest[3]
dust0 = resbest[2][5]
χ2   = [restable[i][1] for i ∈ 1:length(restable)]
ind_good =  [χ2 .< 1.5*χ20]
Mass = [restable[i][3] for i ∈ 1:length(restable)]
dust = [restable[i][2][5] for i ∈ 1:length(restable)]

scatter(dataΔv,dataflx)
for i ∈ 1:length(restable)
    plot!(dataΔv,map(v->OneCompFit(v,restable[i][2]...),dataΔv),linewidth=1,label=:none)
end
plot!(dataΔv,map(v->OneCompFit(v,resbest[2]...),dataΔv),linewidth=2)
plot!(frame=:box)

scatter((Mass), χ2)



scatter(log10.(Mass), χ2)
scatter!([log10.(Mass0)], [χ20])

scatter(10 .^ log10.(Mass), dust)
scatter!(10.0.^[log10.(Mass0)], [dust])

#-----------------------
# Best with with no dust
#-----------------------

resbest = OneCompFitDriver(dataΔv,dataflx; SA=false, degraded = false,
                           fixeddust = true, initguess=[1.0,-0.04,0.2,22.8,-6.5])
restable=[]
#restable =  @time @threads
@time @threads for i ∈ 1:200
    push!(restable,OneCompFitDriver(dataΔv,dataflx; SA=false, degraded = true,
                    fixeddust = true, initguess=resbest[2]))
          end

χ20   = resbest[1]
Mass0 = resbest[3]
dust0 = resbest[2][5]
χ2   = [restable[i][1] for i ∈ 1:length(restable)]
Mass = [restable[i][3] for i ∈ 1:length(restable)]
dust = [restable[i][2][5] for i ∈ 1:length(restable)]

scatter(dataΔv,dataflx)
for i ∈ 1:length(restable)
    plot!(dataΔv,map(v->OneCompFit(v,restable[i][2]...),dataΔv),linewidth=1,label=:none)
end
plot!(dataΔv,map(v->OneCompFit(v,resbest[2]...),dataΔv),linewidth=2)
plot!(frame=:box)

scatter(log10.(Mass), χ2)
scatter!([log10.(Mass0)], [χ20])

scatter(10 .^ log10.(Mass), χ2)


# Best fit two components with dust


#startguesses = [restable[i][2] for i ∈ 1:length(restable)]

resbest = TwoCompFitDriver(dataΔv,dataflx; SA=false, degraded = false,
                           fixeddust = true, initguess=[2.0,-0.04,1.5,21.5,22.9,0.5,-0.5])

resbest[2]

restable=[];
@showprogress for i ∈ 1:20
                 push!(restable,TwoCompFitDriver(dataΔv,dataflx; SA=false, degraded = true,
                                                 fixeddust = true,initguess=resbest[2] )) #startguesses[rand(100)]))
               end

resbest[2][4]

χ20   = resbest[1]
Mass0 = resbest[3]
dust0 = resbest[2][7]
χ2   = [restable[i][1] for i ∈ 1:length(restable)]
Mass = [restable[i][3] for i ∈ 1:length(restable)]
Mass1 = [Mass1FromTwoComp(restable[i][2]) for i ∈ 1:length(restable)]
Mass2 = [Mass2FromTwoComp(restable[i][2]) for i ∈ 1:length(restable)]
Σ1 = [restable[i][2][4] for i ∈ 1:length(restable)]
Σ2 = [restable[i][2][5] for i ∈ 1:length(restable)]
dust = [restable[i][2][7] for i ∈ 1:length(restable)]
fs = [restable[i][2][6] for i ∈ 1:length(restable)]

scatter(dataΔv,dataflx)
for i ∈ 1:length(restable)
    plot!(dataΔv,map(v->TwoCompFit(v,restable[i][2]...),dataΔv),linewidth=1,label=:none)
end
plot!(dataΔv,map(v->TwoCompFit(v,resbest[2]...),dataΔv),linewidth=2)
plot!(frame=:box)

scatter(log10.(Mass), χ2)
scatter!([log10.(Mass0)], [χ20])

scatter(pairmax(Σ1,Σ2),pairmin(Σ1,Σ2))

scatter( log10.(Mass), 10.0 .^dust)


pairmax(v1,v2) = map((x,y)->max(x,y),v1,v2)
pairmin(v1,v2) = map((x,y)->min(x,y),v1,v2)
pairf(v1,v2,f) = map((x,y,f)-> x>y ? f : (1.0-f) ,v1,v2,f)

scatter( log10.(pairmax(Mass1,Mass2)),log10.(pairmin(Mass1,Mass2)))

scatter( log10.(Mass1),log10.(Mass2))



scatter( log10.(pairmax(Mass1,Mass2))-log10.(pairmin(Mass1,Mass2)), pairf(Mass1,Mass2,fs))





















# Two component fit

funfit2(v,Aback,Sback,Aly,pdf,pΣ,pΣ2,f1) =
      Aback + (Sback*v/10000) +
      f1*Aly*passthroughcombined(LyβCS(v),pΣ ,LyβBr(v);df=10.0^pdf) +
      (1.0-f1)*Aly*passthroughcombined(LyβCS(v),pΣ2,LyβBr(v);df=10.0^pdf)

# example solution:
plot(dataΔv,map(v-> funfit2(v,1,-0.03,1.5,-1.4,22.0,23.0,0.5),dataΔv))

function fitfun2(x,Σs,f1,datavel,dataflux)
   pred = map(v->funfit2(v,x...,Σs...,f1),datavel)
   sum( (pred-dataflux).^2 )/length(pred)/errornorm #- 10.0*min(x[3],x[4],0) #+ 10.0*max(Σs[1]-24,Σs[2]-24,0)
end

#-------------------

function bestfitWithDegrade(dataΔv,dataflx;degrade=false,SA=true,initguess=[0.995,-0.0347,1.7,0.7,-1.4,22.2,21.0])
    if degrade==false
        dΔv=dataΔv
        dflx=dataflx
    else
        selection = rand(length(dataΔv)) .< 1-1/exp(1)
        dΔv=dataΔv[selection]
        dflx=dataflx[selection]
    end
    if SA==true
            resSA = optimize(x->fitfun2(x[1:5],x[6:7],dΔv,dflx),
            initguess,
            SimulatedAnnealing(),
            Optim.Options(iterations = 2000))
            newinit = Optim.minimizer(resSA)
        else
            newinit = initguess
    end
    res  = optimize(x->fitfun2(x[1:5],x[6:7],dΔv,dflx),
            newinit,
          #NelderMead(initial_simplex = MySimplexer()),
#            ConjugateGradient(),
            Optim.Options(iterations = 3000))
    sol=Optim.minimizer(res)
    totM = (sol[3]*10.0^sol[6]+sol[4]*10.0^sol[7])/(sol[3]+sol[4])*4*π*(15.0*3.2e7*500.0*1e5)^2 / (0.5 * 2.0e33 * 6.0e23)
    (sol[6:7],Optim.minimum(res)/length(dataΔv)/errornorm,sol,totM)
end

Mass1FromRes(x,f) = (x[3]*10.0^x[5]*f)*4*π*(10.0*3.2e7*500.0*1e5)^2 / (0.5 * 2.0e33 * 6.0e23)
@showprogress [(
            tmp=optimize(x->fitfun2(x[1:4],x[5:6],.5,dataΔv,dataflx),
                [0.97,-0.03,1.7,-3.0,22.5,23.5],
  #             ConjugateGradient(),
                Optim.Options(iterations = 4000));
               (Optim.minimizer(tmp), Optim.minimum(tmp),Mass1FromRes(Optim.minimizer(tmp),0.5) )
               ) for i in 1:2]





restable= @showprogress [(bestfitWithDegrade(dataΔv,dataflx;SA=false,
          initguess=[0.962,-0.0469,2.18,1.41,-1.4,22.702,23.77])) for i ∈ 1:10]
χ2 = [restable[i][2] for i ∈ 1:length(restable)]
Mass = [restable[i][4] for i ∈ 1:length(restable)]
dust = [restable[i][3][5] for i ∈ 1:length(restable)]


scatter(dataΔv,dataflx)
sol = restable[1][3]
plot!(dataΔv,map(v->funfit2(v,sol...),dataΔv),linewidth=5)


scatter(Mass,χ2)
scatter(dust,χ2)


const c = 2.99792458e8
const hp = 6.62607015e-34
const kB = 1.380649e-23
const σSB= 5.67e-8
B(λ,T) = 1.0e-6*2*hp*c^2/((λ*1.0e-6)^5*(exp(hp*c/((λ*1.0e-6)*kB*T))-1))


TforRatio(ratio) = find_zero(T -> B(0.1015,T)/B(0.1015*32/5,T)/(32/5)^3 - ratio,(10000,200000))
TforRatio(2)

