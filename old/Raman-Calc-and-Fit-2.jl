using DataFrames, CSV, Plots
using Quadmath #Provides Float64
using BenchmarkTools
using Optim, ForwardDiff
using Interpolations
using ProgressMeter

data = CSV.read("data/Gemini-etaCarinae-2014-11-fit.csv",DataFrame, header = ["Δv","flux"]);offset=-1
data = CSV.read("data/Magelan-etaCarinae-2015-01-fit.csv",DataFrame, header = ["Δv","flux"]);offset=0

# We cut out the data near zero velocity.
dataΔv = 1000 .* data.Δv[abs.(data.Δv) .>1]
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

function passthrough(pσ, pΣ, BR; df=0.01)
    a = 0.5
    Σ = 10.0^pΣ
    σtot = 10.0^pσ
    κH = df*(1-a)*3.7e-22
    σH = df*a*3.7e-22+6.6e-25
    κLd = df*(1-a)*2.2e-21
    σLd = df*a*0.5*2.2e-21+6.6e-25
    σRe = (1-BR)*σtot
    κR = BR*σtot
    κT = (κR + κLd)
    σL = σLd + σRe

    χLeff = sqrt(κT*(κT+σL))
    χHeff = sqrt(κH*(κH+σH))
    χLtot = κT + σL
    χHtot = κH + σH

    enum=2*κR*χHtot*χLeff*(κH*(χHtot + χLtot)*cosh(Σ*χHeff) -
       κH*(χHtot + χLtot)*cosh(Σ*χLeff) + χHeff*(κH + χLtot)*sinh(Σ*χHeff) -
       (κH*(κT + χHtot)*sqrt(χLtot)*sinh(Σ*χLeff))/sqrt(κT))

    denom=sqrt(κH)*sqrt(χHtot)*(κH*χHtot - κT*χLtot)*(2*χHeff*cosh(Σ*χHeff) +
       (2*κH + σH)*sinh(Σ*χHeff))*(2*χLeff*cosh(Σ*χLeff) + (2*κT + σL)*sinh(Σ*χLeff))

    enum/denom
end

fixNaN(num) = isnan(num) ? 0 : num

function passthroughcombined(pσ, pΣ, BR; df=0.01)
    if pσ + pΣ < 3.2 return passthrough(pσ, pΣ, BR; df=df) |> fixNaN end
    if pσ + pΣ < 4.4 return Float64(passthrough(Float128(pσ), pΣ, BR; df=df)) |> fixNaN end
    Float64(passthrough(Float128(4.4-pΣ), pΣ, BR; df=df)) |> fixNaN
end

#contourf(-23:0.125:-18,20:0.125:26,(x,y)->passthroughcombined((x),y,0.2;df=0.1))

@time contourf(-23:0.125:-18,20:0.125:26,(x,y)->passthroughcombined((x),y,0.2;df=0.01))
@btime [passthroughcombined((x),y,0.2;df=0.001) for x=-23:0.125:-18,y=20:0.125:23]

# simplest fit

funfit(v,Aback,Sback,Aly,pΣ,pdf) =
    Aback + (Sback*v/10000) +
    Aly*passthroughcombined(LyβCS(v),pΣ,LyβBr(v);df=10.0^pdf)

function fitfun(x,datavel,dataflux)
    pred = map(v->funfit(v,x...),datavel)
    sum(map( x->max(min(x,0.8),-0.8),(pred-dataflux).^2))
end

initguess =   [1,-0.05,2.0,23.5,-2.0]
resSA = optimize(x->fitfun(x,dataΔv,dataflx),initguess,SimulatedAnnealing())
res   = optimize(x->fitfun(x,dataΔv,dataflx),Optim.minimizer(resSA))

#=
f(x) = fitfun(x,dataΔv,dataflx)
f(initguess)
g(x) = ForwardDiff.gradient(f,x)
g(initguess)
=#

res = optimize(x->fitfun(x,dataΔv,dataflx),initguess,LBFGS())


sol=Optim.minimizer(res)
scatter(dataΔv,dataflx)
plot!(dataΔv,map(v->funfit(v,sol...),dataΔv),linewidth=5)

# Two component fit

funfit2(v,Aback,Sback,Aly,Aly2,pdf,pΣ,pΣ2) =
      Aback + (Sback*v/10000) +
      Aly*passthroughcombined(LyβCS(v),pΣ ,LyβBr(v);df=10.0^pdf) +
      Aly2*passthroughcombined(LyβCS(v),pΣ2,LyβBr(v);df=10.0^pdf)


plot(dataΔv,map(v-> funfit2(v,initguess...,26.0,23.0),dataΔv))

function fitfun2(x,Σs,datavel,dataflux)
   pred = map(v->funfit2(v,x...,Σs...),datavel)
   sum( (pred-dataflux).^2 ) - 10.0*min(x[3],x[4],0) #+ 10.0*max(Σs[1]-24,Σs[2]-24,0)
end

plot(dataΔv)

initguess
fitfun2(initguess,[24.0,22.0],dataΔv,dataflx)

function fitforΣs(Σs)
    initguess =   [0.905,-0.0347,0.5,0.5,-0.4]
    #resSA = optimize(x->fitfun2(x,Σs,dataΔv,dataflx),
    #            initguess,
    #            SimulatedAnnealing(),
#                Optim.Options(iterations = 5000))
    res = optimize(x->fitfun2(x,Σs,dataΔv,dataflx),
#                Optim.minimizer(resSA),
                initguess,
                Optim.Options(iterations = 40000))
    sol=Optim.minimizer(res)
    totM = (sol[3]*10.0^Σs[1]+sol[4]*10.0^Σs[2])/(sol[3]+sol[4])*4*π*(15.0*3.2e7*500.0*1e5)^2 / (0.5 * 2.0e33 * 6.0e23)
    (Σs,Optim.minimum(res)/length(dataΔv)/errornorm,sol,totM)
end

rngΖhigh = 22.0:0.25:24.0;
rngpratio = 0.5:0.25:2.0;
restable = [fitforΣs([Σhigh,Σhigh-pratio]) for Σhigh = rngΖhigh, pratio = rngpratio]

χ2 = [restable[i][2] for i ∈ 1:length(restable)]
contour!(rngΖhigh,rngpratio,(χ2.- 1.8).* 10)

Mass = [restable[i][4] for i ∈ 1:length(restable)]
contourf(rngΖhigh,rngpratio,log10.(Mass))

dust = [restable[i][3][5] for i ∈ 1:length(restable)]
contourf(rngΖhigh,rngpratio, 10 .^dust .+ 0.001)

#fitforΣs([24.0,22.0])

scatter((Mass),χ2)


#----------------

scatter(dataΔv,dataflx)
plot!(dataΔv,map(v->funfit2(v,sol...,Σs...),dataΔv),linewidth=5)
plot!(dataΔv,map(v->funfit2(v,(sol.*[1,1,0,0,1])...,Σs...),dataΔv),linewidth=5)
plot!(ylims=(0,1.5))

plot!(dataΔv,map(v->funfit2(v,initguess.*[1,1,1,1,2,1,1]...),dataΔv),linewidth=5)


#-----------

##### Find best fit

function simpxvecs(i,n)
  if i ≤ n
      tmp = zeros(n)
      tmp[i] = 0.3
      return(tmp)
  else
      return(0.3 .* ones(n))
  end
end

struct MySimplexer <: Optim.Simplexer end
Optim.simplexer(S::MySimplexer, initial_x) = [simpxvecs(i,length(initial_x)) for i = 1:length(initial_x)+1]


function bestfitWithDegrade(dataΔv,dataflx;degrade=false,SA=true,initguess=[0.995,-0.0347,0.7,0.7,-1.4,22.2,21.0])
    if degrade==true
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
#            NelderMead(initial_simplex = MySimplexer()),
            ConjugateGradient(),
            Optim.Options(iterations = 1000))
    sol=Optim.minimizer(res)
    totM = (sol[3]*10.0^sol[6]+sol[4]*10.0^sol[7])/(sol[3]+sol[4])*4*π*(15.0*3.2e7*500.0*1e5)^2 / (0.5 * 2.0e33 * 6.0e23)
    (sol[6:7],Optim.minimum(res)/length(dataΔv)/errornorm,sol,totM)
end


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




@showprogress [(sleep(1.0);sin(i)) for i in 1:10]


.


0.75e22+0.25e23

10^22.7

6e22*(4*π*(3.2e8*500e5)^2)/(6e23*2e33)

a
#--------------------------------------------- Junk


# with unsimplified solution
function passthrough(pσ, pΣ, BR; df=0.01)
    df = 0.01
    a = 0.5
    Σ = 10.0^pΣ
    σtot = 10.0^pσ
    κH = df*(1-a)*3.7e-22
    σH = df*a*3.7e-22
    κLd = df*(1-a)*2.2e-21
    σLd = df*a*0.5*2.2e-21
    σRe = (1-BR)*σtot
    κR = BR*σtot
    d=Σ
    κT = (κR + κLd)
    σL = σLd + σRe
    (2*κR*(κH + σH)*sqrt(κT*(κT + σL))*(κH*(κH + κT + σH + σL)*cosh(d*sqrt(κH)*sqrt(κH + σH)) -
           (κH*cosh(d*sqrt(κH*(κH + σH)))^2*(sqrt(κT)*(κH + κT + σH + σL)*cosh(d*sqrt(κT)*sqrt(κT + σL)) +
           (κH + κT + σH)*sqrt(κT + σL)*sinh(d*sqrt(κT)*sqrt(κT + σL))))/sqrt(κT) +
           sinh(d*sqrt(κH*(κH + σH)))*(sqrt(κH*(κH + σH))*(κH + κT + σL) +
           (κH*sinh(d*sqrt(κH*(κH + σH)))*(sqrt(κT)*(κH + κT + σH + σL)*cosh(d*sqrt(κT)*sqrt(κT + σL)) +
           (κH + κT + σH)*sqrt(κT + σL)*sinh(d*sqrt(κT)*sqrt(κT + σL))))/sqrt(κT))))/
           (sqrt(κH*(κH + σH))*(κH*(κH + σH) - κT*(κT + σL))*(2*sqrt(κH*(κH + σH))*cosh(d*sqrt(κH*(κH + σH))) + (2*κH + σH)*sinh(d*sqrt(κH*(κH + σH))))*
           (2*sqrt(κT*(κT + σL))*cosh(d*sqrt(κT*(κT + σL))) + (2*κT + σL)*sinh(d*sqrt(κT*(κT + σL)))))
end

500/6500*300000

(1-1/3^2)/(1/2^2-1/3^2)
