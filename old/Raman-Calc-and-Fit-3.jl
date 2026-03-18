using DataFrames, CSV, Plots
using Quadmath #Provides Float64
using BenchmarkTools
using Optim, ForwardDiff
using Interpolations
using ProgressMeter
using Roots
using Base.Threads


# -------------- First lets define the functions we need

# This calculates the Hα flux passed through a slab
function passthrough(pσ, pΣ, BR; df=0.01)
    a = 0.5
    Σ = 10.0^pΣ
    σtot = 10.0^pσ
    κH = df*(1-a)*3.7e-22+6.6e-25/10
    σH = df*a*3.7e-22+6.6e-25
    κLd = df*(1-a)*2.2e-21+6.6e-25/10
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

# This is the actual function which ensures that for large absorption we use higher precision and also fix NaN's
function passthroughcombined(pσ, pΣ, BR; df=0.01)
    if pσ + pΣ < 3.2 return passthrough(pσ, pΣ, BR; df=df) |> fixNaN end
    if pσ + pΣ < 4.4 return Float64(passthrough(Float128(pσ), pΣ, BR; df=df)) |> fixNaN end
    Float64(passthrough(Float128(4.4-pΣ), pΣ, BR; df=df)) |> fixNaN
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


MassFromOneComp(x) = (x[3]*10.0^x[4])*4*π*(10.0*3.2e7*500.0*1e5)^2 / (0.5 * 2.0e33 * 6.0e23)

MassFromTwoComp(x) = x[3]*(10.0^x[4]*x[6]+10.0^x[5]*(1-x[6]))*4*π*(10.0*3.2e7*500.0*1e5)^2 / (0.5 * 2.0e33 * 6.0e23)
Mass1FromTwoComp(x) = x[3]*(10.0^x[4]*x[6])*4*π*(10.0*3.2e7*500.0*1e5)^2 / (0.5 * 2.0e33 * 6.0e23)
Mass2FromTwoComp(x) = x[3]*(10.0^x[5]*(1-x[6]))*4*π*(10.0*3.2e7*500.0*1e5)^2 / (0.5 * 2.0e33 * 6.0e23)


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






#------------

# read the data

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


# Plot the pass through function

@time contourf(-23:0.125:-18,20:0.125:26,(x,y)->passthroughcombined((x),y,0.2;df=0.01))

#@btime [passthroughcombined((x),y,0.2;df=0.001) for x=-23:0.125:-18,y=20:0.125:23]

#------------------------------------

# Best fit with dust

resbest = OneCompFitDriver(dataΔv,dataflx; SA=false, degraded = false,
                           fixeddust = false, initguess=[1.0,-0.04,1.5,22.8,-1.5])
restable =  @showprogress [
          OneCompFitDriver(dataΔv,dataflx; SA=false, degraded = true,
                           fixeddust = false, initguess=[1.0,-0.04,1.5,22.8,-1.5])
       for i ∈ 1:100]

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

scatter(10 .^ log10.(Mass), dust)


# Best with with no dust


resbest = OneCompFitDriver(dataΔv,dataflx; SA=false, degraded = false,
                           fixeddust = true, initguess=[1.0,-0.04,0.2,22.8,-6.5])
restable=[]
#restable =  @time @threads
@time @threads for i ∈ 1:120
    push!(restable,OneCompFitDriver(dataΔv,dataflx; SA=false, degraded = true,
                    fixeddust = true, initguess=[1.0,-0.04,0.2,22.8,-6.5]))
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

scatter(10 .^ log10.(Mass), dust)


# Best fit two components with dust


startguesses = [restable[i][2] for i ∈ 1:length(restable)]

resbest = TwoCompFitDriver(dataΔv,dataflx; SA=false, degraded = false,
                           fixeddust = false, initguess=[1.0,-0.04,1.5,21.5,22.9,0.5,-1.5])

restable=[];
@showprogress for i ∈ 1:200
                 push!(restable,TwoCompFitDriver(dataΔv,dataflx; SA=false, degraded = true,
                                                 fixeddust = false,initguess=[1.0,-0.04,1.5,21.5,22.9,0.5,-1.5] )) #startguesses[rand(100)]))
               end

χ20   = resbest[1]
Mass0 = resbest[3]
dust0 = resbest[2][7]
χ2   = [restable[i][1] for i ∈ 1:length(restable)]
Mass = [restable[i][3] for i ∈ 1:length(restable)]
Mass1 = [Mass1FromTwoComp(restable[i][2]) for i ∈ 1:length(restable)]
Mass2 = [Mass2FromTwoComp(restable[i][2]) for i ∈ 1:length(restable)]
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

scatter( log10.(Mass), 10.0 .^dust)


scatter( log10.(pairmax(Mass1,Mass2)),log10.(pairmin(Mass1,Mass2)))

scatter( log10.(Mass1),log10.(Mass2))


pairmax(v1,v2) = map((x,y)->max(x,y),v1,v2)
pairmin(v1,v2) = map((x,y)->min(x,y),v1,v2)
pairf(v1,v2,f) = map((x,y,f)-> x>y ? f : (1.0-f) ,v1,v2,f)

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
