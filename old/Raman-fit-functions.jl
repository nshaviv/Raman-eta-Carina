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
