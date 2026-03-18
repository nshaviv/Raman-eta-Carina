using Plots               # Standard plotting packages
using Plots.PlotMeasures
using LaTeXStrings
using ColorSchemes


function ramancalc(pσ, pΣ, BR; df=0.01, reflect=true)
    aB = 0.75  # scattering fraction (albedo) for the Balmer lines
    aL = 0.3   # scattering fraction for the Lyman lines
    Σ = 10.0^pΣ
    
    κB         =  df*(1-aB)*3.7e-22          # Balmer frequency absorption cross-section 
    σB         =  df*  aB  *3.7e-22+6.6e-25  # Balmer frequnecy scattering cross-section
    
    σ_raman_tot = 10.0^pσ                 # The Raman scattering cross-section for Lyman photons
    κL_dust    =  df*(1-aL)*2.2e-21          # Lyman frequency absorption cross-section
    σL_dust_th =  df*  aL  *2.2e-21+6.6e-25  # Lyman frequency scattering cross-section
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

    denom = -(sqκB*κL*(sτB*(2*κB + σB) + 2*cτB*χB)*(χB - χL)*(χB + χL)*(cτL*(κL + σL) + sτL*χL))
    enum  = reflect == true ? 
       -(sqχB*σ_raman*(κL + σL)*(κB*κL*(-κB + κL - σB + σL) + sτB*χB*(-(cτL*κL*(-κB + κL + σL)) +
            sτL*(κB - κL)*χL) + cτB*κB*(cτL*κL*(κB - κL + σB - σL) + sτL*(κB - κL + σB)*χL))) :
       -(sqκB*sqχB*σ_raman*(κL + σL)*(sqχB*sτB*κL*(κB + κL + σL) + cτB*sqκB*κL*(κB + κL + σB + σL) 
           - cτL*sqκB*κL*(κB + κL + σB + σL) - sqκB*sτL*(κB + κL + σB)*χL))

    enum/denom
end

function plot_passthrough(dust_fraction; reflect=true)
    contourf(-25:0.125:-18,20:0.125:26,
        (x,y)->passthroughcombined((x),y,0.2;df=dust_fraction, reflect=reflect),
        levels=[0.0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1],
        legend=:none,
        xlabel = L"\sigma (cm^2/g)",ylabel=L"\Sigma ~(g/cm^2)",
        seriescolor=cgrad(:heat,[0.0,1.1,0.925]),
        guidefont=font(24), left_margin=35px, bottom_margin=20px, 
        size=(1000,800),frame=:box, xtickfontsize=18, ytickfontsize=18)
    contour!(-25:0.125:-18,20:0.125:26,
        (x,y)->passthroughcombined((x),y,0.2;df=dust_fraction, reflect=reflect),
        levels=[0.0,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.1],
        legend=:none,
        contour_labels = true,
        seriescolor=:black)
    annotate!(-24.8, 20.25, text("Dust/ISM = $dust_fraction", :blue, :left, 18))
    plotstyle()
end 


# Remove NaNs from results:
fixNaN(num) = isnan(num) ? 0.5 : num

# The expression for the pass through flux has very large exponents. The simplest solution is to 
# do the calculation with higher precision in those regions of phase space that require it.
# We use Julia's "native" arbitrary float precision option (There are packages, such as QuadMath which doesn't 
# run on ARM (new macs) or DoubleFloats, which doesn't have enough exponent bits). 
setprecision(256);
function passthroughcombined(pσ, pΣin, BR; df=0.01, reflect=true)
    if reflect == true
        pΣ = min(pΣin,24-log10(df))
#        if pσ + pΣ > 1.8 && pσ >-23 return ramancalc(-18.8, 20.0, BR; df=df, reflect=reflect) end 
        if pσ + pΣ < 0 return ramancalc(pσ, pΣ, BR; df=df, reflect=reflect)  end
        Float64(ramancalc(BigFloat(pσ), pΣ, BR; df=df, reflect=reflect)) |> fixNaN 
    else
        pΣ = min(pΣin,23-log10(df))
        if pσ + pΣ < 2.8 return ramancalc(pσ, pΣ, BR; df=df, reflect=reflect) |> fixNaN end
        if pσ + pΣ < 4.4 return Float64(ramancalc(BigFloat(pσ), pΣ, BR; df=df, reflect=reflect)) |> fixNaN end
        return(Float64(ramancalc(BigFloat(4.4-pΣ), pΣ, BR; df=df, reflect=reflect)) |> fixNaN)
    end
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

plot_passthrough(0.0001, reflect=true)