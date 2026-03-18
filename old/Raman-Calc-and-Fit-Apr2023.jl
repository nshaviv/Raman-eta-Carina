# Define the functions and set up things
begin
    @info "Initializing"
    include("Raman-Scattering-Functions.jl")
end

# We shall work in velocity space relative to the Balmer alpha frequency (Hα)
# The files have the BR and Cross-section in units of the Lyman β doppler.
# The ratio between the two is 32/5

using CubicSplines

begin
    @info "Reading cross-section data"
    LyβBRRaw = CSV.read("data/Lyβ-BranchingRatio-corrected.csv", DataFrame, header=["Δv0", "br"])
    LyβBr = CubicSpline(LyβBRRaw.Δv0 .* 32 / 5, LyβBRRaw.br)

    LyβCSRaw = CSV.read("data/Lyβ-TotalCrossSection.csv", DataFrame, header=["Δv0", "cs"])
    LyβCS = CubicSpline(LyβCSRaw.Δv0 .* 32 / 5, LyβCSRaw.cs)
end;

Rayleigh_Larger_raw = CSV.read("data/Rayleigh-Larger-Range.csv", DataFrame, header=["λ0", "σ"])
Rayleigh_Larger = CubicSpline( (Rayleigh_Larger_raw.λ0 .- 1026)  .* (300000/1026 * (32 / 5)), Rayleigh_Larger_raw.σ)

Raman_Larger_raw = CSV.read("data/Raman-Larger-Range.csv", DataFrame, header=["λ0", "σ"])
Raman_Larger = CubicSpline( (Raman_Larger_raw.λ0 .- 1026)  .* (300000/1026 * (32 / 5)), Raman_Larger_raw.σ)

vs=-70000:1000:80000
plot(vs, Rayleigh_Larger.(vs))
#plot!(vs, Raman_Larger.(vs))
vs = -40000:1000:63000
#plot!(vs, LyβCS.(vs) .+ log10.(LyβBr.(vs)))
plot!(vs, LyβCS.(vs) .+ log10.(1 .- LyβBr.(vs)))
plot!(vs, LyβCS.(vs) .+ log10.(LyβBr.(vs)))

begin
    xTickLabels(vals) = [vals[i] % 20000 == 0 ? string.(Int(round(vals[i]))) : "" for i in 1:length(vals)]
    x_ticks = [T for T ∈ -40000:5000:63000]
    x_lbls = xTickLabels(x_ticks)
    yTickLabels(vals) = [vals[i] % 1 == 0 ? string.(Int(round(vals[i]))) : "" for i in 1:length(vals)]
    y_ticks = [T for T ∈ -24:0.2:-20]
    y_lbls = yTickLabels(y_ticks)

    rng = -40000:100:63000
    plot(rng, LyβCS.(rng), label="Total Cross Section", legend=:topleft, ylabel=L"\log_{10}(\sigma_i / (cm^2/g))",
        c=RGB(0.8, 0.3, 0.0), ylims=(-25, -20), xlims=(-30000, 30000), xticks=(x_ticks, x_lbls), xlabel=L"v (km/s)")
    plot!(rng, LyβCS.(rng) .+ log10.(LyβBr.(rng)), label="Partial Cross Section", legend=:topleft,
        c=RGB(0.0, 0.3, 0.8), ylims=(-25, -20), xlims=(-30000, 30000), xticks=(x_ticks, x_lbls),
             yticks=(y_ticks, y_lbls))
    plot!(twinx(), rng, LyβBr.(rng), color=:red, xticks=:none, 
         legend=:topright, label=L"\mathrm{B}\mathrm{r}(Ly\beta \rightarrow H\alpha)")
    plotstyle()
end
savefig("figures/Cross-section-vs-v.pdf")

#-------------------------------------------------------------------------- 
## First, plot the pass through function (that we derived from Mathematica)
#-------------------------------------------------------------------------- 


begin
    @info "generating pass through graphs"
    reflect = true
    plot_passthrough(1)
    savefig("figures/Reflect-ISM-dust.pdf")
    plot_passthrough(0.1)
    savefig("figures/Reflect-through-0.1-ISM-dust.pdf")
    plot_passthrough(0.01)
    savefig("figures/Reflect-through-0.01-ISM-dust.pdf")
    plot_passthrough(0.001)
    savefig("figures/Reflect-through-0.001-ISM-dust.pdf")
    plot_passthrough(0.0000001)
    savefig("figures/Reflect-through-No-ISM-dust.pdf")

    reflect = false
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
end

#-------------------------------------------------------------------------- 
## Next, plot solutions
#-------------------------------------------------------------------------- 

begin
    @info "generating model graphs"
    reflect = true
    plot_model_results(-7.0)
    savefig("figures/Model_Results_reflect_no-dust.pdf")
    plot_model_results(-3.0)

    savefig("figures/Model_Results_reflect_dust-pdf-m3.0.pdf")
    plot_model_results(-2.0)
    savefig("figures/Model_Results_reflect_dust-pdf-m2.0.pdf")
    plot_model_results(-1.0)
    savefig("figures/Model_Results_reflect_dust-pdf-m1.0.pdf")
    plot_model_results(-0.0)
    savefig("figures/Model_Results_reflect_dust-pdf-m0.0.pdf")

    reflect = false
    plot_model_results(-7.0)
    savefig("figures/Model_Results_pass_no-dust.pdf")
    plot_model_results(-3.0)
    savefig("figures/Model_Results_pass_dust-pdf-m3.0.pdf")
    plot_model_results(-2.0)
    savefig("figures/Model_Results_pass_dust-pdf-m2.0.pdf")
    plot_model_results(-1.0)
    savefig("figures/Model_Results_pass_dust-pdf-m1.0.pdf")
    plot_model_results(-0.0)
    savefig("figures/Model_Results_pass_dust-pdf-m0.0.pdf")
end

#-------------------------------------------------------------------------- 
#-------------------------------------------------------------------------- 
#-------------------------------------------------------------------------- 
# Now do the fits
#-------------------------------------------------------------------------- 

n_bootstrap = 20

reflect = false
results = [[-3, "irrelevant line to ensure results are of type Any"],
    ["Dust", "Σ1", "Σ2", "f", "M", "M1", "M2", "χ2", "data"]]

# Uncomment one line, to decide which data to fit
data_to_read = "Gemini-1"
#data_to_read = "Magelan-1"

begin

    ## First Read the data
    begin
        if data_to_read == "Gemini-1"
            data_plot = CSV.read("data/Gemini-etaCarinae-2014-11-fit.csv",
                DataFrame, header=["Δv", "flux"])
            offset = -1
            data = CSV.read("data/Gemini-etaCarinae-2014-11-short.csv",
                DataFrame, header=["Δv", "flux"])
            offset = -1
            data_name = "Gemini 11 2014"
        end
        if data_to_read == "Magelan-1"
            data_plot = CSV.read("data/Magelan-etaCarinae-2015-01-fit.csv",
                DataFrame, header=["Δv", "flux"])
            offset = 0
            data = CSV.read("data/Magelan-etaCarinae-2015-01-short.csv",
                DataFrame, header=["Δv", "flux"])
            offset = 0
            data_name = "Magellan 01 2015"
        end
        ## We cut out the data near zero velocity, 
        ## and also remove very large velocities (because there are other lines)
        dataΔv = 0 .+ 1000 .* data.Δv
        dataflx = data.flux
        dataΔv_plot = 0 .+ 1000 .* data_plot.Δv
        dataflx_plot = 1.0 .+ data_plot.flux .+ offset
        display(scatter(dataΔv_plot, dataflx_plot, title=data_to_read))
        errornorm = mean((dataflx[2:end] - dataflx[1:end-1]) .^ 2) / 2
    end


    #-----------------------
    # Best with with no dust
    #-----------------------

    resbest = OneCompFitDriver(dataΔv, dataflx; SA=false, degraded=false,
        fixeddust=true, initguess=[1.0, -0.04, 0.2, 22.8, -6.5])
    restable = []
    #restable =  @time @threads
    @time @threads for i ∈ 1:n_bootstrap
        push!(restable, OneCompFitDriver(dataΔv, dataflx; SA=false, degraded=true,
            fixeddust=true, initguess=resbest[2]))
    end

    begin
        χ20 = resbest[1]
        Mass0 = resbest[3]
        dust0 = resbest[2][5]
        Σ0 = resbest[2][4]
        χ2 = [restable[i][1] for i ∈ 1:length(restable)]
        Mass = [restable[i][3] for i ∈ 1:length(restable)]
        dust = [restable[i][2][5] for i ∈ 1:length(restable)]
        push!(results, ["-", Σ0, "-", "-", Mass0, "-", "-", χ20, data_to_read])
        pretty_table(permutedims(hcat(results[3:end]...)); header=results[2])
    end

    plot_fit(; plot_data=true, label="1 Σ, no dust", color=:red, lw=3, linestyle=:solid)

    #------------------------------------
    # Best fit one component with dust
    #------------------------------------

    resbest = OneCompFitDriver(dataΔv, dataflx; SA=false, degraded=false,
        fixeddust=false, initguess=[1.0, -0.04, 1.5, 22.8, -1.5])
    restable = []
    # 0.5 second per solution 
    @time @threads for i ∈ 1:n_bootstrap
        push!(restable, OneCompFitDriver(dataΔv, dataflx; SA=false, degraded=true,
            fixeddust=false, initguess=resbest[2]))
    end

    begin
        χ20 = resbest[1]
        Mass0 = resbest[3]
        dust0 = resbest[2][5]
        Σ0 = resbest[2][4]
        χ2 = [restable[i][1] for i ∈ 1:length(restable)]
        ind_good = [χ2 .< 1.5 * χ20]
        Mass = [restable[i][3] for i ∈ 1:length(restable)]
        dust = [restable[i][2][5] for i ∈ 1:length(restable)]
    end

    #              ["Dust","Σ1","Σ2","f", "M",  "M1", "M2", "χ2"]                
    push!(results, [dust0, Σ0, "-", "-", Mass0, "-", "-", χ20, data_to_read])

    plot_fit(; plot_data=false, label="1 Σ, with dust", color=:blue, lw=2, linestyle=:dash)

    #------------------------------------
    # Best fit two components with no dust
    #------------------------------------

    resbest = TwoCompFitDriver(dataΔv, dataflx; SA=false, degraded=false,
        fixeddust=true, initguess=[2.0, -0.04, 1.5, 21.5, 22.9, 0.5, -6.0, 0.0])

    restable = []
    @showprogress for i ∈ 1:20
        push!(restable, TwoCompFitDriver(dataΔv, dataflx; SA=false, degraded=true,
            fixeddust=true, initguess=resbest[2])) #startguesses[rand(100)]))
    end

    begin
        χ20 = resbest[1]
        Mass0 = resbest[3]
        Mass01 = Mass1FromTwoComp(resbest[2])
        Mass02 = Mass2FromTwoComp(resbest[2])
        dust0 = resbest[2][7]
        Σ01 = resbest[2][4]
        Σ02 = resbest[2][5]
        fbest = frac(resbest[2][6])
        χ2 = [restable[i][1] for i ∈ 1:length(restable)]
        Mass = [restable[i][3] for i ∈ 1:length(restable)]
        Mass1 = [Mass1FromTwoComp(restable[i][2]) for i ∈ 1:length(restable)]
        Mass2 = [Mass2FromTwoComp(restable[i][2]) for i ∈ 1:length(restable)]
        Σ1 = [restable[i][2][4] for i ∈ 1:length(restable)]
        Σ2 = [restable[i][2][5] for i ∈ 1:length(restable)]
        dust = [restable[i][2][7] for i ∈ 1:length(restable)]
        fs = [restable[i][2][6] for i ∈ 1:length(restable)]
        push!(results, ["-", Σ01, Σ02, fbest, Mass0, Mass01, Mass02, χ20, data_to_read])
    end

    plot_fit(; two_comp=true, plot_data=false, label="2 Σs, no dust", color=:green, lw=3, linestyle=:solid)


    #------------------------------------
    # Best fit two components with dust
    #------------------------------------

    resbest = TwoCompFitDriver(dataΔv, dataflx; SA=false, degraded=false,
        fixeddust=false,
        initguess=[2.0, -0.04, 1.5, 23.5, 20.9, -1.0, -2.5, 0.0])

    restable = []
    @showprogress for i ∈ 1:n_bootstrap
        push!(restable, TwoCompFitDriver(dataΔv, dataflx; SA=false, degraded=true,
            fixeddust=false, initguess=resbest[2])) #startguesses[rand(100)]))
    end

    begin
        χ20 = resbest[1]
        Mass0 = resbest[3]
        Mass01 = Mass1FromTwoComp(resbest[2])
        Mass02 = Mass2FromTwoComp(resbest[2])
        dust0 = resbest[2][7]
        Σ01 = resbest[2][4]
        Σ02 = resbest[2][5]
        fbest = frac(resbest[2][6])
        χ2 = [restable[i][1] for i ∈ 1:length(restable)]
        Mass = [restable[i][3] for i ∈ 1:length(restable)]
        Mass1 = [Mass1FromTwoComp(restable[i][2]) for i ∈ 1:length(restable)]
        Mass2 = [Mass2FromTwoComp(restable[i][2]) for i ∈ 1:length(restable)]
        Σ1 = [restable[i][2][4] for i ∈ 1:length(restable)]
        Σ2 = [restable[i][2][5] for i ∈ 1:length(restable)]
        dust = [restable[i][2][7] for i ∈ 1:length(restable)]
        fs = [restable[i][2][6] for i ∈ 1:length(restable)]

        push!(results, [dust0, Σ01, Σ02, fbest, Mass0, Mass01, Mass02, χ20, data_to_read])
        pretty_table(permutedims(hcat(results[3:end]...)); header=results[2])
    end

    plot_fit(; two_comp=true, plot_data=false, label="2 Σs, with dust", color=:orange, lw=2, linestyle=:dash)
    display(plot!(bottom_margin=8cm))

    plot!(xlabel=L"v (km/s)", ylabel=L" \textrm{N}\textrm{ormalized}~\textrm{F}\textrm{lux}")

end

pretty_table(permutedims(hcat(results[3:end]...)); header=results[2])

savefig("figures/Spectral_fits_to_$(replace(data_name," "=>"_")).pdf")
pretty_table(permutedims(hcat(results[3:end]...)); header=results[2], backend=Val(:latex))



