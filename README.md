# Raman-eta-Carinae

Julia calculations and figure-generation scripts for the paper _On the conditions for wide Raman-scattering "emission" lines and the wide spectral features in η-Carinae's light echoes_ by Shlomi Pistinner and Nir J. Shaviv.

The repository implements the analytic Raman-scattering model used to study broad Hα wings produced by UV continuum photons scattering in neutral hydrogen, including the effects of dust. In the paper context, the calculations are applied to η Carinae light-echo spectra as an alternative to interpreting the very broad wings as extremely fast ejecta.

## Repository contents

- `Raman-Calc-and-Fit-July2024.jl`: main driver script. It reads the Raman/Rayleigh cross-section tables, generates the model figures, and fits the observed Gemini and Magellan spectra.
- `Raman-Scattering-Functions.jl`: model, fitting, plotting, and utility functions used by the main script.
- `data/`: cross-section tables and the processed observational spectra used in the fits.
- `figures/`: output figures produced by the scripts, including model grids and spectral-fit plots.
- `Velocity_distribution/`: auxiliary calculations for angular/velocity-distribution effects.
- `Mathematica/`: derivations and symbolic work notebooks. These are marked as documentation in `.gitattributes` so they do not dominate GitHub language statistics.

## Requirements

The code is written in Julia and currently assumes the needed packages are installed in the active environment.

Main packages used by the Raman calculations:

- `CSV`
- `DataFrames`
- `Plots`
- `Optim`
- `ForwardDiff`
- `Interpolations`
- `ProgressMeter`
- `LaTeXStrings`
- `ColorSchemes`
- `PrettyTables`
- `Measurements`
- `Statistics`
- `CubicSplines`
- `Printf`
- `BenchmarkTools`

The velocity-distribution script also uses:

- `QuadGK`

If you are setting up a fresh Julia environment, install the packages from the Julia REPL with:

```julia
using Pkg
Pkg.add([
    "CSV",
    "DataFrames",
    "Plots",
    "Optim",
    "ForwardDiff",
    "Interpolations",
    "ProgressMeter",
    "LaTeXStrings",
    "ColorSchemes",
    "PrettyTables",
    "Measurements",
    "Statistics",
    "CubicSplines",
    "BenchmarkTools",
    "QuadGK",
])
```

## Running the main calculations

From the repository root:

```bash
julia Raman-Calc-and-Fit-July2024.jl
```

This script:

- loads the Lyβ and Lyγ branching-ratio and cross-section tables from `data/`
- generates the cross-section and model figures written to `figures/`
- fits the Raman-scattering model to the Gemini 2014-11 and Magellan 2015-01 spectra
- prints fit summaries and bootstrap-based parameter estimates to the terminal

To run the auxiliary velocity-distribution calculation:

```bash
cd Velocity_distribution
julia velocity_epsilon.jl
```

## Notes

- The scripts are currently organized as research code rather than as a packaged Julia project.
- The main workflow expects to be run from the repository root so that relative paths to `data/` and `figures/` resolve correctly.
- The `old/` directory is intentionally excluded from version control.
