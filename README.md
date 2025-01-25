# Moist_Dycore_for_zenodo
Strong recommand visit:
https://github.com/CliMA/IdealizedSpectralGCM.jl/tree/master

## Install JGCM
first install Julia, then 
please make sure your dirtionary is at 
`/Moist_Dycore_for_zenodo/IdealizeSpetral.jl/`
And, type
```
julia> ]
pkg> dev .
```
Then load model
```
julia> using JGCM
```
When necessary, yoGCMu can delete the package
```
julia> ]
pkg> rm JGCM 
```

# Run Moist Dycore
Before you run the model for test, please make sure your dirtionary is at `/Moist_Dycore_for_zenodo/IdealizeSpetral.jl/exp/HSt42/`
```
julia For_test_Run_HS.jl
```
Then, if you want to run long period experiment
```
sh PR0_long_run.sh
```
This script is designed for long-duration simulations, producing output results at intervals specified by 'space_day'. 
It helps prevent excessive memory usage by periodically saving results instead of waiting until the end of the simulation.

And, L represent latent heat release efficiency, you can adjust L=0 to 50 (meaing L=0 ~ 0.5 to avoid point in code).
