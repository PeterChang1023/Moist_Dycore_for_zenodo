# Moist_Dycore_for_zenodo
Associated Ref:
https://github.com/CliMA/IdealizedSpectralGCM.jl/tree/master

## Install JGCM
(1) Install Julia under the path 
`/Moist_Dycore_for_zenodo/IdealizeSpetral.jl/`
(2) Switch to developer mode
```
julia> ]
pkg> dev .
```
(3) Load model
```
julia> using JGCM
```
(4) If it's necessary, one can delete the package
```
julia> ]
pkg> rm JGCM 
```

# Run Experiments
Before you run the model for test, please make sure you're at `/Moist_Dycore_for_zenodo/IdealizeSpetral.jl/exp/HSt42/`
```
julia For_test_Run_HS.jl
```
For extended-period experiment
```
sh PR0_long_run.sh
```
This script is designed for long-duration simulations, producing output results at intervals specified by 'space_day'. 
It helps prevent excessive memory usage by periodically saving results instead of waiting until the end of the simulation.

And, L represent latent heat release efficiency, you can adjust L=0 to 50 (meaing L=0 ~ 0.5 to avoid point in code).
