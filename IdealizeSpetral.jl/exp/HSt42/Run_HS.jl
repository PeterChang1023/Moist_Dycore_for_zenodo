using JGCM
import Dates 

include("HS.jl")

start_time = Dates.now() 
### latent heat release ###
touch("Latent_heat.txt")
L_file        = open("Latent_heat.txt", "r")
L_int         = read(L_file, String)
L_float       = parse(Int64, L_int)
L             = L_float / 100
@info "L:", L
close(L_file)
path = "HSt42_" * L_int * "/"

### load interval_day ###
touch(path * "day_interval.txt")
file       = open(path * "day_interval.txt", "r")
data       = read(file, String)
end_day    = parse(Int64, data)
spinup_day = 0
close(file)
### load first day file name ###
touch(path * "firstday_file.txt")
firstday_file        = open(path * "firstday_file.txt", "r")
firstday             = read(firstday_file, String)
@info "filename:", firstday
warm_start_file_name =  firstday 
initial_day          =  end_day   # In this version, it would start at the final day of the warmstart.dat
close(firstday_file)


physics_params = Dict{String,Float64}("σ_b"=>0.7, "k_f" => 1.0, "k_a" => 1.0/40.0, "k_s" => 1.0/4.0, "ΔT_y" => 60.0, "Δθ_z" => 10.0) 
op_man = Atmos_Spectral_Dynamics_Main(physics_params, end_day, spinup_day, L)


Finalize_Output!(op_man, path * "warmstart_$(L_float).dat", path * "all_L"*L_int*".dat") 


### time ###
final_time = Dates.now() 
all_time = final_time - start_time
@info Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(all_time))) 
############



