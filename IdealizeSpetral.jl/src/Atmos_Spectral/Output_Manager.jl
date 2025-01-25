export Output_Manager, Update_Output!, Finalize_Output!
export Lat_Lon_Pcolormesh, Zonal_Mean, Sigma_Zonal_Mean_Pcolormesh, Sigma_Zonal_Mean_Contourf

mutable struct Output_Manager
    nλ::Int64
    nθ::Int64
    nd::Int64
    n_day::Int64
    
    day_to_sec::Int64
    start_time::Int64
    end_time::Int64
    current_time::Int64
    spinup_day::Int64

    λc::Array{Float64, 1}
    θc::Array{Float64, 1}
    σc::Array{Float64, 1}

    n_daily_mean::Array{Float64, 1}    
    ##########################################################################
    # specral vor 
    spe_vor_c_xyzt::Array{ComplexF64,4}
    spe_vor_p_xyzt::Array{ComplexF64,4}
    
    # specral div
    spe_div_c_xyzt::Array{ComplexF64,4}
    spe_div_p_xyzt::Array{ComplexF64,4}

    # specral height or surface pressure
    spe_lnps_c_xyzt::Array{ComplexF64,4}
    spe_lnps_p_xyzt::Array{ComplexF64,4}

    # specral temperature
    spe_t_c_xyzt::Array{ComplexF64,4}
    spe_t_p_xyzt::Array{ComplexF64,4}

    # specral tracer
    spe_tracers_c_xyzt::Array{ComplexF64,4}
    spe_tracers_p_xyzt::Array{ComplexF64,4}
    ##########################################################################
    # grid w-e velocity
    grid_u_n_xyzt::Array{Float64, 4}
    grid_u_c_xyzt::Array{Float64, 4}
    grid_u_p_xyzt::Array{Float64, 4}
    
    # grid n-s velocity
    grid_v_n_xyzt::Array{Float64, 4}
    grid_v_c_xyzt::Array{Float64, 4}
    grid_v_p_xyzt::Array{Float64, 4}

    # grid surface pressure
    grid_ps_c_xyzt::Array{Float64,4}
    grid_ps_p_xyzt::Array{Float64,4}

    # grid temperature
    grid_t_n_xyzt::Array{Float64, 4}
    grid_t_c_xyzt::Array{Float64, 4}
    grid_t_p_xyzt::Array{Float64, 4}

    # grid tracer
    grid_tracers_n_xyzt::Array{Float64,4}
    grid_tracers_c_xyzt::Array{Float64,4}
    grid_tracers_p_xyzt::Array{Float64,4}

    grid_tracers_diff_xyzt::Array{Float64,4}
    grid_δtracers_xyzt::Array{Float64,4}
    
    ##########################################################################
    factor1_xyzt::Array{Float64,4}
    factor2_xyzt::Array{Float64,4}
    factor3_xyzt::Array{Float64,4}
    factor4_xyzt::Array{Float64,4}

    # pressure 
    grid_p_full_xyzt::Array{Float64,4}
    grid_p_half_xyzt::Array{Float64,4}
    
    # geopotential
    grid_geopots_xyzt::Array{Float64,4}
    
    # Memory contrainer for temporal variables
    # vor
    grid_vor_xyzt::Array{Float64,4}
    
    # div
    grid_div_xyzt::Array{Float64,4}

    # w-e velocity tendency
    grid_δu_xyzt::Array{Float64,4}
    
    # n-s velocity tendency
    grid_δv_xyzt::Array{Float64,4}
    #######################################################################
    # equilibrium temperature in HS_Forcing
    convection_xyzt::Array{Float64,4}

    grid_z_full_xyzt::Array{Float64,4}
    grid_w_full_xyzt::Array{Float64,4}
    #######################################################################
    grid_tracers_c_max_Tiffany_xyzt::Array{Float64,4}
    grid_tracers_c_max_xyzt::Array{Float64,4}
    
    
    

    

end

function Output_Manager(mesh::Spectral_Spherical_Mesh, vert_coord::Vert_Coordinate, start_time::Int64, end_time::Int64, spinup_day::Int64, num_grid_tracters::Int64=1, num_spe_tracters::Int64=1) ### By CJY2
    nλ = mesh.nλ
    nθ = mesh.nθ
    nd = mesh.nd
    
    day_to_sec = 86400
    current_time = start_time

    λc = mesh.λc
    θc = mesh.θc

    #todo definition of sigma coordinate
    bk = vert_coord.bk
    σc = (bk[2:nd+1] + bk[1:nd])/2.0
  
    n_day = Int64((end_time - start_time)/ (day_to_sec/4) )
    n_daily_mean = zeros(Float64, n_day)

    grid_geopots_xyzt = zeros(Float64, nλ, nθ, 1, n_day)
    num_fourier, nθ, nd = 42, 64, 20
    num_spherical = num_fourier + 1
    #########################################################
    # specral vor 
    spe_vor_c_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    spe_vor_p_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)

    # specral div
    spe_div_c_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    spe_div_p_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)

    # specral height or surface pressure
    spe_lnps_c_xyzt = zeros(ComplexF64, num_fourier+1, num_spherical+1, 1, n_day)
    spe_lnps_p_xyzt = zeros(ComplexF64, num_fourier+1, num_spherical+1, 1, n_day)

    # specral temperature
    spe_t_c_xyzt    = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    spe_t_p_xyzt    = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)

    # specral tracer
    spe_tracers_c_xyzt = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    spe_tracers_p_xyzt = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    ##########################################################
    # grid w-e velocity
    grid_u_n_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)
    grid_u_c_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)
    grid_u_p_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)

    # grid n-s velocity
    grid_v_n_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)
    grid_v_c_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)
    grid_v_p_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)

    # grid surface pressure
    grid_ps_c_xyzt = zeros(Float64, nλ,  nθ, 1, n_day)
    grid_ps_p_xyzt = zeros(Float64, nλ,  nθ, 1, n_day)

    # grid temperature
    grid_t_n_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_t_c_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_t_p_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day)

    # grid tracer
    grid_tracers_n_xyzt     = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_tracers_c_xyzt     = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_tracers_p_xyzt     = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_tracers_diff_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_δtracers_xyzt      = zeros(Float64, nλ,  nθ, nd, n_day)
    ############################################################
    # factors
    factor1_xyzt = zeros(Float64, nλ,  nθ, nd, n_day) 
    factor2_xyzt = zeros(Float64, nλ,  nθ, nd, n_day) 
    factor3_xyzt = zeros(Float64, nλ,  nθ, nd, n_day) 
    factor4_xyzt = zeros(Float64, nλ,  nθ, nd, n_day) 

    # pressure
    grid_p_full_xyzt = zeros(Float64, nλ,  nθ, nd  , n_day) 
    grid_p_half_xyzt = zeros(Float64, nλ,  nθ, nd+1, n_day) 
    
    # geopotential
    grid_geopots_xyzt = zeros(Float64, nλ,  nθ, 1, n_day)
    ################################################
    # spe_δvor_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    grid_vor_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day)
    # grid_δvor_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
    
    
    # spe_δdiv_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
    grid_div_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day) 
    # grid_δdiv_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
    
    
    # Tendency
    grid_δu_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day) 
    grid_δv_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day) 

    convection_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)

    grid_z_full_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_w_full_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)

    grid_tracers_c_max_Tiffany_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
    grid_tracers_c_max_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
    
    
    
    
    
    
    Output_Manager(nλ, nθ, nd, n_day,
    day_to_sec, start_time, end_time, current_time, spinup_day,
    λc, θc, σc, n_daily_mean, spe_vor_c_xyzt, spe_vor_p_xyzt, spe_div_c_xyzt, spe_div_p_xyzt, spe_lnps_c_xyzt, spe_lnps_p_xyzt, spe_t_c_xyzt, spe_t_p_xyzt, spe_tracers_c_xyzt, spe_tracers_p_xyzt, grid_u_n_xyzt, grid_u_c_xyzt, grid_u_p_xyzt, grid_v_n_xyzt, grid_v_c_xyzt, grid_v_p_xyzt, grid_ps_c_xyzt, grid_ps_p_xyzt, grid_t_n_xyzt, grid_t_c_xyzt, grid_t_p_xyzt, grid_tracers_n_xyzt, grid_tracers_c_xyzt, grid_tracers_p_xyzt, grid_tracers_diff_xyzt, grid_δtracers_xyzt, factor1_xyzt, factor2_xyzt, factor3_xyzt, factor4_xyzt, grid_p_full_xyzt, grid_p_half_xyzt, grid_geopots_xyzt, grid_vor_xyzt, grid_div_xyzt, grid_δu_xyzt, grid_δv_xyzt, convection_xyzt, grid_z_full_xyzt, grid_w_full_xyzt, grid_tracers_c_max_Tiffany_xyzt, grid_tracers_c_max_xyzt)
end

function Update_Output!(output_manager::Output_Manager, dyn_data::Dyn_Data, current_time::Int64)
    @assert(current_time > output_manager.current_time)
    output_manager.current_time = current_time
    day_to_sec, start_time, n_day = output_manager.day_to_sec, output_manager.start_time, output_manager.n_day

    n_daily_mean = output_manager.n_daily_mean
    ############################################################
    # specral vor 
    spe_vor_c_xyzt  = output_manager.spe_vor_c_xyzt
    spe_vor_p_xyzt  = output_manager.spe_vor_p_xyzt
    
    # specral div
    spe_div_c_xyzt  = output_manager.spe_div_c_xyzt
    spe_div_p_xyzt  = output_manager.spe_div_p_xyzt

    # specral height or surface pressure
    spe_lnps_c_xyzt = output_manager.spe_lnps_c_xyzt
    spe_lnps_p_xyzt = output_manager.spe_lnps_p_xyzt
    
    # specral temperature
    spe_t_c_xyzt = output_manager.spe_t_c_xyzt
    spe_t_p_xyzt = output_manager.spe_t_p_xyzt
    
    # specral tracer
    spe_tracers_c_xyzt = output_manager.spe_tracers_c_xyzt
    spe_tracers_p_xyzt = output_manager.spe_tracers_p_xyzt
    ############################################################
    # grid w-e velocity
    grid_u_n_xyzt  = output_manager.grid_u_n_xyzt 
    grid_u_c_xyzt  = output_manager.grid_u_c_xyzt 
    grid_u_p_xyzt  = output_manager.grid_u_p_xyzt 
    
    # grid n-s velocity
    grid_v_n_xyzt  = output_manager.grid_v_n_xyzt
    grid_v_c_xyzt  = output_manager.grid_v_c_xyzt
    grid_v_p_xyzt  = output_manager.grid_v_p_xyzt

    # grid surface pressure
    grid_ps_c_xyzt = output_manager.grid_ps_c_xyzt
    grid_ps_p_xyzt = output_manager.grid_ps_p_xyzt

    # grid temperature
    grid_t_n_xyzt  = output_manager.grid_t_n_xyzt
    grid_t_c_xyzt  = output_manager.grid_t_c_xyzt
    grid_t_p_xyzt  = output_manager.grid_t_p_xyzt

    # grid tracer
    grid_tracers_n_xyzt = output_manager.grid_tracers_n_xyzt
    grid_tracers_c_xyzt = output_manager.grid_tracers_c_xyzt
    grid_tracers_p_xyzt = output_manager.grid_tracers_p_xyzt

    grid_tracers_diff_xyzt = output_manager.grid_tracers_diff_xyzt
    grid_δtracers_xyzt     = output_manager.grid_δtracers_xyzt
    ############################################################
    # factors
    factor1_xyzt = output_manager.factor1_xyzt
    factor2_xyzt = output_manager.factor2_xyzt
    factor3_xyzt = output_manager.factor3_xyzt
    factor4_xyzt = output_manager.factor4_xyzt

    
    # pressure
    grid_p_full_xyzt = output_manager.grid_p_full_xyzt
    grid_p_half_xyzt = output_manager.grid_p_half_xyzt
    
    # geopotential
    grid_geopots_xyzt = output_manager.grid_geopots_xyzt
    ############################################################
    grid_vor_xyzt  = output_manager.grid_vor_xyzt
    
    grid_div_xyzt  = output_manager.grid_div_xyzt
    
    # Tendency
    grid_δu_xyzt   = output_manager.grid_δu_xyzt
    grid_δv_xyzt   = output_manager.grid_δv_xyzt

    convection_xyzt = output_manager.convection_xyzt

    grid_z_full_xyzt = output_manager.grid_z_full_xyzt
    grid_w_full_xyzt = output_manager.grid_w_full_xyzt

    grid_tracers_c_max_Tiffany_xyzt = output_manager.grid_tracers_c_max_Tiffany_xyzt
    grid_tracers_c_max_xyzt = output_manager.grid_tracers_c_max_xyzt
    
    i_day = Int(div(current_time - start_time - 1, day_to_sec/4) + 1)

    if(i_day > n_day)
        @info "Warning: i_day > n_day in Output_Manager:Update!"
        return 
    end
    
    ############################################################
    # specral vor 
    spe_vor_c_xyzt[:,:,:,i_day] .= dyn_data.spe_vor_c[:,:,:]
    spe_vor_p_xyzt[:,:,:,i_day] .= dyn_data.spe_vor_p[:,:,:]
        
    # specral div
    spe_div_c_xyzt[:,:,:,i_day] .= dyn_data.spe_div_c[:,:,:]
    spe_div_p_xyzt[:,:,:,i_day] .= dyn_data.spe_div_p[:,:,:]
    
    # specral height or surface pressure
    spe_lnps_c_xyzt[:,:,:,i_day] .= dyn_data.spe_lnps_c[:,:,:]
    spe_lnps_p_xyzt[:,:,:,i_day] .= dyn_data.spe_lnps_p[:,:,:]
    
    # specral temperature
    spe_t_c_xyzt[:,:,:,i_day] .= dyn_data.spe_t_c[:,:,:]
    spe_t_p_xyzt[:,:,:,i_day] .= dyn_data.spe_t_p[:,:,:]
        
    # specral tracer
    spe_tracers_c_xyzt[:,:,:,i_day] .= dyn_data.spe_tracers_c[:,:,:]
    spe_tracers_p_xyzt[:,:,:,i_day] .= dyn_data.spe_tracers_p[:,:,:]
    ############################################################
    # grid w-e velocity
    grid_u_n_xyzt[:,:,:,i_day]  .= dyn_data.grid_u_n[:,:,:] 
    grid_u_c_xyzt[:,:,:,i_day]  .= dyn_data.grid_u_c[:,:,:] 
    grid_u_p_xyzt[:,:,:,i_day]  .= dyn_data.grid_u_p[:,:,:]
    
    # grid n-s velocity
    grid_v_n_xyzt[:,:,:,i_day]  .= dyn_data.grid_v_n[:,:,:]
    grid_v_c_xyzt[:,:,:,i_day]  .= dyn_data.grid_v_c[:,:,:]
    grid_v_p_xyzt[:,:,:,i_day]  .= dyn_data.grid_v_p[:,:,:]

    # grid surface pressure
    grid_ps_c_xyzt[:,:,:,i_day] .= dyn_data.grid_ps_c[:,:,:]
    grid_ps_p_xyzt[:,:,:,i_day] .= dyn_data.grid_ps_p[:,:,:]

    # grid temperature
    grid_t_n_xyzt[:,:,:,i_day]  .= dyn_data.grid_t_n[:,:,:]
    grid_t_c_xyzt[:,:,:,i_day]  .= dyn_data.grid_t_c[:,:,:]
    grid_t_p_xyzt[:,:,:,i_day]  .= dyn_data.grid_t_p[:,:,:]

    # grid tracer
    grid_tracers_n_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_n[:,:,:]
    grid_tracers_c_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_c[:,:,:]
    grid_tracers_p_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_p[:,:,:]

    grid_tracers_diff_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_diff[:,:,:]
    grid_δtracers_xyzt[:,:,:,i_day]     .= dyn_data.grid_δtracers[:,:,:]
    ############################################################
    # factors
    factor1_xyzt[:,:,:,i_day] .= dyn_data.factor1[:,:,:]
    factor2_xyzt[:,:,:,i_day] .= dyn_data.factor2[:,:,:]
    factor3_xyzt[:,:,:,i_day] .= dyn_data.factor3[:,:,:]
    factor4_xyzt[:,:,:,i_day] .= dyn_data.factor4[:,:,:]

    
    # pressure
    grid_p_full_xyzt[:,:,:,i_day] .= dyn_data.grid_p_full[:,:,:]
    grid_p_half_xyzt[:,:,:,i_day] .= dyn_data.grid_p_half[:,:,:]
    
    # geopotential
    grid_geopots_xyzt[:,:,:,i_day] .= dyn_data.grid_geopots[:,:,:]
    ############################################################
    grid_vor_xyzt[:,:,:,i_day]  .= dyn_data.grid_vor[:,:,:]
    grid_div_xyzt[:,:,:,i_day]  .= dyn_data.grid_div[:,:,:]
    
    # Tendency
    grid_δu_xyzt[:,:,:,i_day]   .= dyn_data.grid_δu[:,:,:]
    grid_δv_xyzt[:,:,:,i_day]   .= dyn_data.grid_δv[:,:,:]

    convection_xyzt[:,:,:,i_day] .= dyn_data.convection[:,:,:]

    grid_z_full_xyzt[:,:,:,i_day] .= dyn_data.grid_z_full[:,:,:]
    grid_w_full_xyzt[:,:,:,i_day] .= dyn_data.grid_w_full[:,:,:]

    grid_tracers_c_max_Tiffany_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_c_max_Tiffany[:,:,:]
    grid_tracers_c_max_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_c_max[:,:,:]
    

    n_daily_mean[i_day] += 1
end

function Finalize_Output!(output_manager::Output_Manager, save_file_name::String = "None", mean_save_file_name::String = "None")

    n_day = output_manager.n_day

    ############################################################    
    # specral vor 
    spe_vor_c_xyzt  = output_manager.spe_vor_c_xyzt
    spe_vor_p_xyzt  = output_manager.spe_vor_p_xyzt
    
    # specral div
    spe_div_c_xyzt  = output_manager.spe_div_c_xyzt
    spe_div_p_xyzt  = output_manager.spe_div_p_xyzt

    # specral height or surface pressure
    spe_lnps_c_xyzt = output_manager.spe_lnps_c_xyzt
    spe_lnps_p_xyzt = output_manager.spe_lnps_p_xyzt
    
    # specral temperature
    spe_t_c_xyzt    = output_manager.spe_t_c_xyzt
    spe_t_p_xyzt    = output_manager.spe_t_p_xyzt
    
    # specral tracer
    spe_tracers_c_xyzt = output_manager.spe_tracers_c_xyzt
    spe_tracers_p_xyzt = output_manager.spe_tracers_p_xyzt
    ############################################################
    # grid w-e velocity
    grid_u_n_xyzt  = output_manager.grid_u_n_xyzt 
    grid_u_c_xyzt  = output_manager.grid_u_c_xyzt 
    grid_u_p_xyzt  = output_manager.grid_u_p_xyzt 
    
    # grid n-s velocity
    grid_v_n_xyzt  = output_manager.grid_v_n_xyzt
    grid_v_c_xyzt  = output_manager.grid_v_c_xyzt
    grid_v_p_xyzt  = output_manager.grid_v_p_xyzt

    # grid surface pressure
    grid_ps_c_xyzt = output_manager.grid_ps_c_xyzt
    grid_ps_p_xyzt = output_manager.grid_ps_p_xyzt

    # grid temperature
    grid_t_n_xyzt  = output_manager.grid_t_n_xyzt
    grid_t_c_xyzt  = output_manager.grid_t_c_xyzt
    grid_t_p_xyzt  = output_manager.grid_t_p_xyzt

    # grid tracer
    grid_tracers_n_xyzt = output_manager.grid_tracers_n_xyzt  
    grid_tracers_c_xyzt = output_manager.grid_tracers_c_xyzt
    grid_tracers_p_xyzt = output_manager.grid_tracers_p_xyzt

    grid_tracers_diff_xyzt = output_manager.grid_tracers_diff_xyzt
    grid_δtracers_xyzt     = output_manager.grid_δtracers_xyzt
    ############################################################
    # factors
    factor1_xyzt = output_manager.factor1_xyzt
    factor2_xyzt = output_manager.factor2_xyzt
    factor3_xyzt = output_manager.factor3_xyzt
    factor4_xyzt = output_manager.factor4_xyzt


    # pressure
    grid_p_full_xyzt = output_manager.grid_p_full_xyzt
    grid_p_half_xyzt = output_manager.grid_p_half_xyzt
    
    # geopotential
    grid_geopots_xyzt = output_manager.grid_geopots_xyzt
    ############################################################
    grid_vor_xyzt  = output_manager.grid_vor_xyzt    
    grid_div_xyzt  = output_manager.grid_div_xyzt
    
    # Tendency
    grid_δu_xyzt   = output_manager.grid_δu_xyzt
    grid_δv_xyzt   = output_manager.grid_δv_xyzt

    convection_xyzt = output_manager.convection_xyzt

    grid_z_full_xyzt = output_manager.grid_z_full_xyzt
    grid_w_full_xyzt = output_manager.grid_w_full_xyzt
    grid_tracers_c_max_Tiffany_xyzt = output_manager.grid_tracers_c_max_Tiffany_xyzt
    grid_tracers_c_max_xyzt = output_manager.grid_tracers_c_max_xyzt
    
    ##############################################################################
    spe_vor_c_final = spe_vor_c_xyzt[:,:,:,end]
    spe_vor_p_final = spe_vor_p_xyzt[:,:,:,end]
    spe_div_c_final = spe_div_c_xyzt[:,:,:,end]
    spe_div_p_final = spe_div_p_xyzt[:,:,:,end]

    spe_lnps_c_final = spe_lnps_c_xyzt[:,:,:,end]
    spe_lnps_p_final = spe_lnps_p_xyzt[:,:,:,end]

    spe_t_c_final = spe_t_c_xyzt[:,:,:,end]
    spe_t_p_final = spe_t_p_xyzt[:,:,:,end]

    spe_tracers_c_final = spe_tracers_c_xyzt[:,:,:,end]
    spe_tracers_p_final = spe_tracers_p_xyzt[:,:,:,end]

    grid_tracers_n_final = grid_tracers_n_xyzt[:,:,:,end]

    grid_u_n_final = grid_u_n_xyzt[:,:,:,end]
    grid_u_c_final = grid_u_c_xyzt[:,:,:,end]
    grid_u_p_final = grid_u_p_xyzt[:,:,:,end]

    grid_v_n_final = grid_v_n_xyzt[:,:,:,end]
    grid_v_c_final = grid_v_c_xyzt[:,:,:,end]
    grid_v_p_final = grid_v_p_xyzt[:,:,:,end]

    grid_ps_c_final = grid_ps_c_xyzt[:,:,:,end]
    grid_ps_p_final = grid_ps_p_xyzt[:,:,:,end]

    grid_t_n_final = grid_t_n_xyzt[:,:,:,end]
    grid_t_c_final = grid_t_c_xyzt[:,:,:,end]
    grid_t_p_final = grid_t_p_xyzt[:,:,:,end]

    grid_tracers_c_final = grid_tracers_c_xyzt[:,:,:,end]
    grid_tracers_p_final = grid_tracers_p_xyzt[:,:,:,end]

    grid_tracers_diff_final = grid_tracers_diff_xyzt[:,:,:,end]
    grid_δtracers_final = grid_δtracers_xyzt[:,:,:,end]

    grid_p_full_final = grid_p_full_xyzt[:,:,:,end]
    grid_p_half_final = grid_p_half_xyzt[:,:,:,end]

    grid_geopots_final = grid_geopots_xyzt[:,:,:,end]

    grid_vor_final = grid_vor_xyzt[:,:,:,end]
    grid_div_final = grid_div_xyzt[:,:,:,end]

    grid_δu_final = grid_δu_xyzt[:,:,:,end]
    grid_δv_final = grid_δv_xyzt[:,:,:,end]

    convection_final = convection_xyzt[:,:,:,end]

    grid_w_full_final = grid_w_full_xyzt[:,:,:,end]
    
    
    
    
    
    
    
    if save_file_name != "None"
        @save save_file_name spe_vor_c_final spe_vor_p_final spe_div_c_final spe_div_p_final spe_lnps_c_final spe_lnps_p_final spe_t_c_final spe_t_p_final spe_tracers_c_final spe_tracers_p_final grid_tracers_n_final grid_u_n_final grid_u_c_final grid_u_p_final grid_v_n_final grid_v_c_final grid_v_p_final grid_ps_c_final grid_ps_p_final grid_t_n_final grid_t_c_final grid_t_p_final grid_tracers_c_final grid_tracers_p_final grid_tracers_diff_final grid_δtracers_final grid_p_full_final grid_p_half_final grid_geopots_final grid_vor_final grid_div_final grid_δu_final grid_δv_final convection_final grid_w_full_final 
    end

    if mean_save_file_name != "None"
        @save mean_save_file_name spe_vor_c_xyzt spe_vor_p_xyzt spe_div_c_xyzt spe_div_p_xyzt spe_lnps_c_xyzt spe_lnps_p_xyzt spe_t_c_xyzt spe_t_p_xyzt spe_tracers_c_xyzt spe_tracers_p_xyzt grid_tracers_n_xyzt grid_u_n_xyzt grid_u_c_xyzt  grid_u_p_xyzt grid_v_n_xyzt grid_v_c_xyzt grid_v_p_xyzt grid_ps_c_xyzt grid_ps_p_xyzt grid_t_n_xyzt grid_t_c_xyzt grid_t_p_xyzt grid_tracers_c_xyzt grid_tracers_p_xyzt grid_tracers_diff_xyzt grid_δtracers_xyzt factor1_xyzt factor2_xyzt factor3_xyzt grid_p_full_xyzt grid_p_half_xyzt grid_geopots_xyzt  grid_vor_xyzt  grid_div_xyzt grid_δu_xyzt grid_δv_xyzt convection_xyzt factor4_xyzt grid_z_full_xyzt grid_w_full_xyzt grid_tracers_c_max_Tiffany_xyzt grid_tracers_c_max_xyzt
    end
end





function Sigma_Zonal_Mean_Contourf(output_manager::Output_Manager, save_file_pref::String)
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)

    t_zonal_mean, t_eq_zonal_mean, u_zonal_mean, v_zonal_mean = output_manager.t_zonal_mean, output_manager.t_eq_zonal_mean, output_manager.u_zonal_mean, output_manager.v_zonal_mean
   
    PyPlot.contourf(X, Y, t_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_T.png")
    PyPlot.close("all")

    PyPlot.contourf(X, Y, t_eq_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_Teq.png")
    PyPlot.close("all")

    PyPlot.contourf(X, Y, u_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_U.png")
    PyPlot.close("all")

    PyPlot.contourf(X, Y, v_zonal_mean, levels = 10)
    PyPlot.gca().invert_yaxis()
    PyPlot.xlabel("Latitude")
    PyPlot.ylabel("σ")
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_V.png")
    PyPlot.close("all")
    
    
end


function Sigma_Zonal_Mean_Pcolormesh(output_manager::Output_Manager, save_file_pref::String)
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)'

    t_zonal_mean, u_zonal_mean, v_zonal_mean = output_manager.t_zonal_mean, output_manager.u_zonal_mean, output_manager.v_zonal_mean
   
    PyPlot.pcolormesh(X, Y, t_zonal_mean, shading= "gouraud", cmap="viridis")
    PyPlot.gca().invert_yaxis()
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_T.png")
    PyPlot.close("all")

    PyPlot.pcolormesh(X, Y, u_zonal_mean, shading= "gouraud", cmap="viridis")
    PyPlot.gca().invert_yaxis()
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_U.png")
    PyPlot.close("all")

    PyPlot.pcolormesh(X, Y, v_zonal_mean, shading= "gouraud", cmap="viridis")
    PyPlot.gca().invert_yaxis()
    PyPlot.colorbar()
    PyPlot.savefig(save_file_pref * "_V.png")
    PyPlot.close("all")
    
    
end


function Lat_Lon_Pcolormesh(output_manager::Output_Manager, grid_dat::Array{Float64,3}, level::Int64, save_file_name::String = "None")
    
    λc, θc = output_manager.λc, output_manager.θc
    nλ, nθ = length(λc), length(θc)
    λc_deg, θc_deg = λc*180/pi, θc*180/pi
    
    
    X,Y = repeat(λc_deg, 1, nθ), repeat(θc_deg, 1, nλ)'
    
    
    PyPlot.pcolormesh(X, Y, grid_dat[:,:,level], shading= "gouraud", cmap="viridis")
    PyPlot.axis("equal")
    PyPlot.colorbar()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end


function Lat_Lon_Pcolormesh(mesh::Spectral_Spherical_Mesh, grid_dat::Array{Float64,3}, level::Int64, save_file_name::String = "None")
    
    λc, θc = mesh.λc, mesh.θc
    nλ, nθ = length(λc), length(θc)
    λc_deg, θc_deg = λc*180/pi, θc*180/pi
    
    X,Y = repeat(λc_deg, 1, nθ), repeat(θc_deg, 1, nλ)'
    
    
    PyPlot.pcolormesh(X, Y, grid_dat[:,:,level], shading= "gouraud", cmap="viridis")
    PyPlot.axis("equal")
    PyPlot.colorbar()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end


function Zonal_Mean(grid_dat::Array{Float64,3})
    
    return dropdims(mean(grid_dat, dims=1), dims=1)
    
end


function Sigma_Zonal_Mean_Pcolormesh(output_manager::Output_Manager,
    zonal_mean_data::Array{Float64,2}, save_file_name::String = "None")
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)'
    
    
    PyPlot.pcolormesh(X, Y, zonal_mean_data, shading= "gouraud", cmap="viridis")
    PyPlot.colorbar()
    PyPlot.gca().invert_yaxis()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end


function Sigma_Zonal_Mean_Contourf(output_manager::Output_Manager, 
    zonal_mean_data::Array{Float64,2}, save_file_name::String = "None")
    
    θc = output_manager.θc
    σc = output_manager.σc
    nθ = length(θc)
    θc_deg = θc*180/pi
    nd = output_manager.nd
    
    X,Y = repeat(θc_deg, 1, nd), repeat(σc, 1, nθ)'
    
    PyPlot.contourf(X, Y, zonal_mean_data, levels = 10)
    PyPlot.gca().invert_yaxis()
    
    if save_file_name != "None"
        PyPlot.savefig(save_file_name)
        PyPlot.close("all")
    end
    
end