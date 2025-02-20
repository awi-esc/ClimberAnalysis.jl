using NCDatasets
using NetCDF
using DimensionalData
using CairoMakie
using GeoMakie
using Statistics
using NaNStatistics
using ColorSchemes
using SkipNan
using CommonDataModel
using DataFrames
using Skipper

include("../src/ClimberPlot.jl")
import .ClimberPlot as plotc


#load climberX data

dir_path="/albedo/work/projects/p_forclima/NaHosMIP/climber_x/";

climber_hos = NCDataset(joinpath(dir_path, "uh03_norestore_run02/ocn.nc"), "r");

#load MOM5 data

dir_path="/albedo/work/projects/p_forclima/NaHosMIP/cm2mc_PISM/";

mom5_data_control= NCDataset(joinpath(dir_path,"concat_results_CM2Mc_PISM_coupling_after_spinup_y1860_backup_control_extended/regridded/ocean-yearly_sub_regridded.nc"), "r"); 
mom5_data_hos= NCDataset(joinpath(dir_path,"concat_results_CM2Mc_PISM_coupling_after_spinup_y1860_norestore_uh03_extended/regridded/ocean-yearly_sub2_regridded.nc"), "r"); 
mom5_data_control_dec= NCDataset(joinpath(dir_path,"concat_results_CM2Mc_PISM_coupling_after_spinup_y1860_backup_control_extended/regridded/oceanRho-decadal_regridded.nc"), "r"); 
mom5_data_hos_dec= NCDataset(joinpath(dir_path,"concat_results_CM2Mc_PISM_coupling_after_spinup_y1860_norestore_uh03_extended/regridded/oceanRho-decadal_regridded.nc"), "r"); 
mom5_data_control_vel= NCDataset(joinpath(dir_path,"concat_results_CM2Mc_PISM_coupling_after_spinup_y1860_backup_control_extended/regridded/oceanVel-decadal_regridded.nc"), "r"); 
mom5_data_hos_vel= NCDataset(joinpath(dir_path,"concat_results_CM2Mc_PISM_coupling_after_spinup_y1860_norestore_uh03_extended/regridded/oceanVel-decadal_regridded.nc"), "r"); 

mom5_data_hos_mld= NCDataset(joinpath(dir_path,"concat_results_CM2Mc_PISM_coupling_after_spinup_y1860_norestore_uh03_extended/regridded/ocean-yearly_max_regridded.nc"), "r"); 


output_path="/albedo/work/user/anhoes003/data_analysis/output_21012025/";



area_control = NCDataset("/albedo/work/projects/p_forclima/NaHosMIP/climber_x/control_run02/ocn.nc", "r");
climber_area = plotc.ClimberAreaAsDimArray(area_control);

clim_hos_temp = plotc.VarAsDimArray(climber_hos, "t", "4D", month=9);
temp_diff= clim_hos_temp .- clim_hos_temp[time=1];

clim_hos_salt = plotc.VarAsDimArray(climber_hos, "s", "4D", month=9);
salt_diff= clim_hos_salt .- clim_hos_salt[time=1];

clim_hos_rho = plotc.VarAsDimArray(climber_hos, "rho", "4D", month=9);
rho_diff= clim_hos_rho .- clim_hos_rho[time=1];

clim_hos_mld = plotc.VarAsDimArray(climber_hos, "mld", "3D", month=9);

clim_hos_u = plotc.VarAsDimArray(climber_hos, "u", "4D", month=9);
u_diff= clim_hos_u .- clim_hos_u[time=1];
clim_hos_v = plotc.VarAsDimArray(climber_hos, "v", "4D", month=9);
v_diff= clim_hos_v .- clim_hos_v[time=1];
clim_hos_w = plotc.VarAsDimArray(climber_hos, "w", "4D", month=9);


mom5_con_temp = plotc.VarAsDimArray(mom5_data_control, "temp", "4D", model="mom5");
mom5_con_temp_ref = dropdims(mapslices(x -> Statistics.mean(skipmissing(x)), mom5_con_temp[:,:,:,1:100], dims=(:time)), dims=(:time)) ;
mom5_hos_temp = plotc.VarAsDimArray(mom5_data_hos, "temp", "4D", model="mom5");

mom_temp_diff= mom5_hos_temp .- mom5_con_temp_ref;

mom5_con_salt = plotc.VarAsDimArray(mom5_data_control, "salt", "4D", model="mom5");

mom5_con_salt_ref = dropdims(mapslices(x -> Statistics.mean(skipmissing(x)), mom5_con_salt[:,:,:,1:100], dims=(:time)), dims=(:time)) ;
mom5_hos_salt = plotc.VarAsDimArray(mom5_data_hos, "salt", "4D", model="mom5");

mom_salt_diff= mom5_hos_salt .- mom5_con_salt_ref;

mom_con_rho = plotc.VarAsDimArray(mom5_data_control_dec, "pot_rho_2k", "4D", model="mom5", decadal=true);
mom_con_rho_ref = dropdims(mapslices(x -> Statistics.mean(skipmissing(x)), mom_con_rho[:,:,:,1:100], dims=(:time)), dims=(:time)) ;
mom_hos_rho = plotc.VarAsDimArray(mom5_data_hos_dec, "pot_rho_2k", "4D", model="mom5", decadal=true);

mom_rho_diff= mom_hos_rho .- mom_con_rho_ref;

mom5_con_u = plotc.VarAsDimArray(mom5_data_control_vel, "u", "4D", model="mom5", decadal=true);
mom5_con_u_ref = dropdims(mapslices(x -> Statistics.mean(skipmissing(x)), mom5_con_u[:,:,:,1:100], dims=(:time)), dims=(:time)) ;
mom_hos_u= plotc.VarAsDimArray(mom5_data_hos_vel, "u", "4D", model="mom5", decadal=true);

mom_u_diff= mom5_con_u .- mom5_con_u_ref;

mom5_con_v = plotc.VarAsDimArray(mom5_data_control_vel, "v", "4D", model="mom5", decadal=true);
mom5_con_v_ref = dropdims(mapslices(x -> Statistics.mean(skipmissing(x)), mom5_con_v[:,:,:,1:100], dims=(:time)), dims=(:time)) ;
mom_hos_v= plotc.VarAsDimArray(mom5_data_hos_vel, "v", "4D", model="mom5");

mom_v_diff= mom5_con_v .- mom5_con_v_ref;

mom_hos_mld= plotc.VarAsDimArray(mom5_data_hos_mld, "mld_max", "3D", model="mom5");

display(lines!(fig, mom_hos_mld[1,5,:]))
#plot climber X variables for cell 70,5 and surrounding

plotc.plot_depth_time(clim_hos_temp[69,4,:,:], "ClimberX_temp_diff; lon= 162.5, lat=-72.5", (-3,3), mld=true, mld_data=clim_hos_mld[69,4,:], save_out=true, output_path=output_path, outfile="C_T_anom_Hov_1.png")
plotc.plot_depth_time(clim_hos_temp[70,4,:,:], "ClimberX_temp_diff; lon= 167.5, lat=-72.5", (-3,3), mld=true, mld_data=clim_hos_mld[70,4,:], save_out=true, output_path=output_path, outfile="C_T_anom_Hov_2.png")
plotc.plot_depth_time(clim_hos_temp[71,4,:,:], "ClimberX_temp_diff; lon= 172.5, lat=-72.5", (-3,3), mld=true, mld_data=clim_hos_mld[71,4,:], save_out=true, output_path=output_path, outfile="C_T_anom_Hov_3.png")

plotc.plot_depth_time(clim_hos_temp[69,5,:,:], "ClimberX_temp_diff; lon= 162.5, lat=-67.5", (-3,3), mld=true, mld_data=clim_hos_mld[69,5,:], save_out=true, output_path=output_path, outfile="C_T_anom_Hov_4.png")
plotc.plot_depth_time(clim_hos_temp[70,5,:,:], "ClimberX_temp_diff; lon= 167.5, lat=-67.5", (-3,3), mld=true, mld_data=clim_hos_mld[70,5,:], save_out=true, output_path=output_path, outfile="C_T_anom_Hov_5.png")
plotc.plot_depth_time(clim_hos_temp[71,5,:,:], "ClimberX_temp_diff; lon= 172.5, lat=-67.5", (-3,3), mld=true, mld_data=clim_hos_mld[71,5,:], save_out=true, output_path=output_path, outfile="C_T_anom_Hov_6.png")

plotc.plot_depth_time(clim_hos_temp[69,6,:,:], "ClimberX_temp_diff; lon= 162.5, lat=-62.5", (-3,3), mld=true, mld_data=clim_hos_mld[69,6,:], save_out=true, output_path=output_path, outfile="C_T_anom_Hov_7.png")
plotc.plot_depth_time(clim_hos_temp[70,6,:,:], "ClimberX_temp_diff; lon= 167.5, lat=-62.5", (-3,3), mld=true, mld_data=clim_hos_mld[70,6,:], save_out=true, output_path=output_path, outfile="C_T_anom_Hov_8.png")
plotc.plot_depth_time(clim_hos_temp[71,6,:,:], "ClimberX_temp_diff; lon= 172.5, lat=-62.5", (-3,3), mld=true, mld_data=clim_hos_mld[71,6,:], save_out=true, output_path=output_path, outfile="C_T_anom_Hov_9.png")



plotc.plot_depth_time(0.8 .*salt_diff[69,4,:,:], "ClimberX_salt_diff; lon= 162.5, lat=-72.5", (-0.8,0.8), mld=true, mld_data=clim_hos_mld[69,4,:], save_out=true, output_path=output_path, outfile="C_S_ru_anom_Hov_1.png")
plotc.plot_depth_time(0.8 .*salt_diff[70,4,:,:], "ClimberX_salt_diff; lon= 167.5, lat=-72.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[70,4,:], save_out=true, output_path=output_path, outfile="C_S_ru_anom_Hov_2.png")
plotc.plot_depth_time(0.8 .*salt_diff[71,4,:,:], "ClimberX_salt_diff; lon= 172.5, lat=-72.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[71,4,:], save_out=true, output_path=output_path, outfile="C_S_ru_anom_Hov_3.png")

plotc.plot_depth_time(0.8 .*salt_diff[69,5,:,:], "ClimberX_salt_diff; lon= 162.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[69,5,:], save_out=true, output_path=output_path, outfile="C_S_ru_anom_Hov_4.png")
plotc.plot_depth_time(0.8 .*salt_diff[70,5,:,:], "ClimberX_salt_diff; lon= 167.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[70,5,:], save_out=true, output_path=output_path, outfile="C_S_ru_anom_Hov_5.png")
plotc.plot_depth_time(0.8 .*salt_diff[71,5,:,:], "ClimberX_salt_diff; lon= 172.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[71,5,:], save_out=true, output_path=output_path, outfile="C_S_ru_anom_Hov_6.png")

plotc.plot_depth_time(0.8 .*salt_diff[69,6,:,:], "ClimberX_salt_diff; lon= 162.5, lat=-62.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[69,6,:], save_out=true, output_path=output_path, outfile="C_S_ru_anom_Hov_7.png")
plotc.plot_depth_time(0.8 .*salt_diff[70,6,:,:], "ClimberX_salt_diff; lon= 167.5, lat=-62.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[70,6,:], save_out=true, output_path=output_path, outfile="C_S_ru_anom_Hov_8.png")
plotc.plot_depth_time(0.8 .*salt_diff[71,6,:,:], "ClimberX_salt_diff; lon= 172.5, lat=-62.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[71,6,:], save_out=true, output_path=output_path, outfile="C_S_ru_anom_Hov_9.png")



plotc.plot_depth_time(rho_diff[69,4,:,:], "ClimberX_rho_diff; lon= 162.5, lat=-72.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[69,4,:], save_out=true, output_path=output_path, outfile="C_r_anom_Hov_1.png")
plotc.plot_depth_time(rho_diff[70,4,:,:], "ClimberX_rho_diff; lon= 167.5, lat=-72.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[70,4,:], save_out=true, output_path=output_path, outfile="C_r_anom_Hov_2.png")
plotc.plot_depth_time(rho_diff[71,4,:,:], "ClimberX_rho_diff; lon= 172.5, lat=-72.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[71,4,:], save_out=true, output_path=output_path, outfile="C_r_anom_Hov_3.png")

plotc.plot_depth_time(rho_diff[69,5,:,:], "ClimberX_rho_diff; lon= 162.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[69,5,:], save_out=true, output_path=output_path, outfile="C_r_anom_Hov_4.png")
plotc.plot_depth_time(rho_diff[70,5,:,:], "ClimberX_rho_diff; lon= 167.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[70,5,:], save_out=true, output_path=output_path, outfile="C_r_anom_Hov_5.png")
plotc.plot_depth_time(rho_diff[71,5,:,:], "ClimberX_rho_diff; lon= 172.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[71,5,:], save_out=true, output_path=output_path, outfile="C_r_anom_Hov_6.png")

plotc.plot_depth_time(rho_diff[69,6,:,:], "ClimberX_rho_diff; lon= 162.5, lat=-62.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[69,6,:], save_out=true, output_path=output_path, outfile="C_r_anom_Hov_7.png")
plotc.plot_depth_time(rho_diff[70,6,:,:], "ClimberX_rho_diff; lon= 167.5, lat=-62.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[70,6,:], save_out=true, output_path=output_path, outfile="C_r_anom_Hov_8.png")
plotc.plot_depth_time(rho_diff[71,6,:,:], "ClimberX_rho_diff; lon= 172.5, lat=-62.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[71,6,:], save_out=true, output_path=output_path, outfile="C_r_anom_Hov_9.png")

alpha=0.052 .+0.012 .*clim_hos_temp[time=1];
beta=0.8;

#plotc.plot_depth_time(-(0.052 .+0.012 .*clim_hos_temp[69,4,:,1]) .*clim_hos_temp[69,4,:,:], "ClimberX_temp_diff; lon= 162.5, lat=-72.5", (-3,3), mld=true, mld_data=clim_hos_mld[69,4,:], save_out=true, output_path=output_path, outfile="C_T_ru_anom_Hov_1.png")
plotc.plot_depth_time(-(0.052 .+0.012 .*clim_hos_temp[70,4,:,1]) .*clim_hos_temp[70,4,:,:], "ClimberX_temp_diff; lon= 167.5, lat=-72.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[70,4,:], save_out=true, output_path=output_path, outfile="C_T_ru_anom_Hov_2.png")
plotc.plot_depth_time(-(0.052 .+0.012 .*clim_hos_temp[71,4,:,1]) .*clim_hos_temp[71,4,:,:], "ClimberX_temp_diff; lon= 172.5, lat=-72.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[71,4,:], save_out=true, output_path=output_path, outfile="C_T_ru_anom_Hov_3.png")

plotc.plot_depth_time(-(0.052 .+0.012 .*clim_hos_temp[69,5,:,1]) .*clim_hos_temp[69,5,:,:], "ClimberX_temp_diff; lon= 162.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[69,5,:], save_out=true, output_path=output_path, outfile="C_T_ru_anom_Hov_4.png")
plotc.plot_depth_time(-(0.052 .+0.012 .*clim_hos_temp[70,5,:,1]) .*clim_hos_temp[70,5,:,:], "ClimberX_temp_diff; lon= 167.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[70,5,:], save_out=true, output_path=output_path, outfile="C_T_ru_anom_Hov_5.png")
plotc.plot_depth_time(-(0.052 .+0.012 .*clim_hos_temp[71,5,:,1]) .*clim_hos_temp[71,5,:,:], "ClimberX_temp_diff; lon= 172.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[71,5,:], save_out=true, output_path=output_path, outfile="C_T_ru_anom_Hov_6.png")

plotc.plot_depth_time(-(0.052 .+0.012 .*clim_hos_temp[69,6,:,1]) .*clim_hos_temp[69,6,:,:], "ClimberX_temp_diff; lon= 162.5, lat=-62.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[69,6,:], save_out=true, output_path=output_path, outfile="C_T_ru_anom_Hov_7.png")
plotc.plot_depth_time(-(0.052 .+0.012 .*clim_hos_temp[70,6,:,1]) .*clim_hos_temp[70,6,:,:], "ClimberX_temp_diff; lon= 167.5, lat=-62.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[70,6,:], save_out=true, output_path=output_path, outfile="C_T_ru_anom_Hov_8.png")
plotc.plot_depth_time(-(0.052 .+0.012 .*clim_hos_temp[71,6,:,1]) .*clim_hos_temp[71,6,:,:], "ClimberX_temp_diff; lon= 172.5, lat=-62.5", (-0.6,0.6), mld=true, mld_data=clim_hos_mld[71,6,:], save_out=true, output_path=output_path, outfile="C_T_ru_anom_Hov_9.png")



##Plot MOM data

plotc.plot_depth_time(mom_temp_diff[72,5,:,:], "MOM_temp_diff; lon= 177.5, lat=-67.5", (-3,3), mld=true, mld_data=mom_hos_mld[72,5,:], save_out=true, output_path=output_path, outfile="M_T_anom_Hov_1.png")
plotc.plot_depth_time(mom_temp_diff[1,5,:,:], "MOM_temp_diff; lon= -177.5, lat=-67.5", (-3,3), mld=true, mld_data=mom_hos_mld[1,5,:], save_out=true, output_path=output_path, outfile="M_T_anom_Hov_2.png")
plotc.plot_depth_time(mom_temp_diff[2,5,:,:], "MOM_temp_diff; lon= -172.5, lat=-67.5", (-3,3), mld=true, mld_data=mom_hos_mld[2,5,:], save_out=true, output_path=output_path, outfile="M_T_anom_Hov_3.png")

plotc.plot_depth_time(mom_temp_diff[72,4,:,:], "MOM_temp_diff; lon= 177.5, lat=-72.5", (-3,3), mld=true, mld_data=mom_hos_mld[72,4,:], save_out=true, output_path=output_path, outfile="M_T_anom_Hov_4.png")
plotc.plot_depth_time(mom_temp_diff[1,4,:,:], "MOM_temp_diff; lon= -177.5, lat=-72.5", (-3,3), mld=true, mld_data=mom_hos_mld[1,4,:], save_out=true, output_path=output_path, outfile="M_T_anom_Hov_5.png")
plotc.plot_depth_time(mom_temp_diff[2,4,:,:], "MOM_temp_diff; lon= -172.5, lat=-72.5", (-3,3), mld=true, mld_data=mom_hos_mld[2,4,:], save_out=true, output_path=output_path, outfile="M_T_anom_Hov_6.png")



plotc.plot_depth_time(mom_salt_diff[72,5,:,:], "MOM_salt_diff; lon= 177.5, lat=-67.5", (-0.8,0.8), mld=true, mld_data=mom_hos_mld[72,5,:], save_out=true, output_path=output_path, outfile="M_S_anom_Hov_1.png")
plotc.plot_depth_time(mom_salt_diff[1,5,:,:], "MOM_salt_diff; lon= -177.5, lat=-67.5", (-0.8,0.8), mld=true, mld_data=mom_hos_mld[1,5,:], save_out=true, output_path=output_path, outfile="M_S_anom_Hov_2.png")
plotc.plot_depth_time(mom_salt_diff[2,5,:,:], "MOM_salt_diff; lon= -172.5, lat=-67.5", (-0.8,0.8), mld=true, mld_data=mom_hos_mld[2,5,:], save_out=true, output_path=output_path, outfile="M_S_anom_Hov_3.png")

plotc.plot_depth_time(mom_salt_diff[72,4,:,:], "MOM_salt_diff; lon= 177.5, lat=-72.5", (-0.8,0.8), mld=true, mld_data=mom_hos_mld[72,4,:], save_out=true, output_path=output_path, outfile="M_S_anom_Hov_4.png")
plotc.plot_depth_time(mom_salt_diff[1,4,:,:], "MOM_salt_diff; lon= -177.5, lat=-72.5", (-0.8,0.8), mld=true, mld_data=mom_hos_mld[1,4,:], save_out=true, output_path=output_path, outfile="M_S_anom_Hov_5.png")
plotc.plot_depth_time(mom_salt_diff[2,4,:,:], "MOM_salt_diff; lon= -172.5, lat=-72.5", (-0.8,0.8), mld=true, mld_data=mom_hos_mld[2,4,:], save_out=true, output_path=output_path, outfile="M_S_anom_Hov_6.png")


plotc.plot_depth_time(mom_rho_diff[72,5,:,:], "MOM_rho_diff; lon= 177.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=mom_hos_mld[72,5,:], save_out=true, output_path=output_path, outfile="M_r_anom_Hov_1.png")
plotc.plot_depth_time(mom_rho_diff[1,5,:,:], "MOM_rho_diff; lon= -177.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=mom_hos_mld[1,5,:], save_out=true, output_path=output_path, outfile="M_r_anom_Hov_2.png")
plotc.plot_depth_time(mom_rho_diff[2,5,:,:], "MOM_rho_diff; lon= -172.5, lat=-67.5", (-0.6,0.6), mld=true, mld_data=mom_hos_mld[2,5,:], save_out=true, output_path=output_path, outfile="M_r_anom_Hov_3.png")

plotc.plot_depth_time(mom_rho_diff[72,4,:,:], "MOM_rho_diff; lon= 177.5, lat=-72.5", (-0.6,0.6), mld=true, mld_data=mom_hos_mld[72,4,:], output_path=output_path, save_out=true, outfile="M_r_anom_Hov_4.png")
plotc.plot_depth_time(mom_rho_diff[1,4,:,:], "MOM_rho_diff; lon= -177.5, lat=-72.5", (-0.6,0.6), mld=true, mld_data=mom_hos_mld[1,4,:], output_path=output_path, save_out=true, outfile="M_r_anom_Hov_5.png")
plotc.plot_depth_time(mom_rho_diff[2,4,:,:], "MOM_rho_diff; lon= -172.5, lat=-72.5", (-0.6,0.6), mld=true, mld_data=mom_hos_mld[2,4,:], output_path=output_path, save_out=true, outfile="M_r_anom_Hov_6.png")







plotc.plot_depth_time(clim_hos_u[69,4,:,:], "ClimberX_u; lon= 162.5, lat=-72.5", (-1,1))
plotc.plot_depth_time(clim_hos_u[70,4,:,:], "ClimberX_u; lon= 167.5, lat=-72.5", (-0.03,0.03))
plotc.plot_depth_time(clim_hos_u[71,4,:,:], "ClimberX_u; lon= 172.5, lat=-72.5", (-0.03,0.03))

plotc.plot_depth_time(clim_hos_u[69,5,:,:], "ClimberX_u; lon= 162.5, lat=-67.5", (-0.03,0.03))
plotc.plot_depth_time(u_diff[70,5,:,:], "ClimberX_u; lon= 167.5, lat=-67.5", (-0.03,0.03))
plotc.plot_depth_time(clim_hos_u[71,5,:,:], "ClimberX_u; lon= 172.5, lat=-67.5", (-0.03,0.03))

plotc.plot_depth_time(mom_u_diff[1,4,:,:], "ClimberX_u; lon= 167.5, lat=-67.5", (-0.01,0.01), mld=true, mld_data=mom_hos_mld[1,4,:])
plotc.plot_depth_time(mom_v_diff[1,4,:,:], "ClimberX_u; lon= 167.5, lat=-67.5", (-0.01,0.01), mld=true, mld_data=mom_hos_mld[1,4,:])


plotc.plot_depth_time(clim_hos_u[69,6,:,:], "ClimberX_u; lon= 162.5, lat=-62.5", (-0.03,0.03))
plotc.plot_depth_time(clim_hos_u[70,6,:,:], "ClimberX_u; lon= 167.5, lat=-62.5", (-0.03,0.03))
plotc.plot_depth_time(clim_hos_u[71,6,:,:], "ClimberX_u; lon= 172.5, lat=-62.5", (-0.03,0.03))


plotc.plot_depth_time(clim_hos_v[69,4,:,:], "ClimberX_v; lon= 162.5, lat=-72.5", (-1,1))
plotc.plot_depth_time(clim_hos_v[70,4,:,:], "ClimberX_v; lon= 167.5, lat=-72.5", (-0.1,0.1))
plotc.plot_depth_time(clim_hos_v[71,4,:,:], "ClimberX_v; lon= 172.5, lat=-72.5", (-0.1,0.1))

plotc.plot_depth_time(clim_hos_v[69,5,:,:], "ClimberX_v; lon= 162.5, lat=-67.5", (-0.03,0.03))
plotc.plot_depth_time(clim_hos_v[70,5,:,:], "ClimberX_v; lon= 167.5, lat=-67.5", (-0.03,0.03))
plotc.plot_depth_time(clim_hos_v[71,5,:,:], "ClimberX_v; lon= 172.5, lat=-67.5", (-0.03,0.03))

plotc.plot_depth_time(clim_hos_v[69,6,:,:], "ClimberX_v; lon= 162.5, lat=-62.5", (-0.03,0.03))
plotc.plot_depth_time(clim_hos_v[70,6,:,:], "ClimberX_v; lon= 167.5, lat=-62.5", (-0.03,0.03))
plotc.plot_depth_time(clim_hos_v[71,6,:,:], "ClimberX_v; lon= 172.5, lat=-62.5", (-0.03,0.03))


plotc.plot_depth_time(clim_hos_w[69,4,:,:], "ClimberX_w; lon= 162.5, lat=-72.5", (-1,1))
plotc.plot_depth_time(clim_hos_w[70,4,:,:], "ClimberX_w; lon= 167.5, lat=-72.5", (-2e-5,2e-5))
plotc.plot_depth_time(clim_hos_w[71,4,:,:], "ClimberX_w; lon= 172.5, lat=-72.5", (-2e-5,2e-5))

plotc.plot_depth_time(clim_hos_w[69,5,:,:], "ClimberX_w; lon= 162.5, lat=-67.5", (-2e-5,2e-5))
plotc.plot_depth_time(clim_hos_w[70,5,:,:], "ClimberX_w; lon= 167.5, lat=-67.5", (-5e-5,5e-5))
plotc.plot_depth_time(clim_hos_w[71,5,:,:], "ClimberX_w; lon= 172.5, lat=-67.5", (-5e-5,5e-5))

plotc.plot_depth_time(clim_hos_w[69,6,:,:], "ClimberX_w; lon= 162.5, lat=-62.5",(-2e-5,2e-5))
plotc.plot_depth_time(clim_hos_w[70,6,:,:], "ClimberX_w; lon= 167.5, lat=-62.5", (-1e-5,1e-5))
plotc.plot_depth_time(clim_hos_w[71,6,:,:], "ClimberX_w; lon= 172.5, lat=-62.5", (-1e-5,1e-5))