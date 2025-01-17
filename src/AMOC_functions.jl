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



#load data


dir_path="/albedo/work/projects/p_forclima/NaHosMIP/cm2mc_PISM/concat_results_CM2Mc_PISM_coupling_after_spinup_y1860_backup_control_extended/regridded/";

Salt = NCDataset(joinpath(dir_path, "oceanSalt-decadal_regridded.nc"), "r");

Velocity = NCDataset(joinpath(dir_path, "oceanVel-decadal_regridded.nc"), "r");

dir_path="/albedo/work/projects/p_forclima/NaHosMIP/climber_x/control_run02";

Salt = NCDataset(joinpath(dir_path, "ocn.nc"), "r");

Velocity = NCDataset(joinpath(dir_path, "ocn.nc"), "r");

#V_atl= Array(Velocity["v"])[(Array(Velocity["lon"]).>(-85)) .& (Array(Velocity["lon"]).<(25)), :,:,:];
#V_atl_new=deepcopy(V_atl[1,:,:,:]);
#Salt_atl=Array(Salt["salt"])[(Array(Salt["lon"]).>(-85)) .& (Array(Salt["lon"]).<(25)), :,:,:];
#Salt_atl_new=deepcopy(Salt_atl[1,:,:,:]);

#for i =1:length(V_atl[1,:,1,1])
 #   for j =1:length(V_atl[1,1,:,1])
  #      V_atl_new[i,j,:].= sum.(skipmissing.(eachcol(V_atl[:,i,j,:])))
   #     Salt_atl_new[i,j,:].= mean.(skipmissing.(eachcol(Salt_atl[:,i,j,:]))) 
    #end
#end

#arg=V_atl_new.*Salt_atl_new;
#Mov=deepcopy(arg[:,1,:]);
#for i =1:length(arg[:,1,1])
 #   Mov[i,:].= sum.(skipnan.(eachcol(arg[i,:,:])))
#end

S0=35;

#plot(((-Mov[(Array(Salt["lat"]).==(-32.5)),:].+Mov[(Array(Salt["lat"]).==(82.5)),:])/S0)[1,:])


lat = Velocity["lat"];
lon = Velocity["lon"];
st_ocean = Velocity["lev"];
#st_ocean = Velocity["st_ocean"];
month= 1:12;
time = 1:10:3000;

#monthly
V_integ=DimArray(Array(Velocity["v"])[:,:,:,1:12,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:st_ocean}(st_ocean),Dim{:month}(month),Dim{:time}(time)));
V_barotrop=DimArray(Array(Velocity["vb"])[:,:,13,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:time}(time)));
S_integ=DimArray(Array(Salt["s"][:,:,:,1:12,:]), (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:st_ocean}(st_ocean),Dim{:month}(month),Dim{:time}(time)));

#yearly
#V_integ=DimArray(Array(Velocity["v"])[:,:,:,13,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:st_ocean}(st_ocean), Dim{:time}(time)));
#V_barotrop=DimArray(Array(Velocity["vb"])[:,:,13,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:time}(time)));
#S_integ=DimArray(Array(Salt["s"][:,:,:,13,:]), (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:st_ocean}(st_ocean), Dim{:time}(time)));

V_barotrop=reshape(V_barotrop,72,36,1,1,300);

V_integ=(V_integ.-V_barotrop).*reshape(Velocity["dxv"][2:end],1,36,1,1,1);

V_dim_atl=V_integ[lon=At(-84:5:25;atol=2.5)];
Salt_dim_atl=S_integ[lon=At(-84:5:25;atol=2.5)];

shifted_salt = Array{Union{Float64, Missing}}(missing, size(Salt_dim_atl));
shifted_salt[:,2:end, :,:, :] .= Array(Salt_dim_atl)[:,1:end-1, :,:, :] ;
Salt_dim_atl=(Salt_dim_atl.+shifted_salt)./2;

V_integ_atl=mapslices(x -> Statistics.sum(skipmissing(x)), V_dim_atl, dims=(:lon));
S_mean_atl=mapslices(x -> Statistics.mean(skipmissing(x)), Salt_dim_atl, dims=(:lon));

Result= mapslices(x -> Statistics.sum(skipnan(x)), V_integ_atl.*S_mean_atl.*reshape(Velocity["dz"],1,1,23,1,1), dims=(:st_ocean))
        
fig = Figure();
ax = Axis(fig[1,1], 
    #title = "Comparison Amoc strength Climberx vs CM2Mc",
    xlabel = "time (years)",
    ylabel = "Fov (Sv)",
    )
scatter!(ax,Array((-Result[lat=At(-32.5)].+Result[lat=At(82.5)])/S0*10^-6)[1,1,9,:])
display(fig)

#end

unique(Array(sum(V_integ,dims=1)));
mask = isfinite.(test) 

test=replace!(Array(Velocity["v"]), missing => NaN);
test=skipmissing(Velocity["v"]);
#sum.(V_integ[mas,dims=1)


mean.(skipmissing.(eachrow(Array(Velocity[:v]))))

valid_mask = .!ismissing(Array(Velocity[:v]))

#function getMov_perLat(data1,data2)
    intg_vel= 
#end

Array(Velocity["v"])
mask
Array(Velocity["v"])[mask]