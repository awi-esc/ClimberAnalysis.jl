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

using GibbsSeaWater

include("../src/ClimberPlot.jl")
import .ClimberPlot as plotc

dir_path="/albedo/work/projects/p_forclima/NaHosMIP/climber_x/";
climber_hos = NCDataset(joinpath(dir_path, "uh03_norestore_run02/ocn.nc"), "r");
mom5_data_hos= NCDataset("/albedo/work/projects/p_forclima/NaHosMIP/cm2mc_PISM/concat_results_CM2Mc_PISM_coupling_after_spinup_y1860_norestore_uh03_extended/regridded/ocean-yearly_sub2_regridded.nc", "r"); 

mom5_hos_temp = plotc.VarAsDimArray(mom5_data_hos, "temp", "4D", model="mom5");
mom5_hos_salt = plotc.VarAsDimArray(mom5_data_hos, "salt", "4D", model="mom5");
mom_temp= plotc.DownscaleTime((mom5_hos_temp ),10);#.- mom5_con_salt_ref
mom_salt= plotc.DownscaleTime((mom5_hos_salt ),10);#.- mom5_con_salt_ref



clim_hos_temp = plotc.VarAsDimArray(climber_hos, "t", "4D");
clim_hos_salt = plotc.VarAsDimArray(climber_hos, "s", "4D");
clim_hos_temp=replace!(mom5_hos_temp, missing => NaN);
clim_hos_salt=replace!(mom5_hos_salt, missing => NaN);


long=Array(clim_hos_salt.dims[1]);
lat=Array(clim_hos_salt.dims[2])[1:18];
lev=Array(clim_hos_salt.dims[3]);
time=Array(clim_hos_salt.dims[4]);


mask = map(x -> isnan(x) ? missing : 1.0, clim_hos_temp[lat=1:18]);
frequen= DimArray(repeat(Array(mask[:,:,2:21,:]), 1, 1, 1,1), (Dim{:lon}(long), Dim{:lat}(lat), Dim{:lev}(lev[2:21]), Dim{:time}(time)));


for t=1:length(time)
    for lo=1:length(long)
        for la=1:length(lat)
            Salt=Array(clim_hos_salt[time=t,lon=lo, lat=la]);
            z=-1 .*Array(lev);
            p= gsw_p_from_z.(z, Vector([lat[la] for i=1:28]));

            abs_salt=gsw_sa_from_sp.(Salt,p,long[lo],lat[la])

            pot_temp=Array(clim_hos_temp[time=t,lon=lo, lat=la]);

            cons_temp=gsw_ct_from_pt.(abs_salt,pot_temp);
            #if any(x<0 for x in gsw_nsquared(abs_salt,cons_temp,p,Vector([lat[la] for i=1:23]))[1])
             #   println( t, lo, la)
            #end

            frequen[lo,la,:,t]=gsw_nsquared(abs_salt[1:21],cons_temp[1:21],p[1:21],Vector([lat[la] for i=1:21]))[1]
        end
    end
end

area_control = NCDataset("/albedo/work/projects/p_forclima/NaHosMIP/climber_x/control_run02/ocn.nc", "r");
climber_area = plotc.ClimberAreaAsDimArray(area_control);

frequen2=replace(frequen, NaN => missing);

frequen_Ross_prof = sqrt.(plotc.get_weighted_Ross_ts(frequen2, climber_area[:,1:18,:,:], "4D", test=false, lat_val=-60)) ./(10^-4);

plotc.plot_depth_time(frequen_Ross_prof, "ClimberX_temp_diff in density units", (0, 8e-6))


fig = Figure();
ax = Axis(fig[1,1], 
        title = "ClimberX: tesr",
        xlabel = "time (years)",
        ylabel = "temp (°C)",
        );
        
lines!(Array(frequen_Ross_prof[time=1]), label="1")

lines!(Array(frequen_Ross_prof[time=10]), label="10")
lines!(Array(frequen_Ross_prof[time=30]), label="30")
lines!(Array(frequen_Ross_prof[time=70]), label="70")
lines!(Array(frequen_Ross_prof[time=80]), label="80")
lines!(Array(frequen_Ross_prof[time=85]), label="85")
lines!(Array(frequen_Ross_prof[time=90]), label="90")
#lines!(Array(frequen_Ross_prof[time=300]), label="220")
axislegend(ax, position = :rb );
fig



Salt=Array(clim_hos_salt[time=1,lon=1, lat=5]);
z=-1 .*Array(lev);
p= gsw_p_from_z.(z, Vector([lat[5] for i=1:23]));

abs_salt=gsw_sa_from_sp.(Salt,p,long[1],lat[5])

pot_temp=Array(clim_hos_temp[time=1,lon=1, lat=5])

cons_temp=gsw_ct_from_pt.(abs_salt,pot_temp)