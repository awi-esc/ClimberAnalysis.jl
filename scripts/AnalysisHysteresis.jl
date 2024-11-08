using NCDatasets
using NetCDF
using DimensionalData
using CairoMakie
using GeoMakie
using Statistics
using ColorSchemes
using SkipNan
using NaturalEarth
using DataFrames

#include("../src/ClimberPlot.jl")
import .ClimberPlot as Cp


# plot data of first Hysteresis test runs with climber x

#Load time series data
dir_path="/albedo/work/user/anhoes003/test_experiments/hysteresis_0.04/";

data_up = NCDataset(joinpath(dir_path, "up/ocn_ts.nc"), "r");

data = NCDataset(joinpath(dir_path, "down/ocn_ts.nc"), "r");

#plot timeseries
test_list=[data[:time], data[:omaxa], data[:sst], data[:t_deep_so]]
test_labels = ["AMOC (Sv)", "SST (°C)", "deep SO temp (°C)"]
timeseries(test_list, test_labels)

#plot hysteresis curve
hysteresis_plot(data_up, data, "15000 year hysteresis run with hosing_ini=-0.3 hosing_trend=0.04")


#load data of restart test for down branch
dir_path="/albedo/work/user/anhoes003/test_experiments/hysteresis_0.04_test_rest/down/";

data_r1 = NCDataset(joinpath(dir_path, "ocn_ts.nc"), "r");

dir_path="/albedo/work/user/anhoes003/test_experiments/hysteresis_0.04_test_rest/down2/";

data_r2 = NCDataset(joinpath(dir_path, "ocn_ts.nc"), "r");

#test if restart shows same hysteresis plot as full 
fig = Figure();
ax = Axis(fig[1,1], 
    title = "restarted",
    xlabel = "Freshwater hosing (Sv)",
    ylabel = "max AMOC (Sv)",
    );
    scatter!(data[:hosing], data[:omaxa])
    scatter!(data_r1[:hosing], data_r1[:omaxa])
    scatter!(data_r2[:hosing], data_r2[:omaxa])
fig


#plot timeseries of down branch and restarted version to be sure it is in line
f = Figure()
ax1 = Axis(f[1, 1], ylabel= test_labels[1])
ax2 = Axis(f[2, 1], ylabel= test_labels[2])
ax3 = Axis(f[3, 1], xtickformat = "{:.0f}", xlabel="time (years)", ylabel= test_labels[3])
# hideydecorations!(ax1, grid = false)
linkxaxes!(ax1, ax2, ax3)
rowgap!(f.layout, 0.2)
hidexdecorations!(ax1, ticks=false, grid=false)
hidexdecorations!(ax2, ticks=false, grid=false)

lines!(ax1, data[:time], data[:omaxa], label="full")
lines!(ax2, data[:time], data[:sst], label="full")
lines!(ax3, data[:time], data[:t_deep_so], label="full")

lines!(ax1, data_r1[:time], data_r1[:omaxa], label="before_rest")
lines!(ax2, data_r1[:time], data_r1[:sst], label="before_rest")
lines!(ax3, data_r1[:time], data_r1[:t_deep_so], label="before_rest")

lines!(ax1, data_r2[:time].+5000, data_r2[:omaxa], label="after_rest")
lines!(ax2, data_r2[:time].+5000, data_r2[:sst], label="after_rest")
lines!(ax3, data_r2[:time].+5000, data_r2[:t_deep_so], label="after_rest")
axislegend(ax3, position = :lt )
xlims!(ax3,(7000,9000))

f


# load spatial ocean data 

dir_path="/albedo/work/user/anhoes003/test_experiments/hysteresis_0.04/up/";

data = NCDataset(joinpath(dir_path, "ocn.nc"), "r");

data

SST = getAveragedVar(data,"sst");
SST_anom = getAnomVar(data, "sst", 15, ref_year=5);

labels=["SST", "SST anomaly", "MLD", "MLD anomaly"]

MLD= getAveragedVar(data, "mld");
MLD_anom= getAnomVar(data, "mld", 15, ref_year=5);
list=[SST, SST_anom, MLD, MLD_anom]

#test to plot variables in panel ###TODO
map_panel(list, labels)

# plot global map
global_map(SST, "SST", colors=:heat)
global_map(MLD, "mld")

#prepare data lists for amoc panel plot
amoc_list = [getDataDepthProfile(data,"opsi_a", year_ind = 8), getDataDepthProfile(data,"opsi_a", year_ind =9), getDataDepthProfile(data,"opsi_a", year_ind =10), getDataDepthProfile(data,"opsi_a", year_ind =11)]
mld_list = [getDataYearlyMax(data, "mld", year_ind=8), getDataYearlyMax(data, "mld", year_ind=9),getDataYearlyMax(data, "mld", year_ind=10), getDataYearlyMax(data, "mld", year_ind=11)]
labels2= ["max MLD yr8000","max MLD yr9000","max MLD yr10000","max MLD yr11000","Atlantic overturning yr8000", "Atlantic overturning yr9000", "Atlantic overturning yr10000", "Atlantic overturning yr11000"]

#load sea ice data
data_si = NCDataset(joinpath(dir_path, "sic.nc"), "r");
SI_frac=getAveragedVar(data_si,"fsic")


#plot amoc panel
amoc_map_panel(mld_list, amoc_list, SI_frac,labels2)



