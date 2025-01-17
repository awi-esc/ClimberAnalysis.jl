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
area_control = NCDataset("/albedo/work/projects/p_forclima/NaHosMIP/climber_x/control_run02/ocn.nc", "r");

climber_area = plotc.ClimberAreaAsDimArray(area_control);


dir_path="/albedo/work/projects/p_forclima/NaHosMIP/climber_x/";
climber_hos = NCDataset(joinpath(dir_path, "uh03_norestore_run02/ocn.nc"), "r");


lon = climber_hos["lon"];
lat = climber_hos["lat"];
lev= climber_hos["lev"];
levw= climber_hos["levw"];
month = 1:12;
time =climber_hos["time"];

#TEMPERATURE

t= DimArray(Array(climber_hos["t"][:,:,:,3,:]), (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:lev}(lev), Dim{:time}(time)));
s= DimArray(Array(climber_hos["s"][:,:,:,3,:]), (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:lev}(lev), Dim{:time}(time)));


a0 = [999.842594,6.793952e-2,-9.095290e-3,1.001685e-4,-1.120083e-6,6.536336e-9];
a1 = [8.24493e-1,-4.0899e-3,7.6438e-5,-8.2467e-7,5.3875e-9];

b1 = [-5.72466e-3,1.0227e-4,-1.6546e-6];
c1 = 4.8314e-4;
c2 = [19659.33, 144.4304, 52.848, 0.3101089, 3.186519, 2.212276e-2, 6.704388e-3];

# EOS80 equation of state, Millero and Poisson 1981, but only limited number of terms!
# parameters for potential instead of in-situ temperature (after Jackett and McDougall, 1995)

dimens = (72, 36, 300, 23);

mask = map(x -> ismissing(x) ? missing : 1.0, t);
#mask=permutedims(mask, (1,2,4,3));
#rho_array= DimArray(repeat(Array(mask), 1, 1, 1,1, length(lev)), (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:lev}(lev), Dim{:time}(time), Dim{:lev2}(lev)));
rho_array= DimArray(repeat(Array(mask), 1, 1, 1,1), (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:lev}(lev), Dim{:time}(time)));

# Create an empty 4D DimensionalArray for rho
#rho_array = DimArray( mask, (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:time}(time), Dim{:lev}(lev)));

cell_area= DimArray(repeat(Array(climber_hos["area"][:,:,:]), 1, 1, 1, length(lev)), (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:time}(time), Dim{:lev}(lev)));
cell_area=permutedims(cell_area, (1,2,4,3));

for z =1:(length(levw)-1)
    println(levw[z+1])

    p = 0.1*(-levw[z+1])  # pressure in bar

    #t2 = t.^2
    #t3 = t.^3
    #t4 = t.^4
    t2_a = t[lev=z].^2
    t3_a = t[lev=z].^3
    t4_a = t[lev=z].^4
    
    # potential density -

    #r00 = a0[1].*mask .+ a0[2].*t .+ a0[3].*t2 .+ a0[4].*t3 .+ a0[5].*t4 .+ a0[6].*t.^5
    #A = a1[1].*mask .+ a1[2].*t .+ a1[3].*t2 .+ a1[4].*t3 .+ a1[5].*t4
    #B = b1[1].*mask .+ b1[2].*t .+ b1[3].*t2
    r00_a = a0[1] .+ a0[2].*t[lev=z] .+ a0[3].*t2_a .+ a0[4].*t3_a .+ a0[5].*t4_a .+ a0[6].*t[lev=z].^5
    A_a = a1[1] .+ a1[2].*t[lev=z] .+ a1[3].*t2_a .+ a1[4].*t3_a .+ a1[5].*t4_a
    B_a = b1[1] .+ b1[2].*t[lev=z] .+ b1[3].*t2_a
    C_a = c1
    s_sqrt_a=map(x -> ismissing(x) ? missing : sqrt(x),s[lev=z])
    #r0 = r00 .+ A.*s .+ B.*s.*s_sqrt .+ C.*s.^2
    r0_a = r00_a .+ A_a.*s[lev=z] .+ B_a.*s[lev=z].*s_sqrt_a .+ C_a.*s[lev=z].^2

    # bulk secant modulus
    #pk = c2[1].*mask .+ c2[2].*t .+ s.*(c2[3].*mask .- c2[4].*t) .+ p.*mask .*(c2[5].*mask .+ c2[6].*t) .+ p.*mask.*s.*c2[7]
    pk_a = c2[1] .+ c2[2].*t[lev=z] .+ s[lev=z].*(c2[3] .- c2[4].*t[lev=z]) .+ p .*(c2[5] .+ c2[6].*t[lev=z]) .+ p .*s[lev=z].*c2[7]

    rho_a = r0_a./(1 .- p./pk_a)

    #t2 = t.^2
    #t3 = t.^3
    #t4 = t.^4
    t2_b = t[lev=z+1].^2
    t3_b = t[lev=z+1].^3
    t4_b = t[lev=z+1].^4
    
  #  mask = map(x -> ismissing(x) ? missing : 1, t[lev=z+1]);

    # potential density -

    #r00 = a0[1].*mask .+ a0[2].*t .+ a0[3].*t2 .+ a0[4].*t3 .+ a0[5].*t4 .+ a0[6].*t.^5
    #A = a1[1].*mask .+ a1[2].*t .+ a1[3].*t2 .+ a1[4].*t3 .+ a1[5].*t4
    #B = b1[1].*mask .+ b1[2].*t .+ b1[3].*t2
    r00_b = a0[1] .+ a0[2].*t[lev=z+1] .+ a0[3].*t2_b .+ a0[4].*t3_b .+ a0[5].*t4_b .+ a0[6].*t[lev=z+1].^5
    A_b = a1[1] .+ a1[2].*t[lev=z+1] .+ a1[3].*t2_b .+ a1[4].*t3_b .+ a1[5].*t4_b
    B_b = b1[1] .+ b1[2].*t[lev=z+1] .+ b1[3].*t2_b
    C_b = c1
    s_sqrt_b=map(x -> ismissing(x) ? missing : sqrt(x),s[lev=z])
    #r0 = r00 .+ A.*s .+ B.*s.*s_sqrt .+ C.*s.^2
    r0_b = r00_b .+ A_b.*s[lev=z+1] .+ B_b.*s[lev=z+1].*s_sqrt_b .+ C_b.*s[lev=z+1].^2

    # bulk secant modulus
    #pk = c2[1].*mask .+ c2[2].*t .+ s.*(c2[3].*mask .- c2[4].*t) .+ p.*mask .*(c2[5].*mask .+ c2[6].*t) .+ p.*mask.*s.*c2[7]
    pk_b = c2[1] .+ c2[2].*t[lev=z+1] .+ s[lev=z+1].*(c2[3] .- c2[4].*t[lev=z+1]) .+ p .*(c2[5] .+ c2[6].*t[lev=z+1]) .+ p .*s[lev=z+1].*c2[7]

    rho_b = r0_b./(1 .- p./pk_b)

    rho_array[:, :, z,:] .= rho_b.-rho_a
end


function get_weighted_Ross_ts(VarData, AreaData, dimension; test=false, lon_val1=0, lon_val2=0, lat_val=-60)
    if (dimension=="4D")
        cell_area= DimArray(permutedims(repeat(AreaData[time=1],1,1,length(VarData.dims[4]),length(VarData.dims[3])), (1,2,4,3)), VarData.dims)[time=1:length(VarData.dims[4])]
    else
        cell_area= AreaData[time=1:length(VarData.dims[3])]
    end

    ross_area= cell_area[lon=Where(x -> x < -130 || x > 160), lat=(Where(<=(lat_val)))];
    var_ross = VarData[lon=Where(x -> x < -130 || x > 160), lat=(Where(<=(lat_val)))];
    mask = map(x -> ismissing(x) ? missing : 1, var_ross);
    ross_area_frac= ross_area.*mask ./ (mapslices(x -> Statistics.sum(skipmissing(x)),ross_area.*mask, dims=(:lon, :lat)));

    if (test)
   # if ((type=="ocn") & test)
        ross_area= cell_area[lon=Where(x -> x < lon_val1 || x > lon_val2), lat=(Where(<=(lat_val)))];
        var_ross = VarData[lon=Where(x -> x < lon_val1 || x > lon_val2), lat=(Where(<=(lat_val)))];
        mask = map(x -> ismissing(x) ? missing : 1, var_ross);
#        display(var_ross)
        ross_area_frac= ross_area.*mask ./ (mapslices(x -> Statistics.sum(skipmissing(x)),ross_area.*mask, dims=(:lon, :lat)));

        var_ross_wgt = var_ross.*ross_area_frac;
        ross_con_mean= dropdims(mapslices(x -> Statistics.sum(skipmissing(x)),var_ross_wgt, dims=(:lon, :lat)), dims=(:lon, :lat)) ;
       
    #elseif (type=="ocn")
       

    else
        var_ross_wgt = var_ross.*ross_area_frac;
        ross_con_mean= dropdims(mapslices(x -> Statistics.sum(skipmissing(x)),var_ross_wgt, dims=(:lon, :lat)), dims=(:lon, :lat)) ;
        println("HERE")
        # Ross Sea: (averaged over 160° −230°E and 60° −90°S) Beadling et al 2022
        
     #   var_ross = var_DimArr[lon=Where(x -> x < -130 || x > 160), lat=(Where(<=(-60)))];
      #  var_ross_wgt = var_ross.*ross_area_frac;
       # ross_con_mean= dropdims(mapslices(x -> Statistics.mean(skipmissing(x)),var_DimArr[lon=Where(x -> x < -130 || x > 160), lat=(Where(<=(-60)))], dims=(:lon, :lat)), dims=(:lon, :lat)) ;

    end

    return ross_con_mean
end


Ross_rho= get_weighted_Ross_ts(rho_array, climber_area, "4D")
function plot_depth_time(variable, name, crange; ax_limits=(nothing, (0, 4500)))
    fig = Figure();
    ax = Axis(fig[1,1], 
            title = name,
            xlabel = "time (years)",
            ylabel = "depth (m)",
            limits = ax_limits #(1500,1800)
            );
    
        colors=reverse(ColorSchemes.RdBu.colors);
        hm2 = heatmap!(ax, transpose(variable), colorrange=crange, colormap=colors)
        cb = Colorbar(fig[1, 2], hm2)
        ax.yreversed=true
    fig
end

plot_depth_time(rho_array[1,6,:,:], "lon=-172.5, lat=-62.5, month=3",(-0.15,0.15))

for z =1:length(lev)
    rho1=rho_array[:,:,:,:,z];
    #rho_array = permutedims(rho_array, (1,2,4,3));
    ross_area= cell_area[lon=Where(x -> x < -130 || x > 160), lat=(Where(<=(-60)))];
    var_ross = rho1[lon=Where(x -> x < -130 || x > 160), lat=(Where(<=(-60)))];
    mask = map(x -> ismissing(x) ? missing : 1, var_ross);
    #display(mask)
    ross_area_frac= ross_area.*mask ./ (mapslices(x -> Statistics.sum(skipmissing(x)),ross_area.*mask, dims=(:lon, :lat)));
    var_ross_wgt = var_ross.*ross_area_frac;
    ross_con_mean= dropdims(mapslices(x -> Statistics.sum(skipmissing(x)),var_ross_wgt, dims=(:lon, :lat)), dims=(:lon, :lat)) ;

    fig = Figure();
    ax = Axis(fig[1,1], 
            title = "ClimberX",
            xlabel = "time (years)",
            ylabel = "depth (m)",
            limits = (nothing, (0, 4500))
            );
    
        colors=reverse(ColorSchemes.RdBu.colors);
        hm = heatmap!(ax, transpose(ross_con_mean[1:22,:]), colormap=colors)
        cb = Colorbar(fig[1, 2], hm)
        ax.yreversed=true

    display(fig)
end