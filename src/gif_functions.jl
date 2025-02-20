using NCDatasets
using NetCDF
using DimensionalData
using CairoMakie
using GeoMakie
using ColorSchemes
using SkipNan
using Animations

include("plot_functions.jl")


function creategif(data, type, depth, output_path; colors_opt= :RdBu, crange=(-5,5), SO=false, cont_anom=false )
    if (type=="S")
        factor=0.8;
        name="Salinity";
    elseif (type=="T")
        factor=0.052 .+0.012 .*data[lev= Near(depth),time=1];
        name="Temperature";
    elseif (type=="R")
        factor=1;
        name="Density";
    else
        println("PROBLEM")
    end

    if cont_anom
        vert_anom= (data[lev= 1] - data[lev= Near(depth)]) .-(data[lev= 1,time=1] - data[lev= Near(depth),time=1]);
    else
        vert_anom= data[lev= 1] - data[lev= Near(depth)];
       # vert_anom =(data[lev=1, time=1] .- data[lev=Near(depth), time=1])

    end
    ti=Observable(1)
    
    fig = Figure(size = (600, 600));
    ax = Axis(fig[1,1], 
        title = (@lift(name * " anomaly, 0 - " * string(depth) * "m, year $($ti*10)")),
        xlabel = "Longitude",
        ylabel = "Latitude",
        );
    dat_to_plot= @lift(factor .* vert_anom[time=$ti])
    p= global_map(dat_to_plot,fig,ax, 1 ,contours=false, crange_opt=true, crange=crange,colors=colors_opt, num_cont=13)
    record(p, output_path, 1:length(vert_anom.dims[3]),framerate = 10) do t
        ti[]=t
    end
end

#creategif(rho_hos,"R", 1000, "/albedo/work/user/anhoes003/data_analysis/gifs/maps_density_SO_con/Anom_rho_1000_sept_fr10.mp4")

function creategif_panel(data_list, data_mld, depth, output_path; colors_opt= :lipari, crange=(-5,0), cont_anom=false )
    types = ["T", "S", "R"]  # Specify the desired types
    names = ["Temperature", "Salinity", "Density"]
    
    ti=Observable(1)
    fig = Figure(size = (1300, 600));
    a=0
    for i in 1:length(types)
        name = names[i];
        if cont_anom
            vert_anom = (data_list[i][lev=1] - data_list[i][lev=Near(depth)]) .- (data_list[i][lev=1, time=1] - data_list[i][lev=Near(depth), time=1])
        else
            vert_anom = data_list[i][lev=1] .- data_list[i][lev=Near(depth)];
        end
    

        ax = GeoAxis(fig[1,i+a],dest="+proj=ortho +lon_0=0 +lat_0=-90", xgridwidth=.5, ygridwidth=.5, limits = ((-180, 180), (-90, -40)), width=300, xticklabelsvisible=false , yticklabelsvisible=false,  title = name) #xgridvisible=false,ygridvisible=false, xticklabelsvisible=false , yticklabelsvisible=false,

        dat_to_plot = @lift(vert_anom[time=$ti])
        global_map(dat_to_plot,fig,ax,contours=false, crange_opt=true, crange=crange,colors=colors_opt, num_cont=13)
        mld_data=@lift(data_mld[time=$ti]);
        contour!(ax, mld_data, labels = true, levels=[0+250*b for b=0:13], color=:black, labelsize=7);
        #contour!(ax, mld_data, labels = true, levels=[0+0.1*b for b=0:11], color=:black, labelsize=7);
        a=a+1
    end
    Label(fig[begin-1, 1:length(types)*2],(@lift("Anomalies, 0 - " * string(depth) * "m, year $($ti*10)")),font = "Nimbus Sans Bold", padding = (0, 0, 0, 0))
    record(fig, output_path, 1:length(data_list[1].dims[4]),framerate = 10) do t
        ti[]=t
    end
end

function creategif_panel_ext(data_list, data_mld, depth, output_path; colors_opt= :lipari, crange=(-5,0), cont_anom=false, sur_data=NaN )
    types = ["T", "S", "R"]  # Specify the desired types
    names = ["Temperature", "Salinity", "Density"]
    
    ti=Observable(1)
    fig = Figure(size = (1300, 800));
    a=0
    for i in 1:length(types)
        name = names[i];
        if cont_anom
            vert_anom = (data_list[i][lev=1] - data_list[i][lev=Near(depth)]) .- (data_list[i][lev=1, time=1] - data_list[i][lev=Near(depth), time=1])
            sur_anom = (data_list[i][lev=1] .- data_list[i][lev=1, time=1])
        else
            vert_anom = data_list[i][lev=1] .- data_list[i][lev=Near(depth)];
            sur_anom= sur_data[i];
        end
    

        ax1 = GeoAxis(fig[1,i],dest="+proj=ortho +lon_0=0 +lat_0=-90", xgridwidth=.5, ygridwidth=.5, limits = ((-180, 180), (-90, -40)), width=300, xticklabelsvisible=false , yticklabelsvisible=false,  title = "Surface" * name) #xgridvisible=false,ygridvisible=false, xticklabelsvisible=false , yticklabelsvisible=false,
        ax2 = GeoAxis(fig[2,i],dest="+proj=ortho +lon_0=0 +lat_0=-90", xgridwidth=.5, ygridwidth=.5, limits = ((-180, 180), (-90, -40)), width=300, xticklabelsvisible=false , yticklabelsvisible=false,  title = name) #xgridvisible=false,ygridvisible=false, xticklabelsvisible=false , yticklabelsvisible=false,

        dat_sur = @lift(sur_anom[time=$ti])
        global_map(dat_sur ,fig,ax1,contours=false, crange_opt=true, crange=crange,colors=colors_opt, num_cont=13)
        mld_data=@lift(data_mld[time=$ti]);
        contour!(ax1, mld_data, labels = true, levels=[0+250*b for b=0:13], color=:black, labelsize=7);
        dat_to_plot = @lift(vert_anom[time=$ti])
        global_map(dat_to_plot,fig,ax2,contours=false, crange_opt=true, crange=crange,colors=colors_opt, num_cont=13)

        a=a+1
    end
    Label(fig[begin-1, 1:length(types)*2],(@lift("Anomalies, 0 - " * string(depth) * "m, year $($ti*10)")),font = "Nimbus Sans Bold", padding = (0, 0, 0, 0))
    record(fig, output_path, 1:length(data_list[1].dims[4]),framerate = 10) do t
        ti[]=t
    end
end

