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

function getColorbarTicks(data::DimArray, num_cont::Number)
    """
        calculates range of colorbar: range from floor(min) to ceil(max) 
        for setting Colorbarticks, step size so that num_cont+1 ticks fit between min and max    
    """
    max = ceil(maximum(skipmissing(data)))
    min = floor(minimum(skipmissing(data)))
    step = Int(floor((max-min)/num_cont))
    
    return min, max, step
end


function global_map(data::DimArray, title::String; crange_opt=false, crange=(0,100), colors::Union{Symbol, Nothing}=nothing, contours::Bool=true, num_cont::Number= 6, save_out=false, outfile="Unknown_output", output_path=" ")
    """
        function to plot a single valiable on worldmap
        #Arguments
        - data: Target Variable as DimArray Use eg getAveragedVar function to get right format of data
        - title: Set plot title

        optional:
        - colors: name of ColorScheme
        - contours: Set false to deactivate
        - num_con: set number of contourlines (#=num_con+1)
    """

    fig = Figure(size = (800, 600));
    ax = Axis(fig[1,1], 
    title = title,
    xlabel = "Longitude",
    ylabel = "Latitude",
    );

    if isnothing(colors)
        colors=reverse(ColorSchemes.RdBu.colors);
    end

    if crange_opt
        min=crange[1];
        max=crange[2];
        step= num_cont;
    else
        max = ceil(maximum(skipmissing(data)))
        min = floor(minimum(skipmissing(data)))
        step = Int(floor((max-min)/num_cont))
        println(min, max, step)
        min, max, step = getColorbarTicks(data, num_cont)
    end

    #plot data as map
    hm = heatmap!(ax, data, colormap= colors, colorrange= (min, max), nan_color= :lightgray);
    
    #plot coastlines
    lines!(ax, GeoMakie.coastlines(); color=:grey);

    #plot contours
    if contours
        contour!(ax, data, levels= min:step:max, labels = false, color=:black);
        Colorbar(fig[1,2], hm, ticks = min:step:max);
    else
        Colorbar(fig[1,2], hm)#, ticks = min:step:max);
    end
    
    if save_out
        save(joinpath(output_path, outfile), fig)
    end
    
    display(fig) 
end

function map_panel(data_list , titles::Vector{String}; SO::Bool=false, colors::Union{Vector{Symbol}, Nothing}=nothing)
    """ 
        plot 2x2 plots 
        
        #Argumenents:
        - data_list: list of 4 DimArrays; one for each variable 
        - titles: list of 4 titles
        - SO: Boolean; if true all plots will have Southern Ocean projection; default: false
        - colors: list of 4 color symbols; default: RdBu
    """
    fig = Figure(size = (1000, 700))
    letters = reshape(collect('a':'d'), (2, 2))
    k_mapping = Dict((1, 1) => 1, (1, 2) => 2,(2, 1) => 3, (2, 2) => 4, (3, 1) => 3, (3, 2) => 4)

    # define 8 axis (4 of them for colorbars)
    if SO
        gas = [GeoAxis(fig[i,j],dest="+proj=ortho +lon_0=0 +lat_0=-90", xgridwidth=.5, ygridwidth=.5, limits = ((-180, 180), (-90, -60)), width=400,  title = titles[k_mapping[(i,j)]]) for i = 1:2:3, j = 1:2] #xgridvisible=false,ygridvisible=false, xticklabelsvisible=false , yticklabelsvisible=false,
    else
        gas = [GeoAxis(fig[i,j],dest="+proj=eqc", xgridwidth=.5, ygridwidth=.5, width=400,  title = titles[k_mapping[(i,j)]]) for i = 1:2:3, j = 1:2] #xgridvisible=false,ygridvisible=false, xticklabelsvisible=false , yticklabelsvisible=false,
    end

    #set default colors to RdBu
    if isnothing(colors)
        colors=[reverse(ColorSchemes.RdBu.colors) for i in 1:4];
    end

    if SO
        hms = [heatmap!(gas[i,j],data_list[k_mapping[(i, j)]],colormap= colors[k_mapping[(i, j)]], nan_color= :lightgray, colorrange = (-1, 1) .* maximum(abs, skipnan(Array(data_list[k_mapping[(i, j)]][lat=At(-90:5:-60;atol=2.5)])))) for i =1:2, j = 1:2]
    else
        hms = [heatmap!(gas[i,j],data_list[k_mapping[(i, j)]],colormap= colors[k_mapping[(i, j)]], nan_color= :lightgray, colorrange = (-1, 1) .* maximum(abs, skipnan(Array(data_list[k_mapping[(i, j)]])))) for i =1:2, j = 1:2]
    end
    
    con= [contour!(gas[i,j],data_list[k_mapping[(i, j)]], levels= ((-1).* maximum(abs, skipnan(Array(data_list[k_mapping[(i, j)]])))):(0.2.*maximum(abs, skipnan(Array(data_list[k_mapping[(i, j)]])))):(maximum(abs, skipnan(Array(data_list[k_mapping[(i, j)]])))), labels = false, color=:black) for i =1:2, j = 1:2];

    lns= [lines!(gas[i,j], GeoMakie.coastlines(), color=:grey) for i =1:2, j = 1:2];

    Colorbar(fig[2, 1], hms[1],  vertical= false, flipaxis = false, width =200, scale = .5)
    Colorbar(fig[2, 2], hms[3],  vertical= false, flipaxis = false, width= 200, scale = .5)
    Colorbar(fig[4, 1], hms[2], vertical= false, flipaxis = false, width =200)
    Colorbar(fig[4, 2], hms[4], vertical= false, flipaxis = false, width= 200)

    [Label(fig[i, j, TopLeft()], "($(letters[k_mapping[(i,j)]]))", fontsize=16,
        padding=(-2, 0, -20, 0)) for i = 1:2:3, j = 1:2]

    fig 
end



function amoc_map_panel(map_list, opsi_list , data_si, titles::Vector{String})#colors::Vector{Union{String, Nothing}}=[nothing], contours::Bool=false,
    fig = Figure(size = (1000, 600))
    letters = collect('a':'h')
    gas = [GeoAxis(fig[1,j], dest = "+proj=ortho +lat_0= 60 +lon_0= -40",xgridvisible=false,ygridvisible=false, xticklabelsvisible=false , yticklabelsvisible=false, limits = ((-90, 90), (20, 690)), height=300, xlabelpadding = 18, title = titles[j]) for j = 1:4]

    axs = [Axis(fig[2, j], height = 150, title = titles[4+j]) for j = 1:4]

    colors = reverse(ColorSchemes.magma.colors);

    hms = [heatmap!(gas[j], map_list[j][:,:,1], colormap = colors, nan_color= :white) for j = 1:4]
    con= [contour!(gas[j], data_si, levels= [0,1], labels = false, color=:black) for j =1:4];

    lns= [lines!(gas[j], GeoMakie.coastlines(), color=:grey,) for j = 1:4];


    Colorbar(fig[1, 5], hms[1], label="m", height =150)

    for ax in axs
        xlims!(ax, (-30,90))
        ax.yreversed=true
    end


    
    extremas = map(extrema, filter(!isnan,opsi_list[4]))
    global_min = minimum(t->first(t), extremas)
    extremas = map(extrema, filter(!isnan,opsi_list[1]))
    global_max = maximum(t->last(t), extremas)
    # these limits have to be shared by the maps and the colorbar
    clims = (-global_max, global_max)
    println(clims)
    prfs = [heatmap!(axs[j], opsi_list[j], colormap= cgrad(reverse(ColorSchemes.RdBu.colors),26, categorical=true), nan_color= :gray; colorrange=clims) for j = 1:4]
    
    Colorbar(fig[2, 5], prfs[1], label="Sv", height= 200)

    [Label(fig[i, j, TopLeft()], "($(letters[(i-1)*4+j]))", fontsize=16,
        padding=(-2, 0, -20, 0)) for i = 1:2, j = 1:4]

    rowgap!(fig.layout, 0.4)

    return fig 
end


function hysteresis_plot(data_up::Dataset, data_down::Dataset, title::String)

    fig = Figure();
    ax = Axis(fig[1,1], 
        title = title,
        xlabel = "Freshwater hosing (Sv)",
        ylabel = "max AMOC (Sv)",
        );
        lines!(ax, data_up[:hosing], data_up[:omaxa]);
        lines!(ax, data_down[:hosing], data_down[:omaxa]);

    fig

end

function timeseries(data::Vector, ylabels::Vector{String}; mov_average::Bool=true, wind::Number=100)
    f = Figure(size = (800, 800))
    list_ax =[]
    for i in 1:(length(data)-1)
        ax = Axis(f[i,1], xtickformat = "{:.0f}", xlabel="time (years)", ylabel = ylabels[i])
        if i!= (length(data)-1)
            hidexdecorations!(ax, ticks=false, grid=false)
        end
        lines!(ax, data[1], data[i+1])
        vlines!(ax, [1630], color=:gray, alpha=0.7)
        if mov_average
            lines!(ax, data[1], movmean(data[i+1], wind))
        end

        push!(list_ax,ax)
    end
    #linkxaxes!(list_ax)
    rowgap!(f.layout, 0.2)

    f
end


#function to create Hovmöller plot

function plot_depth_time(variable, name, crange; ax_limits=(nothing, (0, 4500)), save_out=false, outfile="Unknown_output")
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
    if save_out
        save(joinpath(output_path, outfile), fig)
    end
    fig
end



function plot_Ross_temp_lat_ts(VarData,type)
    var_ross = VarData[lon=Where(x -> x < -130 || x > 160)];
    display(var_ross)

    if (type=="mom5")
        depths=[18+l for l=1:8];
    else
        depths=[13+l for l=1:8];
    end

    for i in depths
        ross_mean= dropdims(mapslices(x -> Statistics.mean(skipmissing(x)),var_ross, dims=(:lon)), dims=(:lon))[lev=i] ;
        fig = Figure();
        ax = Axis(fig[1,1], 
                title = "ClimberX",
                xlabel = "time (years)",
                ylabel = "depth (m)",
                limits = (nothing, (-90, 90))
                );
        
            colors=reverse(ColorSchemes.RdBu.colors);
            hm = heatmap!(ax, transpose(ross_mean), colorrange = (0,8), colormap=colors)
            cb = Colorbar(fig[1, 2], hm)
            ax.yreversed=true
    
        display(fig)
    end
end
    
function plotSI_RossTS(variables, names; save_out=false, outfile="Unknown_output")

    f = Figure(title = "Hosing 0.3 Sv (no restore): Ross Sea Mean",size = (800, 800) )
    list_ax =[]
    for i=1:length(names)
        data1= get_weighted_Ross_ts(variables[i], climber_area, "3D")
        data2= get_weighted_Ross_ts(variables[i], climber_area, "3D", lat_val=-30)
        data3= get_weighted_SO_ts(variables[i], climber_area, "3D", lat_val=-30)
        if (i%2==1)
            ax = Axis(f[i,1],  xtickformat = "{:.0f}", xlabel="time (years)",yticklabelpad = 64.0, ylabel = names[i])
        else
            ax = Axis(f[i,1],  xtickformat = "{:.0f}", xlabel="time (years)", ylabel = names[i])
        end

        if i!= (length(names))
            hidexdecorations!(ax, ticks=false, grid=false)
        end
        lines!(ax, data1, label= "Ross")
        lines!(ax, data2, label= "Ross extended")
        lines!(ax, data3, label= "Southern Ocean")

        vlines!(ax, [1630], color=:gray, alpha=0.7)
        push!(list_ax,ax)

        if (i==1)
            axislegend(ax, position = :lc );
        end
    end
    rowgap!(f.layout, 0.2)

    if save_out
        save(joinpath(output_path, outfile), f)
    end
    f
end
