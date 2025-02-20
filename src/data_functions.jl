
using NCDatasets
using NetCDF
using DimensionalData
using CommonDataModel


function getAveragedVar(data, varname::String; mean::Bool=true, year_ind=nothing, month_ind=9)
# Simple function to convert climberX data into a DimArray. Here: get annual mean of a 3D variable for a specific year. 
# Output: 2D DimArray with dims lat, lon

    # default year should be the end of the simulation
    if year_ind === nothing
        year_ind = length(data[varname]["time"]) 
    end

    if mean
        var= Array(data[varname][:,:,13,year_ind])
    else
        var= Array(data[varname][:,:, month_ind,year_ind])
    end

    lon = data["lon"];
    lat = data["lat"];
    outdata = DimArray(Array(var), (Dim{:lon}(lon), Dim{:lat}(lat)))

    return outdata
end


function getAnomVar(data, varname::String, ind_year::Number; ref_year::Number=1)
# Simple function to calculate anomaly of climberX data into a DimArray. Here: get difference of annual means of two 3D variables for specific years. 
# Output: 2D DimArray with dims lat, lon
    data1 = getAveragedVar(data, varname, year_ind=ind_year)
    data2 = getAveragedVar(data, varname, year_ind=ref_year)
    return data1-data2
end


function getDataYearlyMax(data, varname::String; year_ind::Union{Number,Nothing}=nothing)
# Simple function to get maximum climberX variable in a DimArray. Here: annual max of a 3D variable for a specific year. 
# Output: 2D DimArray with dims lat, lon

    if isnothing(year_ind)
        var= maximum(Array(data[varname][:,:,1:12,end]),dims=3)[:,:,1]
    else
        var= maximum(Array(data[varname][:,:,1:12,year_ind]),dims=3)[:,:,1]
    end

    lon = data["lon"];
    lat = data["lat"];
    outdata = DimArray(Array(var), (Dim{:lon}(lon), Dim{:lat}(lat)))

    return outdata

end

function getDataDepthProfile(data, varname::String; year_ind::Union{Number,Nothing}=nothing)

    if isnothing(year_ind)
        var= Array(data[varname][:,:,13,end])
    else
        var= Array(data[varname][:,:,13,year_ind])
    end

    lat = data["latv1"];
    levw = data["levw"];
    var = ifelse.(iszero.(var), NaN, var)
    outdata = DimArray(Array(var), (Dim{:latv1}(lat),Dim{:levw}(levw)))

    return outdata

end

function getYearlyMean(data::CommonDataModel.CFVariable)
    data_yearly=data[:,:,13,:]
    return data_yearly
end


function getSeasonalMean(data::CommonDataModel.CFVariable, season::String; yearindex::Union{Number, Nothing} = nothing)
    if isnothing(yearindex)
        if season == "DJF"
            data_out=mean(data[:,:,[1,2,12],end],dims=3)
        elseif season == "MAM"
            data_out=mean(data[:,:,[3,4,5],end],dims=3)
        elseif season =="JJA"
            data_out=mean(data[:,:,[6,7,8],end],dims=3)
        elseif season == "SON"
            data_out=mean(data[:,:,[9,10,11],end],dims=3)
        else
            println("Season not defined. Choose from: DJF, MAM, JJA, SON")#
        end
    else
        if season == "DJF"
            data_out=mean(data[:,:,[1,2,12],yearindex],dims=3)
        elseif season == "MAM"
            data_out=mean(data[:,:,[3,4,5],yearindex],dims=3)
        elseif season =="JJA"
            data_out=mean(data[:,:,[6,7,8],yearindex],dims=3)
        elseif season == "SON"
            data_out=mean(data[:,:,[9,10,11],yearindex],dims=3)
        else
            println("Season not defined. Choose from: DJF, MAM, JJA, SON")#
        end
    end    
    return data_out
end



# get DimArray with yearly values for diffent models and dimensions

function VarAsDimArray(data, varname::String, dimension::String; model::String="climberX", month::Number=13, decadal=false)
    lon = data["lon"];
    lat = data["lat"];
    
    if model=="climberX" #has additional dimension for monthly data
        time = data["time"];
        if dimension=="3D"
            outdata = DimArray(data[varname][:,:,month,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:time}(time)));
        else
            lev = data["lev"];
            outdata = DimArray(data[varname][:,:,:,month,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:lev}(lev), Dim{:time}(time)))
        end

    elseif model=="mom5"
        time = 1:length(data["time"]);
        if decadal
            time = 1:10:length(data["time"])*10;
            println(time)
        end
        if dimension=="3D"
            outdata = DimArray(data[varname][:,:,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:time}(time)));
        else
            lev = data["st_ocean"];
            outdata = DimArray(data[varname][:,:,:,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:lev}(lev), Dim{:time}(time)))
        end
    end

    return outdata
end

# get DimArray with yearly values for diffent models and dimensions

function ClimberAreaAsDimArray(data)
    lon = data["lon"];
    lat = data["lat"];
    time = data["time"];
    outdata = DimArray(data["area"][:,:,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:time}(time)));

    return outdata
end


function DownscaleTime(data, timestep; confunc="mean")
#function to change 3D Array with dim 3 = time to a lower time resolution, eg. timestep= 10 allows to go from annual to decadal data
    if (confunc=="max")
        test=mapslices(x -> Statistics.maximum(x),data[:,:,1:timestep], dims=(:time))
        for i in timestep:timestep:(length(data.dims[3])-timestep)
            test=cat(test,mapslices(x -> Statistics.maximum(x),data[:,:,1+i:timestep+i], dims=(:time)),dims=3)
        end

    else
        test=mapslices(x -> Statistics.mean(skipmissing(x)),data[:,:,:,1:timestep], dims=(:time))
        for i in timestep:timestep:(length(data.dims[4])-timestep)
            test=cat(test,mapslices(x -> Statistics.mean(skipmissing(x)),data[:,:,:,1+i:timestep+i], dims=(:time)),dims=4)
        end
    end

    return test
end



function conv_to_rho_unit(data, ref_data, type)
# function for temperature (type=T) or salinity (type=S) data to convert them to density units following Appendix D5 in https://cp.copernicus.org/articles/20/2719/2024/cp-20-2719-2024.pdf 
    if type=="T"
        alpha=-(0.052 .+0.012 .* ref_data);
        data_new= alpha .* data;
    elseif type=="S"
        beta=0.8;
        data_new= beta .* data;
    else
        println("Type not correct! Choose either 'T' or 'S'");
    end
    return data_new
end
