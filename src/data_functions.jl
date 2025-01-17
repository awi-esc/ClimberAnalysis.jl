
using NCDatasets
using NetCDF
using DimensionalData
using CommonDataModel


function getAveragedVar(data, varname::String; mean::Bool=true)
    if mean
        var= Array(data[varname][:,:,13,end])
    else
        var=nothing
    end

    lon = data["lon"];
    lat = data["lat"];
    outdata = DimArray(Array(var), (Dim{:lon}(lon), Dim{:lat}(lat)))

    return outdata
end

function getAnomVar(data, varname::String, ind_year::Number; ref_year::Number=1)
    val1= Array(data[varname][:,:,13,ind_year])
    val2= Array(data[varname][:,:,13,ref_year])
    lon = data["lon"];
    lat = data["lat"];
    data1 = DimArray(Array(val1), (Dim{:lon}(lon), Dim{:lat}(lat)))
    data2 = DimArray(Array(val2), (Dim{:lon}(lon), Dim{:lat}(lat)))
    return data1-data2
end


function getDataYearlyMax(data, varname::String; year_ind::Union{Number,Nothing}=nothing)

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

function VarAsDimArray(data, varname::String, dimension::String; model::String="climberX")
    lon = data["lon"];
    lat = data["lat"];
    
    if model=="climberX" #has additional dimension for monthly data
        time = data["time"];
        if dimension=="3D"
            outdata = DimArray(data[varname][:,:,13,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:time}(time)));
        else
            lev = data["lev"]
            outdata = DimArray(data[varname][:,:,:,13,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:lev}(lev), Dim{:time}(time)))
        end

    elseif model=="mom5"
        time = 1:length(data["time"]);
        if dimension=="3D"
            outdata = DimArray(data[varname][:,:,:], (Dim{:lon}(lon), Dim{:lat}(lat), Dim{:time}(time)));
        else
            lev = data["st_ocean"]
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

mom5_data= NCDataset("/albedo/work/projects/p_forclima/NaHosMIP/cm2mc_PISM/concat_results_CM2Mc_PISM_coupling_after_spinup_y1860_norestore_uh03_extended/regridded/ocean-yearly_sub2_regridded.nc", "r"); 
VarAsDimArray(mom5_data,"temp", "4D", model="mom5")[time=1000:1500]
