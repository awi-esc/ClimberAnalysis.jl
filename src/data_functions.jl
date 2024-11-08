
using NCDatasets
using NetCDF
using DimensionalData
using CommonDataModel


function getAveragedVar(data, varname::String, mean::Bool=true)
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

function getYearlyMean(data::CommonDataModel.CFVariable)
    data_yearly=data[:,:,13,:]
    return data_yearly
end


function getSeasonalMean(data::CommonDataModel.CFVariable, season::String, yearindex::Union{Number, Nothing} = nothing)
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
