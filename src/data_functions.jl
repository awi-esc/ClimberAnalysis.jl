
using NCDatasets
using NetCDF
using DimensionalData
using CommonDataModel

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


