using NCDatasets
using NetCDF
using DimensionalData
using CommonDataModel




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
        ross_area= cell_area[lon=Where(x -> x < lon_val1 || x > lon_val2), lat=(Where(<=(lat_val)))];
        var_ross = VarData[lon=Where(x -> x < lon_val1 || x > lon_val2), lat=(Where(<=(lat_val)))];
        mask = map(x -> ismissing(x) ? missing : 1, var_ross);
        ross_area_frac= ross_area.*mask ./ (mapslices(x -> Statistics.sum(skipmissing(x)),ross_area.*mask, dims=(:lon, :lat)));

        var_ross_wgt = var_ross.*ross_area_frac;
        ross_con_mean= dropdims(mapslices(x -> Statistics.sum(skipmissing(x)),var_ross_wgt, dims=(:lon, :lat)), dims=(:lon, :lat)) ;     
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



function get_weighted_SO_ts(VarData, AreaData, dimension; lat_val=-45)
    if (dimension=="4D")
        cell_area= DimArray(permutedims(repeat(AreaData,1,1,1,length(VarData.dims[3])), (1,2,4,3)), VarData.dims)[time=1:length(VarData.dims[4])]
    else
        cell_area= AreaData[time=1:length(VarData.dims[3])]
    end

    ross_area= cell_area[lat=(Where(<=(lat_val)))];
    var_ross = VarData[lat=(Where(<=(lat_val)))];
    mask = map(x -> ismissing(x) ? missing : 1, var_ross);
    ross_area_frac= ross_area.*mask ./ (mapslices(x -> Statistics.sum(skipmissing(x)),ross_area.*mask, dims=(:lon, :lat)));

    var_ross_wgt = var_ross.*ross_area_frac;
    ross_con_mean= dropdims(mapslices(x -> Statistics.sum(skipmissing(x)),var_ross_wgt, dims=(:lon, :lat)), dims=(:lon, :lat)) ;

    return ross_con_mean
end