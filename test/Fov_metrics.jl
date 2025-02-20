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


mask= NCDataset("/albedo/work/user/anhoes003/climber-x/input/basin_mask_5x5.nc", "r");

mask=mask["basin_mask"];

dir_path="/albedo/work/projects/p_forclima/NaHosMIP/climber_x/control_run02";

ocean = NCDataset(joinpath(dir_path, "ocn.nc"), "r");


# freshwater fluxes by the MOC and gyre into the Atlantic at its southern border (Liu 2017, eq 1)
j = 12 ;   # index of south atlantic border

fovs = zeros(300);
fovs_mean = zeros(300);

for m in 1:12
  for k=1:length(ocean["v"][1,1,:,1,1])
    vzab = zeros(300);
    sza1 = zeros(300);
    sza2 = zeros(300);
    nxa1  = 0;
    nxa2  = 0;
    println(k);

    salt=replace!(Array(ocean["s"]), missing => 0);

    for i=1:length(ocean["v"][:,1,1,1,1])
      if ((mask[i,j]==1) &(mask[i,j+1]==1))
        #zonal integral of baroclinic meridional velocity
        if (ocean["v"][i,j,k,m,:]!=zeros(300))
          vzab = vzab .+ ((ocean["v"][i,j,k,m,:].-ocean["vb"][i,j,m,:]).*ocean["dxv"][j]); # m2/s
        else
          vzab = vzab .+ ((ocean["v"][i,j,k,m,:]).*ocean["dxv"][j]); # m2/s
        end
       # println(unique(vzab));

        #zonal mean salinity

        sza1 = sza1 .+ salt[i,j,k,m,:];
        sza2 = sza2 .+ salt[i,j+1,k,m,:];
        
        if (first(unique(salt[i,j,k,m,:]))!=0)
          nxa1 = nxa1+1;
        end
        if (first(unique(salt[i,j+1,k,m,:]))!=0)
          nxa2 = nxa2+1;
        end
      end
    end

    if (nxa1>0)
      sza1 = sza1./nxa1;
    end
    if (nxa2>0) 
      sza2 = sza2./(nxa2);
    end

    #integrate vertically
    fovs = fovs .- vzab.*0.5 .*(sza1.+sza2).* ocean["dz"][k] ;   #m3/s * psu
  end
  fovs_mean= fovs_mean .+fovs
end

#fovs_mean=fovs_mean./12
fovs = fovs./(34.7*12);  #m3/s

#freshwater fluxes by the MOC and gyre out of the Atlantic into the Arctic (Liu 2017, eq 1)
fovn = zeros(300);

for m in 1:12
  for k=1:length(ocean["v"][1,1,:,1,1])
    vzab = zeros(300);
    sza1 = zeros(300);
    sza2 = zeros(300);
    nxa1  = 0;
    nxa2  = 0;
    println(k);

    salt=replace!(Array(ocean["s"]), missing => 0);

    for i=1:length(ocean["v"][:,1,1,1,1])
      j = 33
      if ((10<i<56) & (mask[i,j]==1) &(mask[i,j+1]==1))
        #zonal integral of baroclinic meridional velocity
        if (ocean["v"][i,j,k,m,:]!=zeros(300))
          vzab = vzab .+ ((ocean["v"][i,j,k,m,:].-ocean["vb"][i,j,m,:]).*ocean["dxv"][j]); # m2/s
        else
          vzab = vzab .+ ((ocean["v"][i,j,k,m,:]).*ocean["dxv"][j]); # m2/s
        end
        # zonal mean salinity

        sza1 = sza1 .+ salt[i,j,k,m,:];
        sza2 = sza2 .+ salt[i,j+1,k,m,:];
        
        if (first(unique(salt[i,j,k,m,:]))!=0)
          nxa1 = nxa1+1;
        end
        if (first(unique(salt[i,j+1,k,m,:]))!=0)
          nxa2 = nxa2+1;
        end
      end
    end

    if (nxa1>0)
      sza1 = sza1./nxa1;
    end
    if (nxa2>0) 
      sza2 = sza2./(nxa2);
    end

    # integrate vertically
    fovn = fovn .- vzab.*0.5 .*(sza1.+sza2).* ocean["dz"][k] ;   #m3/s * psu

  end

end


fovn = fovn./(34.7*12) ;  #m3/s

fov = fovs.-fovn;

plot(fov*10^-6)