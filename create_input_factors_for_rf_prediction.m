%% Create input factors for random forest predictions
lon_modis = -179.975:0.05:179.975;
lat_modis = -89.975:0.05:89.975;
lat_modis = flip(lat_modis);

load('E:\work_for_wildfires\injecton_height\MERRA2_PBLH\lonlat_au_merra2.mat');
for year = 2021
    landuse = hdfread(['E:\work_for_wildfires\random_forest\MCD12C1.A',num2str(year-1),'001.006.hdf'],'Majority_Land_Cover_Type_1');
    for mon = 1
        if mon < 10
            filename = ['E:\work_for_wildfires\biomass_burning_emis\GFED\GFEDv4s_',num2str(year),'_0',num2str(mon),'_OC.nc'];
        else
            filename = ['E:\work_for_wildfires\biomass_burning_emis\GFED\GFEDv4s_',num2str(year),'_',num2str(mon),'_OC.nc'];
        end
        emis_var = ncread(filename,'OC');
        lon_gfed = ncread(filename,'lon');
        lat_gfed = ncread(filename,'lat');
        xdim     = size(emis_var,1);
        ydim     = size(emis_var,2);
        tdim     = size(emis_var,3);
        
        filename1 = ['E:\work_for_wildfires\injecton_height\MERRA2_PBLH\merra2_',num2str(year),num2str(mon),'_pblh_au.mat'];
        load(filename1);
        var_pblh = zeros(size(var_au,1),size(var_au,2),size(var_au,3)/24);
        for t = 1:size(var_au,3)/24
            var_pblh(:,:,t) = mean(var_au(:,:,(t-1)*24+1:24*t),3);
        end
        filename2 =  ['E:\work_for_wildfires\injecton_height\MERRA2_PBLH\merra2_',num2str(year),num2str(mon),'_t2m_au.mat'];
        load(filename2);
        var_t2m = zeros(size(var_au,1),size(var_au,2),size(var_au,3)/24);
        for t = 1:size(var_au,3)/24
            var_t2m(:,:,t) = mean(var_au(:,:,(t-1)*24+1:24*t),3);
        end
        filename3 =  ['E:\work_for_wildfires\injecton_height\MERRA2_PBLH\merra2_',num2str(year),num2str(mon),'_rh_au.mat'];
        load(filename3);
        var_rh = zeros(size(var_au,1),size(var_au,2),size(var_au,3)/8);
        for t = 1:size(var_au,3)/8
            var_rh(:,:,t) = mean(var_au(:,:,(t-1)*8+1:8*t),3);
        end
        filename4 =  ['E:\work_for_wildfires\injecton_height\MERRA2_PBLH\merra2_',num2str(year),num2str(mon),'_u10m_au.mat'];
        load(filename4);
        var_u10m = zeros(size(var_au,1),size(var_au,2),size(var_au,3)/24);
        for t = 1:size(var_au,3)/24
            var_u10m(:,:,t) = mean(var_au(:,:,(t-1)*24+1:24*t),3);
        end
        filename5 =  ['E:\work_for_wildfires\injecton_height\MERRA2_PBLH\merra2_',num2str(year),num2str(mon),'_v10m_au.mat'];
        load(filename5);
        var_v10m = zeros(size(var_au,1),size(var_au,2),size(var_au,3)/24);
        for t = 1:size(var_au,3)/24
            var_v10m(:,:,t) = mean(var_au(:,:,(t-1)*24+1:24*t),3);
        end
        filename6 =  ['E:\work_for_wildfires\injecton_height\MERRA2_PBLH\merra2_',num2str(year),num2str(mon),'_prectot_au.mat'];
        load(filename6);
        var_prectot = zeros(size(var_au,1),size(var_au,2),size(var_au,3)/24);
        for t = 1:size(var_au,3)/24
            var_prectot(:,:,t) = mean(var_au(:,:,(t-1)*24+1:24*t),3);
        end
        
        if mon >= 10
            filename = ['E:\work_for_wildfires\MODIS_C6\MCD14ML\MCD14ML.',num2str(year),num2str(mon),'.006.03.mat'];
        else
            filename = ['E:\work_for_wildfires\MODIS_C6\MCD14ML\MCD14ML.',num2str(year),'0',num2str(mon),'.006.03.mat'];
        end
        load(filename);
        fire_tmp = [YYYYMMDD,HHMM,lat,lon,T21,T31,FRP,conf];
       %% prepare the input factors for machine learning model
        input_factors = zeros(xdim,ydim,tdim,12);
        for ii = 1:xdim
            for jj = 1:ydim
                lon_emis = lon_gfed(ii);
                lat_emis = lat_gfed(jj);
                for tt = 1:tdim
                    if emis_var(ii,jj,tt) > 0 && lon_emis < 155 && lon_emis > 110 && lat_emis > -45 && lat_emis < -9 % for emission grids 
                        % landuse
                        aa1 = find( abs(lon_emis - lon_modis) == min(abs(lon_emis - lon_modis)));
                        bb1 = find( abs(lat_emis - lat_modis) == min(abs(lat_emis - lat_modis)));
                        input_factors(ii,jj,tt,1)  = landuse(bb1(1),aa1(1));
                        % pblh
                        aa = find( abs(lon_emis - lon_au) == min(abs(lon_emis - lon_au)));
                        bb = find( abs(lat_emis - lat_au) == min(abs(lat_emis - lat_au)));
                        input_factors(ii,jj,tt,2) = var_pblh(aa(1),bb(1),tt);
                        % t2m
                        input_factors(ii,jj,tt,3) = var_t2m(aa(1),bb(1),tt);
                        % rh
                        input_factors(ii,jj,tt,4) = var_rh(aa(1),bb(1),tt);
                        % u10m
                        input_factors(ii,jj,tt,5) = var_u10m(aa(1),bb(1),tt);
                        % v10m
                        input_factors(ii,jj,tt,6) = var_v10m(aa(1),bb(1),tt);
                        % prectot
                        input_factors(ii,jj,tt,7) = var_prectot(aa(1),bb(1),tt);
                        % biomass burning emission
                        input_factors(ii,jj,tt,8)    = emis_var(ii,jj,tt);
                        % mon
                        input_factors(ii,jj,tt,9)    = mon;
                        % lon
                        input_factors(ii,jj,tt,10)   = lon_emis;
                        % lat
                        input_factors(ii,jj,tt,11)   = lat_emis;
                        % fire radiative power
                        day_yyyymmdd = year*10000+mon*100+tt;
                        xx = find(fire_tmp(:,1) == day_yyyymmdd & fire_tmp(:,8) > 10 & abs(lon_emis - fire_tmp(:,4)) < 0.5 & abs(lat_emis - fire_tmp(:,3)) < 0.5);
                        if ~isempty(xx)
                            frp_tmp = fire_tmp(min(xx):max(xx),7);
                            input_factors(ii,jj,tt,12) = max(frp_tmp);
                        else
                            xx = find(fire_tmp(:,1) == day_yyyymmdd & fire_tmp(:,8) > 10 & abs(lon_emis - fire_tmp(:,4)) < 1 & abs(lat_emis - fire_tmp(:,3)) < 1);
                            if ~isempty(xx)
                                frp_tmp = fire_tmp(min(xx):max(xx),7);
                                input_factors(ii,jj,tt,12) = max(frp_tmp);
                            else
                                input_factors(ii,jj,tt,12) = 0;
                            end
                        end
                
                    end
                end
            end
        end
        filename_output = ['E:\work_for_wildfires\random_forest\new\input_factors_12_',num2str(year),'_',num2str(mon),'.mat'];
        save(filename_output,'input_factors','-v7.3');
    end
end
        