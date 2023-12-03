%% This script calculates the smoke exposure
%% set some parameters
for years       = 2021
    clearvars -except years
start_mon       = 1;
end_mon         = 1;
start_date      = 1;
end_date        = 31;
backward_hours  = 120;
rundays         = datenum(years,end_mon,end_date) - datenum(years,start_mon,start_date) + 1;
receptors       = 11;

if start_mon == 1 && start_date == 1
    [aemit_daily_1b,aemit_daily_2b,aemit_daily_3b] = emis_dec_before_yr(years-1);
end

if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
    ndays = 366;
    ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
else
    ndays = 365;
    ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
end

if start_mon == 1
    tt_indx1 = 1;
    tt_indx2 = ndays_mon(end_mon);
else
    tt_indx1 = sum(ndays_mon(1:start_mon-1))+1;
    tt_indx2 = sum(ndays_mon(1:end_mon));
end
% AU SE
% 1. Footscray; 2. Albury; 3. Florey; 4. Wollongong; 
% 5. Prospect; 6. Newcastle; 7. Mountain Creek; 8. Springwood; 9. Alphington; 10. Liverpool; 11. Wallsend; 
loc_str = {'144.8728027_-37.80487823','146.93986_-36.05182','149.043539_-35.220606','150.88733_-34.41706','150.91417_-33.79424','151.75965_-32.9312','153.1038_-26.6917','153.1356_-27.6125','145.0306_-37.7784','150.9058_-33.9328','151.6692_-32.8961'};
% AU NORTH
% 1. Palmerston; 2. West Mackay; 3. Coast guard; 4. South Gladstone; 
% loc_str = {'130.94853_-12.50779','146.8257_-19.2542','149.1549_-21.1595','151.2704_-23.8627'};
% loc_lon = [130.94853,146.8257,149.1549,151.2704];
% loc_lat = [-12.50779,-19.2542,-21.1595,-23.8627];
%% Step 1: calculate the biomass burning emissions
disp('Step 1');
location = 'E:\work_for_wildfires\biomass_burning_emis\GFED\';
prefix = 'GFEDv4s_';
% prefix = 'QFEDv2p5r1_';
% prefix = 'FINNv1p5_';
postfix = '.nc';

spc_name = {'BC','OC','CO','PM25'};

dirfiles = dir([location,prefix,'*']);
filename1 = dirfiles.name;
lon_global_gfed = ncread([location,filename1],'lon');
lat_global_gfed = ncread([location,filename1],'lat');

for years = years
    % leap year
    if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
        ndays = 366;
        ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
    else
        ndays = 365;
        ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
    end
    % target emission unit: kg m^-^2 s^-^1
    aemit_daily = zeros(length(lon_global_gfed),length(lat_global_gfed),ndays,2); % 1: BC; 2: OC
    for ispc = 1:2
        
        for imon = 1:12
            if imon < 10
                filename = [location,prefix,num2str(years),'_0',num2str(imon),'_',char(spc_name(ispc)),postfix];
            else
                filename = [location,prefix,num2str(years),'_',num2str(imon),'_',char(spc_name(ispc)),postfix];
            end
            
            if imon == 1
                startday = 1;
                endday = sum(ndays_mon(imon));
            else
                startday = sum(ndays_mon(1:imon-1))+1;
                endday = sum(ndays_mon(1:imon));
            end
            
            disp(filename);
            aemit_daily(:,:,startday:endday,ispc) = ncread(filename,char(spc_name(ispc)));
        end
    end
end
% 
aemit_daily_1 = aemit_daily; % for scale emission mami
aemit_daily_2 = aemit_daily; % no_scale
aemit_daily_3 = aemit_daily; % for scale emission mami & fix non-injection height

filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'01_gfs_nohnf_receptor_se\'];
diroutput = dir(fullfile(filepath,'*_foot.nc'));
filename = {diroutput.name};
filename_1 = char(filename(1));
lon = ncread([filepath,filename_1],'lon');
lat = ncread([filepath,filename_1],'lat');

lon_indx1 = find( abs(min(lon) - lon_global_gfed) == min( abs(min(lon) - lon_global_gfed) ));
lon_indx2 = find( abs(max(lon) - lon_global_gfed) == min( abs(max(lon) - lon_global_gfed) ));
lat_indx1 = find( abs(min(lat) - lat_global_gfed) == min( abs(min(lat) - lat_global_gfed) ));
lat_indx2 = find( abs(max(lat) - lat_global_gfed) == min( abs(max(lat) - lat_global_gfed) ));

gfed_var = aemit_daily(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),tt_indx1:tt_indx2,:);
lon_gfed_au = lon_global_gfed(lon_indx1(1):lon_indx2(1));
lat_gfed_au = lat_global_gfed(lat_indx1(1):lat_indx2(1));
%% Step 2: use injection height to scale the biomass burning emissions
disp('Step 2 - emis 1');
clearvars gfed_all inheight_all pblh_all
% Turn on vertical scale
is_injection_height_scale = 1;

inum = 0;
if is_injection_height_scale
    % combine the injection height in the specified year
    filename_gfas = ['E:\work_for_wildfires\injecton_height\',num2str(years),num2str(1),'_gfas_injection_height.nc']; % unit:m
    if years ~= 2020
        lon_global_gfas = ncread(filename_gfas,'longitude');
        lat_global_gfas = flip(ncread(filename_gfas,'latitude'));
    else
      lon_global_gfas = ncread(filename_gfas,'lon'); % for 2020
      lat_global_gfas = ncread(filename_gfas,'lat'); % for 2020
    end

    for imon = 1%1:12
        filename_gfas = ['E:\work_for_wildfires\injecton_height\',num2str(years),num2str(imon),'_gfas_injection_height.nc'];
        injh_tmp = ncread(filename_gfas,'mami');
        injh_tmp = flip(injh_tmp,2); % reverse latitude
        
        lon_indx1 = find( abs(min(lon) - lon_global_gfas) == min( abs(min(lon) - lon_global_gfas) ));
        lon_indx2 = find( abs(max(lon) - lon_global_gfas) == min( abs(max(lon) - lon_global_gfas) ));
        lat_indx1 = find( abs(min(lat) - lat_global_gfas) == min( abs(min(lat) - lat_global_gfas) ));
        lat_indx2 = find( abs(max(lat) - lat_global_gfas) == min( abs(max(lat) - lat_global_gfas) ));
        injh_tmp  = injh_tmp(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);
        
        if imon == 1
            injh_year = injh_tmp;
        else
            injh_year = cat(3,injh_year,injh_tmp);
        end
        clearvars injh_tmp
    end
    lon_gfas_au = lon_global_gfas(lon_indx1(1):lon_indx2(1));
    lat_gfas_au = lat_global_gfas(lat_indx1(1):lat_indx2(1));
    inheight_var = injh_year(:,:,tt_indx1:tt_indx2);
    
    % combine the pblh in the specified year
    for imon = 1%1:12
        filename = ['E:\work_for_wildfires\injecton_height\MERRA2_PBLH\merra2_',num2str(years),num2str(imon),'_pblh_au.mat'];
        load(filename);
        load('E:\work_for_wildfires\injecton_height\MERRA2_PBLH\lonlat_au_merra2.mat');
        if imon == 1
            pblh_year = var_au;
        else
            pblh_year = cat(3,pblh_year,var_au);
        end
        
        pblh_var = zeros(size(pblh_year,1),size(pblh_year,2),size(pblh_year,3)/24);
        for tt = 1:size(pblh_year,3)/24
            pblh_var(:,:,tt) = mean(pblh_year(:,:,(tt-1)*24+1:24*tt),3);
        end
        lon_pblh_au = lon_au;
        lat_pblh_au = lat_au;
    end
    pblh_var = pblh_var(:,:,tt_indx1:tt_indx2);
    clearvars var_au
    
    % load fire injection height profile GMD or IS4FIRE
    filename = 'E:\work_for_wildfires\injecton_height\emis_mean_normalised_wholeday_smooth.nc';
    lon_prof = ncread(filename,'longitude');
    lat_prof = ncread(filename,'latitude');
    alt_prof = ncread(filename,'z'); % unit: m
    fire_prof = ncread(filename,'ems_norm_wholeday_smooth');
    
    % compare the injection height and pblh and scale the biomass burning emissions
    gfed_var_new = gfed_var;
    ratio_all_1 = zeros(size(gfed_var,1),size(gfed_var,2),size(gfed_var,3));
    
    for ispc = 1:2
        for tt = 1:size(gfed_var,3)
            for jj = 1:size(gfed_var,2)
                for ii = 1:size(gfed_var,1)
                    if gfed_var(ii,jj,tt,ispc) > 0
                        aa = find( abs(lon_gfed_au(ii) - lon_pblh_au) == min(abs(lon_gfed_au(ii) - lon_pblh_au)) );
                        bb = find( abs(lat_gfed_au(jj) - lat_pblh_au) == min(abs(lat_gfed_au(jj) - lat_pblh_au)) );
                        xx = find( abs(lon_gfed_au(ii) - lon_gfas_au) < 1);
                        yy = find( abs(lat_gfed_au(jj) - lat_gfas_au) < 1);
                        % Find the largest inheight
                        injh_area = squeeze(inheight_var(min(xx):max(xx),min(yy):max(yy),tt));
                        inheight_tmp = max(max(injh_area));
                        
                        % Fix non-injection height grid
    %                   inheight_max = 1e4;
    %                   if inheight_tmp < 1e-5
    %                       inheight_tmp = inheight_max;
    %                   end
                        % Find the largest pblh
                        if length(aa) > 1 && length(bb) == 1
                            [m,p] = max([pblh_var(aa(1),bb,tt),pblh_var(aa(2),bb,tt)]);
                            pblh_tmp = pblh_var(aa(p),bb,tt);
                        elseif length(aa) == 1 && length(bb) > 1
                            [m,p] = max([pblh_var(aa,bb(1),tt),pblh_var(aa,bb(2),tt)]);
                            pblh_tmp = pblh_var(aa,bb(p),tt);
                        elseif length(aa) > 1 && length(bb) > 1
                            pblh_values = [pblh_var(aa(1),bb(1),tt),pblh_var(aa(1),bb(2),tt),pblh_var(aa(2),bb(1),tt),pblh_var(aa(2),bb(2),tt)];
                            [m,p] = max(pblh_values);
                            pblh_tmp = pblh_values(p);
                        else
                            pblh_tmp = pblh_var(aa,bb,tt);
                        end

                        inum = inum + 1;
                        inheight_all(inum) = inheight_tmp;
                        pblh_all(inum) = pblh_tmp;
                        gfed_all(inum) = gfed_var(ii,jj,tt,ispc);

                        if pblh_tmp < inheight_tmp
                            cc = find( abs(lon_gfed_au(ii) - lon_prof) == min(abs(lon_gfed_au(ii) - lon_prof)) );
                            dd = find( abs(lat_gfed_au(jj) - lat_prof) == min(abs(lat_gfed_au(jj) - lat_prof)) );
                            tt_days = tt + tt_indx1 - 1;
                            for imon = 1:12
                                if sum(ndays_mon(1:imon)) >= tt_days
                                    month_indx = imon;
                                    break;
                                end
                            end
                            prof_tmp     = squeeze(fire_prof(cc(1),dd(1),:,month_indx));              
                            zz = find( alt_prof > pblh_tmp, 1);
                            gfed_var_new(ii,jj,tt,ispc) = gfed_var(ii,jj,tt,ispc) * ( 1. - sum(prof_tmp(zz:end)) );
                            ratio_all_1(ii,jj,tt) = 1. - sum(prof_tmp(zz:end));
                         end
                    end
                end
            end
        end
    end
lon_indx1 = find( abs(min(lon) - lon_global_gfed) == min( abs(min(lon) - lon_global_gfed) ));
lon_indx2 = find( abs(max(lon) - lon_global_gfed) == min( abs(max(lon) - lon_global_gfed) ));
lat_indx1 = find( abs(min(lat) - lat_global_gfed) == min( abs(min(lat) - lat_global_gfed) ));
lat_indx2 = find( abs(max(lat) - lat_global_gfed) == min( abs(max(lat) - lat_global_gfed) ));
aemit_daily_1(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),tt_indx1:tt_indx2,:) = gfed_var_new;
end    
lon_global = lon_global_gfed;
lat_global = lat_global_gfed;
%
disp('Step 2 - emis 2');
clearvars gfed_all inheight_all pblh_all
% Turn on vertical scale
is_injection_height_scale = 0;
%
disp('Step 2 - emis 3');
clearvars gfed_all inheight_all pblh_all
% Turn on vertical scale
is_injection_height_scale = 1;

inum = 0;
if is_injection_height_scale
    % combine the injection height in the specified year
    filename_gfas = ['E:\work_for_wildfires\injecton_height\',num2str(years),num2str(1),'_gfas_injection_height.nc']; % unit:m
    if years ~= 2020
        lon_global_gfas = ncread(filename_gfas,'longitude');
        lat_global_gfas = flip(ncread(filename_gfas,'latitude'));
    else
      lon_global_gfas = ncread(filename_gfas,'lon'); % for 2020
      lat_global_gfas = ncread(filename_gfas,'lat'); % for 2020
    end

    for imon = 1%1:12
        filename_gfas = ['E:\work_for_wildfires\injecton_height\',num2str(years),num2str(imon),'_gfas_injection_height.nc'];
        injh_tmp = ncread(filename_gfas,'mami');
        injh_tmp = flip(injh_tmp,2); % reverse latitude
        
        lon_indx1 = find( abs(min(lon) - lon_global_gfas) == min( abs(min(lon) - lon_global_gfas) ));
        lon_indx2 = find( abs(max(lon) - lon_global_gfas) == min( abs(max(lon) - lon_global_gfas) ));
        lat_indx1 = find( abs(min(lat) - lat_global_gfas) == min( abs(min(lat) - lat_global_gfas) ));
        lat_indx2 = find( abs(max(lat) - lat_global_gfas) == min( abs(max(lat) - lat_global_gfas) ));
        injh_tmp  = injh_tmp(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);
        
        if imon == 1
            injh_year = injh_tmp;
        else
            injh_year = cat(3,injh_year,injh_tmp);
        end
        clearvars injh_tmp
    end
    lon_gfas_au = lon_global_gfas(lon_indx1(1):lon_indx2(1));
    lat_gfas_au = lat_global_gfas(lat_indx1(1):lat_indx2(1));
    inheight_var = injh_year(:,:,tt_indx1:tt_indx2);
    
    % combine the pblh in the specified year
    for imon = 1%1:12
        filename = ['E:\work_for_wildfires\injecton_height\MERRA2_PBLH\merra2_',num2str(years),num2str(imon),'_pblh_au.mat'];
        load(filename);
        load('E:\work_for_wildfires\injecton_height\MERRA2_PBLH\lonlat_au_merra2.mat');
        if imon == 1
            pblh_year = var_au;
        else
            pblh_year = cat(3,pblh_year,var_au);
        end
        
        pblh_var = zeros(size(pblh_year,1),size(pblh_year,2),size(pblh_year,3)/24);
        for tt = 1:size(pblh_year,3)/24
            pblh_var(:,:,tt) = mean(pblh_year(:,:,(tt-1)*24+1:24*tt),3);
        end
        lon_pblh_au = lon_au;
        lat_pblh_au = lat_au;
    end
    pblh_var = pblh_var(:,:,tt_indx1:tt_indx2);
    clearvars var_au
    
    % load fire injection height profile GMD or IS4FIRE
    filename  = 'E:\work_for_wildfires\injecton_height\emis_mean_normalised_wholeday_smooth.nc';
    lon_prof  = ncread(filename,'longitude');
    lat_prof  = ncread(filename,'latitude');
    alt_prof  = ncread(filename,'z'); % unit: m
    fire_prof = ncread(filename,'ems_norm_wholeday_smooth');
    
    % compare the injection height and pblh and scale the biomass burning emissions
    gfed_var_new = gfed_var;
    ratio_all_3 = zeros(size(gfed_var,1),size(gfed_var,2),size(gfed_var,3));
    
    for ispc = 1:2
        for tt = 1:size(gfed_var,3)
            for jj = 1:size(gfed_var,2)
                for ii = 1:size(gfed_var,1)
                    if gfed_var(ii,jj,tt,ispc) > 0
                        aa = find( abs(lon_gfed_au(ii) - lon_pblh_au) == min(abs(lon_gfed_au(ii) - lon_pblh_au)) );
                        bb = find( abs(lat_gfed_au(jj) - lat_pblh_au) == min(abs(lat_gfed_au(jj) - lat_pblh_au)) );
                        xx = find( abs(lon_gfed_au(ii) - lon_gfas_au) < 1);
                        yy = find( abs(lat_gfed_au(jj) - lat_gfas_au) < 1);
                        % Find the largest inheight
                        injh_area = squeeze(inheight_var(min(xx):max(xx),min(yy):max(yy),tt));
                        inheight_tmp = max(max(injh_area));
                        
                        % Fix non-injection height grid
                        inheight_max = 1e4;
                        if inheight_tmp < 1e-5
                             inheight_tmp = inheight_max;
                        end
                        
                        % Find the largest pblh
                        if length(aa) > 1 && length(bb) == 1
                            [m,p] = max([pblh_var(aa(1),bb,tt),pblh_var(aa(2),bb,tt)]);
                            pblh_tmp = pblh_var(aa(p),bb,tt);
                        elseif length(aa) == 1 && length(bb) > 1
                            [m,p] = max([pblh_var(aa,bb(1),tt),pblh_var(aa,bb(2),tt)]);
                            pblh_tmp = pblh_var(aa,bb(p),tt);
                        elseif length(aa) > 1 && length(bb) > 1
                            pblh_values = [pblh_var(aa(1),bb(1),tt),pblh_var(aa(1),bb(2),tt),pblh_var(aa(2),bb(1),tt),pblh_var(aa(2),bb(2),tt)];
                            [m,p] = max(pblh_values);
                            pblh_tmp = pblh_values(p);
                        else
                            pblh_tmp = pblh_var(aa,bb,tt);
                        end

                        inum = inum + 1;
                        inheight_all(inum) = inheight_tmp;
                        pblh_all(inum) = pblh_tmp;
                        gfed_all(inum) = gfed_var(ii,jj,tt,ispc);

                        if pblh_tmp < inheight_tmp
                            cc = find( abs(lon_gfed_au(ii) - lon_prof) == min(abs(lon_gfed_au(ii) - lon_prof)) );
                            dd = find( abs(lat_gfed_au(jj) - lat_prof) == min(abs(lat_gfed_au(jj) - lat_prof)) );
                            tt_days = tt + tt_indx1 - 1;
                            for imon = 1:12
                                if sum(ndays_mon(1:imon)) >= tt_days
                                    month_indx = imon;
                                    break;
                                end
                            end
                            prof_tmp     = squeeze(fire_prof(cc(1),dd(1),:,month_indx));              
                            zz = find( alt_prof > pblh_tmp, 1);
                            gfed_var_new(ii,jj,tt,ispc) = gfed_var(ii,jj,tt,ispc) * ( 1. - sum(prof_tmp(zz:end)) );
                            ratio_all_3(ii,jj,tt) = 1. - sum(prof_tmp(zz:end));
                         end
                    end
                end
            end
        end
    end
lon_indx1 = find( abs(min(lon) - lon_global_gfed) == min( abs(min(lon) - lon_global_gfed) ));
lon_indx2 = find( abs(max(lon) - lon_global_gfed) == min( abs(max(lon) - lon_global_gfed) ));
lat_indx1 = find( abs(min(lat) - lat_global_gfed) == min( abs(min(lat) - lat_global_gfed) ));
lat_indx2 = find( abs(max(lat) - lat_global_gfed) == min( abs(max(lat) - lat_global_gfed) ));
aemit_daily_3(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),tt_indx1:tt_indx2,:) = gfed_var_new;
end

clearvars gfed_var_new gfed_var inheight_* injh_* 

%% Step 3: Read and calculate the daily footprints
% STILT output files
disp('Step 3');
filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'01_gfs_nohnf_receptor_se\'];
diroutput = dir(fullfile(filepath,'*_foot.nc'));
filename = {diroutput.name};

filename_1 = char(filename(1));
lon = ncread([filepath,filename_1],'lon');
lat = ncread([filepath,filename_1],'lat');
xdim = length(lon);
ydim = length(lat);
    
smoke_exp1    = zeros(xdim,ydim,rundays,receptors,2);
smoke_exp2    = zeros(xdim,ydim,rundays,receptors,2);
smoke_exp3    = zeros(xdim,ydim,rundays,receptors,2);
% loop for each file
for ii = 1:length(filename)
    filename_tmp = char(filename(ii));
    location_tmp = filename_tmp(14:end-10);
    disp(filename_tmp);
    % find which receptors
    for rr = 1:receptors
        if strcmp(char(loc_str(rr)),location_tmp) 
            irec = rr;
            break
        end
    end
    % find which days
    file_year = str2double(filename_tmp(1:4));
    file_mon  = str2double(filename_tmp(5:6));
    file_day  = str2double(filename_tmp(7:8));
    iday      = datenum(file_year,file_mon,file_day) - datenum(years,start_mon,start_date) + 1;

    % for the filename output
    for aa = 1:length(location_tmp)
        if location_tmp(aa) == '_'
            location_tmp(aa) = ',';
        end
    end

    foot = ncread([filepath,filename_tmp],'foot');
    time = ncread([filepath,filename_tmp],'time');

    % convert the format of time
    t0 = datenum(1970,01,01,00,00,00);
    t1 = zeros(length(time),1);
    for tt = 1:length(time)
        t1(tt) = addtodate(t0,time(tt),'second');
    end
    time_str = datestr(t1,'yyyy-mm-dd HH:MM:SS');
    time_vec = datevec(time_str);

    dd = datenum(time_vec(end,1:3)) - datenum(time_vec(end,1),1,1) + 1;

    % calculate the daily backward footprints
    output_days = floor(length(time)/24);
    final_days = mod(length(time),24);

    if final_days == 0
       footprint_daily = zeros(xdim,ydim,output_days);
       for tt = output_days:-1:1
           footprint_daily(:,:,tt) = sum(foot(:,:,(tt-1)*24+1:24*tt),3);
       end
       tt_indx1 = dd-output_days + 1;
       tt_indx2 = dd;
    elseif final_days > 0
        footprint_daily = zeros(xdim,ydim,output_days+1);
        footprint_daily(:,:,1) = sum(foot(:,:,1:final_days),3);
        for tt = 2:output_days + 1
            footprint_daily(:,:,tt) = sum(foot(:,:,final_days + (tt-2)*24+1:final_days + 24*(tt-1)),3);
        end
        tt_indx1 = dd - output_days;
        tt_indx2 = dd;
    end

    % find the domian range in GFED emissions
    lon_indx1 = find( abs(min(lon) - lon_global) == min( abs(min(lon) - lon_global) ));
    lon_indx2 = find( abs(max(lon) - lon_global) == min( abs(max(lon) - lon_global) ));
    lat_indx1 = find( abs(min(lat) - lat_global) == min( abs(min(lat) - lat_global) ));
    lat_indx2 = find( abs(max(lat) - lat_global) == min( abs(max(lat) - lat_global) ));
    
    if tt_indx1 > 0
        for ispc = 1:2
        % read the fire emissions & exchange the dimensions [unit: umol/m2/s]
        fireemis_var1 = aemit_daily_1(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),tt_indx1:tt_indx2,ispc);
        fireemis_var2 = aemit_daily_2(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),tt_indx1:tt_indx2,ispc);
        fireemis_var3 = aemit_daily_3(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),tt_indx1:tt_indx2,ispc);

        % unit conversion for trace gases and aerosols
        % For trace gases: kg m^-2 s^-1 --> umol m^-2 s^-1
        % mw_spc         = 28;
        kg2g             = 1e3;
        % convert_factor_gas = kg2g/mw_spc*1e6;

        % For aerosols: ppm/(kg m^-2 s^-1) --> ug m^-3/(umol m^-2 s^-1)
        mw_air         = 29;
        rho_air        = 1.29; % dry air density at STP [kg m^-3]
        convert_factor = kg2g * 1e6 * rho_air/mw_air * 1e3;
        smoke_exp_tmp1 = fireemis_var1 .* footprint_daily * convert_factor;
        smoke_exp_tmp2 = fireemis_var2 .* footprint_daily * convert_factor;
        smoke_exp_tmp3 = fireemis_var3 .* footprint_daily * convert_factor;

        smoke_exp1(:,:,iday,irec,ispc) = nansum(smoke_exp_tmp1,3);
        smoke_exp2(:,:,iday,irec,ispc) = nansum(smoke_exp_tmp2,3);
        smoke_exp3(:,:,iday,irec,ispc) = nansum(smoke_exp_tmp3,3);
        end
    else
        for ispc = 1:2
        % read the fire emissions & exchange the dimensions [unit: umol/m2/s]
        fireemis_var1b = aemit_daily_1b(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),31+tt_indx1:31,ispc);
        fireemis_var2b = aemit_daily_2b(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),31+tt_indx1:31,ispc);
        fireemis_var3b = aemit_daily_3b(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),31+tt_indx1:31,ispc);
        
        fireemis_var1 = aemit_daily_1(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),1:tt_indx2,ispc);
        fireemis_var2 = aemit_daily_2(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),1:tt_indx2,ispc);
        fireemis_var3 = aemit_daily_3(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),1:tt_indx2,ispc);

        % unit conversion for trace gases and aerosols
        % For trace gases: kg m^-2 s^-1 --> umol m^-2 s^-1
        % mw_spc         = 28;
        kg2g             = 1e3;
        % convert_factor_gas = kg2g/mw_spc*1e6;

        % For aerosols: ppm/(kg m^-2 s^-1) --> ug m^-3/(umol m^-2 s^-1)
        mw_air         = 29;
        rho_air        = 1.29; % dry air density at STP [kg m^-3]
        convert_factor = kg2g * 1e6 * rho_air/mw_air * 1e3;
        smoke_exp_tmp1(:,:,1:size(fireemis_var1b,3))     = fireemis_var1b .* footprint_daily(:,:,1:size(fireemis_var1b,3)) * convert_factor;
        smoke_exp_tmp1(:,:,size(fireemis_var1b,3)+1:size(footprint_daily,3)) = fireemis_var1 .* footprint_daily(:,:,size(fireemis_var1b,3)+1:end) * convert_factor;
        smoke_exp_tmp2(:,:,1:size(fireemis_var1b,3))     = fireemis_var2b .* footprint_daily(:,:,1:size(fireemis_var1b,3)) * convert_factor;
        smoke_exp_tmp2(:,:,size(fireemis_var1b,3)+1:size(footprint_daily,3)) = fireemis_var2 .* footprint_daily(:,:,size(fireemis_var1b,3)+1:end) * convert_factor;
        smoke_exp_tmp3(:,:,1:size(fireemis_var1b,3))     = fireemis_var3b .* footprint_daily(:,:,1:size(fireemis_var1b,3)) * convert_factor;
        smoke_exp_tmp3(:,:,size(fireemis_var1b,3)+1:size(footprint_daily,3)) = fireemis_var3 .* footprint_daily(:,:,size(fireemis_var1b,3)+1:end) * convert_factor;

        smoke_exp1(:,:,iday,irec,ispc) = nansum(smoke_exp_tmp1,3);
        smoke_exp2(:,:,iday,irec,ispc) = nansum(smoke_exp_tmp2,3);
        smoke_exp3(:,:,iday,irec,ispc) = nansum(smoke_exp_tmp3,3);
        end
    end
end
smoke_exp = squeeze(smoke_exp1(:,:,:,:,1));
filename_output1 = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'01_gfs_nohnf_receptor_se\smoke_exp_output\smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat'];
save(filename_output1,'smoke_exp');

clearvars smoke_exp
smoke_exp = squeeze(smoke_exp1(:,:,:,:,2));
filename_output1 = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'01_gfs_nohnf_receptor_se\smoke_exp_output\smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat'];
save(filename_output1,'smoke_exp');

clearvars smoke_exp
smoke_exp = squeeze(smoke_exp2(:,:,:,:,1));
filename_output2 = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'01_gfs_nohnf_receptor_se\smoke_exp_output\smoke_exp_bc_gfed_noscale.mat'];
save(filename_output2,'smoke_exp');

clearvars smoke_exp
smoke_exp = squeeze(smoke_exp2(:,:,:,:,2));
filename_output2 = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'01_gfs_nohnf_receptor_se\smoke_exp_output\smoke_exp_oc_gfed_noscale.mat'];
save(filename_output2,'smoke_exp');

clearvars smoke_exp
smoke_exp = squeeze(smoke_exp3(:,:,:,:,1));
filename_output3 = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'01_gfs_nohnf_receptor_se\smoke_exp_output\smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_fixnon_injh_mami.mat'];
save(filename_output3,'smoke_exp');

clearvars smoke_exp
smoke_exp = squeeze(smoke_exp3(:,:,:,:,2));
filename_output3 = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'01_gfs_nohnf_receptor_se\smoke_exp_output\smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_fixnon_injh_mami.mat'];
save(filename_output3,'smoke_exp');
end
