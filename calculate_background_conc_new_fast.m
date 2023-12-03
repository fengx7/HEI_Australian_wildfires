%% Calcualte the background concentrations for PM2.5
% analyze the time series of fire emissions for each fire season, summed
% over all grids in the burning regions
clear;
location = 'E:\work_for_wildfires\biomass_burning_emis\GFED\';
prefix = 'GFEDv4s_';
% prefix = 'QFEDv2p5r1_';
% prefix = 'FINNv1p5_';
postfix = '.nc';

spc_name   = {'BC','CO','OC','PM25'};
dirfiles   = dir([location,prefix,'*']);
filename1  = dirfiles.name;
lon_global = ncread([location,filename1],'lon');
lat_global = ncread([location,filename1],'lat');
nx         = length(lon_global);
ny         = length(lat_global);
% pm25base_years = zeros(9,1);
site_name = {'palmerston','westmackay','','southgladstone'};
% site_name = {'footscray','albury','florey','wollongong','prospect','newcastle','mountaincreek','springwood','alphington','liverpool','wallsend'};

for iyear = 2009:2020
for ispc = 3
    for years = iyear
        % leap year
        if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
            ndays = 366;
            ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
        else
            ndays = 365;
            ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
        end
 
        emis_daily = zeros(nx,ny,ndays,'single');
        
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
            emis_daily(:,:,startday:endday) = ncread(filename,char(spc_name(ispc)));
        end
    end
end

%% select regions where mean footprints greater than 1e-4
filepath = ['F:\stilt_output\2019\footprints_2019_gfs_nohnf_receptor_north\'];

prefix_foot = [num2str(2019),'*'];
diroutput = dir(fullfile(filepath,prefix_foot)); % *_foot.nc 
filename = {diroutput.name};
filename_1 = char(filename(1));
lon = ncread([filepath,filename_1],'lon');
lat = ncread([filepath,filename_1],'lat');
xdim = length(lon);
ydim = length(lat);

filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_north\'];
load([filepath,'footprints_daily.mat']);
mean_footprint = squeeze(nanmean(footprints_exp,3));
footprint_base = 1e-4;
irec           = 4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change
site_name_tmp  = char(site_name{irec});

is_region = zeros(size(mean_footprint,1),size(mean_footprint,2));
for ii = 1:size(mean_footprint,1)
    for jj = 1:size(mean_footprint,2)
        if mean_footprint(ii,jj,irec) > footprint_base
            is_region(ii,jj) = 1;
        end
    end
end

lon_indx1 = find( abs(min(lon) - lon_global) == min(abs(min(lon) - lon_global)));
lon_indx2 = find( abs(max(lon) - lon_global) == min(abs(max(lon) - lon_global)));
lat_indx1 = find( abs(min(lat) - lat_global) == min(abs(min(lat) - lat_global)));
lat_indx2 = find( abs(max(lat) - lat_global) == min(abs(max(lat) - lat_global)));
emis_daily_au = emis_daily(lon_indx1:lon_indx2,lat_indx1:lat_indx2,:);

for tt = 1:size(emis_daily_au,3)
    emis_daily_au(:,:,tt) = emis_daily_au(:,:,tt) .* is_region;
end
%% select time period and calculate the time series
startday   = 91; % 213; 91
endday     = ndays;

ts_var     = squeeze(sum(sum(emis_daily_au(:,:,startday:endday),2),1));

max_ts     = max(ts_var);
min_ts     = min(ts_var);
edges      = min_ts:(max_ts-min_ts)/1000:max_ts;

[N,edges]  = histcounts(ts_var,edges);
% [N,edges1]  = histcounts(ts_var,1000);

percent_N  = N/length(ts_var);

for ii = 1:length(percent_N)
    sum_percent = sum(percent_N(1:ii));
    if sum_percent > 0.20  % 10% percentile
        index_edge = ii + 1;
        break;
    end
end

emis_threshold = edges(index_edge);

% Mark the Mth day if emissions < emis_threshold; M = 1,2,3,4,5

mark_day = zeros(length(ts_var),1);
% M = 1;
for ii = 1:length(ts_var)
    if ts_var(ii) < emis_threshold
        mark_day(ii) = 1;
    end
end

for ii = 1:length(mark_day)
    if ii == 1
        if mark_day(ii) + mark_day(ii+1) < 2
            mark_day(ii) = 0;
        end
    elseif ii == length(mark_day)
        if mark_day(ii-1) + mark_day(ii) < 2
            mark_day(ii) = 0;
        end
    else
        if mark_day(ii-1) == 0 && mark_day(ii+1) == 0
            mark_day(ii) = 0;
        end
    end
end
%% load observation PM2.5 concentrations
load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
eval(['pm25obs = pm25_daily_',site_name_tmp,'(startday:endday);']);

i = 0;
clearvars pm25baseline
for ii = 1:length(mark_day)
    if mark_day(ii) == 1
      i = i + 1;
      pm25baseline(i) = pm25obs(ii);
    end
end

pm25bg = nanmean(pm25baseline);
pm25base_years(iyear-2008) = pm25bg;

% figure1 = figure;
% plot(pm25_daily_palmerston(91:end));
% plot(pm25obs);
% hold on;
% bar(mark_day*100);
end
filename = ['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat'];
save(filename,'pm25base_years');
