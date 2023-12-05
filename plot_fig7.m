%% Darwin
receptors = 4;
is_omoc   = 0;
loc_lon = [130.94853,146.8257,149.1549,151.2704];
loc_lat = [-12.50779,-19.2542,-21.1595,-23.8627];
site_name = {'palmerston','coastguard','westmackay','southgladstone'};

% Seasonal OM/OC ratios
omoc_site = zeros(receptors,4);
for xx = 1:size(omoc_site,1)
    for yy = 1:size(omoc_site,2)
        omoc_site(xx,yy) = 2.1;
    end
end

yr_num  = 1;
sim_rf  = zeros(274,yr_num,2); % INJ-RF
sim_cl  = zeros(274,yr_num,2); % INJ-CLIM
sim_no  = zeros(274,yr_num,2); % INJ-CTL
obs_all = zeros(274,yr_num,2); 

for iyear = 2019
years = iyear;
% leap year?
if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
    ndays = 366;
    ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
    
else
    ndays = 365;
    ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
end

for irec  = 1
site_name_tmp = char(site_name{irec});
load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
eval(['obs_dy = pm25_daily_',site_name_tmp,';']);

% Remove the day with missing observations more than 8
for ii = 1:length(obs_hr)/24
    obs_tmp = obs_hr((ii-1)*24+1:ii*24);
    nan_num = length(find(isnan(obs_tmp)));
    if nan_num > 8
        obs_dy(ii) = NaN;
    end
end

% Read smoke exposure data
% OC
filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_north\smoke_exp_output\'];
load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
ts_oc1 = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_oc2 = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
ts_oc3 = squeeze(nansum(nansum(smoke_exp,1),2));

% BC
load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
ts_bc1 = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_bc2 = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
ts_bc3 = squeeze(nansum(nansum(smoke_exp,1),2));

% Read background PM2.5
load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
pm25bg = pm25base_years(iyear-2008);

% figure
    ii = irec; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change receptors
    
    % Calculate smoke PM2.5 (OC* OM/OC + BC); We can use seasonal OM/OC ratios 
    ts_1(1:61)    = ts_oc1(1:61,ii)*omoc_site(ii,1) + ts_bc1(1:61,ii);
    ts_1(62:153)  = ts_oc1(62:153,ii)*omoc_site(ii,2) + ts_bc1(62:153,ii);
    ts_1(154:245) = ts_oc1(154:245,ii)*omoc_site(ii,3) + ts_bc1(154:245,ii);
    ts_1(246:275) = ts_oc1(246:275,ii)*omoc_site(ii,4) + ts_bc1(246:275,ii);
    ts_2(1:61)    = ts_oc2(1:61,ii)*omoc_site(ii,1) + ts_bc2(1:61,ii);
    ts_2(62:153)  = ts_oc2(62:153,ii)*omoc_site(ii,2) + ts_bc2(62:153,ii);
    ts_2(154:245) = ts_oc2(154:245,ii)*omoc_site(ii,3) + ts_bc2(154:245,ii);
    ts_2(246:275) = ts_oc2(246:275,ii)*omoc_site(ii,4) + ts_bc2(246:275,ii);
    ts_3(1:61)    = ts_oc3(1:61,ii)*omoc_site(ii,1) + ts_bc3(1:61,ii);
    ts_3(62:153)  = ts_oc3(62:153,ii)*omoc_site(ii,2) + ts_bc3(62:153,ii);
    ts_3(154:245) = ts_oc3(154:245,ii)*omoc_site(ii,3) + ts_bc3(154:245,ii);
    ts_3(246:275) = ts_oc3(246:275,ii)*omoc_site(ii,4) + ts_bc3(246:275,ii);

    start_day  = 1;
    end_day    = 274;
    t          = 1:1:ndays;

    % 10-day moving averaging
    mv_day     = 10;
    sim1 = movmean(ts_1(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim2 = movmean(ts_2(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim3 = movmean(ts_3(start_day:end_day)+pm25bg,mv_day,'omitnan');
    obs  = movmean(obs_dy(ndays-end_day+1:end),mv_day,'omitnan');
    
    yearbase = 2018;
    sim_rf(:,iyear-yearbase,irec) = sim1;
    sim_cl(:,iyear-yearbase,irec) = sim2;
    sim_no(:,iyear-yearbase,irec) = sim3;
    obs_all(:,iyear-yearbase,irec)= obs;
end

end
% Some statistics
% Mean value, NMB, R, and RMSE
mean_sim_rf = nanmean(sim_rf(:,:,irec),2);
mean_sim_cl = nanmean(sim_cl(:,:,irec),2);
mean_sim_no = nanmean(sim_no(:,:,irec),2);
mean_obs_all = nanmean(obs_all(:,:,irec),2);

var_all1 = [mean_sim_rf,mean_obs_all];
var_all1 = rmmissing(var_all1);
nmb1      = (nanmean(mean_sim_rf) - nanmean(mean_obs_all))/nanmean(mean_obs_all);
rr1       = corrcoef(var_all1(:,1),var_all1(:,2));
rmse1     = sqrt(mean((var_all1(:,1) - var_all1(:,2)).^2));
mean1     = mean(mean_sim_rf);

var_all2 = [mean_sim_cl,mean_obs_all];
var_all2 = rmmissing(var_all2);
nmb2      = (nanmean(mean_sim_cl) - nanmean(mean_obs_all))/nanmean(mean_obs_all);
rr2       = corrcoef(var_all2(:,1),var_all2(:,2));
rmse2     = sqrt(mean((var_all2(:,1) - var_all2(:,2)).^2));
mean2     = mean(mean_sim_cl);

var_all3 = [mean_sim_no,mean_obs_all];
var_all3 = rmmissing(var_all3);
nmb3      = (nanmean(mean_sim_no) - nanmean(mean_obs_all))/nanmean(mean_obs_all);
rr3       = corrcoef(var_all3(:,1),var_all3(:,2));
rmse3     = sqrt(mean((var_all3(:,1) - var_all3(:,2)).^2));
mean3     = mean(mean_sim_no);

std_rf    = std(squeeze(sim_rf(:,:,irec)),1,2);
std_cl    = std(squeeze(sim_cl(:,:,irec)),1,2);
std_no    = std(squeeze(sim_no(:,:,irec)),1,2);
std_all    = std(squeeze(obs_all(:,:,irec)),1,2);
mean_obs   = nanmean(mean_obs_all);

%%
figure1 = figure;
subplot1 = subplot(3,2,1,'BoxStyle','full');
hold(subplot1,'on');

plot1 = plot(t(ndays-end_day+1:end),mean_sim_rf','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot1);
set(plot1,'DisplayName','INJ-RF{    }',...
    'Color',[0.851 0.325490206480026 0.0980392172932625]);

plot2 = plot(t(ndays-end_day+1:end),mean_sim_cl','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot1);
set(plot2,'DisplayName','INJ-CLIM{    }',... 
    'Color',[0 0.498039215803146 0]);

plot3 = plot(t(ndays-end_day+1:end),mean_sim_no','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot1);
set(plot3,'DisplayName','CTL{    }','Color',[0 0.447058826684952 0.74117648601532]);

plot5 = plot(t(ndays-end_day+1:end),mean_obs_all','LineWidth',2.5,'Parent',subplot1);
set(plot5,'DisplayName','OBS{    }','Color',[0 0 0]);

% ylabel
% ylabel('PM_2_._5 concentrations (\mug m^-^3)');

% title
title('(a) Darwin (2019)');

xlim(subplot1,[91 365]);
ylim(subplot1,[0 60]);

box(subplot1,'on');
set(subplot1,'FontSize',20,'XTick',[1 32 60 91 121 152 182 213 244 274 305 335 365],...
    'XTickLabel',...
    {'Jan','Feb','Mar','Apr','May','Jun','July','Aug','Sep','Oct','Nov','Dec',''});
% legend
legend1 = legend(subplot1,'show','Location','NorthEast');
set(legend1,'EdgeColor',[0 0 0]);
set(gca,'Fontsize',24);
%% Gladstone
yr_num  = 1;
sim_rf  = zeros(274,yr_num,2);
sim_cl  = zeros(274,yr_num,2);
sim_no  = zeros(274,yr_num,2);
obs_all = zeros(274,yr_num,2);

for iyear = 2019
years = iyear;
% leap year
if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
    ndays = 366;
    ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
    
else
    ndays = 365;
    ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
end

for irec  = 4
site_name_tmp = char(site_name{irec});
load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
eval(['obs_dy = pm25_daily_',site_name_tmp,';']);

for ii = 1:length(obs_hr)/24
    obs_tmp = obs_hr((ii-1)*24+1:ii*24);
    nan_num = length(find(isnan(obs_tmp)));
    if nan_num > 8
        obs_dy(ii) = NaN;
    end
end

filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_north\smoke_exp_output\'];

load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
ts_oc1 = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_oc2 = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
ts_oc3 = squeeze(nansum(nansum(smoke_exp,1),2));

load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
ts_bc1 = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_bc2 = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
ts_bc3 = squeeze(nansum(nansum(smoke_exp,1),2));

load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
pm25bg = pm25base_years(iyear-2008);


% figure
    ii = irec; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change receptors
    
    ts_1(1:61)    = ts_oc1(1:61,ii)*omoc_site(ii,1) + ts_bc1(1:61,ii);
    ts_1(62:153)  = ts_oc1(62:153,ii)*omoc_site(ii,2) + ts_bc1(62:153,ii);
    ts_1(154:245) = ts_oc1(154:245,ii)*omoc_site(ii,3) + ts_bc1(154:245,ii);
    ts_1(246:275) = ts_oc1(246:275,ii)*omoc_site(ii,4) + ts_bc1(246:275,ii);
    ts_2(1:61)    = ts_oc2(1:61,ii)*omoc_site(ii,1) + ts_bc2(1:61,ii);
    ts_2(62:153)  = ts_oc2(62:153,ii)*omoc_site(ii,2) + ts_bc2(62:153,ii);
    ts_2(154:245) = ts_oc2(154:245,ii)*omoc_site(ii,3) + ts_bc2(154:245,ii);
    ts_2(246:275) = ts_oc2(246:275,ii)*omoc_site(ii,4) + ts_bc2(246:275,ii);
    ts_3(1:61)    = ts_oc3(1:61,ii)*omoc_site(ii,1) + ts_bc3(1:61,ii);
    ts_3(62:153)  = ts_oc3(62:153,ii)*omoc_site(ii,2) + ts_bc3(62:153,ii);
    ts_3(154:245) = ts_oc3(154:245,ii)*omoc_site(ii,3) + ts_bc3(154:245,ii);
    ts_3(246:275) = ts_oc3(246:275,ii)*omoc_site(ii,4) + ts_bc3(246:275,ii);

    start_day  = 1;
    end_day    = 274;
    t          = 1:1:ndays;

    sim1 = movmean(ts_1(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim2 = movmean(ts_2(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim3 = movmean(ts_3(start_day:end_day)+pm25bg,mv_day,'omitnan');
    obs  = movmean(obs_dy(ndays-end_day+1:end),mv_day,'omitnan');
     
    yearbase = 2018;
    sim_rf(:,iyear-yearbase,irec) = sim1;
    sim_cl(:,iyear-yearbase,irec) = sim2;
    sim_no(:,iyear-yearbase,irec) = sim3;
    obs_all(:,iyear-yearbase,irec)= obs;  
end

end
%
mean_sim_rf = nanmean(sim_rf(:,:,irec),2);
mean_sim_cl = nanmean(sim_cl(:,:,irec),2);
mean_sim_no = nanmean(sim_no(:,:,irec),2);
mean_obs_all = nanmean(obs_all(:,:,irec),2);

var_all1 = [mean_sim_rf,mean_obs_all];
var_all1 = rmmissing(var_all1);
nmb1      = (nanmean(mean_sim_rf) - nanmean(mean_obs_all))/nanmean(mean_obs_all);
rr1       = corrcoef(var_all1(:,1),var_all1(:,2));
rmse1     = sqrt(mean((var_all1(:,1) - var_all1(:,2)).^2));
mean1     = mean(mean_sim_rf);

var_all2 = [mean_sim_cl,mean_obs_all];
var_all2 = rmmissing(var_all2);
nmb2      = (nanmean(mean_sim_cl) - nanmean(mean_obs_all))/nanmean(mean_obs_all);
rr2       = corrcoef(var_all2(:,1),var_all2(:,2));
rmse2     = sqrt(mean((var_all2(:,1) - var_all2(:,2)).^2));
mean2     = mean(mean_sim_cl);

var_all3 = [mean_sim_no,mean_obs_all];
var_all3 = rmmissing(var_all3);
nmb3      = (nanmean(mean_sim_no) - nanmean(mean_obs_all))/nanmean(mean_obs_all);
rr3       = corrcoef(var_all3(:,1),var_all3(:,2));
rmse3     = sqrt(mean((var_all3(:,1) - var_all3(:,2)).^2));
mean3     = mean(mean_sim_no);

std_rf    = std(squeeze(sim_rf(:,:,irec)),1,2);
std_cl    = std(squeeze(sim_cl(:,:,irec)),1,2);
std_no    = std(squeeze(sim_no(:,:,irec)),1,2);
std_all    = std(squeeze(obs_all(:,:,irec)),1,2);
mean_obs   = nanmean(mean_obs_all);

%%
subplot2 = subplot(3,2,2,'BoxStyle','full');
hold(subplot2,'on');

plot1 = plot(t(ndays-end_day+1:end),mean_sim_rf','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot2);
set(plot1,'DisplayName','INJ-RF',...
    'Color',[0.851 0.325490206480026 0.0980392172932625]);

plot2 = plot(t(ndays-end_day+1:end),mean_sim_cl','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot2);
set(plot2,'DisplayName','INJ-CLIM',...
    'Color',[0 0.498039215803146 0]);

plot3 = plot(t(ndays-end_day+1:end),mean_sim_no','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot2);
set(plot3,'DisplayName','CTL','Color',[0 0.447058826684952 0.74117648601532]);

plot5 = plot(t(ndays-end_day+1:end),mean_obs_all','LineWidth',2.5,'Parent',subplot2);
set(plot5,'DisplayName','OBS','Color',[0 0 0]);

% ylabel
% ylabel('PM_2_._5 concentrations (\mug m^-^3)');

% title
title('(b) Gladstone (2019) ');

xlim(subplot2,[91 365]);
ylim(subplot2,[0 40]);
box(subplot2,'on');

set(subplot2,'FontSize',20,'XTick',[1 32 60 91 121 152 182 213 244 274 305 335 365],...
    'XTickLabel',...
    {'Jan','Feb','Mar','Apr','May','Jun','July','Aug','Sep','Oct','Nov','Dec',''});
set(gca,'Fontsize',24);
%% Newcastle
receptors = 11;
loc_lon = [144.8728027,146.93986,149.043539,150.88733,150.91417,151.75965,153.1038,153.1356,145.0306,150.9058,151.6692];
loc_lat = [-37.80487823,-36.05182,-35.220606,-34.41706,-33.79424,-32.9312,-26.6917,-27.6125,-37.7784,-33.9328,-32.8961];
site_name = {'footscray','albury','florey','wollongong','prospect','newcastle','mountaincreek;','springwood','alphington','liverpool','wallsend'};

omoc_site = zeros(receptors,4);
for xx = 1:size(omoc_site,1)
    for yy = 1:size(omoc_site,2)
        omoc_site(xx,yy) = 2.1;
    end
end

yr_num  = 1;
aaa     = zeros(5,3,yr_num,receptors);
sim_rf  = zeros(152+31,yr_num,receptors);
sim_cl  = zeros(152+31,yr_num,receptors);
sim_no  = zeros(152+31,yr_num,receptors);
obs_all = zeros(152+31,yr_num,receptors);

for iyear = 2019
years = iyear;
% leap year
if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
    ndays = 366;
    ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
    
else
    ndays = 365;
    ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
end

inum = 0;
for irec  =  11 %[4,5,8,10,11]
    inum = inum + 1;
    site_name_tmp = char(site_name{irec});
    load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
    eval(['obs_hr1 = pm25_hourly_',site_name_tmp,';']);
    eval(['obs_dy1 = pm25_daily_',site_name_tmp,';']);
    
    load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear+1),'0101_',num2str(iyear+1),'1231.mat']);
    eval(['obs_hr2 = pm25_hourly_',site_name_tmp,';']);
    eval(['obs_dy2 = pm25_daily_',site_name_tmp,';']);
    
    obs_hr = cat(1,obs_hr1,obs_hr2);
    obs_dy = cat(1,obs_dy1,obs_dy2);
    
    for ii = 1:length(obs_hr)/24
        obs_tmp = obs_hr((ii-1)*24+1:ii*24);
        nan_num = length(find(isnan(obs_tmp)));
        if nan_num > 8
            obs_dy(ii) = NaN;
        end
    end
     obs_dy = obs_dy(214:396); % start from August 2 to Jan 31;

filepath1 = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_se\smoke_exp_output\'];
filepath2 = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear+1),'01_gfs_nohnf_receptor_se\smoke_exp_output\'];

load([filepath1,'smoke_exp_oc_gfed_rf_new.mat']);
ts_oc1a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_rf_new.mat']);
ts_oc1b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc1  = cat(1,ts_oc1a,ts_oc1b);

load([filepath1,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_oc2a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_oc2b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc2  = cat(1,ts_oc2a,ts_oc2b);

load([filepath1,'smoke_exp_oc_gfed_noscale.mat']);
ts_oc3a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_noscale.mat']);
ts_oc3b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc3  = cat(1,ts_oc3a,ts_oc3b);

load([filepath1,'smoke_exp_bc_gfed_rf_new.mat']);
ts_bc1a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_rf_new.mat']);
ts_bc1b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc1  = cat(1,ts_bc1a,ts_bc1b);

load([filepath1,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_bc2a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_bc2b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc2  = cat(1,ts_bc2a,ts_bc2b);

load([filepath1,'smoke_exp_bc_gfed_noscale.mat']);
ts_bc3a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_noscale.mat']);
ts_bc3b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc3  = cat(1,ts_bc3a,ts_bc3b);

load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
pm25bg = pm25base_years(iyear-2008);

% figure
    ii = irec; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change receptors

    ts_1(1:30) = ts_oc1(1:30,ii)*omoc_site(ii,2) + ts_bc1(1:30,ii);
    ts_1(31:121) = ts_oc1(31:121,ii)*omoc_site(ii,3) + ts_bc1(31:121,ii);
    ts_1(122:153) = ts_oc1(122:153,ii)*omoc_site(ii,3) + ts_bc1(122:153,ii);
    ts_1(154:184) = ts_oc1(154:184,ii)*omoc_site(ii,3) + ts_bc1(154:184,ii);
    
    ts_2(1:30) = ts_oc2(1:30,ii)*omoc_site(ii,2) + ts_bc2(1:30,ii);
    ts_2(31:121) = ts_oc2(31:121,ii)*omoc_site(ii,3) + ts_bc2(31:121,ii);
    ts_2(122:153) = ts_oc2(122:153,ii)*omoc_site(ii,3) + ts_bc2(122:153,ii);
    ts_2(154:184) = ts_oc2(154:184,ii)*omoc_site(ii,3) + ts_bc2(154:184,ii);
    
    ts_3(1:30) = ts_oc3(1:30,ii)*omoc_site(ii,2) + ts_bc3(1:30,ii);
    ts_3(31:121) = ts_oc3(31:121,ii)*omoc_site(ii,3) + ts_bc3(31:121,ii);
    ts_3(122:153) = ts_oc3(122:153,ii)*omoc_site(ii,3) + ts_bc3(122:153,ii);
    ts_3(154:184) = ts_oc3(154:184,ii)*omoc_site(ii,3) + ts_bc3(154:184,ii);

    start_day  = 1;
    end_day    = 152+31;
    t          = 214:1:ndays+31;

    sim1 = movmean(ts_1(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim2 = movmean(ts_2(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim3 = movmean(ts_3(start_day:end_day)+pm25bg,mv_day,'omitnan');
    obs  = movmean(obs_dy,mv_day,'omitnan');
     
    yearbase = 2018;
    sim_rf(:,iyear-yearbase,irec) = sim1;
    sim_cl(:,iyear-yearbase,irec) = sim2;
    sim_no(:,iyear-yearbase,irec) = sim3;
    obs_all(:,iyear-yearbase,irec)= obs;
end

end

% Statistics
mean_sim_rf = nanmean(sim_rf(:,:,irec),2);
mean_sim_cl = nanmean(sim_cl(:,:,irec),2);
mean_sim_no = nanmean(sim_no(:,:,irec),2);
mean_obs_all = nanmean(obs_all(:,:,irec),2);

var_all1 = [mean_sim_rf(1:end),mean_obs_all(1:end)];
var_all1 = rmmissing(var_all1);
nmb1      = (nanmean(mean_sim_rf(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr1       = corrcoef(var_all1(:,1),var_all1(:,2));
rmse1     = sqrt(mean((var_all1(:,1) - var_all1(:,2)).^2));
mean1     = mean(mean_sim_rf);

var_all2 = [mean_sim_cl(1:end),mean_obs_all(1:end)];
var_all2 = rmmissing(var_all2);
nmb2      = (nanmean(mean_sim_cl(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr2       = corrcoef(var_all2(:,1),var_all2(:,2));
rmse2     = sqrt(mean((var_all2(:,1) - var_all2(:,2)).^2));
mean2     = mean(mean_sim_cl);

var_all3 = [mean_sim_no(1:end),mean_obs_all(1:end)];
var_all3 = rmmissing(var_all3);
nmb3      = (nanmean(mean_sim_no(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr3       = corrcoef(var_all3(:,1),var_all3(:,2));
rmse3     = sqrt(mean((var_all3(:,1) - var_all3(:,2)).^2));
mean3     = mean(mean_sim_no);
mean_obs   = nanmean(mean_obs_all);

%%
subplot3 = subplot(3,2,3,'BoxStyle','full');
hold(subplot3,'on');

plot1 = plot(t,mean_sim_rf','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot3);
set(plot1,'DisplayName','INJ-RF',...
    'Color',[0.851 0.325490206480026 0.0980392172932625]);

plot2 = plot(t,mean_sim_cl','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot3);
set(plot2,'DisplayName','INJ-CLIM',...
    'Color',[0 0.498039215803146 0]);

plot3 = plot(t,mean_sim_no','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot3);
set(plot3,'DisplayName','CTL','Color',[0 0.447058826684952 0.74117648601532]);

plot5 = plot(t,mean_obs_all','LineWidth',2.5,'Parent',subplot3);
set(plot5,'DisplayName','OBS','Color',[0 0 0]);

% ylabel
ylabel('PM_2_._5 concentrations (\mug m^-^3)');

% title
title('(c) Newcastle (2019) ');

xlim(subplot3,[274 365+31]);
ylim(subplot3,[0 200]);
box(subplot3,'on');

set(subplot3,'FontSize',20,'XTick',[1 32 60 91 121 152 182 213 244 274 305 335 365],...
    'XTickLabel',...
    {'Jan','Feb','Mar','Apr','May','Jun','July','Aug','Sep','Oct','Nov','Dec','Jan',''});
set(gca,'Fontsize',24);
%% Sydney
receptors = 11;
loc_lon = [144.8728027,146.93986,149.043539,150.88733,150.91417,151.75965,153.1038,153.1356,145.0306,150.9058,151.6692];
loc_lat = [-37.80487823,-36.05182,-35.220606,-34.41706,-33.79424,-32.9312,-26.6917,-27.6125,-37.7784,-33.9328,-32.8961];
site_name = {'footscray','albury','florey','wollongong','prospect','newcastle','mountaincreek;','springwood','alphington','liverpool','wallsend'};

omoc_site = zeros(receptors,4);
for xx = 1:size(omoc_site,1)
    for yy = 1:size(omoc_site,2)
        omoc_site(xx,yy) = 2.1;
    end
end

yr_num  = 1;
aaa     = zeros(5,3,yr_num,receptors);
sim_rf  = zeros(152+31,yr_num,receptors);
sim_cl  = zeros(152+31,yr_num,receptors);
sim_no  = zeros(152+31,yr_num,receptors);
obs_all = zeros(152+31,yr_num,receptors);

for iyear = 2019
years = iyear;
% leap year
if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
    ndays = 366;
    ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
    
else
    ndays = 365;
    ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
end

inum = 0;
for irec  =  10 %[4,5,8,10,11]
    inum = inum + 1;
    site_name_tmp = char(site_name{irec});
    load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
    eval(['obs_hr1 = pm25_hourly_',site_name_tmp,';']);
    eval(['obs_dy1 = pm25_daily_',site_name_tmp,';']);
    
    load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear+1),'0101_',num2str(iyear+1),'1231.mat']);
    eval(['obs_hr2 = pm25_hourly_',site_name_tmp,';']);
    eval(['obs_dy2 = pm25_daily_',site_name_tmp,';']);
    
    obs_hr = cat(1,obs_hr1,obs_hr2);
    obs_dy = cat(1,obs_dy1,obs_dy2);
    
    for ii = 1:length(obs_hr)/24
        obs_tmp = obs_hr((ii-1)*24+1:ii*24);
        nan_num = length(find(isnan(obs_tmp)));
        if nan_num > 8
            obs_dy(ii) = NaN;
        end
    end
     obs_dy = obs_dy(214:396); % start from August 2 to Jan 31;

filepath1 = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_se\smoke_exp_output\'];
filepath2 = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear+1),'01_gfs_nohnf_receptor_se\smoke_exp_output\'];

load([filepath1,'smoke_exp_oc_gfed_rf_new.mat']);
ts_oc1a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_rf_new.mat']);
ts_oc1b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc1  = cat(1,ts_oc1a,ts_oc1b);

load([filepath1,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_oc2a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_oc2b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc2  = cat(1,ts_oc2a,ts_oc2b);

load([filepath1,'smoke_exp_oc_gfed_noscale.mat']);
ts_oc3a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_noscale.mat']);
ts_oc3b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc3  = cat(1,ts_oc3a,ts_oc3b);

load([filepath1,'smoke_exp_bc_gfed_rf_new.mat']);
ts_bc1a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_rf_new.mat']);
ts_bc1b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc1  = cat(1,ts_bc1a,ts_bc1b);

load([filepath1,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_bc2a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_bc2b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc2  = cat(1,ts_bc2a,ts_bc2b);

load([filepath1,'smoke_exp_bc_gfed_noscale.mat']);
ts_bc3a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_noscale.mat']);
ts_bc3b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc3  = cat(1,ts_bc3a,ts_bc3b);

load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
pm25bg = pm25base_years(iyear-2008);

% figure
    ii = irec; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change receptors

    ts_1(1:30) = ts_oc1(1:30,ii)*omoc_site(ii,2) + ts_bc1(1:30,ii);
    ts_1(31:121) = ts_oc1(31:121,ii)*omoc_site(ii,3) + ts_bc1(31:121,ii);
    ts_1(122:153) = ts_oc1(122:153,ii)*omoc_site(ii,3) + ts_bc1(122:153,ii);
    ts_1(154:184) = ts_oc1(154:184,ii)*omoc_site(ii,3) + ts_bc1(154:184,ii);
    
    ts_2(1:30) = ts_oc2(1:30,ii)*omoc_site(ii,2) + ts_bc2(1:30,ii);
    ts_2(31:121) = ts_oc2(31:121,ii)*omoc_site(ii,3) + ts_bc2(31:121,ii);
    ts_2(122:153) = ts_oc2(122:153,ii)*omoc_site(ii,3) + ts_bc2(122:153,ii);
    ts_2(154:184) = ts_oc2(154:184,ii)*omoc_site(ii,3) + ts_bc2(154:184,ii);
    
    ts_3(1:30) = ts_oc3(1:30,ii)*omoc_site(ii,2) + ts_bc3(1:30,ii);
    ts_3(31:121) = ts_oc3(31:121,ii)*omoc_site(ii,3) + ts_bc3(31:121,ii);
    ts_3(122:153) = ts_oc3(122:153,ii)*omoc_site(ii,3) + ts_bc3(122:153,ii);
    ts_3(154:184) = ts_oc3(154:184,ii)*omoc_site(ii,3) + ts_bc3(154:184,ii);

    start_day  = 1;
    end_day    = 152+31;
    t          = 214:1:ndays+31;

    % 10-day moving averaging
    sim1 = movmean(ts_1(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim2 = movmean(ts_2(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim3 = movmean(ts_3(start_day:end_day)+pm25bg,mv_day,'omitnan');
    obs  = movmean(obs_dy,mv_day,'omitnan');
    
    yearbase = 2018;
    sim_rf(:,iyear-yearbase,irec) = sim1;
    sim_cl(:,iyear-yearbase,irec) = sim2;
    sim_no(:,iyear-yearbase,irec) = sim3;
    obs_all(:,iyear-yearbase,irec)= obs;
end

end

%
mean_sim_rf = nanmean(sim_rf(:,:,irec),2);
mean_sim_cl = nanmean(sim_cl(:,:,irec),2);
mean_sim_no = nanmean(sim_no(:,:,irec),2);
mean_obs_all = nanmean(obs_all(:,:,irec),2);

var_all1 = [mean_sim_rf(1:end),mean_obs_all(1:end)];
var_all1 = rmmissing(var_all1);
nmb1      = (nanmean(mean_sim_rf(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr1       = corrcoef(var_all1(:,1),var_all1(:,2));
rmse1     = sqrt(mean((var_all1(:,1) - var_all1(:,2)).^2));
mean1     = mean(mean_sim_rf);

var_all2 = [mean_sim_cl(1:end),mean_obs_all(1:end)];
var_all2 = rmmissing(var_all2);
nmb2      = (nanmean(mean_sim_cl(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr2       = corrcoef(var_all2(:,1),var_all2(:,2));
rmse2     = sqrt(mean((var_all2(:,1) - var_all2(:,2)).^2));
mean2     = mean(mean_sim_cl);

var_all3 = [mean_sim_no(1:end),mean_obs_all(1:end)];
var_all3 = rmmissing(var_all3);
nmb3      = (nanmean(mean_sim_no(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr3       = corrcoef(var_all3(:,1),var_all3(:,2));
rmse3     = sqrt(mean((var_all3(:,1) - var_all3(:,2)).^2));
mean3     = mean(mean_sim_no);
mean_obs   = nanmean(mean_obs_all);

%%
subplot4 = subplot(3,2,4,'BoxStyle','full');
hold(subplot4,'on');

plot1 = plot(t,mean_sim_rf','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot4);
set(plot1,'DisplayName','INJ-RF',...
    'Color',[0.851 0.325490206480026 0.0980392172932625]);

plot2 = plot(t,mean_sim_cl','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot4);
set(plot2,'DisplayName','INJ-CLIM',...
    'Color',[0 0.498039215803146 0]);

plot3 = plot(t,mean_sim_no','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot4);
set(plot3,'DisplayName','CTL','Color',[0 0.447058826684952 0.74117648601532]);

plot5 = plot(t,mean_obs_all','LineWidth',2.5,'Parent',subplot4);
set(plot5,'DisplayName','OBS','Color',[0 0 0]);

% ylabel
% ylabel('PM_2_._5 concentrations (\mug m^-^3)');

% title
title('(d) Sydney (2019) ');

xlim(subplot4,[274 365+31]);
ylim(subplot4,[0 200]);
box(subplot4,'on');

set(subplot4,'FontSize',20,'XTick',[1 32 60 91 121 152 182 213 244 274 305 335 365],...
    'XTickLabel',...
    {'Jan','Feb','Mar','Apr','May','Jun','July','Aug','Sep','Oct','Nov','Dec','Jan',''});
set(gca,'Fontsize',24);
%% Canberra
receptors = 11;
loc_lon = [144.8728027,146.93986,149.043539,150.88733,150.91417,151.75965,153.1038,153.1356,145.0306,150.9058,151.6692];
loc_lat = [-37.80487823,-36.05182,-35.220606,-34.41706,-33.79424,-32.9312,-26.6917,-27.6125,-37.7784,-33.9328,-32.8961];
site_name = {'footscray','albury','florey','wollongong','prospect','newcastle','mountaincreek;','springwood','alphington','liverpool','wallsend'};

omoc_site = zeros(receptors,4);
for xx = 1:size(omoc_site,1)
    for yy = 1:size(omoc_site,2)
        omoc_site(xx,yy) = 2.1;
    end
end

yr_num  = 1;
aaa     = zeros(5,3,yr_num,receptors);
sim_rf  = zeros(152+31,yr_num,receptors);
sim_cl  = zeros(152+31,yr_num,receptors);
sim_no  = zeros(152+31,yr_num,receptors);
obs_all = zeros(152+31,yr_num,receptors);

for iyear = 2019
years = iyear;
% leap year
if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
    ndays = 366;
    ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
    
else
    ndays = 365;
    ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
end

inum = 0;
for irec  =  3 %[4,5,8,10,11]
    inum = inum + 1;
    site_name_tmp = char(site_name{irec});
    load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
    eval(['obs_hr1 = pm25_hourly_',site_name_tmp,';']);
    eval(['obs_dy1 = pm25_daily_',site_name_tmp,';']);
    
    load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear+1),'0101_',num2str(iyear+1),'1231.mat']);
    eval(['obs_hr2 = pm25_hourly_',site_name_tmp,';']);
    eval(['obs_dy2 = pm25_daily_',site_name_tmp,';']);
    
    obs_hr = cat(1,obs_hr1,obs_hr2);
    obs_dy = cat(1,obs_dy1,obs_dy2);
    
    for ii = 1:length(obs_hr)/24
        obs_tmp = obs_hr((ii-1)*24+1:ii*24);
        nan_num = length(find(isnan(obs_tmp)));
        if nan_num > 8
            obs_dy(ii) = NaN;
        end
    end
     obs_dy = obs_dy(214:396); % start from August 2 to Jan 31;

filepath1 = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_se\smoke_exp_output\'];
filepath2 = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear+1),'01_gfs_nohnf_receptor_se\smoke_exp_output\'];

load([filepath1,'smoke_exp_oc_gfed_rf_new.mat']);
ts_oc1a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_rf_new.mat']);
ts_oc1b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc1  = cat(1,ts_oc1a,ts_oc1b);

load([filepath1,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_oc2a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_oc2b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc2  = cat(1,ts_oc2a,ts_oc2b);

load([filepath1,'smoke_exp_oc_gfed_noscale.mat']);
ts_oc3a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_noscale.mat']);
ts_oc3b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc3  = cat(1,ts_oc3a,ts_oc3b);

load([filepath1,'smoke_exp_bc_gfed_rf_new.mat']);
ts_bc1a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_rf_new.mat']);
ts_bc1b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc1  = cat(1,ts_bc1a,ts_bc1b);

load([filepath1,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_bc2a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_bc2b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc2  = cat(1,ts_bc2a,ts_bc2b);

load([filepath1,'smoke_exp_bc_gfed_noscale.mat']);
ts_bc3a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_noscale.mat']);
ts_bc3b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc3  = cat(1,ts_bc3a,ts_bc3b);

load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
pm25bg = pm25base_years(iyear-2008);

% figure
    ii = irec; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change receptors

    ts_1(1:30) = ts_oc1(1:30,ii)*omoc_site(ii,2) + ts_bc1(1:30,ii);
    ts_1(31:121) = ts_oc1(31:121,ii)*omoc_site(ii,3) + ts_bc1(31:121,ii);
    ts_1(122:153) = ts_oc1(122:153,ii)*omoc_site(ii,3) + ts_bc1(122:153,ii);
    ts_1(154:184) = ts_oc1(154:184,ii)*omoc_site(ii,3) + ts_bc1(154:184,ii);
    
    ts_2(1:30) = ts_oc2(1:30,ii)*omoc_site(ii,2) + ts_bc2(1:30,ii);
    ts_2(31:121) = ts_oc2(31:121,ii)*omoc_site(ii,3) + ts_bc2(31:121,ii);
    ts_2(122:153) = ts_oc2(122:153,ii)*omoc_site(ii,3) + ts_bc2(122:153,ii);
    ts_2(154:184) = ts_oc2(154:184,ii)*omoc_site(ii,3) + ts_bc2(154:184,ii);
    
    ts_3(1:30) = ts_oc3(1:30,ii)*omoc_site(ii,2) + ts_bc3(1:30,ii);
    ts_3(31:121) = ts_oc3(31:121,ii)*omoc_site(ii,3) + ts_bc3(31:121,ii);
    ts_3(122:153) = ts_oc3(122:153,ii)*omoc_site(ii,3) + ts_bc3(122:153,ii);
    ts_3(154:184) = ts_oc3(154:184,ii)*omoc_site(ii,3) + ts_bc3(154:184,ii);

    start_day  = 1;
    end_day    = 152+31;
    t          = 214:1:ndays+31;

    sim1 = movmean(ts_1(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim2 = movmean(ts_2(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim3 = movmean(ts_3(start_day:end_day)+pm25bg,mv_day,'omitnan');
    obs  = movmean(obs_dy,mv_day,'omitnan');
    
    yearbase = 2018;
    sim_rf(:,iyear-yearbase,irec) = sim1;
    sim_cl(:,iyear-yearbase,irec) = sim2;
    sim_no(:,iyear-yearbase,irec) = sim3;
    obs_all(:,iyear-yearbase,irec)= obs;

end

end
%
mean_sim_rf = nanmean(sim_rf(:,:,irec),2);
mean_sim_cl = nanmean(sim_cl(:,:,irec),2);
mean_sim_no = nanmean(sim_no(:,:,irec),2);
mean_obs_all = nanmean(obs_all(:,:,irec),2);

var_all1 = [mean_sim_rf(1:end),mean_obs_all(1:end)];
var_all1 = rmmissing(var_all1);
nmb1      = (nanmean(mean_sim_rf(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr1       = corrcoef(var_all1(:,1),var_all1(:,2));
rmse1     = sqrt(mean((var_all1(:,1) - var_all1(:,2)).^2));
mean1     = mean(mean_sim_rf);

var_all2 = [mean_sim_cl(1:end),mean_obs_all(1:end)];
var_all2 = rmmissing(var_all2);
nmb2      = (nanmean(mean_sim_cl(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr2       = corrcoef(var_all2(:,1),var_all2(:,2));
rmse2     = sqrt(mean((var_all2(:,1) - var_all2(:,2)).^2));
mean2     = mean(mean_sim_cl);

var_all3 = [mean_sim_no(1:end),mean_obs_all(1:end)];
var_all3 = rmmissing(var_all3);
nmb3      = (nanmean(mean_sim_no(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr3       = corrcoef(var_all3(:,1),var_all3(:,2));
rmse3     = sqrt(mean((var_all3(:,1) - var_all3(:,2)).^2));
mean3     = mean(mean_sim_no);
mean_obs   = nanmean(mean_obs_all);

%%
subplot5 = subplot(3,2,5,'BoxStyle','full');
hold(subplot5,'on');
 
plot1 = plot(t,mean_sim_rf','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot5);
set(plot1,'DisplayName','INJ-RF',...
    'Color',[0.851 0.325490206480026 0.0980392172932625]);

plot2 = plot(t,mean_sim_cl','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot5);
set(plot2,'DisplayName','INJ-CLIM',...
    'Color',[0 0.498039215803146 0]);

plot3 = plot(t,mean_sim_no','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot5);
set(plot3,'DisplayName','CTL','Color',[0 0.447058826684952 0.74117648601532]);

plot5 = plot(t,mean_obs_all','LineWidth',2.5,'Parent',subplot5);
set(plot5,'DisplayName','OBS','Color',[0 0 0]);

% ylabel
% ylabel('PM_2_._5 concentrations (\mug m^-^3)');

% title
title('(e) Canberra (2019) ');

xlim(subplot5,[274 365+31]);
ylim(subplot5,[0 500]);
box(subplot5,'on');

set(subplot5,'FontSize',20,'XTick',[1 32 60 91 121 152 182 213 244 274 305 335 365],...
    'XTickLabel',...
    {'Jan','Feb','Mar','Apr','May','Jun','July','Aug','Sep','Oct','Nov','Dec','Jan',''});
set(gca,'Fontsize',24);

%% Alphington
receptors = 11;
loc_lon = [144.8728027,146.93986,149.043539,150.88733,150.91417,151.75965,153.1038,153.1356,145.0306,150.9058,151.6692];
loc_lat = [-37.80487823,-36.05182,-35.220606,-34.41706,-33.79424,-32.9312,-26.6917,-27.6125,-37.7784,-33.9328,-32.8961];
site_name = {'footscray','albury','florey','wollongong','prospect','newcastle','mountaincreek;','springwood','alphington','liverpool','wallsend'};

omoc_site = zeros(receptors,4);
for xx = 1:size(omoc_site,1)
    for yy = 1:size(omoc_site,2)
        omoc_site(xx,yy) = 2.1;
    end
end

yr_num  = 1;
aaa     = zeros(5,3,yr_num,receptors);
sim_rf  = zeros(152+31,yr_num,receptors);
sim_cl  = zeros(152+31,yr_num,receptors);
sim_no  = zeros(152+31,yr_num,receptors);
obs_all = zeros(152+31,yr_num,receptors);

for iyear = 2019
years = iyear;
% leap year
if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
    ndays = 366;
    ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
    
else
    ndays = 365;
    ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
end

inum = 0;
for irec  =  1 %[4,5,8,10,11]
    inum = inum + 1;
    site_name_tmp = char(site_name{irec});
    load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
%     eval(['obs_hr1 = pm25_hourly_',site_name_tmp,';']);
    eval(['obs_dy1 = pm25_daily_',site_name_tmp,';']);
    
    load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear+1),'0101_',num2str(iyear+1),'1231.mat']);
%     eval(['obs_hr2 = pm25_hourly_',site_name_tmp,';']);
    eval(['obs_dy2 = pm25_daily_',site_name_tmp,';']);
    
    obs_hr = cat(1,obs_hr1,obs_hr2);
    obs_dy = cat(1,obs_dy1,obs_dy2);
    
    for ii = 1:length(obs_hr)/24
        obs_tmp = obs_hr((ii-1)*24+1:ii*24);
        nan_num = length(find(isnan(obs_tmp)));
        if nan_num > 8
            obs_dy(ii) = NaN;
        end
    end
     obs_dy = obs_dy(214:396); % start from August 2 to Jan 31;

filepath1 = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_se\smoke_exp_output\'];
filepath2 = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear+1),'01_gfs_nohnf_receptor_se\smoke_exp_output\'];

load([filepath1,'smoke_exp_oc_gfed_rf_new.mat']);
ts_oc1a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_rf_new.mat']);
ts_oc1b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc1  = cat(1,ts_oc1a,ts_oc1b);

load([filepath1,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_oc2a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_oc2b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc2  = cat(1,ts_oc2a,ts_oc2b);

load([filepath1,'smoke_exp_oc_gfed_noscale.mat']);
ts_oc3a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_oc_gfed_noscale.mat']);
ts_oc3b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_oc3  = cat(1,ts_oc3a,ts_oc3b);

load([filepath1,'smoke_exp_bc_gfed_rf_new.mat']);
ts_bc1a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_rf_new.mat']);
ts_bc1b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc1  = cat(1,ts_bc1a,ts_bc1b);

load([filepath1,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_bc2a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
ts_bc2b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc2  = cat(1,ts_bc2a,ts_bc2b);

load([filepath1,'smoke_exp_bc_gfed_noscale.mat']);
ts_bc3a = squeeze(nansum(nansum(smoke_exp,1),2));
load([filepath2,'smoke_exp_bc_gfed_noscale.mat']);
ts_bc3b = squeeze(nansum(nansum(smoke_exp,1),2));
ts_bc3  = cat(1,ts_bc3a,ts_bc3b);

load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
pm25bg = pm25base_years(iyear-2008);

% figure
    ii = irec; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change receptors

    ts_1(1:30) = ts_oc1(1:30,ii)*omoc_site(ii,2) + ts_bc1(1:30,ii);
    ts_1(31:121) = ts_oc1(31:121,ii)*omoc_site(ii,3) + ts_bc1(31:121,ii);
    ts_1(122:153) = ts_oc1(122:153,ii)*omoc_site(ii,3) + ts_bc1(122:153,ii);
    ts_1(154:184) = ts_oc1(154:184,ii)*omoc_site(ii,3) + ts_bc1(154:184,ii);
    
    ts_2(1:30) = ts_oc2(1:30,ii)*omoc_site(ii,2) + ts_bc2(1:30,ii);
    ts_2(31:121) = ts_oc2(31:121,ii)*omoc_site(ii,3) + ts_bc2(31:121,ii);
    ts_2(122:153) = ts_oc2(122:153,ii)*omoc_site(ii,3) + ts_bc2(122:153,ii);
    ts_2(154:184) = ts_oc2(154:184,ii)*omoc_site(ii,3) + ts_bc2(154:184,ii);
    
    ts_3(1:30) = ts_oc3(1:30,ii)*omoc_site(ii,2) + ts_bc3(1:30,ii);
    ts_3(31:121) = ts_oc3(31:121,ii)*omoc_site(ii,3) + ts_bc3(31:121,ii);
    ts_3(122:153) = ts_oc3(122:153,ii)*omoc_site(ii,3) + ts_bc3(122:153,ii);
    ts_3(154:184) = ts_oc3(154:184,ii)*omoc_site(ii,3) + ts_bc3(154:184,ii);

    start_day  = 1;
    end_day    = 152+31;
    t          = 214:1:ndays+31;

    sim1 = movmean(ts_1(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim2 = movmean(ts_2(start_day:end_day)+pm25bg,mv_day,'omitnan');
    sim3 = movmean(ts_3(start_day:end_day)+pm25bg,mv_day,'omitnan');
    obs  = movmean(obs_dy,mv_day,'omitnan');
    
    yearbase = 2018;
    sim_rf(:,iyear-yearbase,irec) = sim1;
    sim_cl(:,iyear-yearbase,irec) = sim2;
    sim_no(:,iyear-yearbase,irec) = sim3;
    obs_all(:,iyear-yearbase,irec)= obs;

    mean_obs = nanmean(obs);
    mean_sim1 = nanmean(sim1);
    mean_sim2 = nanmean(sim2);
    mean_sim3 = nanmean(sim3);

    var_all1 = [sim1',obs];
    var_all1 = rmmissing(var_all1);
    nmb      = nanmean(mean_sim1 - mean_obs)/mean_obs;
    rr      = corrcoef(var_all1(:,1),var_all1(:,2));
    rmse    = sqrt(mean((var_all1(:,1) - var_all1(:,2)).^2));
    aaa(1,1,iyear-yearbase) = mean_obs;
    aaa(2,1,iyear-yearbase) = mean_sim1;
    aaa(3,1,iyear-yearbase) = rr(1,2);
    aaa(4,1,iyear-yearbase) = nmb;
    aaa(5,1,iyear-yearbase) = rmse;

    var_all2 = [sim2',obs]; %¦Ìg m-3

    var_all2 = rmmissing(var_all2);
    nmb      = nanmean(mean_sim2 - mean_obs)/mean_obs;
    rr      = corrcoef(var_all2(:,1),var_all2(:,2));
    rmse    = sqrt(mean((var_all2(:,1) - var_all2(:,2)).^2));

    aaa(1,2,iyear-yearbase) = mean_obs;
    aaa(2,2,iyear-yearbase) = mean_sim2;
    aaa(3,2,iyear-yearbase) = rr(1,2);
    aaa(4,2,iyear-yearbase) = nmb;
    aaa(5,2,iyear-yearbase) = rmse;

    var_all3 = [sim3',obs];
    var_all3 = rmmissing(var_all3);
    nmb      = nanmean(mean_sim3 - mean_obs)/mean_obs;
    rr       = corrcoef(var_all3(:,1),var_all3(:,2));
    rmse     = sqrt(mean((var_all3(:,1) - var_all3(:,2)).^2));

    aaa(1,3,iyear-yearbase) = mean_obs;
    aaa(2,3,iyear-yearbase) = mean_sim3;
    aaa(3,3,iyear-yearbase) = rr(1,2);
    aaa(4,3,iyear-yearbase) = nmb;
    aaa(5,3,iyear-yearbase) = rmse;
end

end
%
mean_sim_rf = nanmean(sim_rf(:,:,irec),2);
mean_sim_cl = nanmean(sim_cl(:,:,irec),2);
mean_sim_no = nanmean(sim_no(:,:,irec),2);
mean_obs_all = nanmean(obs_all(:,:,irec),2);

var_all1 = [mean_sim_rf(1:end),mean_obs_all(1:end)];
var_all1 = rmmissing(var_all1);
nmb1      = (nanmean(mean_sim_rf(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr1       = corrcoef(var_all1(:,1),var_all1(:,2));
rmse1     = sqrt(mean((var_all1(:,1) - var_all1(:,2)).^2));
mean1     = mean(mean_sim_rf);

var_all2 = [mean_sim_cl(1:end),mean_obs_all(1:end)];
var_all2 = rmmissing(var_all2);
nmb2      = (nanmean(mean_sim_cl(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr2       = corrcoef(var_all2(:,1),var_all2(:,2));
rmse2     = sqrt(mean((var_all2(:,1) - var_all2(:,2)).^2));
mean2     = mean(mean_sim_cl);

var_all3 = [mean_sim_no(1:end),mean_obs_all(1:end)];
var_all3 = rmmissing(var_all3);
nmb3      = (nanmean(mean_sim_no(1:end)) - nanmean(mean_obs_all(1:end)))/nanmean(mean_obs_all(1:end));
rr3       = corrcoef(var_all3(:,1),var_all3(:,2));
rmse3     = sqrt(mean((var_all3(:,1) - var_all3(:,2)).^2));
mean3     = mean(mean_sim_no);
mean_obs   = nanmean(mean_obs_all);

%%
subplot6 = subplot(3,2,6,'BoxStyle','full');
hold(subplot6,'on');

plot1 = plot(t,mean_sim_rf','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot6);
set(plot1,'DisplayName','INJ-RF',...
    'Color',[0.851 0.325490206480026 0.0980392172932625]);

plot2 = plot(t,mean_sim_cl','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot6);
set(plot2,'DisplayName','INJ-CLIM',...
    'Color',[0 0.498039215803146 0]);

plot3 = plot(t,mean_sim_no','MarkerSize',1,'Marker','.','LineWidth',2,'Parent',subplot6);
set(plot3,'DisplayName','CTL','Color',[0 0.447058826684952 0.74117648601532]);

plot5 = plot(t,mean_obs_all','LineWidth',2.5,'Parent',subplot6);
set(plot5,'DisplayName','OBS','Color',[0 0 0]);

% ylabel
% ylabel('PM_2_._5 concentrations (\mug m^-^3)');

% title
title('(f) Melbourne (2019) ');

xlim(subplot6,[274 365+31]);
ylim(subplot6,[0 60]);
box(subplot6,'on');

set(subplot6,'FontSize',20,'XTick',[1 32 60 91 121 152 182 213 244 274 305 335 365],...
    'XTickLabel',...
    {'Jan','Feb','Mar','Apr','May','Jun','July','Aug','Sep','Oct','Nov','Dec','Jan',''});
set(gca,'Fontsize',24);

