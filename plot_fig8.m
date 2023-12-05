%% Plotting script for Figure 8
% 6 sites annual mean PM2.5 comparison
% 1. Darwin
receptors = 4;
loc_lon = [130.94853,146.8257,149.1549,151.2704];
loc_lat = [-12.50779,-19.2542,-21.1595,-23.8627];
site_name = {'palmerston','westmackay','','southgladstone'};
omoc_site = zeros(receptors,4);
for xx = 1:size(omoc_site,1)
    for yy = 1:size(omoc_site,2)
        omoc_site(xx,yy) = 2.1;
    end
end

mon_num = 12;
yr_num  = 12;

sim_yearmean1 = zeros(yr_num,1);
sim_yearmean2 = zeros(yr_num,1);
sim_yearmean3 = zeros(yr_num,1);
obs_yearmean = zeros(yr_num,1);

mon_name = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
start_mon = 4;
end_mon   = 12;

for iyear = 2011:2020
    years = iyear;
    % leap year£¿
    if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
        ndays = 366;
        mon_days = [31,29,31,30,31,30,31,31,30,31,30,31];
    else
        ndays = 365;
        mon_days = [31,28,31,30,31,30,31,31,30,31,30,31];
    end

    isite = 1;
    site_name_tmp = char(site_name{isite});
    if iyear >= 2011
          load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
          eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
          eval(['obs_dy = pm25_daily_',site_name_tmp,';']);
%           for ii = 1:length(obs_hr)/24
%             obs_tmp = obs_hr((ii-1)*24+1:ii*24);
%             nan_num = length(find(isnan(obs_tmp)));
%             if nan_num > 8
%                 obs_dy(ii) = NaN;
%             end
%           end
    end

    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_north\smoke_exp_output\'];
    load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1 = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1 = test_tmp1 + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
    test_tmp2 = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
    test_tmp2 = test_tmp2 + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
    test_tmp3 = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
    test_tmp3 = test_tmp3 + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    filepath = 'F:\stilt_output\GDAS05\footprints_2018_gfs_nohnf_receptor_north\';
    diroutput = dir(fullfile(filepath,'*_foot.nc'));
    filename = {diroutput.name};

    filename_1 = char(filename(1));
    lon = ncread([filepath,filename_1],'lon');
    lat = ncread([filepath,filename_1],'lat');

    for ii = 1:size(test_tmp1,1)
        for jj = 1:size(test_tmp1,2)
            for tt = 1:size(test_tmp1,3)
                if test_tmp1(ii,jj,tt) == 0
                    test_tmp1(ii,jj,tt) = NaN;
                end
                if test_tmp2(ii,jj,tt) == 0
                    test_tmp2(ii,jj,tt) = NaN;
                end
                if test_tmp3(ii,jj,tt) == 0
                    test_tmp3(ii,jj,tt) = NaN;
                end
            end
        end
    end

   sim_yearmean1(iyear-2008)  =  nanmean(nansum(nansum(test_tmp1,1),2),3);
   sim_yearmean2(iyear-2008)  =  nanmean(nansum(nansum(test_tmp2,1),2),3);
   sim_yearmean3(iyear-2008)  =  nanmean(nansum(nansum(test_tmp3,1),2),3);
   
   if iyear >= 2011
       obs_yearmean(iyear-2008) = nanmean(obs_dy(sum(mon_days(1:start_mon-1))+1:sum(mon_days(1:end_mon))));
   end
end
load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
sim_tot1(1:2)    = NaN;
sim_tot1(3:12)   = pm25base_years(3:end)' + sim_yearmean1(3:end);
sim_tot2(1:2)    = NaN;
sim_tot2(3:12)   = pm25base_years(3:end)' + sim_yearmean2(3:end);
sim_tot3(1:2)    = NaN;
sim_tot3(3:12)   = pm25base_years(3:end)' + sim_yearmean3(3:end);
x1 = 1:1:12;
%%
figure1 = figure;
subplot1 = subplot(3,2,1,'BoxStyle','full');

obs_all = [obs_yearmean, obs_yearmean, obs_yearmean];
sim_all = [sim_yearmean3,sim_yearmean1,sim_yearmean2];
tot_all = [sim_tot3', sim_tot1', sim_tot2'];

bar1 = bar(x1,obs_all,'FaceColor','none','Parent',subplot1);
set(bar1(1),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(2),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(3),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
hold on
bar2 = bar(x1,sim_all,'FaceColor','none','Parent',subplot1);
set(bar2(3),'DisplayName','INJ-RF smoke PM_2_._5',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'EdgeColor','none');
set(bar2(2),'DisplayName','INJ-CLIM smoke PM_2_._5',...
    'FaceColor',[0 0.498039215803146 0],...
    'EdgeColor','none');
set(bar2(1),'DisplayName','CTL smoke PM_2_._5',...
    'FaceColor',[0 0.447058826684952 0.74117648601532],...
    'EdgeColor','none');
hold on
bar3 = bar(x1,tot_all,'FaceColor','none','Parent',subplot1);
set(bar3(3),'DisplayName','INJ-RF total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(bar3(2),'DisplayName','INJ-CLIM total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.498039215803146 0]);
set(bar3(1),'DisplayName','CTL total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.447058826684952 0.74117648601532]);

% ylabel('Concentrations [\mug m^-^3]');

xlim(subplot1,[2.5 12.5]);
ylim(subplot1,[0 30]);
set(subplot1,'YColor',[0 0 0]);
title(subplot1,'(a) Darwin');

box(subplot1,'on');
set(subplot1,'FontSize',18,'LineStyleOrderIndex',2,'XTick',...
    [1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'},...
    'YColor',[0 0 0],'YTick',[0 10 20 30 40]);

% legend
% legend1 = legend(subplot1,'show');
% set(legend1,'FontSize',14,'EdgeColor',[1 1 1]);

%%
% 2. South Gladstone
sim_yearmean1 = zeros(yr_num,1);
sim_yearmean2 = zeros(yr_num,1);
sim_yearmean3 = zeros(yr_num,1);
obs_yearmean  = zeros(yr_num,1);

start_mon = 4;
end_mon   = 12;

for iyear = 2009:2020
    years = iyear;
    % leap year
    if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
        ndays = 366;
        mon_days = [31,29,31,30,31,30,31,31,30,31,30,31];
    else
        ndays = 365;
        mon_days = [31,28,31,30,31,30,31,31,30,31,30,31];
    end

    isite = 4;
    site_name_tmp = char(site_name{isite});
    if iyear >= 2009
          load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
          eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
          eval(['obs_dy = pm25_daily_',site_name_tmp,';']);
%           for ii = 1:length(obs_hr)/24
%             obs_tmp = obs_hr((ii-1)*24+1:ii*24);
%             nan_num = length(find(isnan(obs_tmp)));
%             if nan_num > 8
%                 obs_dy(ii) = NaN;
%             end
%           end
    end

    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_north\smoke_exp_output\'];
    load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1 = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1 = test_tmp1 + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
    test_tmp2 = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
    test_tmp2 = test_tmp2 + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
    test_tmp3 = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
    test_tmp3 = test_tmp3 + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    filepath = 'F:\stilt_output\GDAS05\footprints_2018_gfs_nohnf_receptor_north\';
    diroutput = dir(fullfile(filepath,'*_foot.nc'));
    filename = {diroutput.name};

    filename_1 = char(filename(1));
    lon = ncread([filepath,filename_1],'lon');
    lat = ncread([filepath,filename_1],'lat');

    for ii = 1:size(test_tmp1,1)
        for jj = 1:size(test_tmp1,2)
            for tt = 1:size(test_tmp1,3)
                if test_tmp1(ii,jj,tt) == 0
                    test_tmp1(ii,jj,tt) = NaN;
                end
                if test_tmp2(ii,jj,tt) == 0
                    test_tmp2(ii,jj,tt) = NaN;
                end
                if test_tmp3(ii,jj,tt) == 0
                    test_tmp3(ii,jj,tt) = NaN;
                end
            end
        end
    end

   sim_yearmean1(iyear-2008)  =  nanmean(nansum(nansum(test_tmp1,1),2),3);
   sim_yearmean2(iyear-2008)  =  nanmean(nansum(nansum(test_tmp2,1),2),3);
   sim_yearmean3(iyear-2008)  =  nanmean(nansum(nansum(test_tmp3,1),2),3);
   if iyear >= 2009
       obs_yearmean(iyear-2008) = nanmean(obs_dy(sum(mon_days(1:start_mon-1))+1:sum(mon_days(1:end_mon))));
   end
end
load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
sim_tot1  = pm25base_years' + sim_yearmean1;
sim_tot2  = pm25base_years' + sim_yearmean2;
sim_tot3  = pm25base_years' + sim_yearmean3;
x1 = 1:1:12;
%%
subplot2 = subplot(3,2,2,'BoxStyle','full');

obs_all = [obs_yearmean, obs_yearmean, obs_yearmean];
sim_all = [sim_yearmean3, sim_yearmean1, sim_yearmean2];
tot_all = [sim_tot3, sim_tot1, sim_tot2];

bar1 = bar(x1,obs_all,'FaceColor','none','Parent',subplot2);
set(bar1(1),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(2),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(3),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
hold on
bar2 = bar(x1,sim_all,'FaceColor','none','Parent',subplot2);
set(bar2(3),'DisplayName','INJ-RF smoke PM_2_._5',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'EdgeColor','none');
set(bar2(2),'DisplayName','INJ-CLIM smoke PM_2_._5',...
    'FaceColor',[0 0.498039215803146 0],...
    'EdgeColor','none');
set(bar2(1),'DisplayName','CTL smoke PM_2_._5',...
    'FaceColor',[0 0.447058826684952 0.74117648601532],...
    'EdgeColor','none');
hold on
bar3 = bar(x1,tot_all,'FaceColor','none','Parent',subplot2);
set(bar3(3),'DisplayName','INJ-RF total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(bar3(2),'DisplayName','INJ-CLIM total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.498039215803146 0]);
set(bar3(1),'DisplayName','CTL total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.447058826684952 0.74117648601532]);

% ylabel('Concentrations [\mug m^-^3]');

xlim(subplot2,[0.5 12.5]);
ylim(subplot2,[0 30]);
set(subplot2,'YColor',[0 0 0]);
title(subplot2,'(b) South Gladstone');

box(subplot2,'on');
set(subplot2,'FontSize',18,'LineStyleOrderIndex',2,'XTick',...
    [1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'},...
    'YColor',[0 0 0],'YTick',[0 10 20 30 40]);
%% 
% 3. Logan city
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

mon_num   = 13;

sim_yearmean1 = zeros(yr_num,1);
sim_yearmean2 = zeros(yr_num,1);
sim_yearmean3 = zeros(yr_num,1);
obs_yearmean = zeros(yr_num,1);

mon_name = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
start_mon = 8;
end_mon   = 12;

for iyear = 2009:2020
    years = iyear;
    % leap year
    if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
        ndays = 366;
        mon_days = [31,29,31,30,31,30,31,31,30,31,30,31];

    else
        ndays = 365;
        mon_days = [31,28,31,30,31,30,31,31,30,31,30,31];
    end
    
    isite = 8; %%%%%% change receptors %%%%%%%%%
    site_name_tmp = char(site_name{isite});
    if iyear >= 2009
      load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
      eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
      eval(['obs_dy1 = pm25_daily_',site_name_tmp,';']);
%       for ii = 1:length(obs_hr)/24
%         obs_tmp = obs_hr((ii-1)*24+1:ii*24);
%         nan_num = length(find(isnan(obs_tmp)));
%         if nan_num > 8
%             obs_dy1(ii) = NaN;
%         end
%       end
      
      load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear+1),'0101_',num2str(iyear+1),'1231.mat']);
      eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
      eval(['obs_dy2 = pm25_daily_',site_name_tmp,';']);
%       for ii = 1:length(obs_hr)/24
%         obs_tmp = obs_hr((ii-1)*24+1:ii*24);
%         nan_num = length(find(isnan(obs_tmp)));
%         if nan_num > 8
%             obs_dy2(ii) = NaN;
%         end
%       end     
    end
    
    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_se\smoke_exp_output\'];
    load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1a = test_tmp1a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
    test_tmp2a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
    test_tmp2a = test_tmp2a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
    test_tmp3a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
    test_tmp3a = test_tmp3a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear+1),'01_gfs_nohnf_receptor_se\smoke_exp_output\'];
    load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1b = test_tmp1b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;    

    load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
    test_tmp2b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
    test_tmp2b = test_tmp2b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
    test_tmp3b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
    test_tmp3b = test_tmp3b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;
    
    test_tmp1 = cat(3,test_tmp1a,test_tmp1b);
    test_tmp2 = cat(3,test_tmp2a,test_tmp2b);
    test_tmp3 = cat(3,test_tmp3a,test_tmp3b);
    
    filepath = 'F:\stilt_output\GDAS05\footprints_2018_gfs_nohnf_receptor_north\';
    diroutput = dir(fullfile(filepath,'*_foot.nc'));
    filename = {diroutput.name};

    filename_1 = char(filename(1));
    lon = ncread([filepath,filename_1],'lon');
    lat = ncread([filepath,filename_1],'lat');

    for ii = 1:size(test_tmp1,1)
        for jj = 1:size(test_tmp1,2)
            for tt = 1:size(test_tmp1,3)
                if test_tmp1(ii,jj,tt) == 0
                    test_tmp1(ii,jj,tt) = NaN;
                end
                if test_tmp2(ii,jj,tt) == 0
                    test_tmp2(ii,jj,tt) = NaN;
                end
                if test_tmp3(ii,jj,tt) == 0
                    test_tmp3(ii,jj,tt) = NaN;
                end
            end
        end
    end
    
    mon_nextyear = 13;
 
   sim_yearmean1(iyear-2008)  =  nanmean(nansum(nansum(test_tmp1,1),2),3);
   sim_yearmean2(iyear-2008)  =  nanmean(nansum(nansum(test_tmp2,1),2),3);
   sim_yearmean3(iyear-2008)  =  nanmean(nansum(nansum(test_tmp3,1),2),3);
    if iyear >= 2009
       a1 = obs_dy1(sum(mon_days(1:start_mon-1))+1:sum(mon_days(1:end_mon)));
       a2 = obs_dy2(1:31);
       obs_dy = cat(1,a1,a2);
       obs_yearmean(iyear-2008) = nanmean(obs_dy);
    end
end
load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
sim_tot1   = pm25base_years' + sim_yearmean1;
sim_tot2   = pm25base_years' + sim_yearmean2;
sim_tot3   = pm25base_years' + sim_yearmean3;
x1 = 1:1:12;

%%
subplot3 = subplot(3,2,3,'BoxStyle','full');

obs_all = [obs_yearmean, obs_yearmean, obs_yearmean];
sim_all = [sim_yearmean3, sim_yearmean1, sim_yearmean2];
tot_all = [sim_tot3, sim_tot1, sim_tot2];

rtmp     = corrcoef(sim_tot2,obs_yearmean);
rr(3,1)  = rtmp(1,2);
rtmp     = corrcoef(sim_tot1,obs_yearmean);
rr(3,2)  = rtmp(1,2);
rtmp     = corrcoef(sim_tot3,obs_yearmean);
rr(3,3)  = rtmp(1,2);
nmb(3,1) =  sum(sim_tot2 - obs_yearmean)/sum(obs_yearmean);
nmb(3,2) =  sum(sim_tot1 - obs_yearmean)/sum(obs_yearmean);
nmb(3,3) =  sum(sim_tot3 - obs_yearmean)/sum(obs_yearmean);
rmse1(3,1) = rmse(sim_tot2,obs_yearmean);
rmse1(3,2) = rmse(sim_tot1,obs_yearmean);
rmse1(3,3) = rmse(sim_tot3,obs_yearmean);

nmb_year2(:,1,2) = (sim_tot2 - obs_yearmean)./obs_yearmean;
nmb_year2(:,2,2) = (sim_tot1 - obs_yearmean)./obs_yearmean;
nmb_year2(:,3,2) = (sim_tot3 - obs_yearmean)./obs_yearmean;

bar1 = bar(x1,obs_all,'FaceColor','none','Parent',subplot3);
set(bar1(1),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(2),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(3),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
hold on
bar2 = bar(x1,sim_all,'FaceColor','none','Parent',subplot3);
set(bar2(3),'DisplayName','INJ-RF smoke PM_2_._5',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'EdgeColor','none');
set(bar2(2),'DisplayName','INJ-CLIM smoke PM_2_._5',...
    'FaceColor',[0 0.498039215803146 0],...
    'EdgeColor','none');
set(bar2(1),'DisplayName','CTL smoke PM_2_._5',...
    'FaceColor',[0 0.447058826684952 0.74117648601532],...
    'EdgeColor','none');
hold on
bar3 = bar(x1,tot_all,'FaceColor','none','Parent',subplot3);
set(bar3(3),'DisplayName','INJ-RF total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(bar3(2),'DisplayName','INJ-CLIM total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.498039215803146 0]);
set(bar3(1),'DisplayName','CTL total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.447058826684952 0.74117648601532]);

ylabel('PM_2_._5 concentrations [\mug m^-^3]');

xlim(subplot3,[0.5 12.5]);
ylim(subplot3,[0 30]);
set(subplot3,'YColor',[0 0 0]);
title(subplot3,'(c) Brisbane');

box(subplot3,'on');
set(subplot3,'FontSize',18,'LineStyleOrderIndex',2,'XTick',...
    [1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'},...
    'YColor',[0 0 0],'YTick',[0 10 20 30 40]);
%%
% 4. Newcastle
mon_num   = 13;

sim_yearmean1 = zeros(yr_num,1);
sim_yearmean2 = zeros(yr_num,1);
sim_yearmean3 = zeros(yr_num,1);
obs_yearmean = zeros(yr_num,1);

mon_name = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
start_mon = 8;
end_mon   = 12;

for iyear = 2009:2020
    years = iyear;
    % leap year
    if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
        ndays = 366;
        mon_days = [31,29,31,30,31,30,31,31,30,31,30,31];

    else
        ndays = 365;
        mon_days = [31,28,31,30,31,30,31,31,30,31,30,31];
    end
    
    isite = 11; %%%%%% change receptors %%%%%%%%%
    site_name_tmp = char(site_name{isite});
    if iyear >= 2009
      load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
      eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
      eval(['obs_dy1 = pm25_daily_',site_name_tmp,';']);
%       for ii = 1:length(obs_hr)/24
%         obs_tmp = obs_hr((ii-1)*24+1:ii*24);
%         nan_num = length(find(isnan(obs_tmp)));
%         if nan_num > 8
%             obs_dy1(ii) = NaN;
%         end
%       end
      
      load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear+1),'0101_',num2str(iyear+1),'1231.mat']);
      eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
      eval(['obs_dy2 = pm25_daily_',site_name_tmp,';']);
%       for ii = 1:length(obs_hr)/24
%         obs_tmp = obs_hr((ii-1)*24+1:ii*24);
%         nan_num = length(find(isnan(obs_tmp)));
%         if nan_num > 8
%             obs_dy2(ii) = NaN;
%         end
%       end     
    end
    
    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_se\smoke_exp_output\'];
    load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1a = test_tmp1a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
    test_tmp2a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
    test_tmp2a = test_tmp2a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
    test_tmp3a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
    test_tmp3a = test_tmp3a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear+1),'01_gfs_nohnf_receptor_se\smoke_exp_output\'];
    load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1b = test_tmp1b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;    

    load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
    test_tmp2b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
    test_tmp2b = test_tmp2b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
    test_tmp3b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
    test_tmp3b = test_tmp3b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;
    
    test_tmp1 = cat(3,test_tmp1a,test_tmp1b);
    test_tmp2 = cat(3,test_tmp2a,test_tmp2b);
    test_tmp3 = cat(3,test_tmp3a,test_tmp3b);
    
    filepath = 'F:\stilt_output\GDAS05\footprints_2018_gfs_nohnf_receptor_north\';
    diroutput = dir(fullfile(filepath,'*_foot.nc'));
    filename = {diroutput.name};

    filename_1 = char(filename(1));
    lon = ncread([filepath,filename_1],'lon');
    lat = ncread([filepath,filename_1],'lat');

    for ii = 1:size(test_tmp1,1)
        for jj = 1:size(test_tmp1,2)
            for tt = 1:size(test_tmp1,3)
                if test_tmp1(ii,jj,tt) == 0
                    test_tmp1(ii,jj,tt) = NaN;
                end
                if test_tmp2(ii,jj,tt) == 0
                    test_tmp2(ii,jj,tt) = NaN;
                end
                if test_tmp3(ii,jj,tt) == 0
                    test_tmp3(ii,jj,tt) = NaN;
                end
            end
        end
    end
    
    mon_nextyear = 13;

   sim_yearmean1(iyear-2008)  =  nanmean(nansum(nansum(test_tmp1,1),2),3);
   sim_yearmean2(iyear-2008)  =  nanmean(nansum(nansum(test_tmp2,1),2),3);
   sim_yearmean3(iyear-2008)  =  nanmean(nansum(nansum(test_tmp3,1),2),3);
    if iyear >= 2009
       a1 = obs_dy1(sum(mon_days(1:start_mon-1))+1:sum(mon_days(1:end_mon)));
       a2 = obs_dy2(1:31);
       obs_dy = cat(1,a1,a2);
       obs_yearmean(iyear-2008) = nanmean(obs_dy);
    end
end

load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
sim_tot1   = pm25base_years' + sim_yearmean1;
sim_tot2   = pm25base_years' + sim_yearmean2;
sim_tot3   = pm25base_years' + sim_yearmean3;
x1 = 1:1:12;

%%
subplot4 = subplot(3,2,4,'BoxStyle','full');

obs_all = [obs_yearmean, obs_yearmean, obs_yearmean];
sim_all = [sim_yearmean3, sim_yearmean1, sim_yearmean2];
tot_all = [sim_tot3, sim_tot1, sim_tot2];

bar1 = bar(x1,obs_all,'FaceColor','none','Parent',subplot4);
set(bar1(1),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(2),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(3),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
hold on
bar2 = bar(x1,sim_all,'FaceColor','none','Parent',subplot4);
set(bar2(3),'DisplayName','INJ-RF smoke PM_2_._5',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'EdgeColor','none');
set(bar2(2),'DisplayName','INJ-CLIM smoke PM_2_._5',...
    'FaceColor',[0 0.498039215803146 0],...
    'EdgeColor','none');
set(bar2(1),'DisplayName','CTL smoke PM_2_._5',...
    'FaceColor',[0 0.447058826684952 0.74117648601532],...
    'EdgeColor','none');
hold on
bar3 = bar(x1,tot_all,'FaceColor','none','Parent',subplot4);
set(bar3(3),'DisplayName','INJ-RF total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(bar3(2),'DisplayName','INJ-CLIM total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.498039215803146 0]);
set(bar3(1),'DisplayName','CTL total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.447058826684952 0.74117648601532]);

xlim(subplot4,[0.5 12.5]);
ylim(subplot4,[0 30]);
set(subplot4,'YColor',[0 0 0]);
title(subplot4,'(d) Newcastle');

box(subplot4,'on');
set(subplot4,'FontSize',18,'LineStyleOrderIndex',2,'XTick',...
    [1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'},...
    'YColor',[0 0 0],'YTick',[0 10 20 30 40]);

%%
mon_num   = 13;

sim_yearmean1 = zeros(yr_num,1);
sim_yearmean2 = zeros(yr_num,1);
sim_yearmean3 = zeros(yr_num,1);
obs_yearmean = zeros(yr_num,1);

mon_name = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
start_mon = 8;
end_mon   = 12;

for iyear = 2009:2020
    years = iyear;
    % leap year
    if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
        ndays = 366;
        mon_days = [31,29,31,30,31,30,31,31,30,31,30,31];

    else
        ndays = 365;
        mon_days = [31,28,31,30,31,30,31,31,30,31,30,31];
    end
    
    isite = 10; %%%%%% change receptors %%%%%%%%%
    site_name_tmp = char(site_name{isite});
    if iyear >= 2009
      load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
      eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
      eval(['obs_dy1 = pm25_daily_',site_name_tmp,';']);
%       for ii = 1:length(obs_hr)/24
%         obs_tmp = obs_hr((ii-1)*24+1:ii*24);
%         nan_num = length(find(isnan(obs_tmp)));
%         if nan_num > 8
%             obs_dy1(ii) = NaN;
%         end
%       end
      
      load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear+1),'0101_',num2str(iyear+1),'1231.mat']);
      eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
      eval(['obs_dy2 = pm25_daily_',site_name_tmp,';']);
%       for ii = 1:length(obs_hr)/24
%         obs_tmp = obs_hr((ii-1)*24+1:ii*24);
%         nan_num = length(find(isnan(obs_tmp)));
%         if nan_num > 8
%             obs_dy2(ii) = NaN;
%         end
%       end     
    end
    
    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_se\smoke_exp_output\'];
    load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1a = test_tmp1a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
    test_tmp2a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
    test_tmp2a = test_tmp2a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
    test_tmp3a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
    test_tmp3a = test_tmp3a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear+1),'01_gfs_nohnf_receptor_se\smoke_exp_output\'];
    load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1b = test_tmp1b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;    

    load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
    test_tmp2b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
    test_tmp2b = test_tmp2b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
    test_tmp3b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
    test_tmp3b = test_tmp3b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;
    
    test_tmp1 = cat(3,test_tmp1a,test_tmp1b);
    test_tmp2 = cat(3,test_tmp2a,test_tmp2b);
    test_tmp3 = cat(3,test_tmp3a,test_tmp3b);
    
    filepath = 'F:\stilt_output\GDAS05\footprints_2018_gfs_nohnf_receptor_north\';
    diroutput = dir(fullfile(filepath,'*_foot.nc'));
    filename = {diroutput.name};

    filename_1 = char(filename(1));
    lon = ncread([filepath,filename_1],'lon');
    lat = ncread([filepath,filename_1],'lat');

    for ii = 1:size(test_tmp1,1)
        for jj = 1:size(test_tmp1,2)
            for tt = 1:size(test_tmp1,3)
                if test_tmp1(ii,jj,tt) == 0
                    test_tmp1(ii,jj,tt) = NaN;
                end
                if test_tmp2(ii,jj,tt) == 0
                    test_tmp2(ii,jj,tt) = NaN;
                end
                if test_tmp3(ii,jj,tt) == 0
                    test_tmp3(ii,jj,tt) = NaN;
                end
            end
        end
    end
    
    mon_nextyear = 13;

   sim_yearmean1(iyear-2008)  =  nanmean(nansum(nansum(test_tmp1,1),2),3);
   sim_yearmean2(iyear-2008)  =  nanmean(nansum(nansum(test_tmp2,1),2),3);
   sim_yearmean3(iyear-2008)  =  nanmean(nansum(nansum(test_tmp3,1),2),3);
   if iyear >= 2009
       a1 = obs_dy1(sum(mon_days(1:start_mon-1))+1:sum(mon_days(1:end_mon)));
       a2 = obs_dy2(1:31);
       obs_dy = cat(1,a1,a2);
       obs_yearmean(iyear-2008) = nanmean(obs_dy);
   end
end
load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
sim_tot1   = pm25base_years' + sim_yearmean1;
sim_tot2   = pm25base_years' + sim_yearmean2;
sim_tot3   = pm25base_years' + sim_yearmean3;

x1 = 1:1:12;
%%
subplot5 = subplot(3,2,5,'BoxStyle','full');

obs_all = [obs_yearmean, obs_yearmean, obs_yearmean];
sim_all = [sim_yearmean3, sim_yearmean1, sim_yearmean2];
tot_all = [sim_tot3, sim_tot1, sim_tot2];

bar1 = bar(x1,obs_all,'FaceColor','none','Parent',subplot5);
set(bar1(1),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(2),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(3),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
hold on
bar2 = bar(x1,sim_all,'FaceColor','none','Parent',subplot5);
set(bar2(3),'DisplayName','INJ-RF smoke PM_2_._5',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'EdgeColor','none');
set(bar2(2),'DisplayName','INJ-CLIM smoke PM_2_._5',...
    'FaceColor',[0 0.498039215803146 0],...
    'EdgeColor','none');
set(bar2(1),'DisplayName','CTL smoke PM_2_._5',...
    'FaceColor',[0 0.447058826684952 0.74117648601532],...
    'EdgeColor','none');
hold on
bar3 = bar(x1,tot_all,'FaceColor','none','Parent',subplot5);
set(bar3(3),'DisplayName','INJ-RF total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(bar3(2),'DisplayName','INJ-CLIM total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.498039215803146 0]);
set(bar3(1),'DisplayName','CTL total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.447058826684952 0.74117648601532]);

xlim(subplot5,[0.5 12.5]);
ylim(subplot5,[0 30]);
set(subplot5,'YColor',[0 0 0]);
title(subplot5,'(e) Sydney');

box(subplot5,'on');
set(subplot5,'FontSize',18,'LineStyleOrderIndex',2,'XTick',...
    [1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'},...
    'YColor',[0 0 0],'YTick',[0 10 20 30 40]);

%%
sim_yearmean1 = zeros(yr_num,1);
sim_yearmean2 = zeros(yr_num,1);
sim_yearmean3 = zeros(yr_num,1);
obs_yearmean = zeros(yr_num,1);

mon_name = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
start_mon = 8;
end_mon   = 12;

for iyear = 2009:2020
    years = iyear;
    % leap year
    if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
        ndays = 366;
        mon_days = [31,29,31,30,31,30,31,31,30,31,30,31];

    else
        ndays = 365;
        mon_days = [31,28,31,30,31,30,31,31,30,31,30,31];
    end
    
    isite = 9; %%%%%% change receptors %%%%%%%%%
    site_name_tmp = char(site_name{isite});
    if iyear >= 2009
      load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear),'0101_',num2str(iyear),'1231.mat']);
%       eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
      eval(['obs_dy1 = pm25_daily_',site_name_tmp,';']);
%       for ii = 1:length(obs_hr)/24
%         obs_tmp = obs_hr((ii-1)*24+1:ii*24);
%         nan_num = length(find(isnan(obs_tmp)));
%         if nan_num > 8
%             obs_dy1(ii) = NaN;
%         end
%       end
      
      load(['E:\work_for_wildfires\data\air quality\mat_data\',site_name_tmp,'_pm25_',num2str(iyear+1),'0101_',num2str(iyear+1),'1231.mat']);
%       eval(['obs_hr = pm25_hourly_',site_name_tmp,';']);
      eval(['obs_dy2 = pm25_daily_',site_name_tmp,';']);
%       for ii = 1:length(obs_hr)/24
%         obs_tmp = obs_hr((ii-1)*24+1:ii*24);
%         nan_num = length(find(isnan(obs_tmp)));
%         if nan_num > 8
%             obs_dy2(ii) = NaN;
%         end
%       end     
    end
    
    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_se\smoke_exp_output\'];
    load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1a = test_tmp1a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
    test_tmp2a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
    test_tmp2a = test_tmp2a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
    test_tmp3a = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
    test_tmp3a = test_tmp3a + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear+1),'01_gfs_nohnf_receptor_se\smoke_exp_output\'];
    load([filepath,'smoke_exp_bc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_is4fire_wholeday_smooth_scale_mami.mat']);
    test_tmp1b = test_tmp1b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;    

    load([filepath,'smoke_exp_bc_gfed_rf_new.mat']);
    test_tmp2b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_rf_new.mat']);
    test_tmp2b = test_tmp2b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;

    load([filepath,'smoke_exp_bc_gfed_noscale.mat']);
    test_tmp3b = squeeze(smoke_exp(:,:,:,isite));
    load([filepath,'smoke_exp_oc_gfed_noscale.mat']);
    test_tmp3b = test_tmp3b + squeeze(smoke_exp(:,:,:,isite)) * 2.1;
    
    test_tmp1 = cat(3,test_tmp1a,test_tmp1b);
    test_tmp2 = cat(3,test_tmp2a,test_tmp2b);
    test_tmp3 = cat(3,test_tmp3a,test_tmp3b);
    
    filepath = 'F:\stilt_output\GDAS05\footprints_2018_gfs_nohnf_receptor_north\';
    diroutput = dir(fullfile(filepath,'*_foot.nc'));
    filename = {diroutput.name};

    filename_1 = char(filename(1));
    lon = ncread([filepath,filename_1],'lon');
    lat = ncread([filepath,filename_1],'lat');

    for ii = 1:size(test_tmp1,1)
        for jj = 1:size(test_tmp1,2)
            for tt = 1:size(test_tmp1,3)
                if test_tmp1(ii,jj,tt) == 0
                    test_tmp1(ii,jj,tt) = NaN;
                end
                if test_tmp2(ii,jj,tt) == 0
                    test_tmp2(ii,jj,tt) = NaN;
                end
                if test_tmp3(ii,jj,tt) == 0
                    test_tmp3(ii,jj,tt) = NaN;
                end
            end
        end
    end
    
    mon_nextyear = 13;
   sim_yearmean1(iyear-2008)  =  nanmean(nansum(nansum(test_tmp1,1),2),3);
   sim_yearmean2(iyear-2008)  =  nanmean(nansum(nansum(test_tmp2,1),2),3);
   sim_yearmean3(iyear-2008)  =  nanmean(nansum(nansum(test_tmp3,1),2),3);
   if iyear >= 2009
       a1 = obs_dy1(sum(mon_days(1:start_mon-1))+1:sum(mon_days(1:end_mon)));
       a2 = obs_dy2(1:31);
       obs_dy = cat(1,a1,a2);
       obs_yearmean(iyear-2008) = nanmean(obs_dy);
   end
end

load(['F:\stilt_output\background_pm25\',site_name_tmp,'_2009_2020_pm25.mat']);
sim_tot1   = pm25base_years' + sim_yearmean1;
sim_tot2   = pm25base_years' + sim_yearmean2;
sim_tot3   = pm25base_years' + sim_yearmean3;

x1 = 1:1:12;

%%
subplot6 = subplot(3,2,6,'BoxStyle','full');

obs_all = [obs_yearmean, obs_yearmean, obs_yearmean];
sim_all = [sim_yearmean3, sim_yearmean1, sim_yearmean2];
tot_all = [sim_tot3, sim_tot1, sim_tot2];

bar1 = bar(x1,obs_all,'FaceColor','none','Parent',subplot6);
set(bar1(1),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(2),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
set(bar1(3),'DisplayName','OBS total PM_2_._5','LineWidth',2.5);
hold on
bar2 = bar(x1,sim_all,'FaceColor','none','Parent',subplot6);
set(bar2(3),'DisplayName','INJ-RF smoke PM_2_._5',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'EdgeColor','none');
set(bar2(2),'DisplayName','INJ-CLIM smoke PM_2_._5',...
    'FaceColor',[0 0.498039215803146 0],...
    'EdgeColor','none');
set(bar2(1),'DisplayName','CTL smoke PM_2_._5',...
    'FaceColor',[0 0.447058826684952 0.74117648601532],...
    'EdgeColor','none');
hold on
bar3 = bar(x1,tot_all,'FaceColor','none','Parent',subplot6);
set(bar3(3),'DisplayName','INJ-RF total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(bar3(2),'DisplayName','INJ-CLIM total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.498039215803146 0]);
set(bar3(1),'DisplayName','CTL total PM_2_._5','LineWidth',2.5,...
    'EdgeColor',[0 0.447058826684952 0.74117648601532]);

xlim(subplot6,[0.5 12.5]);
ylim(subplot6,[0 30]);
set(subplot6,'YColor',[0 0 0]);
title(subplot6,'(f) Melbourne');

box(subplot6,'on');
set(subplot6,'FontSize',18,'LineStyleOrderIndex',2,'XTick',...
    [1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'},...
    'YColor',[0 0 0],'YTick',[0 10 20 30 40]);
