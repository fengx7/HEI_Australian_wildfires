%% Plotting script for Figure 5
load('E:\work_for_wildfires\random_forest\emis_scale_apr_dec_2009_2020_new.mat');
area = ncread('E:\work_for_wildfires\biomass_burning_emis\GFED\area_GFEDv4s.nc','area');
location = 'E:\work_for_wildfires\biomass_burning_emis\GFED\';
prefix = 'GFEDv4s_';
postfix = '.nc';

% Read lon/lat of GFED output
dirfiles = dir([location,prefix,'*']);
filename1 = dirfiles.name;
lon_global_gfed = ncread([location,filename1],'lon');
lat_global_gfed = ncread([location,filename1],'lat');

% Read lon/lat of STILT output (study area)
filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(2018),'_gfs_nohnf_receptor_north\'];
diroutput = dir(fullfile(filepath,'*_foot.nc'));
filename = {diroutput.name};
filename_1 = char(filename(1));
lon = ncread([filepath,filename_1],'lon');
lat = ncread([filepath,filename_1],'lat');

% Match the study area in GFED file
lon_indx1 = find( abs(min(lon) - lon_global_gfed) == min( abs(min(lon) - lon_global_gfed) ));
lon_indx2 = find( abs(max(lon) - lon_global_gfed) == min( abs(max(lon) - lon_global_gfed) ));
lat_indx1 = find( abs(min(lat) - lat_global_gfed) == min( abs(min(lat) - lat_global_gfed) ));
lat_indx2 = find( abs(max(lat) - lat_global_gfed) == min( abs(max(lat) - lat_global_gfed) ));

% Convert unit to Gg month-1
sf   = 3600 * 24 / 1e6 /9;
for ii = 1:size(aemit_mean_1_oc,3)
    oc_1_emis(:,:,ii) = aemit_mean_1_oc(:,:,ii) .* area * sf;  
    oc_2_emis(:,:,ii) = aemit_mean_2_oc(:,:,ii) .* area * sf;
    oc_4_emis(:,:,ii) = aemit_mean_4_oc(:,:,ii) .* area * sf;
end

aemit_mean_1_oc_au  = oc_1_emis(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);
aemit_mean_2_oc_au  = oc_2_emis(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);
aemit_mean_4_oc_au  = oc_4_emis(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);

mean_1_oc_au        = squeeze(mean(aemit_mean_1_oc_au,3));
mean_2_oc_au        = squeeze(mean(aemit_mean_2_oc_au,3));
mean_4_oc_au        = squeeze(mean(aemit_mean_4_oc_au,3));

diff_perc1           = (mean_4_oc_au - mean_1_oc_au)./mean_2_oc_au*100;

% North Australia
a1 = squeeze(sum(sum(aemit_mean_1_oc_au(25:156,107:142,:)))); % for scale emission mami
a2 = squeeze(sum(sum(aemit_mean_2_oc_au(25:156,107:142,:)))); % no scale
a4 = squeeze(sum(sum(aemit_mean_4_oc_au(25:156,107:142,:)))); % random forest

emis_north = [a1,a4,a2];
rd_rate_1 = (a2 - a4)./a2 * 100;
rd_rate_2 = (a2 - a1)./a2 * 100;
rd_north  = [rd_rate_2,rd_rate_1];
%%
load('E:\work_for_wildfires\random_forest\emis_scale_aug_jan_2009_2020.mat');
area = ncread('E:\work_for_wildfires\biomass_burning_emis\GFED\area_GFEDv4s.nc','area');
sf   = 3600 * 24 / 1e6 /6;
for ii = 1:size(aemit_mean_1_oc,3)
    oc_1_emis(:,:,ii) = aemit_mean_1_oc(:,:,ii) .* area * sf;  % unit: Gg month-1
    oc_2_emis(:,:,ii) = aemit_mean_2_oc(:,:,ii) .* area * sf;
    oc_4_emis(:,:,ii) = aemit_mean_4_oc(:,:,ii) .* area * sf;
end
aemit_mean_1_oc_au  = oc_1_emis(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);
aemit_mean_2_oc_au  = oc_2_emis(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);
aemit_mean_4_oc_au  = oc_4_emis(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);

mean_1_oc_au        = squeeze(mean(aemit_mean_1_oc_au,3));
mean_2_oc_au        = squeeze(mean(aemit_mean_2_oc_au,3));
mean_4_oc_au        = squeeze(mean(aemit_mean_4_oc_au,3));
diff_perc2           = (mean_4_oc_au - mean_1_oc_au)./mean_2_oc_au*100;

% Southeastern Australia
a5 = squeeze(sum(sum(aemit_mean_1_oc_au(113:168,7:86,:)))); % for scale emission mami
a6 = squeeze(sum(sum(aemit_mean_2_oc_au(113:168,7:86,:)))); % no scale
a8 = squeeze(sum(sum(aemit_mean_4_oc_au(113:168,7:86,:)))); % random forest

emis_se = [a5,a8,a6];
rd_rate_1 = (a6 - a8)./a6 * 100;
rd_rate_2 = (a6 - a5)./a6 * 100;
rd_se  = [rd_rate_2,rd_rate_1];
%%
figure1 = figure;
pos1 = [0.05 0.6 0.4 0.3];
subplot1 = subplot('Position',pos1);
bar1 = bar(emis_north,'EdgeColor','none','Parent',subplot1);
set(bar1(1),...
    'DisplayName','INJ-CLIM',...
    'FaceColor',[0 0.498039215803146 0]);
set(bar1(2),...
    'DisplayName','INJ-RF',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(bar1(3),'DisplayName','Total',...
     'FaceColor',[0 0.447058826684952 0.74117648601532]);

ylabel('OC fire emissions(Gg month^-^1)');

title('(a) Northern Australia','fontsize',24);

xlim(subplot1,[0 13]);
ylim(subplot1,[1 500]);
box(subplot1,'on');

set(subplot1,'FontSize',24,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'},...
    'YScale','log','YTick',[10 30 50 100 300 500]);
legend1 = legend(subplot1,'show');
set(legend1,'Location','northwest','EdgeColor',[1 1 1],'FontSize',24);
%
pos2 = [0.5 0.6 0.4 0.3];
subplot2 = subplot('Position',pos2);
bar1 = bar(emis_se,'EdgeColor','none','Parent',subplot2);
set(bar1(1),...
    'DisplayName','INJ-CLIM',...
    'FaceColor',[0 0.498039215803146 0]);
set(bar1(2),...
    'DisplayName','INJ-RF',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);
set(bar1(3),'DisplayName','Total',...
     'FaceColor',[0 0.447058826684952 0.74117648601532]);

title('(b) Southeastern Australia','fontsize',24);

xlim(subplot2,[0 13]);
ylim(subplot2,[1 500]);
% ylabel('Fire emissions of OC (Gg month^-^1)');

box(subplot2,'on');

set(subplot2,'FontSize',24,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'},...
    'YScale','log','YTick',[10 30 50 100 300 500]);
legend1 = legend(subplot2,'show');
set(legend1,'Location','northwest','EdgeColor',[1 1 1],'FontSize',24);

%
pos3 = [0.05 0.2 0.4 0.3];
subplot3 = subplot('Position',pos3);
bar1 = bar(rd_north,'EdgeColor','none','Parent',subplot3);
set(bar1(1),'DisplayName','INJ-CLIM',...
    'FaceColor',[0 0.498039215803146 0]);
set(bar1(2),...
    'DisplayName','INJ-RF',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);

% ylabel
ylabel('Injection percentages (%)');

% title
title({'(c) Northern Australia'},'fontsize',20);

% axis
xlim(subplot3,[0 13]);
ylim(subplot3,[0 60]);
box(subplot3,'on');
set(subplot3,'FontSize',24,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'});

% legend
legend1 = legend(subplot3,'show');
set(legend1,'Location','northwest','FontSize',24,'EdgeColor',[1 1 1]);

pos4 = [0.5 0.2 0.4 0.3];
subplot4 = subplot('Position',pos4);
bar1 = bar(rd_se,'EdgeColor','none','Parent',subplot4);
set(bar1(1),'DisplayName','INJ-CLIM',...
    'FaceColor',[0 0.498039215803146 0]);
set(bar1(2),...
    'DisplayName','INJ-RF',...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625]);

% title
title({'(d) Southeastern Australia'},'fontsize',20);

% axis
xlim(subplot4,[0 13]);
ylim(subplot4,[0 60]);
box(subplot4,'on');
set(subplot4,'FontSize',24,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
    {'2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020'});

% legend
legend1 = legend(subplot4,'show');
set(legend1,'Location','northwest','FontSize',24,'EdgeColor',[1 1 1]);