%% Plotting script for Figure 4
load('E:\work_for_wildfires\random_forest\emis_scale_apr_jan_2009_2020.mat');
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
filepath = 'F:\stilt_output\GDAS05\footprints_2018_gfs_nohnf_receptor_north\';
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
lon_gfed_au = lon_global_gfed(lon_indx1(1):lon_indx2(1));
lat_gfed_au = lat_global_gfed(lat_indx1(1):lat_indx2(1));

% Convert unit to Gg month-1
sf   = 3600 * 24 / 1e6 /10;
for ii = 1:size(aemit_mean_1_oc,3)
    oc_1_emis(:,:,ii) = aemit_mean_1_oc(:,:,ii) .* area * sf;  % unit: Gg month-1
    oc_2_emis(:,:,ii) = aemit_mean_2_oc(:,:,ii) .* area * sf;
    oc_3_emis(:,:,ii) = aemit_mean_3_oc(:,:,ii) .* area * sf;
    oc_4_emis(:,:,ii) = aemit_mean_4_oc(:,:,ii) .* area * sf;
end

aemit_mean_1_oc_au  = oc_1_emis(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);
aemit_mean_2_oc_au  = oc_2_emis(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);
aemit_mean_3_oc_au  = oc_3_emis(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);
aemit_mean_4_oc_au  = oc_4_emis(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);

mean_1_oc_au        = squeeze(mean(aemit_mean_1_oc_au,3));
mean_2_oc_au        = squeeze(mean(aemit_mean_2_oc_au,3));
mean_4_oc_au        = squeeze(mean(aemit_mean_4_oc_au,3));
diff_perc1          = (mean_2_oc_au - mean_4_oc_au)./mean_2_oc_au*100;
diff_perc2          = (mean_2_oc_au - mean_1_oc_au)./mean_2_oc_au*100;
diff_perc           = (mean_1_oc_au - mean_4_oc_au)./mean_2_oc_au*100;

% North Australia
a1 = squeeze(sum(sum(aemit_mean_1_oc_au(13:168,107:142,:)))); % for scale emission mami
a2 = squeeze(sum(sum(aemit_mean_2_oc_au(13:168,107:142,:)))); % no scale
a3 = squeeze(sum(sum(aemit_mean_3_oc_au(13:168,107:142,:)))); % for scale emission mami & fix non-injection height
a4 = squeeze(sum(sum(aemit_mean_4_oc_au(13:168,107:142,:)))); % random forest

emis_north = [a4,a1,a2];
rd_rate_1 = (a2 - a4)./a2 * 100;
rd_rate_2 = (a2 - a1)./a2 * 100;
rd_north  = [rd_rate_1,rd_rate_2];
%%
for ii = 1:length(lon_gfed_au)
    for jj = 1:length(lat_gfed_au)
        if mean_4_oc_au(ii,jj) == 0
            mean_4_oc_au(ii,jj) = NaN;
        end
        if mean_1_oc_au(ii,jj) == 0
            mean_1_oc_au(ii,jj) = NaN;
        end
        if mean_2_oc_au(ii,jj) == 0
            mean_2_oc_au(ii,jj) = NaN;
        end
    end
end

% Plot for injection fraction (INJ-RF)
figure1 = figure;
pos1 = [0.3 0.3 0.4 0.4];
subplot1 = subplot('Position',pos1,'Parent',figure1);
m_proj('mercator','lon',[min(lon_gfed_au),max(lon_gfed_au)-10],'lat',[min(lat_gfed_au),max(lat_gfed_au)]);
m_pcolor(lon_gfed_au,lat_gfed_au,diff_perc1');
m_grid('linestyle','none','tickdir','in','fontsize',20,'FontWeight','bold');
m_plotbndry('C:\Program Files\MATLAB\R2017b\m_map\world','color','k','linewidth',1.5);
hold on
shading flat
climax = 60;
climin = 0;
set(gca,'clim',[climin,climax]);
cmap_data = flipud(othercolor('RdBu9',10));
cmap_data(5:6,:) = 1;
colormap(cmap_data(6:end,:));
hcol = colorbar;
set(hcol,'YTick',climin:(climax-climin)/5:climax);
set(hcol,'FontSize',26);
title('(b) Fraction above the PBL (%)','FontSize',26);
hold on

% Add text
annotation(figure1,'textbox',...
    [0.380885118194842 0.2370712469625 0.198 0.114],...
    'String',{'INJ-RF'},...
    'LineStyle','none',...
    'FontSize',26,...
    'FitBoxToText','off');
	
% Plot for injection fraction (INJ-CLIM)
figure1 = figure;
subplot5 = subplot('Position',pos1,'Parent',figure1);
m_proj('mercator','lon',[min(lon_gfed_au),max(lon_gfed_au)-10],'lat',[min(lat_gfed_au),max(lat_gfed_au)]);
m_pcolor(lon_gfed_au,lat_gfed_au,diff_perc2');
m_grid('linestyle','none','tickdir','in','fontsize',20,'FontWeight','bold');
m_plotbndry('C:\Program Files\MATLAB\R2017b\m_map\world','color','k','linewidth',1.5);
hold on
shading flat
climax = 60;
climin = 0;
set(gca,'clim',[climin,climax]);
cmap_data = flipud(othercolor('RdBu9',10));
cmap_data(5:6,:) = 1;
colormap(cmap_data(6:end,:));
hcol = colorbar;
set(hcol,'YTick',climin:(climax-climin)/5:climax);
set(hcol,'FontSize',26);
title('(a) Fraction above the PBL (%)','FontSize',26);

% Add text
annotation(figure1,'textbox',...
    [0.380885118194842 0.2370712469625 0.198 0.114],...
    'String',{'INJ-CLIM'},...
    'LineStyle','none',...
    'FontSize',26,...
    'FitBoxToText','off');
	
% Plot for relative difference between INJ-CLIM and INJ-RF
figure1 = figure;
subplot6 = subplot('Position',pos1,'Parent',figure1);
m_proj('mercator','lon',[min(lon_gfed_au),max(lon_gfed_au)-10],'lat',[min(lat_gfed_au),max(lat_gfed_au)]);
m_pcolor(lon_gfed_au,lat_gfed_au,diff_perc');
m_grid('linestyle','none','tickdir','in','fontsize',20,'FontWeight','bold');
m_plotbndry('C:\Program Files\MATLAB\R2017b\m_map\world','color','k','linewidth',1.5);
hold on
shading flat
climax = 40;
climin = -40;
set(gca,'clim',[climin,climax]);
cmap_data = flipud(othercolor('RdBu9',10));
cmap_data(5:6,:) = 1;
colormap(cmap_data);
hcol = colorbar;
set(hcol,'YTick',climin:(climax-climin)/5:climax);
set(hcol,'FontSize',26);
title('(c) Fraction in PBL, relative difference (%)','FontSize',26);

% Add text
annotation(figure1,'textbox',...
    [0.380885118194842 0.2370712469625 0.198 0.114],...
    'String',{'INJ-CLIM vs. INJ-RF'},...
    'LineStyle','none',...
    'FontSize',26,...
    'FitBoxToText','off');

