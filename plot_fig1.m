%% Plotting script for Figure 1
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
    oc_2_emis(:,:,ii) = aemit_mean_2_oc(:,:,ii) .* area * sf;
end

aemit_mean_2_oc_au  = oc_2_emis(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:);
mean_2_oc_au        = squeeze(mean(aemit_mean_2_oc_au,3));
%%
for ii = 1:length(lon_gfed_au)
    for jj = 1:length(lat_gfed_au)
        if mean_2_oc_au(ii,jj) == 0
            mean_2_oc_au(ii,jj) = NaN;
        end
    end
end

figure1 = figure;
pos1 = [0.3 0.3 0.4 0.4];
subplot1 = subplot('Position',pos1,'Parent',figure1);
m_proj('mercator','lon',[min(lon_gfed_au)-3,max(lon_gfed_au)-10],'lat',[min(lat_gfed_au),max(lat_gfed_au)]);
m_pcolor(lon_gfed_au,lat_gfed_au,log10(mean_2_oc_au'));
hold on
m_line([118.125 118.125],[-18.875, -10.125],'color','k','linewidth',1.5,'linestyle','-.');
m_line([150.875 150.875],[-18.875, -10.125],'color','k','linewidth',1.5,'linestyle','-.');
m_line([118.125 150.875],[-18.875, -18.875],'color','k','linewidth',1.5,'linestyle','-.');
m_line([118.125 150.875],[-10.125, -10.125],'color','k','linewidth',1.5,'linestyle','-.');

m_line([140.125 140.125],[-43.875, -24.125],'color','k','linewidth',1.5,'linestyle','-.');
m_line([153.875 153.875],[-43.875, -24.125],'color','k','linewidth',1.5,'linestyle','-.');
m_line([140.125 153.875],[-43.875, -43.875],'color','k','linewidth',1.5,'linestyle','-.');
m_line([140.125 153.875],[-24.125, -24.125],'color','k','linewidth',1.5,'linestyle','-.');

m_grid('linestyle','none','tickdir','in','fontsize',10,'FontWeight','bold');
m_plotbndry('C:\Program Files\MATLAB\R2017b\m_map\world','color','k','linewidth',1.5);
hold on
shading flat
climax = 1;
climin = -4;
set(gca,'clim',[climin,climax]);
colormap(subplot1,flipud(othercolor('Spectral10',10)));
hcol = colorbar;
set(hcol,'YTick',climin:(climax-climin)/5:climax);
set(hcol,'YTicklabel',{'10^-^4','10^-^3','10^-^2','10^-^1','10^0','10^1'});
set(hcol,'FontSize',24);
title('Mean total OC fire emissions (Gg month^-^1)','FontSize',24);