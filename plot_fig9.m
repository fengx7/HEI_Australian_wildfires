mon_name = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
start_mon = 4;
end_mon   = 12;
receptors = 4;
rundays   = 275;
yr_num    = 12;
loc_str = {'130.94853_-12.50779','146.8257_-19.2542','149.1549_-21.1595','151.2704_-23.8627'};
loc_lon = [130.94853,146.8257,149.1549,151.2704];
loc_lat = [-12.50779,-19.2542,-21.1595,-23.8627];

iyear = 2011;
filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_north\'];
prefix_foot = [num2str(iyear),'*'];
diroutput = dir(fullfile(filepath,prefix_foot)); % *_foot.nc 
filename = {diroutput.name};
filename_1 = char(filename(1));
lon = ncread([filepath,filename_1],'lon');
lat = ncread([filepath,filename_1],'lat');
xdim = length(lon);
ydim = length(lat);

footprints_sum = zeros(xdim,ydim,rundays*yr_num,receptors);
for iyear = 2009:2020
    for years = iyear
        % leap year
        if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
            ndays = 366;
            ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
        else
            ndays = 365;
            ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
        end
    end
    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_north\'];
    load([filepath,'footprints_daily.mat']);
    footprints_sum(:,:,(iyear-2008-1)*rundays+1:(iyear-2008)*rundays,:) = footprints_exp;  
end
footprints_mean = nanmean(footprints_sum,3);
%%
irec = 1;
var1 = squeeze(footprints_mean(:,:,:,irec));

figure1 = figure;
subplot1 = subplot(2,3,1,'BoxStyle','full');
m_proj('mercator','lon',[min(lon),max(lon)],'lat',[min(lat),max(lat)]);
m_pcolor(lon,lat,log10(var1'));
hold on
m_plot(loc_lon(irec),loc_lat(irec),'ko','linewidth',2,'MarkerSize',10);
m_grid('Fontsize',20);
m_plotbndry('C:\Program Files\MATLAB\R2017b\m_map\world','color','k','linewidth',1.5);
shading flat
hold on
climax = 0;
climin = -5;
set(gca,'clim',[climin,climax]);
color_map = othercolor('Reds7',5);
color_map(1,:) = 1;
colormap(color_map);
hcol = colorbar;
set(hcol,'FontSize',24);
set(hcol,'TickLabels',{'10^-^5','10^-^4','10^-^3','10^-^2','10^-^1','10^0'});
title('(a) Darwin','FontSize',24);
set(gca,'FontSize',24);

irec = 4;
var2 = squeeze(footprints_mean(:,:,:,irec));
subplot2 = subplot(2,3,2,'BoxStyle','full');
m_proj('mercator','lon',[min(lon),max(lon)],'lat',[min(lat),max(lat)]);
m_pcolor(lon,lat,log10(var2'));
hold on
m_plot(loc_lon(irec),loc_lat(irec),'ko','linewidth',2,'MarkerSize',10);
m_grid('Fontsize',20);
m_plotbndry('C:\Program Files\MATLAB\R2017b\m_map\world','color','k','linewidth',1.5);
shading flat
hold on
climax = 0;
climin = -5;
set(gca,'clim',[climin,climax]);
color_map = othercolor('Reds7',5);
color_map(1,:) = 1;
colormap(color_map);
hcol = colorbar;
set(hcol,'FontSize',24);
set(hcol,'TickLabels',{'10^-^5','10^-^4','10^-^3','10^-^2','10^-^1','10^0'});
title('(b) Gladstone','FontSize',24);
set(gca,'FontSize',24);

%%
mon_name = {'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'};
start_mon = 8;
end_mon   = 12;
receptors = 11;
rundays   = 153+31;

% 1. Footscray; 2. Albury; 3. Florey; 4. Wollongong; 
% 5. Prospect; 6. Newcastle; 7. Mountain Creek; 8. Springwood; 9. Alphington; 10. Liverpool; 11. Wallsend; 
loc_str = {'144.8728027_-37.80487823','146.93986_-36.05182','149.043539_-35.220606','150.88733_-34.41706','150.91417_-33.79424','151.75965_-32.9312','153.1038_-26.6917','153.1356_-27.6125','145.0306_-37.7784','150.9058_-33.9328','151.6692_-32.8961'};
loc_lon = [144.8728027,146.93986,149.043539,150.88733,150.91417,151.75965,153.1038,153.1356,145.0306,150.9058,151.6692];
loc_lat = [-37.80487823,-36.05182,-35.220606,-34.41706,-33.79424,-32.9312,-26.6917,-27.6125,-37.7784,-33.9328,-32.8961];

iyear = 2011;
filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_se\'];
prefix_foot = [num2str(iyear),'*'];
diroutput = dir(fullfile(filepath,prefix_foot)); % *_foot.nc 
filename = {diroutput.name};
filename_1 = char(filename(1));
lon = ncread([filepath,filename_1],'lon');
lat = ncread([filepath,filename_1],'lat');
xdim = length(lon);
ydim = length(lat);

footprints_sum = zeros(xdim,ydim,rundays*yr_num,receptors);
for iyear = 2009:2020
    for years = iyear
        % leap year
        if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
            ndays = 366;
            ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
        else
            ndays = 365;
            ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
        end
    end
    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear),'_gfs_nohnf_receptor_se\'];
    load([filepath,'footprints_daily.mat']);
    footprints_exp1 = footprints_exp;
    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(iyear+1),'01_gfs_nohnf_receptor_se\'];
    load([filepath,'footprints_daily.mat']);
    footprints_exp2 = footprints_exp;
    footprints_exp = cat(3,footprints_exp1,footprints_exp2);
    footprints_sum(:,:,(iyear-2008-1)*rundays+1:(iyear-2008)*rundays,:) = footprints_exp;  
end
footprints_mean = nanmean(footprints_sum,3);

irec = 8;
var3 = squeeze(footprints_mean(:,:,:,irec));

subplot3 = subplot(2,3,3,'BoxStyle','full');
m_proj('mercator','lon',[min(lon),max(lon)],'lat',[min(lat),max(lat)]);
m_pcolor(lon,lat,log10(var3'));
hold on
m_plot(loc_lon(irec),loc_lat(irec),'ko','linewidth',2,'MarkerSize',10);
m_grid('Fontsize',20);
m_plotbndry('C:\Program Files\MATLAB\R2017b\m_map\world','color','k','linewidth',1.5);
shading flat
hold on
climax = 0;
climin = -5;
set(gca,'clim',[climin,climax]);
color_map = othercolor('Reds7',5);
color_map(1,:) = 1;
colormap(color_map);
hcol = colorbar;
set(hcol,'FontSize',24);
set(hcol,'TickLabels',{'10^-^5','10^-^4','10^-^3','10^-^2','10^-^1','10^0'});
title('(c) Brisbane','FontSize',24);
set(gca,'FontSize',24);

irec = 11;
var4 = squeeze(footprints_mean(:,:,:,irec));
subplot4 = subplot(2,3,4,'BoxStyle','full');
m_proj('mercator','lon',[min(lon),max(lon)],'lat',[min(lat),max(lat)]);
m_pcolor(lon,lat,log10(var4'));
hold on
m_plot(loc_lon(irec),loc_lat(irec),'ko','linewidth',2,'MarkerSize',10);
m_grid('Fontsize',20);
m_plotbndry('C:\Program Files\MATLAB\R2017b\m_map\world','color','k','linewidth',1.5);
shading flat
hold on
climax = 0;
climin = -5;
set(gca,'clim',[climin,climax]);
color_map = othercolor('Reds7',5);
color_map(1,:) = 1;
colormap(color_map);
hcol = colorbar;
set(hcol,'FontSize',24);
set(hcol,'TickLabels',{'10^-^5','10^-^4','10^-^3','10^-^2','10^-^1','10^0'});
title('(d) Newcastle','FontSize',24);
set(gca,'FontSize',24);

irec = 10;
var5 = squeeze(footprints_mean(:,:,:,irec));
subplot5 = subplot(2,3,5,'BoxStyle','full');
m_proj('mercator','lon',[min(lon),max(lon)],'lat',[min(lat),max(lat)]);
m_pcolor(lon,lat,log10(var5'));
hold on
m_plot(loc_lon(irec),loc_lat(irec),'ko','linewidth',2,'MarkerSize',10);
m_grid('Fontsize',20);
m_plotbndry('C:\Program Files\MATLAB\R2017b\m_map\world','color','k','linewidth',1.5);
shading flat
hold on
climax = 0;
climin = -5;
set(gca,'clim',[climin,climax]);
color_map = othercolor('Reds7',5);
color_map(1,:) = 1;
colormap(color_map);
hcol = colorbar;
set(hcol,'FontSize',24);
set(hcol,'TickLabels',{'10^-^5','10^-^4','10^-^3','10^-^2','10^-^1','10^0'});
title('(e) Sydney','FontSize',24);
set(gca,'FontSize',24);

irec = 9;
var6 = squeeze(footprints_mean(:,:,:,irec));
subplot6 = subplot(2,3,6,'BoxStyle','full');
m_proj('mercator','lon',[min(lon),max(lon)],'lat',[min(lat),max(lat)]);
m_pcolor(lon,lat,log10(var6'));
hold on
m_plot(loc_lon(irec),loc_lat(irec),'ko','linewidth',2,'MarkerSize',10);
m_grid('Fontsize',20);
m_plotbndry('C:\Program Files\MATLAB\R2017b\m_map\world','color','k','linewidth',1.5);
shading flat
hold on
climax = 0;
climin = -5;
set(gca,'clim',[climin,climax]);
color_map = othercolor('Reds7',5);
color_map(1,:) = 1;
colormap(color_map);
hcol = colorbar;
set(hcol,'FontSize',24);
set(hcol,'TickLabels',{'10^-^5','10^-^4','10^-^3','10^-^2','10^-^1','10^0'});
title('(f) Melbourne','FontSize',24);
set(gca,'FontSize',24);