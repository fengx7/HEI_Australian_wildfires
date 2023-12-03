% Set up some parameters
path = 'E:\work_for_wildfires\injecton_height\MERRA2\';
prefix_name = 'MERRA2.';
postfix_name = '.A1.05x0625.nc4';

startyr        = 2021;
endyr          = 2021;
startmon       = 1;
endmon         = 1;
startdate      = 1; 
enddate        = 31;
num_days       = datenum(endyr,endmon,enddate) - datenum(startyr,startmon,startdate) + 1;
datestd        = zeros(num_days,6);
i = 0;
for mm = datenum(startyr,startmon,startdate):datenum(endyr,endmon,enddate)
    i = i + 1;
        datestd(i,:) = datevec(mm);
end
formatOut2     = 'yyyymmdd';      % For create 'date' variable
date_str1      = datestr(datestd,formatOut2);

filename1 = [path,prefix_name,date_str1(1,:),postfix_name];

spc_name = 'V10M';
var1 = ncread(filename1,spc_name);
lon  = ncread(filename1,'lon');
lat  = ncread(filename1,'lat');

xdim = size(var1,1);
ydim = size(var1,2);
tdim = size(var1,3)*31;
var = zeros(xdim,ydim,tdim);

for tt = 1:num_days
    filename_tmp = [path,prefix_name,date_str1(tt,:),postfix_name];
    disp(filename_tmp);
    var(:,:,(tt-1)*24+1:24*tt) = ncread(filename_tmp,spc_name);
end

lon_indx1 = find(lon==113.1250);
lon_indx2 = find(lon==153.7500);

lat_indx1 = find(lat==-44.0000);
lat_indx2 = find(lat==-10.0000);

var_au = var(lon_indx1:lon_indx2,lat_indx1:lat_indx2,:);
lon_au = lon(lon_indx1:lon_indx2);
lat_au = lat(lat_indx1:lat_indx2);

save(['E:\work_for_wildfires\injecton_height\MERRA2_PBLH\merra2_202101_v10m_au.mat'],'var_au','-v7.3');
%%
path = 'E:\work_for_wildfires\injecton_height\MERRA2\';
prefix_name = 'MERRA2.';
postfix_name = '.A3dyn.05x0625.nc4';

startyr        = 2021;
endyr          = 2021;
startmon       = 1;
endmon         = 1;
startdate      = 1; 
enddate        = 31;
num_days       = datenum(endyr,endmon,enddate) - datenum(startyr,startmon,startdate) + 1;
datestd        = zeros(num_days,6);
i = 0;
for mm = datenum(startyr,startmon,startdate):datenum(endyr,endmon,enddate)
    i = i + 1;
        datestd(i,:) = datevec(mm);
end
formatOut2     = 'yyyymmdd';      % For create 'date' variable
date_str1      = datestr(datestd,formatOut2);

filename1 = [path,prefix_name,date_str1(1,:),postfix_name];

spc_name = 'RH';
var1 = ncread(filename1,spc_name);
lon  = ncread(filename1,'lon');
lat  = ncread(filename1,'lat');

xdim = size(var1,1);
ydim = size(var1,2);
tdim = size(var1,4)*31;
var = zeros(xdim,ydim,tdim);

for tt = 1:num_days
    filename_tmp = [path,prefix_name,date_str1(tt,:),postfix_name];
    disp(filename_tmp);
    tmp = ncread(filename_tmp,spc_name);
    var(:,:,(tt-1)*8+1:8*tt) = squeeze(tmp(:,:,1,:));
end

lon_indx1 = find(lon==113.1250);
lon_indx2 = find(lon==153.7500);

lat_indx1 = find(lat==-44.0000);
lat_indx2 = find(lat==-10.0000);

var_au = var(lon_indx1:lon_indx2,lat_indx1:lat_indx2,:);
lon_au = lon(lon_indx1:lon_indx2);
lat_au = lat(lat_indx1:lat_indx2);

save(['E:\work_for_wildfires\injecton_height\MERRA2_PBLH\merra2_20211_rh_au.mat'],'var_au','-v7.3');