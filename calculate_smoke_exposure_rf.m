%% set some parameters
for years       = 2019
    start_mon       = 4;
    end_mon         = 12;
    start_date      = 1;
    end_date        = 30;
    backward_hours  = 120;
    rundays         = datenum(years,end_mon,end_date) - datenum(years,start_mon,start_date) + 1;
    receptors       = 4;

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
    postfix = '.nc';

    spc_name = {'BC','OC'};

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
    aemit_daily_4 = aemit_daily; % for scale emission random forest
    lon_global = lon_global_gfed;
    lat_global = lat_global_gfed;

    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'_gfs_nohnf_receptor_north\'];
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

    %% Step 2: use predicted injection percentages from RF to scale the biomass burning emissions
    disp('Step 2 - emis 1');
    clearvars gfed_all inheight_all pblh_all
    load(['E:\work_for_wildfires\random_forest\predicted_injection_percentages_rf\new\predict_output_au_',num2str(years),'.mat']);
    eval(['ratio_all','=','predict_output_au_',num2str(years),';']);
    for ispc = 1:2
         gfed_var_new(:,:,:,ispc) = (1 - ratio_all(:,:,tt_indx1:tt_indx2)) .* gfed_var(:,:,:,ispc);
    end
    aemit_daily_4(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),tt_indx1:tt_indx2,:) = gfed_var_new;
    %% Step 3: Read and calculate the daily footprints
    % STILT output files
    disp('Step 3');
    filepath = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'_gfs_nohnf_receptor_north\'];
    diroutput = dir(fullfile(filepath,'*_foot.nc'));
    filename = {diroutput.name};

    filename_1 = char(filename(1));
    lon = ncread([filepath,filename_1],'lon');
    lat = ncread([filepath,filename_1],'lat');
    xdim = length(lon);
    ydim = length(lat);

    smoke_exp4    = zeros(xdim,ydim,rundays,receptors,2);
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

        for ispc = 1:2
        % read the fire emissions & exchange the dimensions [unit: umol/m2/s]
        fireemis_var4 = aemit_daily_4(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),tt_indx1:tt_indx2,ispc);

        % unit conversion for trace gases and aerosols
        % For trace gases: kg m^-2 s^-1 --> umol m^-2 s^-1
        % mw_spc         = 28;
        kg2g             = 1e3;
        % convert_factor_gas = kg2g/mw_spc*1e6;

        % For aerosols: ppm/(kg m^-2 s^-1) --> ug m^-3/(umol m^-2 s^-1)
        mw_air         = 29;
        rho_air        = 1.29; % dry air density at STP [kg m^-3]
        convert_factor = kg2g * 1e6 * rho_air/mw_air * 1e3;
        smoke_exp_tmp4 = fireemis_var4 .* footprint_daily * convert_factor;

        smoke_exp4(:,:,iday,irec,ispc) = nansum(smoke_exp_tmp4,3);
        end
    end
    smoke_exp = squeeze(smoke_exp4(:,:,:,:,1));
%     filename_output1 = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'_gfs_nohnf_receptor_north\smoke_exp_output\smoke_exp_bc_gfed_rf_new.mat'];
%     save(filename_output1,'smoke_exp');

    clearvars smoke_exp
    smoke_exp = squeeze(smoke_exp4(:,:,:,:,2));
%     filename_output1 = ['F:\stilt_output\GDAS05\footprints_',num2str(years),'_gfs_nohnf_receptor_north\smoke_exp_output\smoke_exp_oc_gfed_rf_new.mat'];
%     save(filename_output1,'smoke_exp');
end