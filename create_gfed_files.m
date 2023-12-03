%% load emission factors for each species and fire type
prefix='GFED4.1s_';
% postfix='_beta.hdf5';
postfix='.hdf5';
% path='F:\Work_for_wildfires\GFED4\';
path='E:\work_for_wildfires\biomass_burning_emis\GFED\ancill\';
for years=2003:2016
    nyears=length(years);
    % load table of emission factors: unit: kg species / kg DM burned !!!!!
    input        = importdata([path 'GFED4_Emission_Factors.txt']);
    EF           = input.data; % unit: g species / kg DM burned
    EF           = EF*1e-3;    % convert unit; 
    [nspec,ncat] = size(EF);
    if ncat~=6; disp('Error: There should be six columns for different categories of emissions. Returning ...'); return; end
    % get names of chemical species in emissions table:
    species=input.textdata(size(input.textdata,1)-nspec+1:end,1);

    mask25=h5read([path prefix num2str(years) postfix],'/ancill/basis_regions');
    area25=h5read([path prefix num2str(years) postfix],'/ancill/grid_cell_area').*(mask25>0);
    [ny25,nx25]=size(mask25);

    GC_gfed_species = {'NOx','CO','C3H6O','MEK','C2H4O','C3H8','CH2O','C2H6','SO2','NH3','BC','OC','C10H16','C6H6','C7H8','C8H10','C2H5OH','CH3OH','Higher_Alkanes','PRPE'};
    GC_gfed_output  = {'NO','CO','ACET','MEK','ALD2','C3H8','CH2O','C2H6','SO2','NH3','BC','OC','MTPA','BENZ','TOLU','XYLE','EOH','MOH','ALK4','PRPE'};
    %%
    % set up arrays to store targetted arrays:
    emit=zeros(ny25,nx25,12); % total biomass burning emissions climatology[kg DM / m^2 / month]
    aemit_month=zeros(ny25,nx25,12,length(GC_gfed_species)); % mean annual total biomass burning emissions per scpecies [kg / m^2 / month]

    for year = years
        fname=sprintf('%s%4.4i%s',prefix,year,postfix);
        fprintf('Loading from %s ...\n',fname);
        % leap year
        if(( rem(year,100)~= 0 && rem(year,4) == 0) || (rem(year,100) == 0 && rem(year,400) == 0))
            ndays = 366;
            ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
        else
            ndays = 365;
            ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
        end
        % daily fraction & lon & lat
        daily_fraction = zeros(ny25,nx25,ndays);
        lon            = h5read([path fname],'/lon');
        lat            = h5read([path fname],'/lat');
        lon            = lon(:,1);
        lat            = flip(lat(1,:));
        % loop for months
        for m = 1:12
            DM     = h5read([path fname],sprintf('/emissions/%2.2i/DM',m));
            f_BORF = h5read([path fname],sprintf('/emissions/%2.2i/partitioning/DM_BORF',m));
            f_SAVA = h5read([path fname],sprintf('/emissions/%2.2i/partitioning/DM_SAVA',m));
            f_TEMF = h5read([path fname],sprintf('/emissions/%2.2i/partitioning/DM_TEMF',m));
            f_DEFO = h5read([path fname],sprintf('/emissions/%2.2i/partitioning/DM_DEFO',m));
            f_PEAT = h5read([path fname],sprintf('/emissions/%2.2i/partitioning/DM_PEAT',m));
            f_AGRI = h5read([path fname],sprintf('/emissions/%2.2i/partitioning/DM_AGRI',m));
            emit(:,:,m)         = DM;

            % loop for species
            for k = 1:length(GC_gfed_species)-2
                spc_name_inp = char(GC_gfed_species(k));
                disp(spc_name_inp);
                for ispc = 1:length(species)
                    if strcmp(spc_name_inp,char(species(ispc)))
                        index = ispc;
                    end
                end
                disp(['index in EF: ',num2str(index)]);
                % Calculate the monthly emissions
                aemit_month(:,:,m,k) = DM.*(EF(index,2)*f_BORF + EF(index,1)*f_SAVA + EF(index,3)*f_TEMF + EF(index,4)*f_DEFO + EF(index,5)*f_PEAT + EF(index,6)*f_AGRI);
            end
            % for special species ALK4 (High_Alkanes); Convert unit of EF from kgC to kg
            k = 19;
            convert_factor  = 58.12/( 4.3*12.);
            spc_name_inp = char(GC_gfed_species(k));
            % find index in EF matrix
            for ispc = 1:length(species)
                if strcmp(spc_name_inp,char(species(ispc)))
                    index = ispc;
                end
            end
            aemit_month(:,:,m,k) = DM.*(EF(index,2)*f_BORF + EF(index,1)*f_SAVA + EF(index,3)*f_TEMF + EF(index,4)*f_DEFO + EF(index,5)*f_PEAT + EF(index,6)*f_AGRI) * convert_factor;
            % for special species PRPE (C3H6 + Higher_Alkenes): Convert unit of EF (Higher_Alkenes) from kgC to kg
            k = 20;
            convert_factor       = 42.08/( 3.0*12.);
            ef_prpe              = EF(22,:) + EF(29,:)*convert_factor;
            aemit_month(:,:,m,k) = DM.*(ef_prpe(2)*f_BORF + ef_prpe(1)*f_SAVA + ef_prpe(3)*f_TEMF + ef_prpe(4)*f_DEFO + ef_prpe(5)*f_PEAT + ef_prpe(6)*f_AGRI);
        end
        % loop for species to create var
        for k = 1:length(GC_gfed_species)
            GCouput_name = char(GC_gfed_output(k));
            aemit_daily = zeros(ny25,nx25,ndays);
            for m = 1:12
                if m == 1
                    startday = 1;
                    endday = ndays_mon(m);
                else
                    startday = sum(ndays_mon(1:m-1))+1;
                    endday = sum(ndays_mon(1:m));
                end

                % load daily_fraction
                day2s  = 1/24/3600;
                for tt =startday:endday
                    daily_fraction(:,:,tt) = h5read([path fname],[sprintf('/emissions/%2.2i/daily_fraction/',m),'day_',num2str(tt-startday+1)]);
                    aemit_daily(:,:,tt) = aemit_month(:,:,m,k).* daily_fraction(:,:,tt)*day2s; % unit: kg/m2/s
                end
            end
            aemit_daily = flip(aemit_daily,2);
            eval([GCouput_name,'=','aemit_daily',';']);
            save_filename = [GCouput_name,'.mat'];
            save(['E:\work_for_wildfires\wrfgc\GFED4s_daily\',num2str(years),'\',save_filename],GCouput_name,'-v7.3');
            eval(['clearvars ',GCouput_name,';']);
            disp(['Finish ',GCouput_name]);
        end
    end
    %% Generate GFED emission files during the whole year
    location        = ['E:\work_for_wildfires\wrfgc\GFED4s_daily\',num2str(years),'\'];
    GC_gfed_output  = {'NO','CO','ACET','MEK','ALD2','C3H8','CH2O','C2H6','SO2','NH3','BC','OC','MTPA','BENZ','TOLU','XYLE','EOH','MOH','ALK4','PRPE'};

    year = years;
    fname=sprintf('%s%4.4i%s',prefix,year,postfix);
    fprintf('Loading from %s ...\n',fname);
    % leap year
    if(( rem(year,100)~= 0 && rem(year,4) == 0) || (rem(year,100) == 0 && rem(year,400) == 0))
        ndays = 366;
        ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
    else
        ndays = 365;
        ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
    end

    % lon & lat variables
    lon            = h5read([path fname],'/lon');
    lat            = h5read([path fname],'/lat');
    lon            = lon(:,1);
    lat            = flip(lat(1,:));
    % time variable
    time = zeros(ndays,1);
    iday = 0;
    for m = 1:12
        for dd = 1:ndays_mon(m)
            iday = iday + 1;
            t0 = datenum(1970,1,1,0,0,0);
            tm = datenum(years,m,dd,0,0,0);
            time(iday) = (tm - t0)*24;
        end
    end       
    % for m = 1:12
    %     if m < 10
    %         filename = ['GFAS_20190',num2str(m),'.nc'];
    %     else
    %         filename = ['GFAS_2019',num2str(m),'.nc'];
    %     end
    %     time_tmp = ncread(filename,'time');
    %     if m == 1
    %         startday = 1;
    %         endday = ndays_mon(m);
    %     else
    %         startday = sum(ndays_mon(1:m-1))+1;
    %         endday = sum(ndays_mon(1:m));
    %     end
    %     time(startday:endday) = time_tmp;
    % end
    time = int32(time);
    %%
    for m = 1:12
        if m < 10
            filename_output      =  ['E:\work_for_wildfires\wrfgc\GFED4s_daily\GFED4_daily_025x025.',num2str(years),'0',num2str(m),'.nc'];
        else
            filename_output      =  ['E:\work_for_wildfires\wrfgc\GFED4s_daily\GFED4_daily_025x025.',num2str(years),num2str(m),'.nc'];
        end

        if m == 1
            startday  = 1;
            endday    = ndays_mon(m);
        else
            startday  = sum(ndays_mon(1:m-1))+1;
            endday    = sum(ndays_mon(1:m));
        end
        clearvars time_mon;
        time_mon      = time(startday:endday);

        gfedSchema.Filename  =  filename_output;
        gfedSchema.Name      = '/';
        % Dimensions
        gfedSchema.Dimensions(1).Name      = 'lon';
        gfedSchema.Dimensions(1).Length    = length(lon);
        gfedSchema.Dimensions(1).Unlimited = false;
        gfedSchema.Dimensions(2).Name      = 'lat';
        gfedSchema.Dimensions(2).Length    = length(lat);
        gfedSchema.Dimensions(2).Unlimited = false;
        gfedSchema.Dimensions(3).Name      = 'time';
        gfedSchema.Dimensions(3).Length    = length(time_mon);
        gfedSchema.Dimensions(3).Unlimited = false;
        % Attributes (based on GFED4_dailyfrac_gen.025x025.201901.nc)
        gfedSchema.Attributes(1).Name      = '_NCProperties';
        gfedSchema.Attributes(1).Value     = 'version=1|netcdflibversion=4.6.1|hdf5libversion=1.10.2';
        gfedSchema.Attributes(2).Name      = 'history';
        gfedSchema.Attributes(2).Value     = 'Thu Apr 11 21:00:00 2022; GFED4_daily_025x025.2019.nc';
        gfedSchema.Attributes(3).Name      = 'NCO';
        gfedSchema.Attributes(3).Value     = 'netCDF Operators version 4.7.9 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)';
        % Groups
        gfedSchema.Groups                  = [];
        % Format                           
        gfedSchema.Format                  = 'netcdf4';
        ncwriteschema(filename_output,gfedSchema);

        % Variables -- time
        timeschema.Name                               = 'time';
        timeschema.Dimensions.Name                    = 'time';
        timeschema.Dimensions.Length                  = length(time_mon);
        timeschema.Dimensions.Unlimited               = false;
        timeschema.Size                               = length(time_mon);
        timeschema.Datatype                           = 'int32';
        timeschema.Attributes(1).Name                 = 'units';
        timeschema.Attributes(1).Value                = 'hours since 1970-01-01 00:00:00';
        timeschema.Attributes(2).Name                 = 'long_name';
        timeschema.Attributes(2).Value                = 'time';
        timeschema.Attributes(3).Name                 = 'calendar';
        timeschema.Attributes(3).Value                = 'gregorian';
        timeschema.ChunkSize                          = length(time_mon);
        timeschema.Shuffle                            = false;
        timeschema.DeflateLevel                       = 1;
        ncwriteschema(filename_output,timeschema);
        ncwrite(filename_output,'time',time_mon);
        % Variables -- lon
        lonschema.Name                               = 'lon';
        lonschema.Dimensions.Name                    = 'lon';
        lonschema.Dimensions.Length                  = length(lon);
        lonschema.Dimensions.Unlimited               = false;
        lonschema.Size                               = length(lon);
        lonschema.Datatype                           = 'single';
        lonschema.Attributes(1).Name                 = 'units';
        lonschema.Attributes(1).Value                = 'degrees_east';
        lonschema.Attributes(2).Name                 = 'long_name';
        lonschema.Attributes(2).Value                = 'longitude';
        lonschema.ChunkSize                          = length(lon);
        lonschema.Shuffle                            = false;
        lonschema.DeflateLevel                       = 1;
        ncwriteschema(filename_output,lonschema);
        ncwrite(filename_output,'lon',lon);
        % Variables -- lat
        latschema.Name                               = 'lat';
        latschema.Dimensions.Name                    = 'lat';
        latschema.Dimensions.Length                  = length(lat);
        latschema.Dimensions.Unlimited               = false;
        latschema.Size                               = length(lat);
        latschema.Datatype                           = 'single';
        latschema.Attributes(1).Name                 = 'units';
        latschema.Attributes(1).Value                = 'degrees_north';
        latschema.Attributes(2).Name                 = 'long_name';
        latschema.Attributes(2).Value                = 'latitude';
        latschema.ChunkSize                          = length(lat);
        latschema.Shuffle                            = false;
        latschema.DeflateLevel                       = 1;
        ncwriteschema(filename_output,latschema);
        ncwrite(filename_output,'lat',lat);

        % Variables -- fire emissions of each species
        for ii = 1:length(GC_gfed_output)
            spc_name   = char(GC_gfed_output(ii));
            load([location,spc_name,'.mat']);
            disp(spc_name);
            eval(['var = ',spc_name,';']);
            eval(['clearvars ',spc_name,';']);
            clearvars var_mon
            var_mon   = single(var(:,:,startday:endday));
            varschema.Name                            = spc_name;
            varschema.Dimensions(1).Name              = 'lon';
            varschema.Dimensions(1).Length            = length(lon);
            varschema.Dimensions(1).Unlimited         = false;
            varschema.Dimensions(2).Name              = 'lat';
            varschema.Dimensions(2).Length            = length(lat);
            varschema.Dimensions(2).Unlimited         = false;
            varschema.Dimensions(3).Name              = 'time';
            varschema.Dimensions(3).Length            = length(time_mon);
            varschema.Dimensions(3).Unlimited         = false;
            varschema.Size                            = [size(var_mon,1),size(var_mon,2),size(var_mon,3)];
            varschema.Datatype                        = 'single';
            varschema.Attributes(1).Name              = 'units';
            varschema.Attributes(1).Value             = 'kg/m2/s';
            varschema.Attributes(2).Name              = 'long_name';
            varschema.Attributes(2).Value             = ['Biomass burning emission flux of ',spc_name];
            varschema.ChunkSize                       = [size(var_mon,1),size(var_mon,2),1];
            varschema.Shuffle                         = false;
            varschema.DeflateLevel                    = 1;
            ncwriteschema(filename_output,varschema);
            ncwrite(filename_output,spc_name,var_mon);
        end
        aa = ncinfo(filename_output);
        ncdisp(filename_output);
    end
end



