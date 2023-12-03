filepath = 'E:\work_for_wildfires\injecton_height\MISR\plumeheights_website\MISR_mat\';
diroutput = dir(fullfile(filepath,'Plumes*'));
filename = {diroutput.name};

i = 0;
j = 0;
for ii = 1:length(diroutput)
    filename_tmp = char(filename(ii));
    if filename_tmp(end-6) == 'B'
        i = i + 1;
        blue_files{i} = filename_tmp;
    else
        j = j + 1;
        red_files{j}  = filename_tmp;
    end
end

%
inum = 0;
jnum = 0;
clearvars fp_all *_summary
altitude = 0:250:8000;
for ii = 1:length(blue_files)
    filename_1 = char(blue_files(ii));
    filename_2 = filename_1;
    filename_2(end-6) = 'R';
    red_exist = ismember(filename_2,red_files);
    
    if red_exist
        load([filepath,filename_1]);
        blue_quality = retrieval_quality;
        load([filepath,filename_2]);
        red_quality  = retrieval_quality;

        if blue_quality == 'Good'
            load([filepath,filename_1]);
            jnum = jnum + 1;
            filename_summary{jnum} = filename_1;
            file_exist = 1;
        elseif blue_quality == 'Fair' && red_quality ~= 'Good'
            load([filepath,filename_1]);
            jnum = jnum + 1;
            filename_summary{jnum} = filename_1;
            file_exist = 1;
        elseif blue_quality == 'Fair' && red_quality == 'Good'
            load([filepath,filename_2])
            jnum = jnum + 1;
            filename_summary{jnum} = filename_2;
            file_exist = 1;
        elseif blue_quality == 'Poor' && red_quality ~= 'Poor'
            load([filepath,filename_2]);
            jnum = jnum + 1;
            filename_summary{jnum} = filename_2;
            file_exist = 1;
        elseif blue_quality == 'Poor' && red_quality == 'Poor'
            file_exist = 0;
        end
        if file_exist
            inum = inum + 1;
            fp_all(inum,1) = str2double(fp_lon);
            fp_all(inum,2) = str2double(fp_lat);
            maxht_summary(inum,1) = str2double(max_ht);
            medianht_summary(inum,1) = str2double(median_ht);
            totalfrp_summary(inum,1) = str2double(total_frp);

            date_tmp = char(retrieval_date);
            date_summary(inum,1)  = str2double(date_tmp(1:4));
            date_summary(inum,2)  = str2double(date_tmp(6:7));
            date_summary(inum,3)  = str2double(date_tmp(9:10));
            
            time_tmp = char(retrieval_time);
            min_tmp  = str2double(time_tmp(4:5));
            hr_tmp   = str2double(time_tmp(1:2));
            if min_tmp > 30 && hr_tmp ~= 23
                time_summary(inum,1) = hr_tmp + 1;
            elseif min_tmp > 30 && hr_tmp == 23
                time_summary(inum,1) = 0;
            elseif min_tmp < 30
                time_summary(inum,1) = hr_tmp;
            end
   
            
            biome  = str2double(retrieval_biome);
            biome_summary(inum)   = biome(~isnan(biome));
            lonlat(inum,1)        = Plumes_data(1,2);
            lonlat(inum,2)        = Plumes_data(1,3);
            nowind_ht_tmp         = Plumes_data(:,10);
            windcor_ht_tmp        = Plumes_data(:,11);
            terr_elev_tmp         = Plumes_data(:,9);
            kmtop1_tmp            = Plumes_data(:,7);

            ik1  =  0;
            ik2  =  0;
            clearvars nowind_ht windcor_ht
            for kk = 1:length(nowind_ht_tmp)
                if nowind_ht_tmp(kk) - terr_elev_tmp(kk) > 0 && nowind_ht_tmp(kk) - terr_elev_tmp(kk) < 8000
                    ik1 = ik1+1;
                    nowind_ht(ik1) = nowind_ht_tmp(kk) - terr_elev_tmp(kk);
                end
                if windcor_ht_tmp(kk) - terr_elev_tmp(kk) > 0 && windcor_ht_tmp(kk) - terr_elev_tmp(kk) < 8000 
                    ik2 = ik2 + 1;
                    windcor_ht(ik2) = windcor_ht_tmp(kk) - terr_elev_tmp(kk) ;
                end
            end
            if ik1 > 0
                valid_pixels_summary(inum,1) = length(nowind_ht);
                max_ht_summary(inum,1) = max(nowind_ht);
                min_ht_summary(inum,1) = min(nowind_ht);
                mean_ht_summary(inum,1) = mean(nowind_ht);
                median_ht_summary(inum,1) = median(nowind_ht);
            else
                valid_pixels_summary(inum,1) = 0;
                max_ht_summary(inum,1) = NaN;
                min_ht_summary(inum,1) = NaN;
                mean_ht_summary(inum,1) = NaN;
                median_ht_summary(inum,1) = NaN;
            end
            if ik2 > 0
                valid_pixels_summary(inum,2) = length(windcor_ht);
                max_ht_summary(inum,2) = max(windcor_ht);
                min_ht_summary(inum,2) = min(windcor_ht);
                mean_ht_summary(inum,2) = mean(windcor_ht);
                median_ht_summary(inum,2) = median(windcor_ht);
            else
                valid_pixels_summary(inum,2) = 0;
                max_ht_summary(inum,2) = NaN;
                min_ht_summary(inum,2) = NaN;
                mean_ht_summary(inum,2) = NaN;
                median_ht_summary(inum,2) = NaN;
            end

            if ik1 > 0
                [n1,~]=histcounts(nowind_ht,altitude);
                [n3,~]=histcounts(nowind_ht,altitude,'Normalization','probability');
                histcounts_inij_nowind(inum,:)  = n1;
                perc_inij_nowind(inum,:)        = n3;
            else
                histcounts_inij_nowind(inum,:)  = NaN;
                perc_inij_nowind(inum,:)        = NaN;
            end

            if ik2 > 0
                [n2,~]=histcounts(windcor_ht,altitude);
                [n4,~]=histcounts(windcor_ht,altitude,'Normalization','probability');
                histcounts_inij_windcor(inum,:)  = n2;
                perc_inij_windcor(inum,:) = n4;
            else
                histcounts_inij_windcor(inum,:)  = NaN;
                perc_inij_windcor(inum,:)        = NaN;
            end

        end
    else 
        load([filepath,filename_1]);
        blue_quality = retrieval_quality;
         if blue_quality == 'Good'|| blue_quality == 'Fair'
            load([filepath,filename_1]);
            jnum = jnum + 1;
            filename_summary{jnum} = filename_1;
            file_exist = 1;
        elseif blue_quality == 'Poor'
            file_exist = 0;
         end
        
        if file_exist
            inum = inum + 1;
            fp_all(inum,1) = str2double(fp_lon);
            fp_all(inum,2) = str2double(fp_lat);
            maxht_summary(inum,1) = str2double(max_ht);
            medianht_summary(inum,1) = str2double(median_ht);
            totalfrp_summary(inum,1) = str2double(total_frp);

            date_tmp = char(retrieval_date);
            date_summary(inum,1)  = str2double(date_tmp(1:4));
            date_summary(inum,2)  = str2double(date_tmp(6:7));
            date_summary(inum,3)  = str2double(date_tmp(9:10));
            biome  = str2double(retrieval_biome);
            biome_summary(inum)   = biome(~isnan(biome));
            lonlat(inum,1)        = Plumes_data(1,2);
            lonlat(inum,2)        = Plumes_data(1,3);
            nowind_ht_tmp         = Plumes_data(:,10);
            windcor_ht_tmp        = Plumes_data(:,11);
            terr_elev_tmp         = Plumes_data(:,9);
            kmtop1_tmp            = Plumes_data(:,7);

            ik1  =  0;
            ik2  =  0;
            clearvars nowind_ht windcor_ht
            for kk = 1:length(nowind_ht_tmp)
                if nowind_ht_tmp(kk) - terr_elev_tmp(kk) > 0 && nowind_ht_tmp(kk) - terr_elev_tmp(kk) < 8000
                    ik1 = ik1+1;
                    nowind_ht(ik1) = nowind_ht_tmp(kk) - terr_elev_tmp(kk);
                end
                if windcor_ht_tmp(kk) - terr_elev_tmp(kk) > 0 && windcor_ht_tmp(kk) - terr_elev_tmp(kk) < 8000 
                    ik2 = ik2 + 1;
                    windcor_ht(ik2) = windcor_ht_tmp(kk) - terr_elev_tmp(kk) ;
                end
            end
            if ik1 > 0
                valid_pixels_summary(inum,1) = length(nowind_ht);
                max_ht_summary(inum,1) = max(nowind_ht);
                min_ht_summary(inum,1) = min(nowind_ht);
                mean_ht_summary(inum,1) = mean(nowind_ht);
                median_ht_summary(inum,1) = median(nowind_ht);
            else
                valid_pixels_summary(inum,1) = 0;
                max_ht_summary(inum,1) = NaN;
                min_ht_summary(inum,1) = NaN;
                mean_ht_summary(inum,1) = NaN;
                median_ht_summary(inum,1) = NaN;
            end
            if ik2 > 0
                valid_pixels_summary(inum,2) = length(windcor_ht);
                max_ht_summary(inum,2) = max(windcor_ht);
                min_ht_summary(inum,2) = min(windcor_ht);
                mean_ht_summary(inum,2) = mean(windcor_ht);
                median_ht_summary(inum,2) = median(windcor_ht);
            else
                valid_pixels_summary(inum,2) = 0;
                max_ht_summary(inum,2) = NaN;
                min_ht_summary(inum,2) = NaN;
                mean_ht_summary(inum,2) = NaN;
                median_ht_summary(inum,2) = NaN;
            end

            if ik1 > 0
                [n1,~]=histcounts(nowind_ht,altitude);
                [n3,~]=histcounts(nowind_ht,altitude,'Normalization','probability');
                histcounts_inij_nowind(inum,:)  = n1;
                perc_inij_nowind(inum,:)        = n3;
            else
                histcounts_inij_nowind(inum,:)  = NaN;
                perc_inij_nowind(inum,:)        = NaN;
            end

            if ik2 > 0
                [n2,~]=histcounts(windcor_ht,altitude);
                [n4,~]=histcounts(windcor_ht,altitude,'Normalization','probability');
                histcounts_inij_windcor(inum,:)  = n2;
                perc_inij_windcor(inum,:) = n4;
            else
                histcounts_inij_windcor(inum,:)  = NaN;
                perc_inij_windcor(inum,:)        = NaN;
            end

        end
    end
end

%%
load('E:\work_for_wildfires\injecton_height\MERRA2_PBLH\lonlat_au_merra2.mat');
plumes_number = length(date_summary);
ratio_uppblh_nowind = zeros(plumes_number,1);
ratio_uppblh_windcor = zeros(plumes_number,1);

alt_prof  = 0:250:8000;

for ii = 1:plumes_number
    yr_tmp    = date_summary(ii,1);
    mon_tmp   = date_summary(ii,2);
    day_tmp   = date_summary(ii,3);
    hr_tmp    = time_summary(ii);
    filename =  ['E:\work_for_wildfires\injecton_height\MERRA2_PBLH\merra2_',num2str(yr_tmp),num2str(mon_tmp),'_pblh_au.mat'];
    load(filename);
    var_mean = zeros(size(var_au,1),size(var_au,2),size(var_au,3)/24);
    for tt = 1:size(var_au,3)/24
        var_mean(:,:,tt) = mean(var_au(:,:,(tt-1)*24+1:24*tt),3);
    end
    
    lon_tmp    = fp_all(ii,1);
    lat_tmp    = fp_all(ii,2);
    aa = find( abs(lon_tmp - lon_au) == min(abs(lon_tmp - lon_au)));
    bb = find( abs(lat_tmp - lat_au) == min(abs(lat_tmp - lat_au)));
%     pblh_tmp   = var_au(aa(1),bb(1),(day_tmp-1)*24 + hr_tmp+1);
    pblh_tmp   = var_mean(aa(1),bb(1),day_tmp);
    zz = find( alt_prof > pblh_tmp, 1);
    prof_nowind_tmp = perc_inij_nowind(ii,:);
    prof_windcor_tmp = perc_inij_windcor(ii,:);
    
    ratio_uppblh_nowind(ii)  = sum(prof_nowind_tmp(zz-1:end));
    ratio_uppblh_windcor(ii) = sum(prof_windcor_tmp(zz-1:end));
end

