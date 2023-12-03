%% Use other input data to predict the plume injection fractions
for years  = 2021
    disp(['year:',num2str(years)]);
    if(( rem(years,100)~= 0 && rem(years,4) == 0) || (rem(years,100) == 0 && rem(years,400) == 0))
        ndays = 366;
        ndays_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
    else
        ndays = 365;
        ndays_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
    end
    
    load('C:\Users\fox\Desktop\demo_HEI\random_forest\predicting_data\input_factors_lonlat.mat');
    filepath = ['C:\Users\fox\Desktop\demo_HEI\STILT_output\2019_se\'];
    diroutput = dir(fullfile(filepath,'*_foot.nc'));
    filename = {diroutput.name};
    filename_1 = char(filename(1));
    lon_stilt = ncread([filepath,filename_1],'lon');
    lat_stilt = ncread([filepath,filename_1],'lat');

    lon_indx1 = find( abs(min(lon_stilt) - lon) == min( abs(min(lon_stilt) - lon) ));
    lon_indx2 = find( abs(max(lon_stilt) - lon) == min( abs(max(lon_stilt) - lon) ));
    lat_indx1 = find( abs(min(lat_stilt) - lat) == min( abs(min(lat_stilt) - lat) ));
    lat_indx2 = find( abs(max(lat_stilt) - lat) == min( abs(max(lat_stilt) - lat) ));

    lon_gfed_au = lon(lon_indx1(1):lon_indx2(1));
    lat_gfed_au = lat(lat_indx1(1):lat_indx2(1));
    xdim        = length(lon_gfed_au);
    ydim        = length(lat_gfed_au);
    eval(['predict_output_au_',num2str(years),'=','zeros(xdim,ydim,ndays)',';']);
    
    % loop for 12 months
    for mon = 1%1:12
        disp(['month:',num2str(mon)]);
        
        if mon == 1
            tt_indx1 = 1;
            tt_indx2 = ndays_mon(mon);
        else
            tt_indx1 = sum(ndays_mon(1:mon-1))+1;
            tt_indx2 = sum(ndays_mon(1:mon));
        end
        
        load('C:\Users\fox\Desktop\demo_HEI\random_forest\output_random_forest\rf_final_regression_ratio_pblh_training_2008_2019_reduced_input9.mat');
        load(['C:\Users\fox\Desktop\demo_HEI\random_forest\predicting_data\input_factors_12_',num2str(years),'_',num2str(mon),'.mat']);

        input_factors_au = input_factors(lon_indx1(1):lon_indx2(1),lat_indx1(1):lat_indx2(1),:,:);
        xdim = size(input_factors_au,1);
        ydim = size(input_factors_au,2);
        tdim = size(input_factors_au,3);
        %
        for ii = 1:xdim
            for jj = 1:ydim
                for tt = tt_indx1%:tt_indx2
                    var_tmp = squeeze(input_factors_au(ii,jj,tt-(tt_indx1-1),2));
                    if var_tmp == 0
                        eval(['predict_output_au_',num2str(years),'(ii,jj,tt)','=','0',';']);
                    end
                    if var_tmp > 0
                        input_test(1:8)  = squeeze(input_factors_au(ii,jj,tt-(tt_indx1-1),1:8));
                        input_test(9:11) = squeeze(input_factors_au(ii,jj,tt-(tt_indx1-1),10:12));
                        output_tmp       = final_b.predict(input_test);

                        if output_tmp > 1
                            output_tmp = 1;
                        end
                        
                        if output_tmp < 0
                            output_tmp = 0;
                        end
                        eval(['predict_output_au_',num2str(years),'(ii,jj,tt)','=','output_tmp',';']);
                    end
                end
            end
        end
    end
% filename_output = ['E:\work_for_wildfires\random_forest\predicted_injection_percentages_rf\new\predict_output_au_',num2str(years),'.mat'];

% save(filename_output,['predict_output_au_',num2str(years)]);
end

