%% Plotting script for trained random forest
load('E:\work_for_wildfires\random_forest\rf_final_regression_ratio_pblh_training_2008_2019_reduced_input9.mat');
X_name = {'LANDUSE','PBLH','T2','RH','U10','V10','PRECIP','EMIS','LON','LAT','FRP'};
%%
figure1 = figure;
pos1 = [0.3 0.3 0.4 0.4];
aa = final_b.OOBPermutedPredictorDeltaError;
[bb,index] = sort(aa,'ascend');
X_name_new = cell(1,11);
for ii = 1:length(index)
    X_name_new{ii} = X_name{index(ii)};
end

subplot2 = subplot('Position',pos1,'Parent',figure1);
hold(subplot2,'on');
axes(subplot2);
xx = 1:11;
barh(xx,bb,'Parent',subplot2,...
    'FaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'EdgeColor','none',...
    'BarWidth',0.6);
ylim([0,6]);
ylim([0,12]);
xlabel('Feature Importance');
set(subplot2,'FontSize',24,'XTick',[0 2 4 6],'XTickLabel',[0 2 4 6],'YTick',[1 2 3 4 5 6 7 8 9 10 11],'YTickLabel',X_name_new);
box(subplot2,'on');
axis(subplot2,'square');

rmse1 = rmse(predict_output,output_test);
rr  = corrcoef(predict_output,output_test);
nmb = mean(predict_output - output_test)/mean(output_test);
r   = rr(1,2);
