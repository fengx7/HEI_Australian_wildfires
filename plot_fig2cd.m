%% Plotting script for Figure 2c-2d
% load trained random forest
load('E:\work_for_wildfires\random_forest\rf_final_regression_ratio_pblh_training_2008_2019_reduced_input9.mat');
% load historical dataset - input variables
load('E:\work_for_wildfires\random_forest\predictors.mat');
datanum1 = length(input_1);
load('E:\work_for_wildfires\random_forest\predictors_201911.mat');
datanum2 = length(input_1);
X        = zeros(datanum1+ datanum2,11);
Y        = zeros(datanum1+ datanum2,1);
X_name = {'LANDUSE','PBLH','T2','RH','U10','V10','PRECIP','EMIS','LON','LAT','FRP'};
load('E:\work_for_wildfires\random_forest\predictors.mat');
X(1:datanum1,:) = [input_1, input_2, input_3, input_4, input_5, input_6, input_7, input_8, input_10,input_11,input_12];
load('E:\work_for_wildfires\random_forest\predictors_201911.mat');
X(datanum1+1:datanum1+datanum2,:) = [input_1,input_2, input_3, input_4, input_5, input_6, input_7, input_8, input_10,input_11,input_12];

% load historical dataset - output variables
load('E:\work_for_wildfires\random_forest\output_variables_ratio.mat');
Y(1:datanum1) = ratio_uppblh_nowind;
load('E:\work_for_wildfires\random_forest\output_variables_ratio_201911.mat');
Y(datanum1+1:datanum1+datanum2) = ratio_uppblh_nowind;

% divide the historical data into training and test datasets
num_test = 200;
num_record = size(X,1);
index = 1:1:num_record;

index_test = 1:10:2000;
index_train = zeros(num_record - num_test,1);
i = 0;
for ii = 1:num_record
    aa = find(ii == index_test);
    if isempty(aa)
        i = i + 1;
        index_train(i) = ii;
    end
end

input = zeros(num_record-num_test,size(X,2));
input_test = zeros(num_test,size(X,2));
output = zeros(num_record-num_test,1);
output_test = zeros(num_test,1);

for ii = 1:num_test
    input_test(ii,:)  = X(index_test(ii),:);
    output_test(ii,:) = Y(index_test(ii),:);
end

for ii = 1:num_record - num_test
    input(ii,:) = X(index_train(ii),:);
    output(ii,:) = Y(index_train(ii),:);
end

% Output target variable
predict_output = final_b.predict(input_test);

for ii = 1:length(predict_output)
    if predict_output(ii) > 1
        predict_output(ii) = 1;
    end 
    
    if predict_output(ii) < 0
        predict_output(ii) = 0;
    end
end
%%
figure1 = figure;
pos1 = [0.3 0.3 0.4 0.4];
subplot1 = subplot('Position',pos1,'Parent',figure1);
plot(subplot1,output_test,predict_output,'Marker','o','LineWidth',2,'LineStyle','none','Color',[1 0 0]);
hold on
ss = 0:0.1:1;
plot(ss,ss,'LineStyle','--','LineWidth',1.5,'color','k');
hold on
ylabel('RF predictions');
% xlabel
xlabel('MISR observations');
title('(c) Plume injection fractions above PBL (%)');

box(subplot1,'on');
axis(subplot1,'square');
set(subplot1,'FontSize',24,'XTick',[0 0.2 0.4 0.6 0.8 1],'XTickLabel',[0 0.2 0.4 0.6 0.8 1],'YTick',...
    [0 0.2 0.4 0.6 0.8 1],'YTickLabel',[0 0.2 0.4 0.6 0.8 1]);
annotation(figure1,'textbox',...
    [0.40 0.58 0.198 0.114],...
    'String',{'R^2 = 0.53','RMSE = 0.22'},...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');

predict_emis = predict_output .* input_test(:,8) * 1e8;
output_emis  = output_test .* input_test(:,8)  * 1e8 ;
ss = 0:1e-4:1;
%%
figure1 = figure;
subplot3 = subplot('Position',pos1,'Parent',figure1);
loglog(subplot3, output_emis,predict_emis,'Marker','o','LineWidth',2,'LineStyle','none','Color',[1 0 0]);
hold on
plot(ss,ss,'LineStyle','--','LineWidth',1.5,'color','k');

ylabel('RF predictions');
% xlabel
xlabel('MISR observations');
title('(d) Fire emissions above PBL (10^-^8 kg m^-^2 s^-^1)');
xlim([1e-4,1]);
ylim([1e-4,1]);
box(subplot3,'on');
axis(subplot3,'square');
set(subplot3,'FontSize',24,'XMinorTick','on','XScale','log','XTick',...
    [0.0001 0.001 0.01 0.1 1],'XTickLabel',{' ','10^-^3','10^-^2','10^-^1','10^0'},'YMinorTick','on','YScale','log','YTick',...
    [0.0001 0.001 0.01 0.1 1],'YTickLabel',{' ','10^-^3','10^-^2','10^-^1','10^0'});
annotation(figure1,'textbox',...
    [0.40 0.58 0.198 0.114],...
    'String',{'R^2 = 0.88','RMSE = 2.1e-10 kg m^-^2 s^-^1'},...
    'FontSize',24,...
    'FitBoxToText','off',...
    'EdgeColor','none');