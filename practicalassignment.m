%% Practical Assignment - Metals Mining - Scripts
clearvars;close all;clc

% Run the whole file to get nice final figures (and include Excel files to 
% current path) :-)


warning('off', 'MATLAB:readtable:ModifiedVarnames');
%% Gold price - prediction with gbm

% Read gold data 
golddata = readmatrix("DataSheet.xlsx", 'Range', 'A2:C240');
prices = golddata(1:end,2);
time = 1:length(prices);
change = golddata(2:end, 3);

% Quick plots
%plot(time, prices, 'r')
%figure
%plot(time(2:end), change, 'g')

% Compute parameters for GBM
drift = mean(change)*12; % for mean all historical data
volatility = std(change)*sqrt(12);
initialPrice = prices(end);

% Create GMB distribution
gbm_pd = gbm(drift, volatility, 'StartState', initialPrice);
nYrs = 20;
nSim = 10000;
simTime = 5;
goldPriceSims = zeros(nYrs, nSim);

% Simulating first 5 years
for i = 1:nSim
    goldPriceSims(1:simTime,i) = gbm_pd.simBySolution(simTime-1);
end

curLast = simTime;

% Simulating 3 more 5 year periods
for i = 1:3
    for j = 1:nSim
        newStart = goldPriceSims(curLast, j);
        gbm_pd = gbm(drift, volatility, 'StartState', newStart);
        curSim = gbm_pd.simBySolution(simTime);
        goldPriceSims(curLast+1:curLast + simTime, j) = curSim(2:end);
    end
    curLast = curLast + simTime;
end




%% Plotting the gold price

% Modifying the prices to yearly values
yearlyprices = [];
for i = 1:12:length(prices)
    yearlyprices = [yearlyprices prices(i)];
end

% 10% and 90% quantiles
Q_10 = quantile(goldPriceSims, 0.1, 2);
Q_90 = quantile(goldPriceSims, 0.9, 2);

% Visualizing historical prices + quantiles and medians of simulated ones
figure
plot(-18:1', yearlyprices, 'r', 'LineWidth', 3)
hold on
grid on
title('Estimated gold price')
plot(1:20', goldPriceSims')
plot(1:20', median(goldPriceSims,2), 'k--', 'LineWidth', 3)
plot(1:20', Q_10, 'k--', 'LineWidth', 3)
plot(1:20', Q_90, 'k--', 'LineWidth', 3)
xline(5, 'k--', 'LineWidth', 2); 
xline(10, 'k--', 'LineWidth', 2); 
xline(15, 'k--', 'LineWidth', 2);
xlim([-20 19]);
xlabel('Year')
ylabel('Gold price [USD/tn]')

%% Mine 4Mt
% Simulating the 4Mt plan

% Reading data
data = readtable("Planning_Data_with_NPV_values.xlsx", 'Sheet', 'Plan_4Mt', 'Range','B2:U22');

% Fixed variables
dataFixed = readmatrix("Planning_Data_with_NPV_values.xlsx", 'Sheet', 'Plan_4Mt', 'Range','D25:D28');
initMiningCost = dataFixed(1);
miningCostInflation = dataFixed(2);
discountRate = dataFixed(4);

% Defining model parameters
model = Simulink.SimulationInput('MineModel.slx');
model = model.setVariable('initialMiningCost', initMiningCost);
model = model.setVariable('miningCostInflationRate', miningCostInflation);
model = model.setVariable('discountRate', discountRate);

model = model.setModelParameter('SolverType', 'Fixed-step');
model = model.setModelParameter('FixedStep', '1');
model = model.setModelParameter('StartTime', '0');
model = model.setModelParameter('StopTime', '19');

model = model.setVariable('period', [data.Period-1 data.Period]);
model = model.setVariable('tonnes', [data.Period-1 data.Tonnes]);
model = model.setVariable('mill', [data.Period-1 data.Mill1]);
model = model.setVariable('millAuGrade', [data.Period-1 data.Mill_Au_GRADE_g_t_]);
model = model.setVariable('waste', [data.Period-1 data.Waste]);
model = model.setVariable('stockpileIn', [data.Period-1 data.Stockpile_t__in]);
model = model.setVariable('stockpileOut', [data.Period-1 data.Stockpile_t__out]);
model = model.setVariable('recoveryRate', [data.Period-1 data.RecoveryRate]);
model = model.setVariable('unitProcessingCost', [data.Period-1 data.UnitProcessingCost_USD_tn_]);
model = model.setVariable('capitalExpenditure', [data.Period-1 data.CapitalExpenditure_USD_]);
model = model.setVariable('taxAndRoyalty', [data.Period-1 data.TaxAndRoyalty___]);
model = model.setVariable('goldPrice', [data.Period-1 goldPriceSims(1:data.Period(end), :)]);

% Simulating the model and storing the resulting NPV's
results = sim(model);
NPVvals = squeeze(results.yout{1}.Values.Data)';

% Computing probability of NPV's being below 0
count = zeros(size(NPVvals,1)-1,1);
for i = 1:size(NPVvals,1)
    count(i) = sum(NPVvals(i,:) < 0);
end
percentages = count ./ length(NPVvals) * 100;

% Visualizing 
figure
subplot(3,3,1)
sgtitle('NPV results for each mine: 4Mt, 6Mt, 8Mt')
plot(1:20, median(NPVvals, 2), 'r')
hold on
plot(1:20, quantile(NPVvals, 0.9, 2), 'k--')
plot(1:20, quantile(NPVvals, 0.1, 2), 'k--')
grid on
xlim([1 20])
xlabel('Year')
ylabel('NPV')
title('NPV values')
legend('Median NPV', '10% quantile', '90% quantile')

subplot(3,3,2)
hold on
grid on
xlim([1 20])
plot(1:20, cumsum(median(NPVvals, 2)), 'g')  
xlabel('Year')
ylabel('NPV')
title('Cumul. median NPV')

subplot(3,3,3)
grid on
hold on
xlim([2 20])
xlabel('Year')
ylabel('%')
title('% of NPV values below threshold')

plot(2:20, percentages(2:end), 'r*')

%% Mine 6Mt, 14 per
% Simulating the 6Mt plans as before


data = readtable("Planning_Data_with_NPV_values.xlsx", 'Sheet', 'Plan_6Mt', 'Range','B2:U16');
dataFixed = readmatrix("Planning_Data_with_NPV_values.xlsx", 'Sheet', 'Plan_6Mt', 'Range','D25:D28');

initMiningCost = dataFixed(1);
miningCostInflation = dataFixed(2);
discountRate = dataFixed(4);

model = Simulink.SimulationInput('MineModel.slx');
model = model.setModelParameter('SolverType', 'Fixed-step');
model = model.setModelParameter('FixedStep', '1');
model = model.setModelParameter('StartTime', '0');
model = model.setModelParameter('StopTime', '13');

model = model.setVariable('period', [data.Period-1 data.Period]);
model = model.setVariable('tonnes', [data.Period-1 data.Tonnes]);
model = model.setVariable('mill', [data.Period-1 data.Mill1]);
model = model.setVariable('millAuGrade', [data.Period-1 data.Mill_Au_GRADE_g_t_]);
model = model.setVariable('waste', [data.Period-1 data.Waste]);
model = model.setVariable('stockpileIn', [data.Period-1 data.Stockpile_t__in]);
model = model.setVariable('stockpileOut', [data.Period-1 data.Stockpile_t__out]);
model = model.setVariable('recoveryRate', [data.Period-1 data.RecoveryRate]);
model = model.setVariable('unitProcessingCost', [data.Period-1 data.UnitProcessingCost_USD_tn_]);
model = model.setVariable('capitalExpenditure', [data.Period-1 data.CapitalExpenditure_USD_]);
model = model.setVariable('taxAndRoyalty', [data.Period-1 data.TaxAndRoyalty___]);

model = model.setVariable('initialMiningCost', initMiningCost);
model = model.setVariable('miningCostInflationRate', miningCostInflation);
model = model.setVariable('discountRate', discountRate);


model = model.setVariable('goldPrice', [data.Period-1 goldPriceSims(1:data.Period(end), :)]);
results = sim(model);
NPVvals = squeeze(results.yout{1}.Values.Data)';

count = zeros(size(NPVvals,1)-1,1);
for i = 1:size(NPVvals,1)
    count(i) = sum(NPVvals(i,:) < 0);
end

percentages = count ./ length(NPVvals) * 100;
subplot(3,3,4)
plot(1:14, median(NPVvals, 2), 'r')
hold on
plot(1:14, quantile(NPVvals, 0.9, 2), 'k--')
plot(1:14, quantile(NPVvals, 0.1, 2), 'k--')
grid on
xlim([1 14])
xlabel('Year')
ylabel('NPV')
title('NPV values')
legend('Median NPV', '10% quantile', '90% quantile')

subplot(3,3,5)
hold on
grid on
xlim([1 14])
plot(1:14, cumsum(median(NPVvals, 2)), 'g')  
xlabel('Year')
ylabel('NPV')
title('Cumul. median NPV')

subplot(3,3,6)
grid on
hold on
xlim([2 14])
xlabel('Year')
ylabel('%')
title('% of NPV values below threshold')

plot(2:14, percentages(2:end), 'r*')


%% Mine 8Mt, 11 per
% Simulating the 8Mt plan as before

data = readtable("Planning_Data_with_NPV_values.xlsx", 'Sheet', 'Plan_8Mt', 'Range','B2:U13');
dataFixed = readmatrix("Planning_Data_with_NPV_values.xlsx", 'Sheet', 'Plan_8Mt', 'Range','D25:D28');

initMiningCost = dataFixed(1);
miningCostInflation = dataFixed(2);
discountRate = dataFixed(4);

model = Simulink.SimulationInput('MineModel.slx');
model = model.setModelParameter('SolverType', 'Fixed-step');
model = model.setModelParameter('FixedStep', '1');
model = model.setModelParameter('StartTime', '0');
model = model.setModelParameter('StopTime', '10');

model = model.setVariable('period', [data.Period-1 data.Period]);
model = model.setVariable('tonnes', [data.Period-1 data.Tonnes]);
model = model.setVariable('mill', [data.Period-1 data.Mill1]);
model = model.setVariable('millAuGrade', [data.Period-1 data.Mill_Au_GRADE_g_t_]);
model = model.setVariable('waste', [data.Period-1 data.Waste]);
model = model.setVariable('stockpileIn', [data.Period-1 data.Stockpile_t__in]);
model = model.setVariable('stockpileOut', [data.Period-1 data.Stockpile_t__out]);
model = model.setVariable('recoveryRate', [data.Period-1 data.RecoveryRate]);
model = model.setVariable('unitProcessingCost', [data.Period-1 data.UnitProcessingCost_USD_tn_]);
model = model.setVariable('capitalExpenditure', [data.Period-1 data.CapitalExpenditure_USD_]);
model = model.setVariable('goldPrice', [data.Period-1 data.GoldPrice_USD_oz_]);
model = model.setVariable('taxAndRoyalty', [data.Period-1 data.TaxAndRoyalty___]);

model = model.setVariable('initialMiningCost', initMiningCost);
model = model.setVariable('miningCostInflationRate', miningCostInflation);
model = model.setVariable('discountRate', discountRate);

model = model.setVariable('goldPrice', [data.Period-1 goldPriceSims(1:data.Period(end), :)]);
results = sim(model);
NPVvals = squeeze(results.yout{1}.Values.Data)';

count = zeros(size(NPVvals,1)-1,1);
for i = 1:size(NPVvals,1)
    count(i) = sum(NPVvals(i,:) < 0);
end

percentages = count ./ length(NPVvals) * 100;

subplot(3,3,7)
plot(1:11, median(NPVvals, 2), 'r')
hold on
plot(1:11, quantile(NPVvals, 0.9, 2), 'k--')
plot(1:11, quantile(NPVvals, 0.1, 2), 'k--')
grid on
xlim([1 11])
xlabel('Year')
ylabel('NPV')
title('NPV values')
legend('Median NPV', '10% quantile', '90% quantile')

subplot(3,3,8)
hold on
grid on
xlim([1 11])
plot(1:11, cumsum(median(NPVvals, 2)), 'g')  
xlabel('Year')
ylabel('NPV')
title('Cumul. median NPV')

subplot(3,3,9)
grid on
hold on
xlim([2 11])
xlabel('Year')
ylabel('%')
title('% of NPV values below threshold')

plot(2:11, percentages(2:end), 'r*')
