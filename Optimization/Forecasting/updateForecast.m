function Forecast = updateForecast(Date,RealData)
%Date is the date number, Time is a vector of times (in hours)
global Plant
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
switch Plant.optimoptions.forecast
    case 'SES'
        %%
    case 'ARIMA'
        %%
    case 'NeuralNet'
        %%
    case 'Surface'
        Forecast =  CreateForecast(Date,Time,RealData);
    case 'Perfect'
        Forecast = GetCurrentData(Date+Time/24);
end

nS = length(Time);
nG = length(Plant.Generator);
Forecast.Renewable = zeros(nS,nG);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
        Forecast.Renewable(:,i) = RenewableOutput(Plant.Generator(i).VariableStruct,Date,Time,'Predict');
    end
end