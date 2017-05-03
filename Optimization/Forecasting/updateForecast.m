function [Forecast, Renewable] = updateForecast(Date,RealData)
%Date is the date number, Time is a vector of times (in hours)
global Plant
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
Forecast =  CreateForecast(Date,Time,RealData);
nS = length(Time);
nG = length(Plant.Generator);
Renewable = zeros(nS,nG);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
        Renewable(:,i) = RenewableOutput(Plant.Generator(i).VariableStruct,Date,Time,'Predict');
    end
end