function [Forecast, Renewable] = updateForecast(Date,Time)
%Date is the date number, Time is a vector of times (in hours)
global Plant
Forecast =  CreateForecast(Date,Time);
nS = length(Time);
nG = length(Plant.Generator);
Renewable = zeros(nS,nG);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
        Renewable(:,i) = RenewableOutput(Plant.Generator(i).VariableStruct,Date,Time,'Predict');
    end
end