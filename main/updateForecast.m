function Forecast = updateForecast(Date,Time)
%Date is the date number, Time is a vector of times (in hours), renew is
%Plant.Generator for any renewable generators
global Last24hour RealTimeData
if ~isempty(Last24hour)
    S = fieldnames(Last24hour);
    Data =[];
else
   S = ['T',fieldnames(RealTimeData.Demand)'];
   Data = RealTimeData;
end
for i = 1:1:length(S) %repeat for electric, cooling, heating, and steam as necessary
    Days = ceil(rem(Date,1)+Time(end)/24);%this says that you may need to use surface fits for multiple days.
    if ~strcmp(S{i},'Timestamp')
        Forecast.(S{i}) =  CreateForecast((S{i}),Date,Time,Days,Last24hour,Data,'single');
    end
end
%%Predict Renewables
[renewGen, renew]= RenewableOutput(Date,Time,'Predict');%Add forecast of renewable power @ resulution equal to the Time vector
if ~isempty(renew)
    Forecast.E = Forecast.E-sum(renewGen,2)';%if there are no renewables then renewGen is a vector of zeros
end