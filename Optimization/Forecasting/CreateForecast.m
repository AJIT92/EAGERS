function Forecast = CreateForecast(Date,Time)
global Plant Last24hour 
Days = ceil(rem(Date,1)+Time(end)/24);%this says that you may need to use surface fits for multiple days.
s12hour = nnz(Time<12); %steps in next 12 hours
A = datevec(Date);
dateDays = floor(Date):ceil(Date+Days);
Weekday = zeros(length(dateDays),1);
for i = 1:1:length(dateDays)
    Weekday(i) = isbusday(dateDays(i),Plant.Data.Holidays);
end
monthName = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec' 'Jan'};

%create forward forecast of temperature
hour = (Date+Time/24-floor(Date))*24;
n = length(Last24hour.Timestamp);
h = round((Last24hour.Timestamp(1) - floor(Last24hour.Timestamp(1)))*length(Last24hour.Timestamp));
if h ==0
    Tlast24 = [Last24hour.Temperature;Last24hour.Temperature(1)];
else Tlast24 = [Last24hour.Temperature(n-h+1:n);Last24hour.Temperature(1:n-h+1)];
end
YestFit = interp1(linspace(0,24,n+1),Tlast24,mod(hour,24));
HistFit = interp1(0:24,[Plant.Data.HistProf.Temperature(A(2),end),Plant.Data.HistProf.Temperature(A(2),:)],mod(hour,24));

W = interp1([Date,Date+3,Date+100],[0.9,0,0],Date+Time/24);%weight between yesterday and historical average
W(isnan(W)) = 0;
Forecast.T = (W.*YestFit + (1-W).*HistFit)';% Balanced between yesterdays T and historical T

S = fieldnames(Last24hour.Demand);
for s = 1:1:length(S) %repeat for electric, cooling, heating, and steam as necessary
    [nS,nD] = size(Last24hour.Demand.(S{s}));
    Forecast.(S{s}) = zeros(nS,nD);
    for k = 1:1:nD
        month = A(2);
        HistFitDem = nan(1,length(Time));
        for j = 1:1:ceil(Days) %% load surfaces
            if A(3)+ (j-1) > (datenum(A(1), A(2)+1, 1)-datenum(A(1), A(2), 1))
                month = A(2)+1;
            end
            list = fieldnames(Plant.Data.HistProf.(S{s})(1));
            if length(list)>1
                if nnz(strcmp(list,strcat(monthName(month),'WeekEnd')))>0 %historical profile is split to weekday/weekend
                    if Weekday(j) == 0
                        Surface = Plant.Data.HistProf.(S{s})(k).(char(strcat(monthName(month),'WeekEnd'))); 
                    else Surface = Plant.Data.HistProf.(S{s})(k).(char(strcat(monthName(month),'WeekDay')));
                    end
                elseif nnz(strcmp(list,char(monthName(month))))>0 %historical profile is broken by month, but not weekend/weekday
                    Surface = Plant.Data.HistProf.(S{s})(k).(monthName{month}); 
                else 
                    Surface = Plant.Data.HistProf.(S{s})(k).(list{1});
                end
            else Surface = Plant.Data.HistProf.(S{s})(k).(list{1});% only 1 historical surface fit 
            end
            index = (hour<=(24*j)).*(hour>(24*(j-1)))>0; %index of Time corresponding to this day
            HistFitDem(index) = Surface(mod(hour(index),24),Forecast.T(index));
        end
        YestPred = Surface(mod(hour,24),Last24hour.Temperature);
        YestError = Last24hour.Demand.(S{s})(:,k) - YestPred';
        Forecast.(S{s})(1:s12hour,k) =   HistFitDem(1:s12hour)+ linspace(1,0,s12hour).*interp1(linspace(0,1,length(Last24hour.Temperature)),YestError,Time(1:s12hour)/24);     
        Forecast.(S{s})(s12hour+1:end,k) = HistFitDem(s12hour+1:end);
        bias = Last24hour.Demand.(S{s})(end,k)/YestPred(end);
        error = Last24hour.Demand.(S{s})(end,k) - YestPred(end);
        % scale small values proportionally (up to double demand at t = 0), and add 2* error to large values
        for t = 1:1:s12hour
            if Forecast.(S{s})(t,k)>2*YestPred(end)
                Forecast.(S{s})(t,k) = Forecast.(S{s})(t,k) + (1-t/s12hour)*error;
            else
                Forecast.(S{s})(t,k) = Forecast.(S{s})(t,k).*(1+(bias-1)*(1-t/s12hour));
            end
        end
    end
end


%% check if the day you are forecasting is during the week (mon through fri) and that it is not a holiday
function workday = isbusday(dateDay, Holidays)
if isempty(Holidays)
    hol = [];
else
    hol = find(Holidays==dateDay);
end
if weekday(dateDay)<=6 && weekday(dateDay)>=2 && isempty(hol)
     workday = 1;
else workday = 0;
end