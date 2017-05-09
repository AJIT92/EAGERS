function Forecast = CreateForecast(Date,Time,RealData)
global Plant Last24hour
Days = ceil(rem(Date,1)+Time(end)/24); %this says that you may need to use surface fits for multiple days.
s12hour = nnz(Time<12); %steps in next 12 hours
A = datevec(Date);
dateDays = floor(Date)-1:ceil(Date+Days);
Weekday = zeros(length(dateDays),1);
for i = 1:1:length(dateDays)
    % Is dateDays(i) a weekday and not a holiday?
    if weekday(dateDays(i))<=6 && weekday(dateDays(i))>=2 && ~any(Plant.Data.Holidays==dateDays(i))
        Weekday(i) = 1;
    else
        Weekday(i) = 0;
    end
    % Is the historical data from a year in the past? (DOES NOT HANDLE
    % HISTORICAL DATA FROM MORE THAN 1 YEAR IN THE PAST.)
    if Weekday(i)==1 && any(Plant.Data.Holidays == dateDays(i)+365) && dateDays(i)<min(Plant.Data.Holidays)
        Weekday(i) = 0;
    end
end
monthName = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec' 'Jan'};

%create forward forecast of temperature
hour = (Date+[0;Time]/24-floor(Date))*24;
YestFit = interp1(linspace(0,24,length(Last24hour.Temperature)+1)',[Last24hour.Temperature;RealData.Temperature],mod([0; Time],24));
HistFit = interp1(0:24,[Plant.Data.HistProf.Temperature(A(2),end),Plant.Data.HistProf.Temperature(A(2),:)],mod(hour,24));

W = interp1([Date,Date+3,Date+100],[0.9,0,0],Date+[0;Time]/24);%weight between yesterday and historical average
W(isnan(W)) = 0;
Tforecast = (W.*YestFit + (1-W).*HistFit);% Balanced between yesterdays T and historical T (includes forecast for current time)
S = fieldnames(Last24hour.Demand);
for s = 1:1:length(S) %repeat for electric, cooling, heating, and steam as necessary
    [~,nD] = size(Last24hour.Demand.(S{s}));
    Forecast.(S{s}) = zeros(length(Time),nD);
    list = fieldnames(Plant.Data.HistProf.(S{s})(1));
    for k = 1:1:nD
        month = A(2);
        HistFitDem = nan(length(Time),1);
        for j = 0:1:ceil(Days)
            if A(3)+ (j-1) > (datenum(A(1), A(2)+1, 1)-datenum(A(1), A(2), 1))
                month = A(2)+1;
            end
            %load suface
            if length(list)>1
                if nnz(strcmp(list,strcat(monthName(month),'WeekEnd')))>0 %historical profile is split to weekday/weekend
                    if Weekday(j+1) == 0
                        Surface = Plant.Data.HistProf.(S{s})(k).(strcat(monthName{month},'WeekEnd'));
                    else
                        Surface = Plant.Data.HistProf.(S{s})(k).(strcat(monthName{month},'WeekDay'));
                    end
                elseif nnz(strcmp(list,char(monthName(month))))>0 %historical profile is broken by month, but not weekend/weekday
                    Surface = Plant.Data.HistProf.(S{s})(k).(monthName{month});
                else 
                    Surface = Plant.Data.HistProf.(S{s})(k).(list{1});
                end
            else
                Surface = Plant.Data.HistProf.(S{s})(k).(list{1});% only 1 historical surface fit
            end
            if j ==0 %surface for yesterday
                YestPred = Surface(mod((Last24hour.Timestamp-floor(Last24hour.Timestamp(1)))*24,24),Last24hour.Temperature);
                BiasPerc = mean((Last24hour.Demand.(S{s})(:,k)- YestPred)./YestPred)*100; %average percent error between hisorical expectations for yesterday, and actual yesterday
                HourlyError = Last24hour.Demand.(S{s})(:,k) - YestPred*(1+0.7*BiasPerc/100);%hourly variation with equal mean (historical and last 24 hours)
            else %surface for current day or subsequent days
                index = (hour<=(24*j)) & (hour>(24*(j-1))); %index of Time corresponding to this day
                HistFitDem(index) = Surface(mod(hour(index),24),Tforecast(index));
            end
        end
        Forecast.(S{s})(1:s12hour,k) =   HistFitDem(2:s12hour+1)*(1+0.7*BiasPerc/100) + linspace(1,0,s12hour)'.*interp1(Last24hour.Timestamp,HourlyError,Date+Time(1:s12hour)/24-1);
        Forecast.(S{s})(s12hour+1:end,k) = HistFitDem(s12hour+2:end)*(1+0.7*BiasPerc/100);
    end
end
Forecast.T = Tforecast(2:end);%exclude forecast for current time