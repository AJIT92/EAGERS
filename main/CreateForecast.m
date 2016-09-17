function pred = CreateForecast(Param,Date,Time,Days,PrevDay,Data,HiLow)
global Holidays HistProf
nS = length(Time); %Time is in hours
dt = Time-[0,Time(1:end-1)];
pred = zeros(3,nS);
s12hour = nnz(Time<12); %steps in next 12 hours
sAfter12 = nnz(Time>=12); %steps after 12 hours
A = datevec(Date);
dateDays = floor(Date):ceil(Date+Days);
Weekday = zeros(length(dateDays),1);
for i = 1:1:length(dateDays)
    Weekday(i) = isbusday(dateDays(i),Holidays);
end
monthName = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec' 'Jan'};
  
if isempty(PrevDay)%% load previous day profile if not specified (only when finding initial condition!!)
    date1 = datevec(Data.Timestamp(1));
    year1 = date1(1);
    interpDate = datenum([year1, A(2:end)])-1;
    if interpDate<Data.Timestamp(1)
        interpDate = interpDate + ceil(Data.Timestamp(1)-interpDate); %move forward a whole # of days
    elseif (interpDate+1)>Data.Timestamp(end);
        interpDate = interpDate - ceil(interpDate + 1 - Data.Timestamp(end)); %move backward a whole # of days
    end
    interpDate = interpDate + (Time - Time(1))/24; %previous 24 hours @ same frquency as dispatch loop timestep
    %% need to make sure interpDate is in the Data range
    PrevDay.T = interp1(Data.Timestamp,Data.Temperature,interpDate);
    if ~strcmp(Param,'T')
        PrevDay.(Param) = interp1(Data.Timestamp,Data.Demand.(Param),interpDate);
    end
    if Date-1<interpDate(1)
        Date=interpDate(1)+1;
        A = datevec(Date);
    end
    pred = PrevDay.(Param);
else
    interpDate = PrevDay.Timestamp; %previous 24 hours in same increment as PrevDay.T

    %create forward forecast of temperature
    YestFit = nan(1,nS);
    HistFit = nan(1,nS);
    for i = 1:1:Days+1
        YestFit(isnan(YestFit)) = interp1(interpDate,PrevDay.T,Date+Time(isnan(YestFit))/24-i);
        HistFit(isnan(HistFit)) = interp1(linspace(datenum(A(1:3))-1,datenum(A(1:3)),length(HistProf.Temperature)+1),[HistProf.Temperature(A(2),end) HistProf.Temperature(A(2),:)],Date+Time(isnan(HistFit))/24-i);
    end

    if nnz(dt>1)>0 && ~strcmp(Param,'T')%if there are any timesteps that are greater than the data sample size save the Yest and Hisfit interps
%         Ttrap = [Time-dt./2,Time(end)+dt(end)/2];%this is the time steps for the trapezoidal shapes
%         dtTrap = Ttrap(2:end)-Ttrap(1:end-1);
        Timeinterp = unique([Time,0:floor(Time(end))]);
    else Timeinterp = Time;
    end
    
    YestFitinterp = nan(1,length(Timeinterp));
    HistFitinterp = nan(1,length(Timeinterp));
    HistProfTemp = [HistProf.Temperature(A(2),end),HistProf.Temperature(A(2),:)];
    for i = 1:1:Days+1 %time average the historical and yesterday fits
        if interpDate(end)<interpDate(1)+1 && Days>1 %make sure there are 24 hours of data to interpolate from
            YestFitinterp(isnan(YestFitinterp)) = interp1([interpDate,interpDate(1)+1],[PrevDay.T,PrevDay.T(1)],Date+Timeinterp(isnan(YestFitinterp))/24-i);
        else YestFitinterp(isnan(YestFitinterp)) = interp1(interpDate,PrevDay.T,Date+Timeinterp(isnan(YestFitinterp))/24-i);
        end
        HistFitinterp(isnan(HistFitinterp)) = interp1(linspace(datenum(A(1:3))-1,datenum(A(1:3)),length(HistProf.Temperature)+1),HistProfTemp,Date+Timeinterp(isnan(HistFitinterp))/24-i);
    end
    % now integrate the temperature and divide by dt if parameter it T
    if strcmp(Param,'T') && nnz(dt>1)>0 %only do the intigration division if your parameter is T and you have a timestep larger than the data timesteps
        TIntFit = (interpDate-Date+1).*24;
        HIntFit = 0:length(HistProfTemp)-1;%this is necessary for the online loops where T is shorter
%         ytrap = YestFitinterp(ismember(Timeinterp,Ttrap));
%         htrap = HistFitinterp(ismember(Timeinterp,Ttrap));
        lastX = 0;
        ylast = YestFitinterp(end);
        hlast = HistFitinterp(end);
        for i = 1:length(Time)
            if dt(i)>1 %only integrate if the timestep is larger than one timestep of data, otherwise stick with the interpolation
                %centered trapezoidal integration is commented out. this
                %method would be used for instantaneous power data (kW)
%                 xy = [Ttrap(i), TIntFit(and(TIntFit>Ttrap(i),TIntFit<Ttrap(i+1))), Ttrap(i+1)];
%                 y = [ytrap(i),PrevDay.T(and(TIntFit>Ttrap(i),TIntFit<Ttrap(i+1))),ytrap(i+1)];
%                 xh = [Ttrap(i), HIntFit(and(HIntFit>Ttrap(i),HIntFit<Ttrap(i+1))), Ttrap(i+1)];
%                 h = [htrap(i), HistProfTemp(and(HIntFit>Ttrap(i),HIntFit<Ttrap(i+1))), htrap(i+1)];
                %left handed trapezoidal integration is used because it is
                %assumed that data is given in energy use for the past
                %amount of time (kWh)
                xy = [lastX,TIntFit(and(TIntFit>lastX,TIntFit<Time(i))),Time(i)];
                y = [ylast,PrevDay.T(and(TIntFit>lastX,TIntFit<Time(i))),YestFitinterp(i)];
                xh = [lastX, HIntFit(and(HIntFit>lastX,HIntFit<Time(i))),Time(i)];
                h = [hlast,HistProfTemp(and(HIntFit>lastX,HIntFit<Time(i))),HistFitinterp(i)];
                YestFit(i) = trapz(xy,y)/dt(i);
                HistFit(i) = trapz(xh,h)/dt(i);
            end
            lastX = Time(i);
            ylast = YestFitinterp(i);
            hlast = HistFitinterp(i);
        end
    else YestFit = YestFitinterp;
        HistFit = HistFitinterp;
    end

    W = interp1([Date,Date+3,Date+100],[0.9,0,0],Date+Timeinterp/24);%weight between yesterday and historical average
    W(isnan(W)) = 0;
    T = W.*YestFit + (1-W).*HistFit;% Balanced between yesterdays T and historical T
    if Date>interpDate(end) %this occurs when MPCloop has moved the date forward and is predicting, but one dispatch loop has not finished, so the structure Last24hour has not been updated becuase recordFromMPCloop has not been called yet
        Date = Date-1;
    end
    bias = (PrevDay.T(end)-interp1(interpDate,PrevDay.T,Date)); %connect current T to prediction
    T(1,1:s12hour) = T(1,1:s12hour)+ linspace(1,0,s12hour)*bias;

    if strcmp(Param,'T')
        pred(1,:) = T;
        if strcmp(HiLow,'HiLow')%%% High/Low prediction
            error =([linspace(.5,3,s12hour-1) interp1([3,8],[12,72],Time(s12hour:end))]);
            pred(2,:) = pred(1,:)+error;
            pred(3,:) = pred(1,:)-error;
        else  pred = max(0,pred(1,:));
        end
    else%Demand
        month = A(2);
        HistFitDem = nan(1,length(Time));
        tlast = 0;
        for j = 1:1:ceil(Days) %% load surfaces
            if A(3)+ (j-1) > (datenum(A(1), A(2)+1, 1)-datenum(A(1), A(2), 1))
                month = A(2)+1;
            end
            if isstruct(HistProf.(Param))
                list = fieldnames(HistProf.(Param));
                if nnz(strcmp(list,char(strcat(monthName(month),'WeekEnd'))))>0 %historical profile is split to weekday/weekend
                    if Weekday(j) == 0
                        Surface = HistProf.(Param).(char(strcat(monthName(month),'WeekEnd'))); 
                    else Surface = HistProf.(Param).(char(strcat(monthName(month),'WeekDay')));
                    end
                elseif nnz(strcmp(list,char(monthName(month))))>0 %historical profile is broken by month, but not weekend/weekday
                    Surface = HistProf.(Param).(char(monthName(month))); 
                else 
                    Surface = HistProf.(Param).(char(list(1)));
                end
            else Surface = HistProf.(Param);% only 1 historical surface fit 
            end
            index = (((A(4)+A(5)/60+A(6)/3600)+Time)<=(24*j)).*(((A(4)+A(5)/60+A(6)/3600)+Time)>(24*(j-1)))>0; %index of Time corresponding to this day
            tDay = ((A(4)+A(5)/60+A(6)/3600))+Time(index)-(24*(j-1));
            HistFitDem(index) = Surface(tDay,T(index));

    %         HistFitDem(index) = Surface(tDay,T(ismember(Time(index),Timeinterp)));%Timeinterp is the time vector setup when finding T, it includes Time and 0:horizon
            if nnz(dt>1)>0 %if there are stepsizes greater than an hour, then do an integration/time to find demand fits
                TsEnd = (Time(index));
                TsEnd = TsEnd(end);
                HistFitinterp = Surface(Timeinterp(and(Timeinterp>=tlast,Timeinterp<=TsEnd)),T(and(Timeinterp>=tlast,Timeinterp<=TsEnd)));%first create a profile for that day
                timetoday = Timeinterp(and(Timeinterp<=tDay(end)+24*(j-1),Timeinterp>=tlast));
                dttoday = dt(index);
    %             clockNow = A(4)+A(5)/60+A(6)/3600;
                for i = 1:1:length(tDay)
                    if dttoday(i)>1 %only integrate if the step is larger than 1 hour, otherwise, keep interpolated demand value
                        HistFitDem(i+find(index,1)-1) = trapz(Timeinterp(and(Timeinterp<=Time(i+find(index,1)-1),Timeinterp>=tlast)),HistFitinterp(and(timetoday<=Time(i+find(index,1)-1),timetoday>=tlast)))/dttoday(i);
                    end
                    tlast = Time(i+find(index,1)-1);
                end
            end
        end
        B = datevec(interpDate);
        YestPred = Surface((B(:,4)+B(:,5)/60+B(:,6)/3600),PrevDay.T);
        YestError = PrevDay.(Param) - YestPred';
        pred(1,1:s12hour) =   HistFitDem(1:s12hour)+ linspace(1,0,s12hour).*interp1(linspace(0,1,length(PrevDay.T)),YestError,Time(1:s12hour)/24);     
        pred(1,s12hour+1:end) = HistFitDem(s12hour+1:end);
        bias = PrevDay.(Param)(end)/YestPred(end);
        error = PrevDay.(Param)(end) - YestPred(end);
        % scale small values proportionally (up to double demand at t = 0), and add 2* error to large values
        for t = 1:1:s12hour
            if pred(1,t)>2*YestPred(end)
                pred(1,t) = pred(1,t) + (1-t/s12hour)*error;
            else
                pred(1,t) = pred(1,t).*(1+(bias-1)*(1-t/s12hour));
            end
        end

        if strcmp(HiLow,'HiLow')%%% High/Low prediction
            if ~isfield(Data,'Demand')
                Floor = .8*min(PrevDay.(Param));
            else Floor = .8*min(Data.Demand.(Param));
            end
            pred(2,:) = Floor+[linspace(1.1,1.3,s12hour) linspace(1.3,1.3,sAfter12)].*(pred(1,:)-Floor);
            pred(3,:) = Floor+[linspace(.9,.7,s12hour) linspace(.7,.7,sAfter12)].*(pred(1,:)-Floor);
            pred = max(0,pred(1:3,:));
        else  pred = max(0,pred(1,:));
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
