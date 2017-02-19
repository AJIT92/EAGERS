function [DGpower]=ControllerWeekendDip(DemandE,SysSize,RampRate,ElecOnly)

steps = length(DemandE);
Ts = 8760/steps;

%Find holidays/weekends
Year=2006;
startDate = datenum([Year 1 1  0 0 0]);
dayMax = zeros(365,1);
for day = 1:1:365
    dayMax(day) = max(ElecOnly(1+24/Ts*(day-1):24/Ts*day));
end
holidays = [];
for i = 1:52
    lowDays =7*(i-1) + find(dayMax(1+7*(i-1):7*i)<(median(dayMax(1+7*(i-1):7*i))-.5*std(dayMax(1+7*(i-1):7*i))));
    holidays(end+1:end+length(lowDays)) = lowDays;
end
if weekday(startDate+364)==7 || weekday(startDate+364)==1
    holidays(end+1) = 365;
end


%Find Weekday & weekend minimums
CHPsize = sum(SysSize);
WeekdayMin = zeros(53,1);
WeekendMin = zeros(53,1);
week = 1;
RampTs = sum(RampRate.*SysSize)/CHPsize;
WeekendProf = [];
WeekdayProf = [];
for day = 1:1:365
    
    dayOfWeek = weekday(startDate+day-1);
    if nnz(day == holidays)==1% ||  dayOfWeek==1 ||  dayOfWeek==7
        WeekendProf(end+1:end+24/Ts) = DemandE(1+24/Ts*(day-1):24/Ts*day);
    else
        WeekdayProf(end+1:end+24/Ts) = DemandE(1+24/Ts*(day-1):24/Ts*day);
    end
    if dayOfWeek == 6 || day == 365
        if ~isempty(WeekdayProf)
            WeekdayMin(week) = min(min(WeekdayProf(10/Ts:end-12/Ts)),CHPsize);
            RampUp(week) = 1;
            while WeekdayProf(RampUp(week))<WeekdayMin(week) && RampUp(week)+1<12/Ts
                RampUp(week) = RampUp(week)+1;
            end
            RampDown(week) = length(WeekdayProf)-12/Ts;
            while WeekdayProf(RampDown(week)+1)<WeekdayMin(week) && RampDown(week)+1<24/Ts
                RampDown(week) = RampDown(week)+1;
            end
        end
        if ~isempty(WeekendProf)
            WeekendMin(week) = min(min(WeekendProf),CHPsize);
        else WeekendMin(week) = WeekdayMin(week);
        end
        week = week+1;
        WeekendProf = [];
        WeekdayProf = [];
    end
end


%% Meet Electric Demand
DGpower = zeros(steps,1);
week = 1;
for day = 1:1:365
    dayOfWeek = weekday(startDate+day-1);
    MinPow = WeekendMin(week);
    MaxPow = WeekdayMin(week);

    if nnz(day == holidays)==1 %% weekend/holiday
        DayDisp(1:24/Ts) = MinPow;
        if day<365 && nnz(day+1 == holidays)==0 %% tomorrow = weekday
            if week<=52 && dayOfWeek == 6
                ramp = WeekdayMin(week+1)-WeekendMin(week);
            else ramp = WeekdayMin(week)-WeekendMin(week);
            end
            Hz = max(ceil(abs(ramp/(CHPsize*RampTs))),1);
            I = 1+24/Ts*(day-1);
            I = max(find(DemandE(I:I+36/Ts-1)<MaxPow,1,'last'),Hz);
            if I<24/Ts
                RampUp = I;
                while min(DemandE(RampUp-Hz+1+24/Ts*(day-1):RampUp+24/Ts*(day-1))- (MinPow +(ramp)/Hz*(1:Hz))')<0 && RampUp<I+Hz
                    RampUp = RampUp+1;
                end
                DayDisp(RampUp-Hz+1:RampUp) = MinPow +(ramp)/Hz*(1:Hz);
                if RampUp<24/Ts
                    DayDisp(RampUp:24/Ts) = MinPow +ramp;
                end
            end
        end
        if day>1 && nnz(day-1 == holidays)==0 %% yesterday = weekday
            if week>1 && dayOfWeek == 7
                ramp = WeekendMin(week)-WeekdayMin(week-1);
            else ramp = WeekendMin(week)-WeekdayMin(week);
            end
            Hz = max(ceil(abs(ramp/(CHPsize*RampTs))),1);
            I = 1+24/Ts*(day-2);
            I = find(DemandE(I:I+48/Ts-1)<MinPow-(ramp),1);
            if isempty(I)
                 DayDisp(1:24/Ts) = MinPow -ramp;
                 holidays(find(max(0,holidays-day+1),1)) = 1;
            elseif I>24/Ts
                RampDown = I-24/Ts-1;
                while min(DemandE(RampDown+1+24/Ts*(day-1):RampDown+Hz+24/Ts*(day-1))- (MinPow -(ramp)/Hz*(Hz+1-(1:Hz)))')<0 && RampDown>(I-24/Ts-1)-Hz
                    RampDown = RampDown-1;
                end
                if RampDown>=0
                    DayDisp(1:RampDown) = MinPow -ramp;
                    DayDisp(RampDown+1:RampDown+Hz) = MinPow -ramp/Hz*(Hz+1-(1:Hz));
                else DayDisp(1:RampDown+Hz) = MinPow -ramp/Hz*((Hz+RampDown+1)-1:Hz+RampDown);
                end
            end
        end
    else DayDisp(1:24/Ts) = MaxPow;
        if day>1 && nnz(day-1 == holidays)==1 %% yesterday == weekend
            if week>1 && dayOfWeek == 7
                ramp = WeekdayMin(week)-WeekendMin(week-1);
            else ramp = WeekdayMin(week)-WeekendMin(week);
            end
            Hz = max(ceil(abs(ramp/(CHPsize*RampTs))),1);
            I = 1+24/Ts*(day-1);
            I = find(DemandE(I:I+12/Ts-1)<MaxPow,1,'last');
            if ~isempty(I)
                RampUp = I;
                while min(DemandE(RampUp-Hz+1+24/Ts*(day-1):RampUp+24/Ts*(day-1))- (MinPow +(ramp)/Hz*(1:Hz))')<0 && RampUp<I+Hz
                    RampUp = RampUp+1;
                end
                if RampUp-Hz>=0
                    DayDisp(1:RampUp-Hz) = MinPow;
                    DayDisp(RampUp-Hz+1:RampUp) = MinPow +(ramp)/Hz*(1:Hz);
                else DayDisp(1:RampUp) = MinPow +(ramp)/Hz*((1-RampUp+Hz):Hz);
                end
            end
        end
        if day>1 && nnz(day+1 == holidays)==1 %% tomorrow = weekend
            if week<=52 && dayOfWeek == 6
                ramp = WeekendMin(week+1)-WeekdayMin(week);
            else ramp = WeekendMin(week)-WeekdayMin(week);
            end
            Hz = max(ceil(abs(ramp/(CHPsize*RampTs))),1);
            I = 1+12/Ts+24/Ts*(day-1);
            I = find(DemandE(I:I+24/Ts-1)<MinPow-(ramp),1);
            if I<12/Ts
                RampDown = I+12/Ts-1;
                while min(DemandE(RampDown+1+24/Ts*(day-1):RampDown+Hz+24/Ts*(day-1))- (MinPow -(ramp)/Hz*(Hz+1-(1:Hz)))')<0 && RampDown>(I+12/Ts-1)-Hz
                    RampDown = RampDown-1;
                end
                DayDisp(RampDown+1:RampDown+Hz) = DayDisp(RampDown) +ramp/Hz*(1:Hz);
                if RampDown+Hz<24/Ts
                    DayDisp(RampDown+Hz+1:24/Ts) = MinPow;
                end
            end
        end
    end
    DGpower(1+24/Ts*(day-1):24/Ts*day) =  min(CHPsize,DayDisp(1:24/Ts)); 
   
    if dayOfWeek == 6
        week = week+1;
    end
end