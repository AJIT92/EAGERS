function OptimalFCsizing(I, sizeChange)
% determines OptSize from SysSiz & vice versa for different control cases
global Project
demand=Project.Building; %the building load structure
systems=Project.System;
controlType=Project.Control.Name;
grid=Project.Utilities.Grid;

if isfield(Project,'Renewable')
    renew = Project.Renewable;
else renew = 1;
end
[DemandE, DemandH,out] = nonCHPdispatch(systems,demand,grid,renew);

steps = length(DemandE);
Ts = 8760/steps;    
[gridRate, sellBack] = ElecRate(grid);
%% look up electricity costs
SumStart = datenum(2013,grid.summerStartDate(1), grid.summerStartDate(2))-datenum(2013,1,1);
SumEnd = datenum(2013,grid.summerEndDate(1), grid.summerEndDate(2))-datenum(2013,1,1);

%% Load Generator Features
nSys = length(systems.CHP);
SysSize = zeros(nSys,1);
RampRate = zeros(nSys,1);
for g= 1:nSys
    SysSize(g)              = systems.CHP(g).SysSize(1,1);
    RampRate(g)             = systems.CHP(g).SysRamp/100*3600*Ts;
end
RampTs = sum(RampRate.*SysSize);

%Find holidays/weekends
ElecOnly = demand.DemandE-demand.CoolingElectricalLoad;
MidDayAvg = zeros(365,1);
for day = 1:1:365
    MidDayAvg(day) = mean(ElecOnly(1+11/Ts+24/Ts*(day-1):14/Ts+24/Ts*(day-1)));
end
holidays = [];
for i = 1:52
    lowDays =7*(i-1) + find(MidDayAvg(1+7*(i-1):7*i)<(median(MidDayAvg(1+7*(i-1):7*i))-.5*std(MidDayAvg(1+7*(i-1):7*i))));
    holidays(end+1:end+length(lowDays)) = lowDays;
end

%Find Minimum demand without weekends/Holidays
SumNonHoliday = zeros(24/Ts,1);
SummerHolidays = nnz((holidays>=SumStart).*(holidays<=SumEnd));
count = 0;
dayProf= zeros((SumEnd-SumStart+1)-SummerHolidays,24/Ts);
PartOfWeek = zeros(SumEnd-SumStart-SummerHolidays,1);
week = 1;
for day = SumStart:1:SumEnd  %%Size system for summer loads
    if nnz(day == holidays)==0;
        if day>=1  && nnz(day-1 == holidays)==1
            week = week+1;
        end
        SumNonHoliday = SumNonHoliday+DemandE(1+24/Ts*(day-1):24/Ts*day);
        count = count+1;
        dayProf(count,1:24/Ts) = DemandE(1+24/Ts*(day-1):24/Ts*day);
        PartOfWeek(count) = week;
    end
end
AvgDay = SumNonHoliday/count;

MinDemandE = min(min(dayProf(1:count,:)));
% %option 1:
% PeakAvg = max(AvgDay); %Peak power on an average summer day

%Option 2:
OrderedDemand = sort(DemandE);
PeakAvg = OrderedDemand(floor(.98*length(DemandE)));% Capacity necessary to meet demand 98% of the time

if strcmp(controlType,'3_DiurnalPeak') % find mid-day power setting for average summer day
    A = zeros(24/Ts,SumEnd-SumStart+1);
    B = zeros(24/Ts,SumEnd-SumStart+1);
    for day = 1:1:SumEnd-SumStart+1
        for j = 1:1:24
            A(4*(j-1)+1:4*j,day) = gridRate.CentskWh(24*(day-1)+j);
            B(4*(j-1)+1:4*j,day) = sellBack.CentskWh(24*(day-1)+j);
        end
    end
    DayRate = max(A,[],2);
    DaySellback = max(B,[],2);
    morningMin = min(AvgDay(1:12/Ts));
    afternoonMin = min(AvgDay(12/Ts+1:24/Ts));
    dayMax = max(AvgDay);
    GenProf = zeros(10,24/Ts);
    ElecValue = zeros(10,1);
    %% Find the best time to ramp up & down to maximize the value of electricity generated on the day
    for repeat = 1:1:3
        if repeat ==1
            guessMidday = linspace(morningMin,dayMax,10); %% midday power setting
        else guessMidday = linspace(guessMidday(max(1,Imax-1)),guessMidday(min(10,Imax+1)),10); %% midday power setting
        end
        RampUpTs = min(44,max(1,ceil(round(abs(guessMidday-morningMin)/RampTs*1e6)/1e6)));
        RampDownTs = min(40,max(1,ceil(round(abs(guessMidday-afternoonMin)/RampTs*1e6)/1e6)));
        for k = 1:1:length(guessMidday)
            tRampUpEnd = min(12/Ts,max(RampUpTs(k)+2,find(AvgDay>=guessMidday(k),1,'first')));
            if isempty(tRampUpEnd)
                tRampUpEnd = 8/Ts;
            end
            tRampDownStart = max(13/Ts,min(24/Ts-RampDownTs(k)-2,find(AvgDay>=guessMidday(k),1,'last')));
            if isempty(tRampDownStart)
                tRampDownStart = 16/Ts;
            end
            RampUpProfile = zeros(1,RampUpTs(k))+morningMin; 
            for j = 2:1:RampUpTs(k)
                RampUpProfile(j) = RampUpProfile(j-1)+ RampTs;
            end
            RampDownProfile = zeros(1,RampDownTs(k))+afternoonMin; 
            for j = RampDownTs(k)-1:-1:1
                RampDownProfile(j) = RampDownProfile(j+1)  + RampTs;
            end
            tRampUpStart = tRampUpEnd-RampUpTs(k);
            tRampDownEnd = tRampDownStart+RampDownTs(k);
            GenProf(k,1:24/Ts) =[zeros(1,tRampUpStart-1)+morningMin RampUpProfile zeros(1,tRampDownStart-tRampUpEnd+1)+guessMidday(k) RampDownProfile zeros(1,24/Ts-tRampDownEnd)+afternoonMin] ;
            GridElec = max(GenProf(k,1:24/Ts),0)';
            SellBack = max(-GenProf(k,1:24/Ts),0)';
            ElecValue(k) = sum(DayRate.*GridElec) - sum(DaySellback.*SellBack);
            if grid.SellBackRate==0
                if guessMidday(k)>min(AvgDay(tRampUpEnd+1:tRampDownStart-1))
                  ElecValue(k) = 0;
                end
            end
        end
        [MaxValue, Imax] = max(ElecValue);
    end
    PeakAvg = guessMidday(Imax);
elseif strcmp(controlType,'2_WeekendDip') %find the highest weekday setting using this dispatch strategy
    % This should be equal to the minimum between noon and noon on
    % consecutive non-holiday days.
    PeakAvg = 0;
    for holiday = 1:1:length(holidays)-1  %%Size system for summer loads
        if holidays(holiday)>SumStart
            day1 = holidays(holiday);
            day2 = holidays(holiday+1);
            if day2-day1>2
                Profile = DemandE(12/Ts+24/Ts*(day1):12/Ts+24/Ts*(day2-2));
                PeakAvg = max(PeakAvg,min(Profile));
            end
        end
    end
end


%check to see if the peak to minimum exceeds turndown ratio
CHPSize = 0;
minSize = 0;
turndown = zeros(length(systems.CHP),1);
for i = 1:1:length(systems.CHP)
    CHPSize = CHPSize+systems.CHP(i).SysSize(1);
    minSize = minSize+systems.CHP(i).SysSize(2);
    turndown(i) = systems.CHP(i).TurnDown;
end
TurnDown = CHPSize/minSize;

if PeakAvg/MinDemandE>TurnDown %if the ratio does exceed that specifed by the CHP eithier increase the minumum or decrase the maximum depending upon the sellback ability
    if Project.Utilities.Grid.SellBackRate>0
        MinDemandE = PeakAvg/TurnDown;
    else
        PeakAvg = MinDemandE*TurnDown;
    end  
end
if strcmp(controlType,'1_BaseLoad') %Reset upper limit if using baseload control
    PeakAvg = MinDemandE;
end

if strcmp(sizeChange,'editMaxPower') %Find OptSize from SysSize
    OptSize = CHPSize/PeakAvg;
    systems.CHP(I).OptSize = OptSize*100;
elseif strcmp(sizeChange,'OptSize') %Find SysSize from OptSize
    OptSize = systems.CHP(I).OptSize/100;
    for i = 1:1:length(systems.CHP)
        systems.CHP(i).SysSize(1) = PeakAvg*OptSize*(systems.CHP(i).SysSize(1)/CHPSize);
        systems.CHP(i).SysSize(2) = systems.CHP(i).SysSize(1)/turndown(i);%set min FC size equal to 10% max
        systems.CHP(i).OptSize = OptSize*100;  
    end
end

totalSize = 0;
for i=1:1:length(systems.CHP)
    totalSize = totalSize + systems.CHP(i).SysSize;
end

Project.System.CHP = systems.CHP;
Project.Economic.LifekWh = CHPSize*8760*Project.Economic.LifeYrs;
