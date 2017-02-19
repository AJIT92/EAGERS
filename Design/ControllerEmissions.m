function [DGpower]=ControllerEmissions(GHGCriteria,grid,DemandE,SysSize,RampRate,TurnDown,State)
% Set max system power for 5% of time with highest emission & min sys power
% for 5% of time with lowest grid emissions

stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};

StateNum = find(strcmp(State,stateName),1);
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};
State = stateAbrev(StateNum);

steps = length(DemandE);
Ts = 8760/steps;

%% import state emission profiles
[CO2, NOx, SO2] = AvgEmissionProfile(State,0,1);
switch GHGCriteria
    case 'CO2'
        Pollutant = CO2;
    case 'SO2'
        Pollutant = SO2;
    case 'NOx'
        Pollutant = NOx;
end
A = sort(Pollutant);
L = length(Pollutant);
MaxGrid = A(floor(.95*L));
MinGrid = A(floor(.05*L));
[gridRate, sellback] = ElecRate(grid);

%% Meet Electric Demand
CHPsize=sum(SysSize);
MinPower = sum(SysSize./TurnDown);

DemandProfile = MinPower+(CHPsize-MinPower)*max(0,min(1,(Pollutant-MinGrid)/(MaxGrid-MinGrid)));
DemandProfile2 = zeros(8760/Ts,1);
for i = 1:1:8760
    DemandProfile2(1+4*(i-1):4*i) = DemandProfile(i);
end
%%I am having it ignore the zero export constraint so that it follows the emission profile only
DemandProfile = DemandProfile2;
% DemandE2 = min(CHPsize,max(MinPower,DemandE));
% DemandProfile = min(DemandProfile2,DemandE2);

% %% the following is like the load following control
% DGpower = zeros(steps+1,1);
% DGpower(1) = DemandProfile(1);
% RampTs = sum(RampRate.*SysSize);
% Hz = min(linspace(8760/Ts-1,0,8760/Ts)',ceil(DemandProfile/RampTs));
% for i = 1:1:steps
%     if DemandProfile(i) > DGpower(i)
%         nextPow = min([CHPsize DGpower(i)+RampTs DemandProfile(i)]);
%     else nextPow = max([0 DGpower(i)-RampTs DemandProfile(i)]);
%     end
%     [Y, I] = min(DemandProfile(i:i+Hz(i))-(DGpower(i)-(0:1:Hz(i))'*RampTs));
%     if Y<0
%         DGpower(i+1) = max(0,DGpower(i)+Y);
%         for j = 1:1:I-1
%             i = i+1;
%             DGpower(i+1) = max(MinPower,DGpower(i)-RampTs);
%         end
%     else DGpower(i+1) = nextPow;
%     end
% 
% end
% DGpower = DGpower(2:end);

%% the following is like the diurnal dispatch control
DGpower = zeros(steps,1);
DayRate = zeros(24/Ts,1);
DaySellback = zeros(24/Ts,1);
RampTs = sum(RampRate.*SysSize);

morningMin = zeros(365,1);
afternoonMin = zeros(365,1);
dayMax = zeros(365,1);
for day = 1:1:365
    morningMin(day) = min(CHPsize,max(MinPower,min(DemandProfile(1+24/Ts*(day-1):12/Ts+24/Ts*(day-1)))));
    afternoonMin(day) = min(CHPsize,max(MinPower,min(DemandProfile(1+12/Ts+24/Ts*(day-1):24/Ts*day))));
    dayMax(day) = min(CHPsize,max(MinPower,max(DemandProfile(1+24/Ts*(day-1):24/Ts*day))));
end
morningMin(2:365,1) = min(morningMin(2:365),afternoonMin(1:364)); %% morning power setting 
afternoonMin(1:364,1) = morningMin(2:365);%% evening power setting
GenProf = zeros(10,24/Ts);
ElecValue = zeros(1,10);
for day = 1:1:365
    dayProf(1:24/Ts) = DemandProfile(1+24/Ts*(day-1):24/Ts*day);
    for j = 1:1:24
        DayRate(4*(j-1)+1:4*j) =  gridRate.CentskWh(j+24*(day-1));
        DaySellback(4*(j-1)+1:4*j) =  sellback.CentskWh(j+24*(day-1));
    end
    %% Find the best time to ramp up & down to maximize the value of electricity generated on the day
    for repeat = 1:1:3
        if repeat ==1
            guessMidday = linspace(morningMin(day),dayMax(day),10); %% midday power setting
        else guessMidday = linspace(guessMidday(max(1,Imax-1)),guessMidday(min(10,Imax+1)),10); %% midday power setting
        end
        RampUpTs = min(44,max(1,ceil(round(abs(guessMidday-morningMin(day))/RampTs*1e6)/1e6)));
        RampDownTs = min(40,max(1,ceil(round(abs(guessMidday-afternoonMin(day))/RampTs*1e6)/1e6)));
        for k = 1:1:length(guessMidday)
            tRampUpEnd = min(12/Ts,max(RampUpTs(k)+2,find(dayProf>=guessMidday(k),1,'first')));
            if isempty(tRampUpEnd)
                tRampUpEnd = 8/Ts;
            end
            RampUpProfile = zeros(1,RampUpTs(k))+morningMin(day); 
            for j = 2:1:RampUpTs(k)
                RampUpProfile(j) = RampUpProfile(j-1)+ RampTs;
            end
            %% now move ramp later if exporting during ramp
            if max(RampUpProfile-dayProf(tRampUpEnd-RampUpTs(k):tRampUpEnd-1))>0
                if grid.SellBackRate>0 % minimize exchange with grid
                    RampEndOptions = tRampUpEnd-RampUpTs(k):1:tRampUpEnd+RampUpTs(k);
                    NetExchange = zeros(length(RampEndOptions),1);
                    for rE = 1:1:length(RampEndOptions)
                        NetExchange(rE) = sum(abs(RampUpProfile-dayProf(RampEndOptions(rE)-RampUpTs(k):RampEndOptions(rE)-1)));
                    end
                    [minVal,Index] = min(NetExchange);
                    tRampUpEnd=min(12/Ts,RampEndOptions(Index));
                else %zero export
                    while max(RampUpProfile-dayProf(tRampUpEnd-RampUpTs(k):tRampUpEnd-1))>0 && tRampUpEnd<12/Ts
                        tRampUpEnd = tRampUpEnd+1;
                    end
                end
            end
            tRampDownStart = max(13/Ts,min(24/Ts-RampDownTs(k)-2,find(dayProf>=guessMidday(k),1,'last')));
            if isempty(tRampDownStart)
                tRampDownStart = 16/Ts;
            end
            RampDownProfile = zeros(1,RampDownTs(k))+afternoonMin(day); 
            for j = RampDownTs(k)-1:-1:1
                RampDownProfile(j) = RampDownProfile(j+1)  + RampTs;
            end
            %% now move ramp earlier if exporting during ramp
            if max(RampDownProfile-dayProf(tRampDownStart+1:tRampDownStart+RampDownTs(k)))>0
                if grid.SellBackRate>0 % minimize exchange with grid
                    RampStartOptions = tRampDownStart-RampDownTs(k):1:tRampDownStart+RampDownTs(k);
                    NetExchange = zeros(length(RampStartOptions),1);
                    for rS = 1:1:length(RampStartOptions)
                        NetExchange(rS) = sum(abs(RampDownProfile-dayProf(RampStartOptions(rS)+1:RampStartOptions(rS)+RampDownTs(k))));
                    end
                    [minVal,Index] = min(NetExchange);
                    tRampDownStart=max(13/Ts,RampStartOptions(Index));
                else %zero export
                    while max(RampDownProfile-dayProf(tRampDownStart+1:tRampDownStart+RampDownTs(k)))>0 && tRampDownStart>13/Ts
                        tRampDownStart = tRampDownStart-1;
                    end
                end
            end
            tRampUpStart = tRampUpEnd-RampUpTs(k);
            tRampDownEnd = tRampDownStart+RampDownTs(k);
            GenProf(k,1:24/Ts) =[zeros(1,tRampUpStart-1)+morningMin(day) RampUpProfile zeros(1,tRampDownStart-tRampUpEnd+1)+guessMidday(k) RampDownProfile zeros(1,24/Ts-tRampDownEnd)+afternoonMin(day)] ;
            GridElec = max(GenProf(k,1:24/Ts),0)';
            SellBack = max(-GenProf(k,1:24/Ts),0)';
            ElecValue(k) = sum(DayRate.*GridElec) - sum(DaySellback.*SellBack);
        end
        [MaxValue, Imax] = max(ElecValue);
    end
    DGpower(1+24/Ts*(day-1):24/Ts*day) = GenProf(Imax,1:24/Ts);
end