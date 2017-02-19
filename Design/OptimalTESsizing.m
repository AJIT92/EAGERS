function OptimalTESsizing(sizeChange)
% determines OptSize from SysSiz & vice versa for different control cases
global Project
demand=Project.Building; %the building load structure
systems=Project.System;
grid=Project.Utilities.Grid;
TES = Project.System.TES;
Chiller = Project.System.Chiller;
DemandH                         = demand.DemandH;
DemandC                         = demand.DemandC;

steps = length(DemandC);
Ts = 8760/steps;


kWh2GalDeg = 3600/4.186*.264; %3600 seconds / 4.186 kJ/kg*K * .264 gal/L

Year=2006;
tmx = [ones(steps,1)*[Year 1 1]  [0:Ts:8760-Ts]' ones(steps,1)*[0 0]];
date = datenum(tmx);
SummerStart = datenum([num2str(grid.summerStartDate(1)),'/',num2str(grid.summerStartDate(2)),'/',num2str(Year)]);
SummerEnd = datenum([num2str(grid.summerEndDate(1)),'/',num2str(grid.summerEndDate(2)),'/',num2str(Year)]);
SumOn = find(grid.summerRateTable(2,:)-1,1,'first')/24;
SumOff = find(grid.summerRateTable(2,:)-1,1,'last')/24;
if isempty(SumOn)
    SumOn = 9/24;
    SumOff = 20/24;
end
summer = date>=SummerStart & date<=(SummerEnd+.99);
t_day = linspace(Ts,24,24/Ts)/24;
SumMidOn = t_day>=SumOn & t_day<=SumOff;
DemandCSummer = summer.*DemandC;

PeakCool = zeros(365,1);
for day = 1:1:365
    PeakCool(day) = sum(DemandCSummer(1+24/Ts*(day-1):24/Ts*day).*SumMidOn')*Ts;
end

%% Load CHP Features
nSys = length(systems.CHP);
CHP_SysSize = zeros(nSys,1);
uiTable = length(systems.CHP(1).SysEff(:,1));
CHP_SysEff = zeros(uiTable,3,nSys);
uiCHP = max(length(systems.CHP(1).SysCHP(:,1)));
SysCHP = zeros(uiCHP,2,nSys);
CHPDerate = zeros(1,nSys);
HeatCaptureG = zeros(nSys,1);
CHP_NominalEff= zeros(nSys,1);
CHPT = min(demand.CHPtemp(:));
for g= 1:nSys
    CHP_SysSize(g)            = systems.CHP(g).SysSize(1,1);
    CHP_SysEff(:,1:3,g)       = systems.CHP(g).SysEff(:,1:3);
    SysCHP(:,:,g)             = systems.CHP(g).SysCHP;
    CHPDerate(g)              = interp1(SysCHP(:,1,g),SysCHP(:,2,g),CHPT);
    CHP_NominalEff(g)         = CHP_SysEff(find(CHP_SysEff(:,1,g)-.999,1,'first'),2,g);
    HeatCaptureG(g)          = CHP_SysSize(g)*(1/CHP_NominalEff(g)-1)*CHPDerate(g);
end
HeatCapture = sum(HeatCaptureG(1:nSys));

%% Determine existing chiller size
nChill = length(Chiller);
ChillSize = zeros(nChill,1);
ChillerCOP = zeros(nChill,1);
for g = 1:1:nChill
    ChillSize(g)   = Chiller(g).SysSize;
    ChillerCOP(g)  = Chiller(g).COP;
end  
name = cellstr(char(Chiller(:).ChillType));
elec = strcmp(name,'electric');
absorp = 1-elec;
ElecChillCap = sum(ChillSize.*elec);
AbsorpChillCap = sum(ChillSize.*absorp);

%% Determine Existing TES capacity
TESsize = zeros(length(TES),1);
OptSize = zeros(length(TES),1);
SizeHours = zeros(length(TES),1);
name = cellstr(char(TES(:).TEStype));
cold = strcmp(name,'cold');
hot = 1-cold;
for i = 1:1:length(TES)
    TESsize(i)   = TES(i).SysSize;
    OptSize(i)   = TES(i).OptSize/100;
    SizeHours(i) = TES(i).SizeHours;
end
coldTEScap = sum(TESsize.*cold);
hotTEScap = sum(TESsize.*hot);


%% 12 hour heat demand
DemandH12 = zeros(steps,1);
for t_i = 1:1:(365*24/Ts)-(12/Ts-1)
    DemandH12(t_i) = sum(DemandH(t_i:t_i+(12/Ts-1)));
end

if strcmp(sizeChange,'editSize') %Find OptSize & sizeHours from SysSize
    for i = 1:1:length(TES)
    if strcmp(TES(i).TEStype,'cold')
        %Optimal size is the largest summer cooling day (sum of on-peak cooling)
        OptSize(i) = coldTEScap/max(PeakCool);%determine % optimal size
        %# of hours depends on chilller capcaity
        SizeHours(i) = coldTEScap/(ElecChillCap+AbsorpChillCap);
    elseif strcmp(TES(i).TEStype,'hot')
        %set TES size equal to 12 hours of max available heat recovery
        OptSize(i) = hotTEScap/(HeatCapture*12);
        %Reduce size if this is bigger than any 12hour heat demand period
        if (HeatCapture*12)>max(DemandH12)
            OptSize(i) = hotTEScap/max(DemandH12);
        end
        % # of hours also depends on heat capture
        SizeHours(i) = hotTEScap/HeatCapture;
    end
    end
    TES(i).OptSize = OptSize(i)*100;
    TES(i).SizeHours = SizeHours(i);
elseif strcmp(sizeChange,'OptSize') %Find SysSize & sizeHours from OptSize
    for i = 1:1:length(TES)
        if strcmp(TES(i).TEStype,'cold')
             %Optimal size is the largest summer cooling day (sum of on-peak cooling)
            TESsize(i) = OptSize(i)*max(PeakCool)/sum(cold);%determine % optimal size
            %# of hours depends on chilller capcaity
            SizeHours(i) = TESsize(i)/(ElecChillCap+AbsorpChillCap);
        elseif strcmp(TES(i).TEStype,'hot')
            TESsize(i) = OptSize(i)*(HeatCapture*12)/sum(hot);
            %Reduce size if this is bigger than any 12hour heat demand period
            if (HeatCapture*12)>max(DemandH12)
                TESsize(i) = OptSize(i)*max(DemandH12)/sum(hot);
            end
            SizeHours(i) = OptSize(i)*12/sum(hot);
        end
        TES(i).SysSize = TESsize(i);
        TES(i).SizeHours = SizeHours(i);
    end
elseif strcmp(sizeChange,'editSizeHours') %Find SysSize & OptSize from sizeHours
    for i = 1:1:length(TES)
        if strcmp(TES(i).TEStype,'cold')
            TESsize(i) = (ElecChillCap+AbsorpChillCap)*SizeHours(i);
            %Optimal size is the largest summer cooling day (sum of on-peak cooling)
            OptSize(i) = coldTEScap/max(PeakCool)*(TESsize(i)/TES(i).SysSize);%determine % optimal size
        elseif strcmp(TES(i).TEStype,'hot')
            %set TES size equal to 12 hours of max available heat recovery
            TESsize(i) = HeatCapture*SizeHours(i);
            %set TES size equal to 12 hours of max available heat recovery
            OptSize(i) = hotTEScap/(HeatCapture*12)*(TESsize(i)/TES(i).SysSize);
            %Reduce size if this is bigger than any 12hour heat demand period
            if (HeatCapture*12)>max(DemandH12)
                OptSize(i) = hotTEScap/max(DemandH12)*(TESsize(i)/TES(i).SysSize);
            end
        end
        TES(i).SysSize = TESsize(i);
        TES(i).OptSize = OptSize(i)*100;
    end
end
for i = 1:1:length(TES)
    Tcold = TES(i).Tcold;
    Thot = TES(i).Thot;
    TES(i).SizeGal = TES(i).SysSize*kWh2GalDeg/(Thot-Tcold);
end
Project.System.TES = TES;