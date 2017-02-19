function OptimalChillerSizing(sizeChange)
% determines OptSize from SysSiz & vice versa for different control cases
global Project
Chiller = Project.System.Chiller;
demand=Project.Building; %the building load structure
systems=Project.System;
grid = Project.Utilities.Grid;
DemandC = demand.DemandC;
DemandE = demand.DemandE;  
DemandH = demand.DemandH; 

%% Load Generator Features
nSys = length(systems.CHP);
CHP_SysSize = zeros(nSys,1);
uiTable = length(systems.CHP(1).SysEff(:,1));
CHP_SysEff = zeros(uiTable,3,nSys);
uiCHP = max(length(systems.CHP(1).SysCHP(:,1)));
SysCHP = zeros(uiCHP,2,nSys);
CHPDerate = zeros(1,nSys);
if strcmp(Project.Control.Name,'1_BaseLoad')
    HeatCaptureG = zeros(1,nSys);
else HeatCaptureG = zeros(length(DemandE),nSys);
    Eff= zeros(length(DemandE),nSys);
end
CHP_NominalEff= zeros(nSys,1);
CHPT = min(demand.CHPtemp(:));
for g= 1:nSys
    CHP_SysSize(g)            = systems.CHP(g).SysSize(1,1);
    CHP_SysEff(:,1:3,g)       = systems.CHP(g).SysEff(:,1:3);
    SysCHP(:,:,g)             = systems.CHP(g).SysCHP;
    CHPDerate(g)              = interp1(SysCHP(:,1,g),SysCHP(:,2,g),CHPT);
    CHP_NominalEff(g)         = CHP_SysEff(find(CHP_SysEff(:,1,g)-.999,1,'first'),2,g);
    if strcmp(Project.Control.Name,'1_BaseLoad')
        HeatCaptureG(1,g)     = CHP_SysSize(g)*(1/CHP_NominalEff(g)-1)*CHPDerate(g);
    else 
        Eff(:,g)              = interp1(CHP_SysEff(:,1,g),CHP_SysEff(:,2,g),DemandE/max(DemandE));
        HeatCaptureG(:,g)     = CHP_SysSize(g)*(DemandE/max(DemandE)).*(1/Eff(:,g)-1)'*CHPDerate(g);
    end
    
end
HeatCapture = sum(HeatCaptureG(:,1:nSys));

%% TES
DemandE = DemandE-demand.CoolingElectricalLoad; %Chiller repaces HVAC power
if isfield(Project,'Renewable')
    renew = Project.Renewable;
    if isfield(renew,'Solar')
        SolPow = zeros(length(DemandE),length(renew.Solar));
        for i = 1:1:length(renew.Solar)
            Irrad = renew.Solar(i).Irrad;
            Az = renew.Solar(i).SunAz;
            Zen = renew.Solar(i).SunZen;
            if strcmp(renew.Solar(i).Tracking,'fixed')
                SolPow(:,i) = renew.Solar(i).Size*Irrad.*cosd(Zen-renew.Solar(i).Tilt).*cosd(Az-renew.Solar(i).Azimuth)*renew.Solar(i).Eff;
            elseif strcmp(renew.Solar(i).Tracking,'1axis')
                SolPow(:,i) = renew.Solar(i).Size*Irrad.*cosd(Zen-renew.Solar(i).Tilt)*renew.Solar(i).Eff;
            else SolPow(:,i) = renew.Solar(i).Size*Irrad*renew.Solar(i).Eff;
            end
        end
        DemandE = DemandE-sum(SolPow,2);
    end
    %% Add wind Power
end
if isfield(systems,'TES')%TES is part of DG system, decouples cooling demand & generation
    TESsize = 0;
    for i = 1:1:length(systems.TES)
        TESsize = TESsize + systems.TES(i).SysSize;
    end
    if TESsize>0
        if strcmp(cellstr(char(systems.TES(:).TEStype)),'cold')
            DemandC = ColdEnergyStorageShift(DemandE,DemandC, systems.Chiller, systems.TES,grid);
        end
    end
end

%% Determine existing chiller size
nChill = length(Chiller);
ChillSize = zeros(nChill,1);
ChillerCOP = zeros(nChill,1);
OptSize  = zeros(nChill,1);
for g = 1:1:nChill
    ChillSize(g)   = Chiller(g).SysSize;
    if ChillSize(g) ==0
        ChillSize(g) =1;
    end
    ChillerCOP(g)  = Chiller(g).COP;
    OptSize(g)     = Chiller(g).OptSize/100;
end  
name = cellstr(char(Chiller(:).ChillType));
elec = strcmp(name,'electric');
absorp = 1-elec;
ElecChillCap = sum(ChillSize.*elec);
AbsorpChillCap = sum(ChillSize.*absorp);
AvgAbsorpCOP = sum(ChillerCOP.*absorp.*ChillSize./max(AbsorpChillCap,.01));
if AvgAbsorpCOP == 0 && AbsorpChillCap == 0
    Heat2driveChill = 0;
else Heat2driveChill = AbsorpChillCap/AvgAbsorpCOP;
end
%base load heat available
Heat2Achill = min(max(HeatCapture-DemandH,0), Heat2driveChill);
%cooling met by absoprtion Chiller
AbsorpChill = min(DemandC, Heat2Achill*AvgAbsorpCOP*max(absorp)); 

if strcmp(sizeChange,'editSize') %Find OptSize from SysSize
    for i = 1:1:length(Chiller)
        if strcmp(Chiller(i).ChillType,'electric') %% Optimally size electric chillers to meet remaining load
            if max(DemandC-AbsorpChill)<=0
                OptSize(i) = ElecChillCap/max(DemandC); %Optimal size meets 50% of peak cooling demand
            else OptSize(i) = ElecChillCap/max(DemandC-AbsorpChill); %Optimal size would meet cooling demand not met by absorption chiller
            end
        elseif strcmp(Chiller(i).ChillType,'absorption') %% Optimally size any absorption chillers to use all the heat
            OptSize(i) = max(Heat2driveChill/max(HeatCapture-DemandH),AbsorpChillCap/max(DemandC));
        end
        Chiller(i).OptSize = OptSize(i)*100;
    end
elseif strcmp(sizeChange,'OptSize') %Find SysSize from OptSize
    OptSizeOld = zeros(length(Chiller),1);
    for i = 1:1:length(Chiller)
        if strcmp(Chiller(i).ChillType,'electric')
            if max(DemandC-AbsorpChill)<=0
                OptSizeOld(i) = ElecChillCap/max(DemandC); %Optimal size meets 50% of peak cooling demand
            else OptSizeOld(i) = ElecChillCap/max(DemandC-AbsorpChill); %Optimal size would meet cooling demand not met by absorption chiller
            end
            Frac = (OptSize(i)-OptSizeOld(i))+((Chiller(i).SysSize/ElecChillCap)*OptSizeOld(i));
            Chiller(i).SysSize = max(DemandC-AbsorpChill)*Frac;%set chiller size equal to max cooling load
        elseif strcmp(Chiller(i).ChillType,'absorption')
            OptSizeOld(i) = max(Heat2driveChill/max(HeatCapture-DemandH),AbsorpChillCap/max(DemandC));
            Frac = (OptSize(i)-OptSizeOld(i))+((Chiller(i).SysSize/AbsorpChillCap)*OptSizeOld(i));
            %set chiller size equal to smaller of max of available heat recovery or building load
            PeakHeatCapSize = max(HeatCapture-DemandH)*ChillerCOP(i)*Frac;%set absoption chiller output size using max available heat recovery
            PeakCoolSize = max(DemandC)*Frac; %Set to building demand
            Chiller(i).SysSize = min(PeakCoolSize,PeakHeatCapSize);
        end
    end
end
Project.System.Chiller = Chiller;