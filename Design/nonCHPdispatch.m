function [DemandE, DemandH,out] = nonCHPdispatch(systems,demand,grid,renew)
%% Load load demands
DemandE                         = demand.DemandE;
DemandH                         = demand.DemandH;
DemandC                         = demand.DemandC;
steps = length(DemandE);
Ts = 8760/steps;
%% Load Generator Features
nSys = length(systems.CHP);
SysSize = zeros(nSys,1);
SysRamp = zeros(nSys,1);
uiCHP = max(length(systems.CHP(1).SysCHP(:,1)));
SysCHP = zeros(uiCHP,2,nSys);
uiTable = length(systems.CHP(1).SysEff(:,1));
SysEff = zeros(uiTable,3,nSys);
CHPDerate = zeros(nSys,1);
CHPT = min(demand.CHPtemp(:));
CHPeff = zeros(nSys);
RampRate = zeros(nSys,1);
HeatCaptureG = zeros(length(DemandE),nSys);
Eff= zeros(length(DemandE),nSys);

for g= 1:nSys
    SysSize(g)              = systems.CHP(g).SysSize(1,1);
    SysRamp(g)              = systems.CHP(g).SysRamp;  % %rated power per second
    SysCHP(:,:,g)           = systems.CHP(g).SysCHP;
    SysEff(:,1:3,g)         = systems.CHP(g).SysEff(:,1:3);
    CHPDerate(g)            = interp1(SysCHP(:,1,g),SysCHP(:,2,g),CHPT);
    CHPeff(g)               = interp1(SysEff(:,1,g),SysEff(:,2,g),1); % 1 because of baseload
    RampRate(g)             = systems.CHP(g).SysRamp/100*3600*Ts;
    Eff(:,g)                = interp1(SysEff(:,1,g),SysEff(:,2,g),DemandE/max(DemandE));
    HeatCaptureG(:,g)       = (DemandE/max(DemandE))*SysSize(g)./Eff(:,g)*CHPDerate(g).*interp1(SysEff(:,2,g),SysEff(:,3,g),Eff(:,g));
end
HeatCapture = sum(HeatCaptureG,2);
%% Apply any renewable power as must-take power
if isfield(renew,'Solar')
    SolPow = zeros(length(DemandE),length(renew.Solar));
    for i = 1:1:length(renew.Solar)
        IrradNormalized = renew.Solar(i).Irrad/1000; %Solar irradiance normalized by standard test condition of 1000W/m2
        Az = renew.Solar(i).SunAz;
        Zen = renew.Solar(i).SunZen;
        if strcmp(renew.Solar(i).Tracking,'fixed')
            SolPow(:,i) = renew.Solar(i).Size*IrradNormalized.*cosd(Zen-renew.Solar(i).Tilt).*cosd(Az-renew.Solar(i).Azimuth)*renew.Solar(i).Eff;
        elseif strcmp(renew.Solar(i).Tracking,'1axis')
            SolPow(:,i) = renew.Solar(i).Size*IrradNormalized.*cosd(Zen-renew.Solar(i).Tilt)*renew.Solar(i).Eff;
        else SolPow(:,i) = renew.Solar(i).Size*IrradNormalized*renew.Solar(i).Eff;
        end
    end
    DemandE = DemandE-sum(SolPow,2);
    out.SolarGen = sum(SolPow,2);
end
if isfield(renew,'Wind')
    WindPow = zeros(length(DemandE),length(renew.Wind));
    for i = 1:1:length(renew.Wind)
        wind = renew.Wind(i);
        rho = 1.1798-1.3793e-4*wind.Elev+5.667e-9*wind.Elev^2;
        Wind = wind.Wind.*(wind.Wind>wind.CutIn).*(wind.Wind<wind.ShutDown);
        WindPow(:,i) = wind.Eff*.5*rho*(3.1416*(wind.Diam)^2/4).*Wind.^3;
    end
    DemandE = DemandE-sum(WindPow,2);
    out.WindGen = sum(WindPow,2);
end
%% Meet cooling demand
if isfield(systems,'Chiller')%chiller is part of DG system
    DemandE = demand.DemandE-demand.CoolingElectricalLoad; %Chiller repaces HVAC power
    if isfield(systems,'TES')%TES is part of DG system, decouples cooling demand & generation
        ColdTES = 0;
        for g = 1:1:length(systems.TES)
            if strcmp(systems.TES(g).TEStype,'cold') && systems.TES(g).SysSize>0
                ColdTES = 1;
            end
        end
        if ColdTES
            DemandC = ColdEnergyStorageShift(DemandE,DemandC, systems.Chiller, systems.TES,grid);
        end
    end
    Chiller = systems.Chiller;
    nChill = length(Chiller);
    ChillSize = zeros(nChill,1);
    ChillerCOP = zeros(nChill,1);
    uiTable2 = length(systems.Chiller(1).COPcurve(:,1));
    COPcurveA = zeros(uiTable2,nChill);
    COPcurveB = zeros(uiTable2,nChill);
    ChillCOP = zeros(steps,nChill);
    for g = 1:1:nChill
        ChillSize(g,1) = Chiller(g).SysSize;
        ChillerCOP(g)  = Chiller(g).COP;
        COPcurveA(:,g) = Chiller(g).COPcurve(:,1);
        COPcurveB(:,g) = Chiller(g).COPcurve(:,2);
    end
   
    name = cellstr(char(Chiller(:).ChillType));
    elec = strcmp(name,'electric');
    absorp = 1-elec;
    ElecChillCap = sum(ChillSize(:,1).*elec);
    AbsorpChillCap = sum(ChillSize(:,1).*absorp);
    AvgAbsorpCOP = sum(ChillerCOP.*absorp.*ChillSize(:,1)./max(AbsorpChillCap,.1));
    if AvgAbsorpCOP == 0 && AbsorpChillCap == 0
        Heat2driveChill = 0;
    else Heat2driveChill = AbsorpChillCap/AvgAbsorpCOP;
    end
    Heat2Achill = min(max(0,HeatCapture-DemandH), Heat2driveChill);
    AbsorpChill = min(DemandC, Heat2Achill*AvgAbsorpCOP*max(absorp)); %cooling met by absoprtion Chiller
    Cooling = max(min((DemandC-AbsorpChill),ElecChillCap),0);
    HVACload = max(((DemandC-AbsorpChill)-Cooling)./DemandC.*demand.CoolingElectricalLoad,0);
    PowerFracChill = AbsorpChill/max(AbsorpChillCap,.01)*absorp' +Cooling./max(ElecChillCap,.01)*elec'; %Operating state of each chiller
    for g = 1:1:nChill
        ChillCOP(:,g) = interp1(COPcurveA(:,g),COPcurveB(:,g),PowerFracChill(:,g)).*ChillerCOP(g);
    end
    ChillHeatCapt = PowerFracChill./ChillCOP*(ChillSize(:,1).*absorp);
    ChillLoadElec =  PowerFracChill./ChillCOP*(ChillSize(:,1).*elec);

    DemandE = DemandE+ChillLoadElec+HVACload;
    DemandH = DemandH+ChillHeatCapt;
end
out.DemandEshift = DemandE;
if isfield(systems,'Battery')
    DemandE = BatteryDemandShift(DemandE,systems.Battery);
end
out.DemandEshift2 = DemandE;

out.ChillerSize(1:2) = 0;
if isfield(systems,'Chiller')%chiller is part of DG system
    for i = 1:1:length(systems.Chiller)
        if strcmp(systems.Chiller(i).ChillType,'electric')
            out.ChillerSize(1) = out.ChillerSize(1)+systems.Chiller(i).SysSize;
        else out.ChillerSize(2) = out.ChillerSize(2)+systems.Chiller(i).SysSize;
        end
    end
end
out.TESsize = 0;
if isfield(systems,'TES')%TES is part of DG system, decouples cooling demand & generation
    for i = 1:1:length(systems.TES)
        if strcmp(systems.TES(i).TEStype,'cold')
            out.TESsize = out.TESsize+systems.TES(i).SysSize;
        end 
    end
end
out.BatterySize = 0;
if isfield(systems,'Battery')
    for i = 1:1:length(systems.Battery)
        out.BatterySize = out.BatterySize  + systems.Battery(i).Size;
    end
end
out.RenewableSize(1:2) = 0;
if isstruct(renew)
    if isfield(renew,'Solar')
        for i = 1:1:length(renew.Solar)
            out.RenewableSize(1) = out.RenewableSize(1) + renew.Solar(i).Size;
        end
    end
    if isfield(renew,'Wind')
        for i = 1:1:length(renew.Wind)
            out.RenewableSize(2) = out.RenewableSize(2) + renew.Wind(i).Size;
        end
    end
end