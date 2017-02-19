function Out=dispatchSystem(systems,demand,grid,renew,control,state)

[DemandE, DemandH,Out.eOut] = nonCHPdispatch(systems,demand,grid,renew);
NonDistH                        = demand.NonDistH;
steps = length(DemandE);
Ts = 8760/steps;

%% Load Generator Features
nSys = length(systems.CHP);
SysSize = zeros(nSys,1);
uiCHP = max(length(systems.CHP(1).SysCHP(:,1)));
SysCHP = zeros(uiCHP,2,nSys);
uiTable = length(systems.CHP(1).SysEff(:,1));
SysEff = zeros(uiTable,3,nSys);
CHPDerate = zeros(nSys,1);
CHPT = min(demand.CHPtemp(:));
RampRate = zeros(nSys,1);
TurnDown = zeros(nSys,1);
for g= 1:nSys
    SysSize(g)              = systems.CHP(g).SysSize(1,1);
    SysCHP(:,:,g)           = systems.CHP(g).SysCHP;
    SysEff(:,1:6,g)         = systems.CHP(g).SysEff(:,1:6);
    CHPDerate(g)            = interp1(SysCHP(:,1,g),SysCHP(:,2,g),CHPT);
    RampRate(g)             = systems.CHP(g).SysRamp/100*3600*Ts;
    TurnDown(g)             = systems.CHP(g).TurnDown;
end

ElecOnly = demand.DemandE-demand.CoolingElectricalLoad;

%Dispatch the CHP systems
if strcmp(control.Name,'1_BaseLoad')
    [DGpower] = ControllerBaseLoad(DemandE,SysSize);
elseif strcmp(control.Name,'2_WeekendDip')
    [DGpower] = ControllerWeekendDip(DemandE,SysSize,RampRate,ElecOnly);
elseif strcmp(control.Name,'3_DiurnalPeak')
    [DGpower] = ControllerDailyPeak(grid,DemandE,SysSize,RampRate,TurnDown);
elseif strcmp(control.Name,'4_LoadFollow')
    [DGpower] = ControllerLoadFollow(DemandE,SysSize,RampRate,TurnDown);
elseif strcmp(control.Name,'5_LoadFollow_NonPredict')
    [DGpower] = ControllerLoadFollow2(DemandE,SysSize,RampRate,TurnDown);
elseif strcmp(control.Name,'6_CO2EmissionControl')
    [DGpower] = ControllerEmissions('CO2',grid,DemandE,SysSize,RampRate,TurnDown,state);
elseif strcmp(control.Name,'7_SO2EmissionControl')
    [DGpower] = ControllerEmissions('SO2',grid,DemandE,SysSize,RampRate,TurnDown,state);
elseif strcmp(control.Name,'8_NOxEmissionControl')
    [DGpower] = ControllerEmissions('NOx',grid,DemandE,SysSize,RampRate,TurnDown,state);
    % elseif strcmp(control.Name,'6_CO2EmissionControl')
    %     [DGpower] = ControllerCO2Emissions(grid,DemandE,SysSize,RampRate,TurnDown,state);
    % elseif strcmp(control.Name,'7_SO2EmissionControl')
    %     [DGpower] = ControllerSO2Emissions(grid,DemandE,SysSize,RampRate,TurnDown,state);
    % elseif strcmp(control.Name,'8_NOxEmissionControl')
    %     [DGpower] = ControllerNOxEmissions(grid,DemandE,SysSize,RampRate,TurnDown,state);
else 
    error(['No control defined for: ' control.Name])
end

%%% Dispatch generators from most to least efficient
PowerFracDG = zeros(steps,nSys);
EProduced = zeros(steps,nSys);
FuelUsedDG = zeros(steps,nSys);
CHPheatRecover = zeros(steps,nSys);
A = SysEff(end,1,:);
[Y, I] = sort(A);
for i = 1:1:steps
    j = 1;
    size = 0;
    while size+SysSize(I(j))<DGpower(i)
        size = size+SysSize(I(j));
        PowerFracDG(i,I(j))=1;
        j = j+1;
        if numel(I)<j
            j = j-1;
            break;
        end
    end
    PowerFracDG(i,I(j))=(DGpower(i)-size)/SysSize(I(j));
    if PowerFracDG(i,I(j))<(1/TurnDown(I(j))); %enforce the turndown constraint
        if j>1
            shed = SysSize(I(j))*(PowerFracDG(i,I(j))-(1/TurnDown(I(j))));
            PowerFracDG(i,I(j-1)) = max(1/TurnDown(I(j-1)),(SysSize(I(j-1))-shed)/SysSize(I(j-1)));
            PowerFracDG(i,I(j))=(1/TurnDown(I(j)));
        else % turn CHP off or leave at min power
            if PowerFracDG(i,I(j))<.5*(1/TurnDown(I(j)))
                PowerFracDG(i,I(j))=0;
            else PowerFracDG(i,I(j))=(1/TurnDown(I(j)));
            end
        end
    end
end

for g = 1:nSys
    EProduced(:,g) = PowerFracDG(:,g).*SysSize(g).*Ts; %Energy produced by each DG system (kWh)  
    FuelUsedDG(:,g)   = EProduced(:,g)./interp1(SysEff(:,1,g),SysEff(:,2,g),PowerFracDG(:,g)); % power/efficiency = fuel use in kWh
    CHPheatRecover(:,g) = FuelUsedDG(:,g).*CHPDerate(g).*interp1(SysEff(:,1,g),SysEff(:,3,g),PowerFracDG(:,g));
    
    CO2rate = interp1(SysEff(:,1,g),SysEff(:,4,g),PowerFracDG(:,g));
    CO2(:,g) = CO2rate.*PowerFracDG(:,g)*SysSize(g)*Ts./1000;%lbsC02/MWhr*kW*.25hr*(1MWhr/1000kWhr) ->lbsCO2
    NOxrate = interp1(SysEff(:,1,g),SysEff(:,5,g),PowerFracDG(:,g));
    NOx(:,g) = NOxrate.*PowerFracDG(:,g)*SysSize(g)*Ts./1000;
    SO2rate = interp1(SysEff(:,1,g),SysEff(:,6,g),PowerFracDG(:,g));
    SO2(:,g) = SO2rate.*PowerFracDG(:,g)*SysSize(g)*Ts./1000;
end
CO2 = sum(CO2,2);
NOx = sum(NOx,2);
SO2 = sum(SO2,2);

CO2Total = sum(CO2);
NOxTotal = sum(NOx);
SO2Total = sum(SO2);

PeakBurnerHeat = max((DemandH-sum(CHPheatRecover,2)),0);     
% Report out key economic values
GridPurchaseHour = zeros(8760,1);
BoilerHeatHour = zeros(8760,1);
NOxhour = zeros(8760,1);
CO2hour = zeros(8760,1);
SO2hour = zeros(8760,1);

for i = 1:1:8760
    GridPurchaseHour(i) = sum(DemandE(1+1/Ts*(i-1):1/Ts*i)*Ts)-sum(sum(EProduced(1+1/Ts*(i-1):1/Ts*i,:),2));
    BoilerHeatHour(i) = sum(PeakBurnerHeat(1+1/Ts*(i-1):1/Ts*i)*Ts)+sum(NonDistH(1+1/Ts*(i-1):1/Ts*i)*Ts);
    CO2hour(i) = sum(CO2(1+1/Ts*(i-1):1/Ts*i)*Ts);
    SO2hour(i) = sum(SO2(1+1/Ts*(i-1):1/Ts*i)*Ts);
    NOxhour(i) = sum(NOx(1+1/Ts*(i-1):1/Ts*i)*Ts);
end
Out.eOut.Grid_purchases_hourly_kWh = GridPurchaseHour;

Out.eOut.GenPower = sum(EProduced,2)/Ts;
Out.eOut.Total_Electricity_produced_kWh = sum(sum(EProduced,2));
Out.eOut.SelfGen = Out.eOut.Total_Electricity_produced_kWh/sum(demand.DemandE*Ts);
Out.eOut.Total_DG_fuel_used_kWh = sum(sum(FuelUsedDG,2));
Out.eOut.Total_peak_burner_heat_out_kWh = sum(PeakBurnerHeat)*Ts + sum(NonDistH)*Ts;
Out.eOut.Fuel = sum(FuelUsedDG,2)+PeakBurnerHeat*Ts + NonDistH*Ts; %fuel used in kWh (no need to multiply by Ts later
Out.eOut.BoilerHeatHour = BoilerHeatHour; %hourly heat not met by CHP
Out.eOut.Total_DG_Heat_out_kWh =  sum(DemandH-PeakBurnerHeat); %DG heat captured & utilized by building  
Out.eOut.SysSize = sum(SysSize);
Out.eOut.CO2Total = CO2Total;%lbs of CO2 produced
Out.eOut.NOxTotal = NOxTotal;
Out.eOut.SO2Total = SO2Total;
Out.eOut.CO2 = CO2hour;%lbs of CO2 produced
Out.eOut.NOx = NOxhour;
Out.eOut.SO2 = SO2hour;