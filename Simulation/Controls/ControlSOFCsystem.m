function Out = ControlSOFCsystem(t,Y, Inlet,block,string1)
% Controls for SOFC system, control blower power, heater bypass, anode recirculation and fuel flow rate
% Four (4) inlets: T hot, T cold, T average, Voltage
% Five (5) outlets:  Heater bypass, blower power, fuel flow rate, current, anode recirculation
% Two (2) states: Heater bypass, blower power
%Need to add state for current back in to avoid fuel starvation during step changes
global Tags F
I_Gain = block.Gain.*block.Scale;
TinletError = ((mean(Inlet.Cold) + .5*block.Target(2) + block.dT_cath_PEN) - block.Target(1))/block.Target(2); %target a fixed inlet temperature
ToutletError = ((mean(Inlet.Hot) - .5*block.Target(2) + block.dT_cath_PEN) - block.Target(1))/block.Target(2); %target a fixed outlet temperature


% deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
% dTerror =(deltaT-block.Target(2))/block.Target(2); %target a fixed temperature gradient
% BlowerPower = Y(2)+(dTerror*block.PropGain(2))*block.Scale(2);

BlowerPower = Y(2)+(ToutletError*block.PropGain(2))*block.Scale(2);
StackPower = PowerDemandLookup(t) + BlowerPower;
% Current = Y(3)+block.PropGain(3)*(StackPower*1000/(Inlet.Voltage*block.Cells) - Y(3));
% Power = Current*Inlet.Voltage*block.Cells/1000;
% PowerError = (StackPower - Power)/StackPower;

Current = StackPower*1000/(Inlet.Voltage*block.Cells);

Bypass = min(1,max(0,Y(1)+TinletError*block.PropGain(1)));
Utilization = min(block.Utilization,block.Utilization + .1*(Y(1)+TinletError*block.PropGain(1))); %lower utilization if additional pre-heating is necessary
FuelFlow = block.Cells*Current/(2*F*Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000;

if block.AnodeRecirc.IC == 0
    Recirculation = 0;
else
    Steam2Carbon = block.Target(3);
    Recirculation = anodeRecircHumidification(block.Fuel,FuelFlow,block.WGSeffective,Steam2Carbon,block.Cells*Current/(2*F*1000),block.AnodeRecirc.IC);
end
    
if strcmp(string1,'Outlet')
    Out.HeaterBypass = Bypass;
    Out.Blower = BlowerPower;
    Out.Current = Current;
    Out.AnodeRecirc = Recirculation;
    Out.FuelFlow = FuelFlow;
    Tags.(block.name).Bypass = Bypass;
    Tags.(block.name).Blower = BlowerPower;
    Tags.(block.name).Recirculation = Recirculation;
    Tags.(block.name).FuelFlow = FuelFlow;
    Tags.(block.name).Current = Current;
    Tags.(block.name).Utilization = Utilization;
elseif strcmp(string1,'dY')
    dY = Y*0;
    if Y(1)+TinletError*block.PropGain(1) < 0 
        dY(1) = 0.02*I_Gain(1)*TinletError; %controlling utilization
    else
        dY(1) = I_Gain(1)*TinletError;
    end
    %% need to deterimine saturation on blower
    if Y(1) <0
        dY(2) = 0;
    else
        dY(2) = I_Gain(2)*ToutletError;
    end

%     dY(3) = I_Gain(3)*PowerError;
    Out = dY;  
end