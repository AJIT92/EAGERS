function Out = ControlSOFCsystem(t,Y, Inlet,block,string1)
% Controls for SOFC system, control blower power, heater bypass, anode recirculation and fuel flow rate
% Four (4) inlets: T hot, T cold, T average, Voltage
% Five (5) outlets:  Heater bypass, blower power, fuel flow rate, current, anode recirculation
% Three (3) states: Heater bypass, blower power, Current 
global Tags F
NetPower = PowerDemandLookup(t);
averageT = mean(Inlet.PEN_Temp);
TavgError = (averageT - block.Target(1))/block.Target(2);
deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
dTerror =(deltaT-block.Target(2))/block.Target(2);
Current = (Y(3)+ block.PropGain(3))/(1+block.Scale(3)*Inlet.Voltage*block.Cells/(1000*NetPower)*block.PropGain(3))*block.Scale(3); 
Power = Current*Inlet.Voltage*block.Cells/1000;
PowerError = (NetPower - Power)/NetPower;
FuelFlow = block.Cells*Current/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000;

Bypass = min(1,max(0,Y(1)+TavgError*block.PropGain(1)));
if block.AnodeRecirc.IC == 0
    Recirculation = 0;
else
    Steam2Carbon = block.Target(3);
    CH4 = block.Fuel.CH4*FuelFlow;
    r = block.AnodeRecirc.IC; %Initial guess of anode recirculation
    dr = 1e-6;
    error = 1;
    % Inlet = (Inlet + generated - consumed)*r  + New, thus inlet = New/(1-r) + (generated - consumed)*r/(1-r)
    while abs(error)>1e-6
        COin = block.Fuel.CO*FuelFlow/(1-r) + (block.Fuel.CH4 - block.WGSeffective*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r/(1-r);
        Inlet.FuelMix.H2O = block.Fuel.H2O*FuelFlow/(1-r) + (block.Cells*Current /(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*block.WGSeffective)*FuelFlow)*r/(1-r);
        S2C = Inlet.FuelMix.H2O/(CH4 + 0.5*COin);
        error = Steam2Carbon - S2C;
        r2 = r+dr;
        COin2 = block.Fuel.CO*FuelFlow/(1-r2) + (block.Fuel.CH4 - block.WGSeffective*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r2/(1-r2);
        H2Oin2 = block.Fuel.H2O*FuelFlow/(1-r2) + (block.Cells*Current/(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*block.WGSeffective)*FuelFlow)*r2/(1-r2);
        S2C2 = H2Oin2/(CH4 + 0.5*COin2);
        dSdr = (S2C2 - S2C)/dr;
        r = r + error/dSdr;
    end
    Recirculation = r;
end
    
if strcmp(string1,'Outlet')
    Out.HeaterBypass = Bypass;
    Out.Blower = (Y(2)+dTerror*block.PropGain(2))*block.Scale(2);
    Out.Current = Current;
    Out.AnodeRecirc = Recirculation;
    Out.FuelFlow = FuelFlow;
elseif strcmp(string1,'dY')
    dY = Y*0;
    if (Bypass<=0 && TavgError<0) || (Bypass>=1 && TavgError>0)
        dY(1) = 0; %saturation, anti-windup
    else
        dY(1) = block.Gain(1)*TavgError;
    end
    %% need to deterimine saturation on blower
    dY(2) = block.Gain(2)*dTerror*block.Scale(2);
    dY(3) = block.Gain(3)*PowerError*block.Scale(3);
    Out = dY;
    Tags.(block.name).Bypass = Bypass;
    Tags.(block.name).Blower = (Y(2)+dTerror*block.PropGain(2))*block.Scale(2);
    Tags.(block.name).Recirculation = Recirculation;
    Tags.(block.name).FuelFlow = FuelFlow;
    Tags.(block.name).Current = Current;
end