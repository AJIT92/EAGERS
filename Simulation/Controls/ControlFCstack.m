function Out = ControlFCstack(t,Y, Inlet,block,string1)
% Controls for Fuel cell stack only, control air flow and inlet temperature and fuel flow rate, net current and anode recirculation
% Four (4) inlets: T hot, T cold, T average, Voltage
% Six (6) outlets: OxidentTemp, oxidant flow rate, fuel temp, fuel flow rate, current, anode recirculation
% Three (3) states: OxidentTemp, oxidant flow rate, Current 
% if OxyFC: % Three (3) states: anode recirculation, Fuel flow rate, net current
global Tags F
NetPower = PowerDemandLookup(t);
averageT = mean(Inlet.PEN_Temp);
TavgError = (block.Target(1)-averageT)/block.Target(2);
deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
dTerror =(deltaT-block.Target(2))/block.Target(2);
    
if isfield(block,'OxyFC')
    Current = (Y(3)+ TavgError*block.PropGain(3))*block.Scale(3); % current controller for managing temperature in closed-end cathode FC (proportional control only)
    % O2flow = block.Cells*Current/(4*F*1000)*32*3600*24/1000; %Ton/day
    Parasitic = 0;%(1.0101*O2flow^-.202)*(O2flow*1000/24); %parasitic in kW
    Power = Current*Inlet.Voltage*block.Cells/1000;
    PowerError = (NetPower - (Power - Parasitic))/NetPower;
    
    FuelFlow = (Y(2)+PowerError*block.PropGain(2))*block.Scale(2);
    CH4_Util = Y(3)*block.Scale(3)*block.Cells/(2000*F)/(FuelFlow*4*block.Fuel.CH4);%hydrogen used vs. ideal H2 available
    R.CH4 = Current/(8000*F*CH4_Util);
    a = 4352.2./block.Target(1) - 3.99;
    R.WGS = R.CH4*exp(a)*block.WGSeffective;% Water gas shift equilibrium constant
    R.H2 = Current/(2000*F); %# of electrochemical reactions
    [h,~] = enthalpy(averageT,{'H2','H2O','O2','CO','CO2','CH4'});
    h_rxn1 = h.CO+3*h.H2-h.CH4-h.H2O;
    h_rxn2 = h.CO2+h.H2-h.CO-h.H2O;
    h_rxn3 = h.H2O - h.H2 -.5*h.O2;
    
    FuelIn = block.Fuel;
    FuelIn.T = mean(Inlet.Cold);
    Cp = SpecHeat(FuelIn);

    BalancedPower = R.H2*(-h_rxn3) - (R.CH4*h_rxn1+R.WGS*h_rxn2) - FuelFlow/block.Cells*Cp*block.Target(2);
    BalancedVoltage = BalancedPower*1000/Current;
    VoltageError = Inlet.Voltage - BalancedVoltage;
else
    Current = (Y(3)+ block.PropGain(3))/(1+block.Scale(3)*Inlet.Voltage*block.Cells/(1000*NetPower)*block.PropGain(3))*block.Scale(3); 
    Power = Current*Inlet.Voltage*block.Cells/1000;
    PowerError = (NetPower - Power)/NetPower;
    FuelFlow = block.Cells*Current/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000;
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
end
if strcmp(string1,'Outlet')
    if isfield(block,'OxyFC')
        Out.OxidantTemp = block.Target(1);
        Out.OxidantFlow = block.Cells*Current/(4000*F*block.Oxidant.O2); % kmol/s
        Out.AnodeRecirc = Y(1)*block.Scale(1)+ dTerror*block.PropGain(1)*block.Target(2);
    else
        Out.OxidantTemp = (Y(1)+TavgError*block.PropGain(1))*block.Scale(1);
        Out.OxidantFlow = (Y(2)+dTerror*block.PropGain(2))*block.Scale(2);
        Out.AnodeRecirc =  Recirculation;
    end
    Out.FuelFlow = FuelFlow;
    Out.Current = Current;
    
    Tags.(block.name).FuelFlow = Out.FuelFlow;
    Tags.(block.name).Current = Out.Current;
    Tags.(block.name).Recirculation = Out.AnodeRecirc;
    Tags.(block.name).OxidantTemp = Out.OxidantTemp;
    Tags.(block.name).OxidantFlow = Out.OxidantFlow;
    Tags.(block.name).Power = Power;
elseif strcmp(string1,'dY')
    dY = Y*0;
    if isfield(block,'OxyFC')
        dY(1) = block.Gain(1)*dTerror; %recirculation changes temperature gradient
        dY(2) = block.Gain(2)*PowerError;%fuel flow
        dY(3) = block.Gain(3)*VoltageError;%current
    else    
        dY(1) = block.Gain(1)*TavgError;
        dY(2) = block.Gain(2)*dTerror;
        dY(3) = block.Gain(3)*PowerError;
    end
    Out = dY;
end