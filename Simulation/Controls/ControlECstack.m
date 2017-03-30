function Out = ControlECstack(t,Y, Inlet,block,string1)
% Controls for electrolyzer stack only, control air flow to anode and inlet temperature and steam flow rate
% Four (4) inlets: T hot, T cold, T average, Voltage
% Five (5) outlets: OxidentTemp, oxidant flow rate, steam temp, steam flow rate, Current
% Two  (2) states: Oxidant flow rate, net current
global Tags F
averageT = mean(Inlet.PEN_Temp);
% averageT = (Inlet.Hot - Inlet.Cold + block.dT_cath_PEN);
[h,~] = enthalpy(block.Target(1),{'H2','H2O','O2'});
h_rxn3 = h.H2+.5*h.O2-h.H2O;
Vbalance = 1./(2*F)*h_rxn3; %voltage that balances heat

if block.HasFlow
    NetPower = PowerDemandLookup(t);
    Current = -NetPower*1000/(Inlet.Voltage*block.Cells);
else
    Current = Y(1);
    VoltError = Vbalance - Inlet.Voltage;
end
SteamFlow = block.Cells*abs(Current)/(2*F*block.Utilization*block.Steam.H2O)/1000;

if block.HasFlow
    Q_cathode = Y(1)*40*block.Target(2);
        
    if ((block.Cells*(Inlet.Voltage - Vbalance)*abs(Current)/1000) - Q_cathode)>0
        Out.OxidantTemp = block.Target(1)-100; %cooling stack
        TavgError = (averageT -block.Target(1))/block.Target(2); % too hot = increase flow
    else
        Out.OxidantTemp = block.Target(1)+100;%heating stack
        TavgError = (block.Target(1)-averageT)/block.Target(2); %too hot = reduce flow
    end
    Out.OxidantFlow = max(0,Y(1)+(TavgError*block.PropGain(1))*block.Scale(1));
else
    Out.OxidantTemp = block.Target(1);
    Out.OxidantFlow = 0;
end

if strcmp(string1,'Outlet')
    Out.SteamTemp = block.SteamTemp.IC;%currently no control of fuel inlet temp
    Out.SteamFlow = SteamFlow;
    Out.Current = Current;
    Tags.(block.name).OxidantTemp = Out.OxidantTemp;
    Tags.(block.name).OxidantFlow = Out.OxidantFlow;
    Tags.(block.name).Tsteam = block.SteamTemp.IC;%currently no control of steam inlet temp
    Tags.(block.name).SteamFlow = SteamFlow;
    Tags.(block.name).Current = Current;
elseif strcmp(string1,'dY')
    dY = Y*0;
    if block.HasFlow
        if Y(1)>0 || TavgError>0 %anti-windup so flow does not go negative
            dY(1) = block.Gain(1)*TavgError;
        end
    else
        dY(1) = block.Gain(2)*VoltError;
    end
    Out = dY;
end