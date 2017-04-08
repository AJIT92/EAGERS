function block = InitializeControlFCstack(varargin)
% Controls for Fuel cell stack only, control air flow and inlet temperature and fuel flow rate, net current and anode recirculation
% Three (3) inlets: T hot, T cold, Voltage
% Six (6) outlets: OxidentTemp, oxidant flow rate, fuel temp, fuel flow rate, current, anode recirculation
% One (1) state: , oxidant flow rate %% %Need to add state for current back in to avoid fuel starvation during step changes
% Current, Fuel Flow, oxidant temp, and Anode Recirculation are calculated directly (feed-forward)
% if OxyFC: % Three (3) states: anode recirculation, Fuel flow rate, Current
global F
block = varargin{1};
if length(varargin)==1 %first initialization
    if isfield(block,'OxyFC')
        block.description = {'Anode Recirculation';'Fuel Flow';'Fuel Cell Current';};
    else
        block.description = {'Oxidant Inlet Temperature';'Oxidant Flow Rate';'Fuel Cell Current';};
    end
    
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(block.Target)
        Target(j) = ComponentProperty(block.Target{j});
    end
    block.Target = Target;
    
    block.Cells = ComponentProperty(block.Cells);
    block.Fuel = ComponentProperty(block.Fuel);

    a = ComponentProperty(block.InitConditions{3});
    Current = sum(a.H2 + a.CO);%net current
    if isfield(block,'OxyFC')
        Recirculation = ComponentProperty(block.InitConditions{1});
        FuelFlow = NetFlow(ComponentProperty(block.InitConditions{2}));
        block.Oxidant = ComponentProperty(block.Oxidant);
        OxidantTemp = block.Target(1);
        block.OxidantUtilization = ComponentProperty(block.OxidantUtilization);
        OxidantFlow = Current*block.Cells/(4000*F*block.Oxidant.O2*block.OxidantUtilization);
    else
        Recirculation = ComponentProperty(block.InitialAnodeRecirc);
        block.Utilization = ComponentProperty(block.Utilization);
        FuelFlow = (block.Cells*Current/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
        OxidantTemp = ComponentProperty(block.InitConditions{1});
        OxidantFlow = NetFlow(ComponentProperty(block.InitConditions{2}));
    end
    
    if Recirculation > 0
        %use reciculation from 1st initialization to set 
        Steam2Carbon = block.Target(3);
        if Steam2Carbon ==2 %steam to carbon ratio is unaffected by WGS effectiveness
            block.WGSeffective = 0.7;
        else
            block.WGSeffective = effectiveWGS(block.Fuel,FuelFlow,0.7,Steam2Carbon,block.Cells*Current/(2*F*1000),Recirculation);
        end
    end
        
    block.InletPorts = {'Hot','Cold','Voltage'};
    block.Hot.IC = block.Target(2)+.5*block.Target(1); 
    block.Cold.IC = block.Target(2)-.5*block.Target(1); 
    block.Voltage.IC = 0.85; %needed for flow outlet
    
    block.OutletPorts = {'OxidantTemp','OxidantFlow','AnodeRecirc','FuelFlow','Current'};
    block.OxidantTemp.IC = OxidantTemp;
    block.OxidantFlow.IC = OxidantFlow; 
    block.AnodeRecirc.IC = Recirculation;
    block.FuelFlow.IC = FuelFlow;
    block.Current.IC = Current; %needed for flow outlet
    

    if isfield(block,'OxyFC')
        block.IC = [Recirculation;1;1;]; % inital condition
        block.Scale = [1;FuelFlow;Current];
    else
        block.IC = [1 1]; % inital condition
        block.Scale =  [OxidantTemp; OxidantFlow;];
    end
    
    block.P_Difference = {};
end

if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};

    NetPower = PowerDemandLookup(0);
    Power = block.Current.IC*Inlet.Voltage*block.Cells/1000;
%     PEN_Temperature = mean(ComponentProperty('FC1.T.Elec'));
    averageT = (mean(Inlet.Hot)+mean(Inlet.Cold))/2; %average temperature of cathode
%     block.dT_cath_PEN = PEN_Temperature - averageT; %temperature differencce between cathode and PEN    
    TavgError = (block.Target(1)-averageT)/block.Target(2);
    deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
    dTerror =(deltaT-block.Target(2))/block.Target(2);
     
    Steam2Carbon = block.Target(3);
%     dTerror
%     TavgError
    if isfield(block,'OxyFC')
        FuelFlow = block.FuelFlow.IC*(1-.02*TavgError);
        O2flow = block.Cells*block.Current.IC/(4*F*1000)*32*3600*24/1000; %Ton/day
        Parasitic = (1.0101*O2flow^-.202)*(O2flow*1000/24); %parasitic in kW
        PowerError = (NetPower + Parasitic - Power)/NetPower;
        block.Current.IC = block.Current.IC*(1 + PowerError);
%         block.Current.IC = (NetPower + Parasitic)*1000/(Inlet.Voltage*block.Cells);
    else
        block.Current.IC = NetPower*1000/(Inlet.Voltage*block.Cells);
        FuelFlow = (block.Cells*block.Current.IC/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
        a = .5;
        block.OxidantFlow.IC = block.OxidantFlow.IC*(1+a*dTerror); 
        block.OxidantTemp.IC = block.OxidantTemp.IC + a*(TavgError + .75*dTerror)*block.Target(2);
    end
    
    block.FuelFlow.IC = FuelFlow;
    if block.AnodeRecirc.IC > 0
        block.AnodeRecirc.IC = anodeRecircHumidification(block.Fuel,FuelFlow,block.WGSeffective,Steam2Carbon,block.Cells*block.Current.IC/(2*F*1000),block.AnodeRecirc.IC);
    end
    if isfield(block,'OxyFC')
        block.InitializeError = 0;%max(abs(PowerError));
%         block.Scale = [1, block.FuelFlow.IC];
%         block.IC = [block.AnodeRecirc.IC-dTerror*block.PropGain(1), 1-PowerError*block.PropGain(2)];%-VoltageError*block.PropGain(3)];
        block.Scale = [1, block.FuelFlow.IC, block.Current.IC];
        block.IC = [block.AnodeRecirc.IC-dTerror*block.PropGain(1), 1-PowerError*block.PropGain(2),1];%-VoltageError*block.PropGain(3)];
    else
        block.InitializeError = max(abs(TavgError),abs(dTerror));
        block.Scale = [block.OxidantFlow.IC];
        block.IC = [1-dTerror*block.PropGain(1)];
    end
end