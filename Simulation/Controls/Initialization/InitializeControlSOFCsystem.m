function block = InitializeControlSOFCsystem(varargin)
% Controls for SOFC system with external reformer, control heat exchanger bypass, blower power, anode recirculation and fuel flow rate
% Four (4) inlets: T hot, T cold, T average, Voltage
% Five (5) outlets:  Heater bypass, blower power, fuel flow rate, current, anode recirculation
% Two (2) states: Heater bypass, blower power
%Need to add state for current back in to avoid fuel starvation during step changes
global F
block = varargin{1};
if length(varargin)==1 %first initialization
    block.description = {'Heater Bypass';'Blower Power';'Fuel Cell Current';};
    
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(Target)
        Target(j) = ComponentProperty(block.Target{j});
    end
    block.Target  = Target;
    
    block.Cells = ComponentProperty(block.Cells);
    block.Fuel = ComponentProperty(block.Fuel);
    block.Utilization = ComponentProperty(block.Utilization);
    Recirculation = ComponentProperty(block.InitialAnodeRecirc);

    Bypass = ComponentProperty(block.InitConditions{1});
    Blower = ComponentProperty(block.InitConditions{2});
    BlowerMass = ComponentProperty('Blower.FlowDesign');
    CathodeMass = NetFlow(ComponentProperty('FC1.Flow2Out.IC'))*28.84;
    ComponentProperty('Blower.FlowDesign',CathodeMass);
    Blower = Blower*CathodeMass/BlowerMass; %scale blower power
    a = ComponentProperty(block.InitConditions{3});
    Current = sum(a.H2 + a.CO);
    
    
    FuelFlow = (block.Cells*Current/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
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
    block.Voltage.IC = 0.85; 
    
    block.OutletPorts = {'HeaterBypass','Blower','Current','AnodeRecirc','FuelFlow'};
    block.HeaterBypass.IC = Bypass;
    block.Blower.IC = Blower; 
    block.Current.IC =  Current;
    block.AnodeRecirc.IC = Recirculation;
    block.FuelFlow.IC = FuelFlow;
    
    block.P_Difference = {};
    
    block.IC = [Bypass;1;]; % inital condition
    block.Scale = [1;Blower;];
%     block.IC = [Bypass;1;1;]; % inital condition
%     block.Scale = [1;Blower;Current];
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    PEN_Temperature = mean(ComponentProperty('FC1.T.Elec'));
    blowerPower = ComponentProperty('Blower.NominalPower');
    StackPower = PowerDemandLookup(0) + blowerPower;
%     PowerError = (StackPower-block.Current.IC*Inlet.Voltage*block.Cells/1000)/StackPower;
%     block.Current.IC = block.Current.IC*(1 + PowerError);
    
    block.Current.IC = StackPower*1000/(Inlet.Voltage*block.Cells);
    averageT = (mean(Inlet.Hot)+mean(Inlet.Cold))/2; %average temperature of cathode
    block.dT_cath_PEN = PEN_Temperature - averageT; %temperature differencce between cathode and PEN
    TinletError = (block.Target(1) - (mean(Inlet.Cold) + .5*block.Target(2) + block.dT_cath_PEN))/block.Target(2);
    
    deltaT = (mean(Inlet.Hot)-mean(Inlet.Cold));
    dTerror =(deltaT-(block.Target(2)))/block.Target(2);

    Steam2Carbon = block.Target(3);
    
    %% adjust blower mass flow to get deltaT correct
    blowerMassFlow = ComponentProperty('Blower.FlowDesign');
    ComponentProperty('Blower.FlowDesign',blowerMassFlow*(1+.5*dTerror));
    block.Blower.IC = blowerPower*(1+.5*dTerror);
%     blower = block.Blower.IC
    %%if too cold, add more fuel, it will burn and more heat will be recovered
%     Utilization = block.Utilization*(1-0.04*TavgError)
%     block.Utilization = Utilization;
    
%     %% adjust heat exchanger effectiveness to control inlet temperature
%     HXtarget = ComponentProperty('HX1.Effectiveness');
%     HXtarget = min(0.96,HXtarget*(1 + .05*TavgError));
%     ComponentProperty('HX1.Target',HXtarget);
%     ComponentProperty('HX1.sizemethod','Effectiveness');
    
    %% adjust heat exchanger cold exit temperature to control FC inlet temperature
    HXtarget = ComponentProperty('HX1.Target');
    HXtarget = HXtarget + TinletError*block.Target(2);
    ComponentProperty('HX1.Target',HXtarget);

    %%adjust exit temperature to avoid energy feedback during iterations
    Terror = (block.Target(1) - (averageT+block.dT_cath_PEN));
    CathodeOut = ComponentProperty('FC1.Flow2Out.IC');
    CathodeOut.T = CathodeOut.T + Terror;
    ComponentProperty('FC1.Flow2Out.IC',CathodeOut);
    AnodeOut = ComponentProperty('FC1.Flow1Out.IC');
    AnodeOut.T = AnodeOut.T + Terror;
    ComponentProperty('FC1.Flow1Out.IC',AnodeOut);

    
    FuelFlow = (block.Cells*block.Current.IC/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
    block.FuelFlow.IC = FuelFlow;
    
    if block.AnodeRecirc.IC > 0
        block.AnodeRecirc.IC = anodeRecircHumidification(block.Fuel,FuelFlow,block.WGSeffective,Steam2Carbon,block.Cells*block.Current.IC/(2*F*1000),block.AnodeRecirc.IC);
    end
    block.InitializeError = max(abs(TinletError),abs(dTerror)); 
%     block.InitializeError = abs(TavgError);
    block.Scale = [1; block.Blower.IC;];
    block.IC = [block.HeaterBypass.IC-TinletError*block.PropGain(1); 1-dTerror*block.PropGain(2);];
%     block.Scale = [1; block.Blower.IC; block.Current.IC];
%     block.IC = [block.HeaterBypass.IC-TinletError*block.PropGain(1); 1-dTerror*block.PropGain(2);1-PowerError*block.PropGain(3);];
end