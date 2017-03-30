function block = InitializeControlSOFCsystem(varargin)
% Controls for SOFC system with external reformer, control heat exchanger bypass, blower power, anode recirculation and fuel flow rate
% Four (4) inlets: T hot, T cold, T average, Voltage
% Five (5) outlets:  Heater bypass, blower power, fuel flow rate, current, anode recirculation
% Three (3) states: Heater bypass, blower power, Current 
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
    Current = sum(ComponentProperty(block.InitConditions{3}));
    
    
    FuelFlow = (block.Cells*Current/(2*F*block.Utilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))/1000);
    if Recirculation > 0
        %use reciculation from 1st initialization to set 
        Steam2Carbon = block.Target(3);
        if Steam2Carbon ==2 %steam to carbon ratio is unaffected by WGS effectiveness
            block.WGSeffective = 0.7;
        else
            eWGS = .7; %Initial guess of WGS efftiveness
            CH4 = block.Fuel.CH4*FuelFlow;
            r = Recirculation; 
            dWGS = 1e-6;
            error = 1;
            % Inlet = (Inlet + generated - consumed)*r  + New, thus inlet = New/(1-r) + (generated - consumed)*r/(1-r)
            while abs(error)>1e-6
                COin = block.Fuel.CO*FuelFlow/(1-r) + (block.Fuel.CH4 - eWGS*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r/(1-r);
                Inlet.FuelMix.H2O = block.Fuel.H2O*FuelFlow/(1-r) + (block.Cells*Current/(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*eWGS)*FuelFlow)*r/(1-r);
                S2C = Inlet.FuelMix.H2O/(CH4 + 0.5*COin);
                error = Steam2Carbon - S2C;
                eWGS2 = eWGS+dWGS;
                COin2 = block.Fuel.CO*FuelFlow/(1-r) + (block.Fuel.CH4 - eWGS2*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r/(1-r);
                H2Oin2 = block.Fuel.H2O*FuelFlow/(1-r) + (block.Cells*Current/(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*eWGS2)*FuelFlow)*r/(1-r);
                S2C2 = H2Oin2/(CH4 + 0.5*COin2);
                dSdWGS = (S2C2 - S2C)/dWGS;
                if dSdWGS ==0
                    error = 0;
                    eWGS = 1;
                else
                    eWGS = max(0,eWGS + error/dSdWGS);
                end
            end
            block.WGSeffective = eWGS;
        end
    end
    
    block.PortNames = {'Hot','Cold','Voltage','HeaterBypass','Blower','Current','AnodeRecirc','FuelFlow'};
    block.Hot.type = 'in';
    block.Hot.IC = block.Target(2)+.5*block.Target(1); 
    block.Cold.type = 'in';
    block.Cold.IC = block.Target(2)-.5*block.Target(1); 
    block.Voltage.type = 'in';
    block.Voltage.IC = 0.85; 
    block.HeaterBypass.type = 'out';
    block.HeaterBypass.IC = Bypass;
    block.Blower.type = 'out';
    block.Blower.IC = Blower; 
    block.Current.type = 'out';
    block.Current.IC =  Current;
    block.AnodeRecirc.type = 'out';
    block.AnodeRecirc.IC = Recirculation;
    block.FuelFlow.type = 'out';
    block.FuelFlow.IC = FuelFlow;
    
    block.P_Difference = {};
    
    block.IC = [Bypass;1;]; % inital condition
    block.Scale = [1;Blower;];
%     block.IC = [Bypass;1;1;]; % inital condition
%     block.Scale = [1;Blower;Current];
    
    for i = 1:1:length(block.PortNames)
        if length(block.connections)<i || isempty(block.connections{i})
            block.(block.PortNames{i}).connected={};
        else
            if ischar(block.connections{i})
                block.(block.PortNames{i}).connected = block.connections(i);
            else
                block.(block.PortNames{i}).IC = block.connections{i};
                block.(block.PortNames{i}).connected={};
            end
        end
    end
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
        r = block.AnodeRecirc.IC; %Initial guess of anode recirculation
        CH4 = block.Fuel.CH4*FuelFlow;
        dr = 1e-6;
        error = 1;
        % Inlet = (Inlet + generated - consumed)*r  + New, thus inlet = New/(1-r) + (generated - consumed)*r/(1-r)
        while abs(error)>1e-6
            COin = block.Fuel.CO*FuelFlow/(1-r) + (block.Fuel.CH4 - block.WGSeffective*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r/(1-r);
            Inlet.FuelMix.H2O = block.Fuel.H2O*FuelFlow/(1-r) + (block.Cells*block.Current.IC /(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*block.WGSeffective)*FuelFlow)*r/(1-r);
            S2C = Inlet.FuelMix.H2O/(CH4 + 0.5*COin);
            error = Steam2Carbon - S2C;
            r2 = r+dr;
            COin2 = block.Fuel.CO*FuelFlow/(1-r2) + (block.Fuel.CH4 - block.WGSeffective*(block.Fuel.CH4+block.Fuel.CO))*FuelFlow*r2/(1-r2);
            H2Oin2 = block.Fuel.H2O*FuelFlow/(1-r2) + (block.Cells*block.Current.IC /(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*block.WGSeffective)*FuelFlow)*r2/(1-r2);
            S2C2 = H2Oin2/(CH4 + 0.5*COin2);
            dSdr = (S2C2 - S2C)/dr;
            r = r + error/dSdr;
        end
        block.AnodeRecirc.IC = r;
    end
    block.InitializeError = max(abs(TinletError),abs(dTerror)); 
%     block.InitializeError = abs(TavgError);
    block.Scale = [1; block.Blower.IC;];
    block.IC = [block.HeaterBypass.IC-TinletError*block.PropGain(1); 1-dTerror*block.PropGain(2);];
%     block.Scale = [1; block.Blower.IC; block.Current.IC];
%     block.IC = [block.HeaterBypass.IC-TinletError*block.PropGain(1); 1-dTerror*block.PropGain(2);1-PowerError*block.PropGain(3);];
end