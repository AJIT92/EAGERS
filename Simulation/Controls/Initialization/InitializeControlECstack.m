function block = InitializeControlECstack(varargin)
% Controls for electrolyzer stack only, control air flow and inlet temperature and fuel flow rate
% Four (4) inlets: T hot, T cold, T average, Voltage
% Five (5) outlets: OxidentTemp, oxidant flow rate, steam temp, steam flow rate, Current
% Two  (2) states: Oxidant flow rate, net current
global F
block = varargin{1};
if length(varargin)==1 %first initialization
    block.description = {'Oxidant Flow Rate';'Fuel Cell Current';};
    
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(Target)
        Target(j) = ComponentProperty(block.Target{j});
    end
    block.Target  = Target;
    
    block.Cells = ComponentProperty(block.Cells);
    block.Steam = ComponentProperty(block.Steam);
    block.Utilization = ComponentProperty(block.Utilization);
    block.SteamTemperature = ComponentProperty(block.SteamTemperature);
    
    OxFlow = NetFlow(ComponentProperty(block.InitConditions{1})); % inital condition: oxidant flow rate, 
    a = ComponentProperty(block.InitConditions{2});
    Current = sum(a.H2 + a.CO);% net current
    
    SteamFlow = (block.Cells*abs(Current)/(2*F*block.Utilization*block.Steam.H2O)/1000);
    
    if OxFlow>0
        block.HasFlow = true;
    else
        block.HasFlow = false;
    end
    
    block.InletPorts = {'Hot','Cold','Voltage','PEN_Temp'};
    block.Hot.IC = block.Target(1)+.5*block.Target(2); 
    block.Cold.IC = block.Target(1)-.5*block.Target(2);    
    block.Voltage.IC = 1.3; 
    block.PEN_Temp.IC = block.Target(1);
    
    block.OutletPorts = {'OxidantTemp','OxidantFlow','SteamTemp','SteamFlow','Current'};
    block.OxidantTemp.IC = block.Target(1);
    block.OxidantFlow.IC = OxFlow; 
    block.SteamTemp.IC = block.SteamTemperature;
    block.SteamFlow.IC = SteamFlow;
    block.Current.IC = Current;
    
    block.P_Difference = {};
    
    
    block.IC = 1; % inital condition
    if block.HasFlow
        block.Scale = [OxFlow;];
    else 
        block.Scale = [Current;];
    end
end

if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    NetPower = PowerDemandLookup(0);
    PowerError = (NetPower-abs(block.Current.IC)*Inlet.Voltage*block.Cells/1000)/NetPower;
    block.Current.IC = block.Current.IC*(1 + PowerError);
    block.SteamFlow.IC = (block.Cells*abs(block.Current.IC)/(2*F*block.Utilization*block.Steam.H2O)/1000);
    block.SteamTemp.IC = block.SteamTemp.IC; %currently not manipulated during initialization
    
    if block.HasFlow
        averageT = mean(Inlet.PEN_Temp);
        [h,~] = enthalpy(block.Target(1),{'H2','H2O','O2'});
        h_rxn3 = h.H2+.5*h.O2-h.H2O;
        Vbalance = 1./(2*F)*h_rxn3; %voltage that balances heat
        Q_cathode = block.OxidantFlow.IC*40*block.Target(2);
        
        if ((block.Cells*(Inlet.Voltage - Vbalance)*abs(block.Current.IC)/1000) - Q_cathode)>0
            block.OxidantTemp.IC = block.Target(1)-100; %cooling stack
            TavgError = (averageT -block.Target(1))/block.Target(2); % too hot = increase flow
        else
            block.OxidantTemp.IC = block.Target(1)+100;%heating stack
            TavgError = (block.Target(1)-averageT)/block.Target(2); %too hot = reduce flow
        end
        block.InitializeError = abs(TavgError); 

        newFlow = block.OxidantFlow.IC*(1+.5*TavgError);
        block.OxidantFlow.IC = newFlow;
        block.Scale = [block.OxidantFlow.IC];
        block.IC = [1-TavgError*block.PropGain(1)]; %Oxidant flow rate
    else
        block.InitializeError = 0;
        block.Scale = [block.Current.IC];
        block.IC = [1-PowerError*block.PropGain(2)]; % inital condition with no anode flow in: net current
    end
end