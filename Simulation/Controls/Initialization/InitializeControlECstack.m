function block = InitializeControlECstack(varargin)
% Controls for electrolyzer stack only, control air flow and inlet temperature and fuel flow rate
% Four (4) inlets: T hot, T cold, T average, Voltage
% Five (5) outlets: OxidentTemp, oxidant flow rate, steam temp, steam flow rate, Current
% Two  (2) states: Oxidant flow rate, net current
global F
block = varargin{1};
if length(varargin)==1 %first initialization
    block.IC = [1 1]; % inital condition: OxidentTemp, oxidant flow rate, net current
    block.Scale = block.IC;
    block.description = {'Oxidant Flow Rate';'Fuel Cell Current';};
    
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(Target)
        Target(j) = lookupVal(block.Target{j});
    end
    block.Target  = Target;
    
    block.Cells = lookupVal(block.Cells);
    block.Steam = lookupVal(block.Steam);
    block.Utilization = lookupVal(block.Utilization);
    block.SteamTemperature = lookupVal(block.SteamTemperature);
    
    for j = 1:1:length(block.Scale)
        block.Scale(j) = lookupVal(block.InitConditions{j});
    end
    Current = block.Scale(2);
    
    SteamFlow = (block.Cells*abs(Current)/(2*F*block.Utilization*block.Steam.H2O)/1000);
    
    if block.Scale(1)>0
        block.HasFlow = true;
    else
        block.HasFlow = false;
    end

    block.PortNames = {'Hot','Cold','PEN_Temp','Voltage','OxidantTemp','OxidantFlow','SteamTemp','SteamFlow','Current'};
    block.Hot.type = 'in';
    block.Hot.IC = block.Target(1)+.5*block.Target(2); 
    block.Cold.type = 'in';
    block.Cold.IC = block.Target(1)-.5*block.Target(2); 
    block.PEN_Temp.type = 'in';
    block.PEN_Temp.IC = block.Target(1);
    block.Voltage.type = 'in';
    block.Voltage.IC = 1.3; 
    block.OxidantTemp.type = 'out';
    block.OxidantTemp.IC = block.Target(1);
    block.OxidantFlow.type = 'out';
    block.OxidantFlow.IC = block.Scale(1); 
    block.SteamTemp.type = 'out';
    block.SteamTemp.IC = block.SteamTemperature;
    block.SteamFlow.type = 'out';
    block.SteamFlow.IC = SteamFlow;
    block.Current.type = 'out';
    block.Current.IC = Current;
    
    block.P_Difference = {};
    
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
        block.Scale = [block.OxidantFlow.IC, block.Current.IC];
        block.IC = [1-TavgError*block.PropGain(1),1-PowerError*block.PropGain(2)]; %Oxidant flow rate, net current
    else
        block.InitializeError = 0;
        block.Scale = [block.SteamFlow.IC, block.Current.IC];
        block.IC = [0,1-PowerError*block.PropGain(2)]; % inital condition with no anode flow in: OxidentTemp, oxidant flow rate, net current
    end
end

function const = lookupVal(initCond)
global modelParam
if ischar(initCond)
    r = strfind(initCond,'.');
    if ~isempty(r)
        r = [r,length(initCond)+1];
        A = modelParam.(initCond(1:r(1)-1));
        for i = 2:1:length(r)
            field = initCond(r(i-1)+1:r(i)-1);
            A = A.(field);
        end
    else
        A = modelParam.(initCond);
    end
    if isstruct(A) && isfield(A,'T')
        const = NetFlow(A);
    elseif length(A)>1
        const  = sum(A);
    else
        const = A;
    end
else const  = initCond;
end