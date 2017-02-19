function block = InitializeControlMCFCstack(varargin)
% Controls for MCFC stack + oxidizer fed with some incoming air, and oxidizer output going into cathode, control air flow and inlet temperature and fuel flow rate
% Five (5) inlets: T hot, T cold, T average, Current, Equivelance ratio
% Five (5) outlets: OxidentTemp, oxidant flow rate, fuel temp, fuel flow rate, Valve position
% Three (3) states: OxidentTemp, oxidant flow rate, valve postion 
global modelParam F
block = varargin{1};
if length(varargin)==1 %first initialization
    block.States=[];
    block.TagInf = {}; %Tags to record at each step
    block.TagFinal = {}; %Tags to record at the final step
    block.IC = [1 1 1]; % inital condition
    block.Scale = block.IC;
    if ischar(block.Cells)
        block.Cells = eval(strcat('modelParam.',block.Cells,';'));
    end
    if ischar(block.Fuel)
        block.Fuel = eval(strcat('modelParam.',block.Fuel,';'));
    end
    if ischar(block.Utilization)
        block.Utilization = eval(strcat('modelParam.',block.Utilization,';'));
    end
    if ischar(block.DesignTemp)
        block.Scale(1) = eval(strcat('modelParam.',block.DesignTemp,';'));
    else block.Scale(1) = block.DesignTemp;
    end
    if ischar(block.DesignFlow)
        block.Scale(2) = eval(strcat('modelParam.',block.DesignFlow,';'));
    else block.Scale(2) = block.DesignFlow;
    end
    block.Scale(3) = block.valveIC;%Initial condition normalized to 1
    Target = zeros(length(block.Target),1);
    for j = 1:1:length(Target)
        if ischar(block.Target{j})
            Target(j) = eval(strcat('modelParam.',block.Target{j},';'));
        else Target(j) = block.Target{j};
        end
    end
    block.Target  = Target;
    block.StackPower = modelParam.('FC1').StackPower;
    
    block.PortNames = {'Hot','Cold','PEN_Temp','Current','Equiv','OxidantTemp','OxidantFlow','FuelTemp','FuelFlow','ValvePerc','PowerStack'};
    block.Hot.type = 'in';
    block.Hot.IC = block.Target(2)+.5*block.Target(1); 
    block.Cold.type = 'in';
    block.Cold.IC = block.Target(2)-.5*block.Target(1); 
    block.PEN_Temp.type = 'in';
    block.PEN_Temp.IC = block.Target(2);
    block.Current.type = 'in';
    block.Current.IC = modelParam.('FC1').Current; %needed for flow outlet
    block.Equiv.type = 'in';
    block.Equiv.IC = 0.8;
    block.OxidantTemp.type = 'out';
    block.OxidantTemp.IC = block.Scale(1);
    block.OxidantFlow.type = 'out';
    block.OxidantFlow.IC = block.Scale(2); 
    block.FuelTemp.type = 'out';
    block.FuelTemp.IC = block.Scale(1);
    block.FuelFlow.type = 'out';
    block.FuelFlow.IC = block.Cells*block.Current.IC/(2*F*block.Utilization*(4*block.Fuel(1)+block.Fuel(2)+block.Fuel(4)))/1000;
    block.ValvePerc.type = 'out';
    block.ValvePerc.IC = block.Scale(3);
    block.PowerStack.type = 'out';
    block.PowerStack.IC =  modelParam.('FC1').StackPower;
    
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

% if length(varargin)==2 %% Have inlets connected, re-initialize
%     Inlet = varargin{2};
%     
% end