function block = InitializeLeakageValve(varargin)
%2 inlets: Inlet flow and valve position
%3 outlets: Leakage and Flow1 and flow 2
block = varargin{1};
if length(varargin) == 1%first initialization
    block.IC = [];
    block.Scale = [];
    
    block.PortNames = {'ValvePos','Intake','Flow1','Flow2','Leakage'};
    
    block.ValvePos.type = 'in';
    block.ValvePos.IC = 0;
    Inlet.ValvePos = block.ValvePos.IC;
    
    block.Intake.type = 'in';
    block.Intake.IC.T = 500;
    block.Intake.IC.N2 = 0.0131;
    block.Intake.IC.O2 = 0.0035;
    
    Inlet.Intake = block.Intake.IC;
    spec = fieldnames(Inlet.Intake);
        for i = 1:length(spec)
            Out.Leakage.(spec{i}) = block.leakVal*Inlet.Intake.(spec{i});
            Out.Flow1.(spec{i}) = Inlet.Intake.(spec{i})*(1-block.leakVal)*(1-Inlet.ValvePos);
            Out.Flow2.(spec{i}) = Inlet.Intake.(spec{i})*(1-block.leakVal)*Inlet.ValvePos;
        end
    Out.Leakage.T = Inlet.Intake.T;
    Out.Flow1.T = Inlet.Intake.T;
    Out.Flow2.T = Inlet.Intake.T;

    block.Flow2.type = 'out';
    block.Flow2.IC = Out.Flow2;
    
    block.Leakage.type = 'out';
    block.Leakage.IC = Out.Leakage;
    
    block.Flow1.type = 'out';
    block.Flow1.IC = Out.Flow1;
    
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
elseif length(varargin) == 2
    Inlet = varargin{2};
    block.ValvePos.IC = Inlet.ValvePos;
    block.Intake.IC = Inlet.Intake;
    spec = fieldnames(Inlet.Intake);
    for i = 1:length(spec)
        Out.Leakage.(spec{i}) = block.leakVal*Inlet.Intake.(spec{i});
        Out.Flow1.(spec{i}) = Inlet.Intake.(spec{i})*(1-block.leakVal)*(1-Inlet.ValvePos);
        Out.Flow2.(spec{i}) = Inlet.Intake.(spec{i})*(1-block.leakVal)*Inlet.ValvePos;
    end
    Out.Leakage.T = Inlet.Intake.T;
    Out.Flow1.T = Inlet.Intake.T;
    Out.Flow2.T = Inlet.Intake.T;
end