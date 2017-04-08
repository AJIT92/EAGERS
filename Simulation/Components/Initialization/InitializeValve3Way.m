function block = InitializeValve3Way(varargin)
% a simple flow splitting block with 1 inlet flow,and 2 outlets
% Two (2) inlets: inlet flow ,  valve postion
% Two (2) outlets:  Flow1 , and Flow2
% Zero (0) states: 
block = varargin{1};
if length(varargin)==1 %first initialization
    block.IC =[];% no states
    block.Scale =[];

    if ischar(block.InitialFlowIn)
        block.InitialFlowIn = ComponentProperty(block.InitialFlowIn);
    end
    spec = fieldnames(block.InitialFlowIn);
    for i = 1:1:length(spec)
        FlowIn.(spec{i}) = block.InitialFlowIn.(spec{i});
        Out1.(spec{i}) = block.InitialFlowIn.(spec{i})*block.PercOpen;
        Out2.(spec{i}) = block.InitialFlowIn.(spec{i})*(1-block.PercOpen);
    end
    FlowIn.T = block.InitialFlowIn.T;
    Out1.T = block.InitialFlowIn.T;
    Out2.T = block.InitialFlowIn.T;
    
    block.InletPorts = {'Inlet','ValvePos'};
    block.Inlet.IC = FlowIn;
    block.ValvePos.IC = block.PercOpen;
    
    block.OutletPorts = {'Out1','Out2'};
    block.Out1.IC = Out1;
    block.Out2.IC = Out2;  
    
    block.P_Difference = {};
end

if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    spec = fieldnames(Inlet.Inlet);
    for i = 1:1:length(spec)
        if ~strcmp(spec{i},'T')
            Out1.(spec{i}) = Inlet.Inlet.(spec{i})*Inlet.ValvePos;
            Out2.(spec{i}) = Inlet.Inlet.(spec{i})*(1-Inlet.ValvePos);
        end
    end
    Out1.T = Inlet.Inlet.T;
    Out2.T = Inlet.Inlet.T;
    block.Out1.IC = Out1;
    block.Out2.IC = Out2;  
end 