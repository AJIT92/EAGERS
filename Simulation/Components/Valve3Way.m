function Out = Valve3Way(t,Y,Inlet,block,string1)
% a simple flow splitting block with 1 inlet flow,and 2 outlets
% Two (2) inlets: inlet flow ,  valve postion
% Two (2) outlets:  Flow1 , and Flow2 
% Zero (0) states:
if strcmp(string1,'Outlet')
    spec = fieldnames(Inlet.Inlet);
    for i = 1:length(spec)
        Out.Out1.(spec{i}) = Inlet.ValvePos*Inlet.Inlet.(spec{i});
        Out.Out2.(spec{i}) = (1-Inlet.ValvePos)*Inlet.Inlet.(spec{i});
    end
    Out.Out1.T = Inlet.Inlet.T;
    Out.Out2.T = Inlet.Inlet.T;
elseif strcmp(string1,'dY')
    %no states
end