function Out = LeakageValve(t,Y,Inlet,block)
%2 inlets: Inlet flow and valve position
%3 outlets: Leakage and Flow1 and flow 2
if strcmp(string1,'Outlet')
    specName = fieldnames(Inlet.Intake);
    for i = 1:length(specName)
        Out.Leakage.(specName{i}) = block.leakVal*Inlet.Intake.(specName{i});
        Out.Flow1.(specName{i}) = Inlet.Intake.(specName{i})*(1-block.leakVal)*(1-Inlet.ValvePos);
        Out.Flow2.(specName{i}) = Inlet.Intake.(specName{i})*(1-block.leakVal)*Inlet.ValvePos;
    end
    Out.Leakage.T = Inlet.Intake.T;
    Out.Flow1.T = Inlet.Intake.T;
    Out.Flow2.T = Inlet.Intake.T;
elseif strcmp(string1,'dY')
    %no states
end