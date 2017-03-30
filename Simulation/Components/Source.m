function Out = Source(t,Y,Inlet,block,string1)
% a source of regular fuel with no states
% the type of fuel (species) must be initialized first
% Two (2) inlets: Temperature and flow rate
% no state (will only ever get called upon to calculate outlet)
global Tags
if strcmp(string1,'Outlet')
    Out.Outlet.T = Inlet.Temperature;
    speciesNames = fieldnames(Inlet.Species);
    for i = 1:1:length(speciesNames)
        Out.Outlet.(speciesNames{i}) = Inlet.Species.(speciesNames{i})*Inlet.Flow;
    end
    Tags.(block.name).Outlet = Out.Outlet;
elseif strcmp(string1,'dY')
    %no states
end