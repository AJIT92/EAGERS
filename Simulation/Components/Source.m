function Out = Source(t,Inlet,block)
% a source of regular fuel with no states
% the type of fuel (species) must be initialized first
% Two (2) inlets: Temperature and flow rate
% no state (will only ever get called upon to calculate outlet)
global Tags
Out.Outlet.T = Inlet.Temperature;
speciesNames = fieldnames(Inlet.Species);
for i = 1:1:length(speciesNames)
    Out.Outlet.(speciesNames{i}) = Inlet.Species.(speciesNames{i})*Inlet.Flow;
end
Tags.(block.name).Outlet = Out.Outlet;