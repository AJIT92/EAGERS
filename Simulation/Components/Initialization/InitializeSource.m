function block = InitializeSource(varargin)
% a source block with no states
% Three (3) inlets: Temperature, species, and flow rate
block = varargin{1};

if length(varargin)==1 %first initialization
    block.Scale =[];%no states exist
    block.IC = ones(length(block.Scale),1);
    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    block.InletPorts = {'Temperature','Species','Flow'};
    block.Temperature.IC = 298; 
    block.Species.IC = block.InitialComposition; 
    block.Flow.IC = 1;  
    
    block.OutletPorts = {'Outlet'};
    block.Outlet.IC.T = block.Temperature.IC;
    speciesNames = fieldnames(block.Species.IC);
    for i = 1:1:length(speciesNames)
        block.Outlet.IC.(speciesNames{i}) = block.Species.IC.(speciesNames{i})*block.Flow.IC;
    end
    
    block.P_Difference = {};
end

if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    block.Temperature.IC = Inlet.Temperature;
    block.Species.IC = Inlet.Species;  
    block.Flow.IC = Inlet.Flow;
    block.Outlet.IC.T = Inlet.Temperature;
    speciesNames = fieldnames(Inlet.Species);
    for i = 1:1:length(speciesNames)
        block.Outlet.IC.(speciesNames{i}) = Inlet.Species.(speciesNames{i})*Inlet.Flow;
    end
end