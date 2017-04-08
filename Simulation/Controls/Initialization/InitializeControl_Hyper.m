function block = InitializeControl_Hyper(varargin)
% Controls for recoperated micro-turbine
% Three (3) inlets: Turbine exit temperature (TET), RPM, Combustor equivelance
% Three (3) outlets: Generator Power, fuel flow rate and byass vale postion
% Four (4) states: desired RPM, generator power, fuel flow and bypass valve position
block = varargin{1};
if length(varargin)==1 %first initialization
    
    block.InitialFuel = block.NominalPower/(block.EstimatedEfficiency*(block.Fuel.CH4*802952.15+block.Fuel.CO*305200+block.Fuel.H2*240424)); %molar flow rate of fuel inlet
    block.InletPorts = {'TET','RPM','WTurbine','WCompressor'};
    block.TET.IC = block.TETset; 
    block.RPM.IC = block.RPMdesign; 
    block.WTurbine.IC = 4*block.NominalPower;
    block.WCompressor.IC = 3*block.NominalPower;
    
    block.OutletPorts = {'GenPower','FuelFlow','ColdBypass','HotBypass','BleedValve'};
    block.GenPower.IC = block.NominalPower/block.GenEfficiency; 
    block.FuelFlow.IC = block.InitialFuel;
    block.InitializeError = 1;
    block.ColdBypass.IC = 0;
    block.HotBypass.IC = 0;
    block.BleedValve.IC = 0;

    block.Scale = [block.NominalPower/block.GenEfficiency; block.InitialFuel; 1; 1; 1];
    block.IC = [1; 1; block.ColdBypass.IC;block.HotBypass.IC; block.BleedValve.IC;];% inital condition 
    block.P_Difference = {};
end

if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    PowerSet = PowerDemandLookup(0);
    block.RPM.IC = Inlet.RPM;
    block.WTurbine.IC = Inlet.WTurbine;
    block.WCompressor.IC = Inlet.WCompressor;
    block.TET.IC = Inlet.TET;
    
    block.GenPower.IC = PowerSet/block.GenEfficiency;
    PowerError =(PowerSet - (Inlet.WTurbine-Inlet.WCompressor)*block.GenEfficiency)/block.NominalPower;
    block.InitializeError = abs(PowerError);
    
    if abs(PowerError)<0.06
        a = .3;
    else a = .15;
    end
    block.FuelFlow.IC = block.FuelFlow.IC + a*PowerError*block.Scale(3); %adjust fuel flow to control power;
    block.ColdBypass.IC = block.IC(3);
    block.HotBypass.IC = block.IC(4);
    block.BleedValve.IC = block.IC(5);
    
    block.IC(1) = block.GenPower.IC/(block.NominalPower/block.GenEfficiency);%adjust initial Power
    block.IC(2) = block.FuelFlow.IC/block.Scale(3); %adjust fuel flow to control TET
end