function block = InitializeControlRecouperatedGasTurbine(varargin)
% Controls for recoperated micro-turbine
% Three (3) inlets: Turbine exit temperature (TET), RPM, Combustor equivelance
% Three (3) outlets: Generator Power, fuel flow rae and byass vale postion
% Four (4) states: desired RPM, generator power, fuel flow and bypass valve position
block = varargin{1};
if length(varargin)==1 %first initialization
    block.description = {'Target RPM';'Power Pulled by the Generator';'Fuel Flow';};
    
    block.InitialFuel = block.NominalPower/(block.EstimatedEfficiency*(block.Fuel.CH4*802952.15+block.Fuel.CO*305200+block.Fuel.H2*240424)); %molar flow rate of fuel inlet
    
    block.InletPorts = {'TET','RPM','WTurbine','WCompressor'};
    block.TET.IC = block.TETset; 
    block.RPM.IC = block.RPMdesign; 
    block.WTurbine.IC = 4*block.NominalPower;
    block.WCompressor.IC = 3*block.NominalPower;
    
    block.OutletPorts = {'GenPower','FuelFlow'};
    block.GenPower.IC = block.NominalPower/block.GenEfficiency; 
    block.FuelFlow.IC = block.InitialFuel;
    block.InitializeError = 1;
    
    block.Scale = [block.RPMdesign; block.NominalPower/block.GenEfficiency; block.InitialFuel;];
    block.IC = [1; 1; 1;]; % inital condition 
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
    TETerror = (block.TETset-Inlet.TET)/400; %error in TET -- should change RPM
    PowerError =(PowerSet - (Inlet.WTurbine-Inlet.WCompressor)*block.GenEfficiency)/block.NominalPower;
    block.InitializeError = abs(PowerError);
    
    if abs(PowerError)<0.06
        a = .3;
    else a = .15;
    end
    block.FuelFlow.IC = block.FuelFlow.IC + a*PowerError*block.Scale(3); %adjust fuel flow to control power;
    block.IC(1) = Inlet.RPM/block.RPMdesign;%adjust initial RPM
    block.IC(2) = block.GenPower.IC/(block.NominalPower/block.GenEfficiency);%adjust initial Power
    block.IC(3) = block.FuelFlow.IC/block.Scale(3); %adjust fuel flow to control TET
end