function Out = Control_Hyper(t,Y, Inlet,block,string1)
% Controls for recoperated micro-turbine
% Three (3) inlets: Turbine exit temperature (TET), RPM, Leakage Valve Interp Value
% Three (3) outlets: Generator Power, fuel flow rate to combustor, fuel flow rate to oxidizer, and leakage valve postion
% Four (4) states: desired RPM, generator power, molar fuel flow and compressor leakage valve position
global Tags 
SteadyPower = Inlet.WTurbine - Inlet.WCompressor;
PowerSet = PowerDemandLookup(t);
P_Gain = block.PropGain.*block.Scale;

PowerError = PowerSet-Y(2);
GenPower =(P_Gain(2)*PowerError + Y(2));
FuelError =(PowerSet - GenPower*block.GenEfficiency)/PowerSet;

if strcmp(string1,'Outlet')
    if t>0
        Out.GenPower = GenPower;
        Out.FuelFlow = (P_Gain(2)*FuelError + Y(2));
        Out.ColdBypass = Y(3);
        Out.HotBypass = Y(4);
        Out.BleedValve = Y(5);
    else
        Out.GenPower = block.GenPower.IC;
        Out.FuelFlow = block.FuelFlow.IC;
        Out.ColdBypass = Y(3);
        Out.HotBypass = Y(4);
        Out.BleedValve = Y(5);
    end
    Tags.(block.name).TET = Inlet.TET;
    Tags.(block.name).GenPower = Out.GenPower*block.GenEfficiency;
    Tags.(block.name).FuelFlow = Out.FuelFlow;
    Tags.(block.name).SteadyPower = SteadyPower;
    Tags.(block.name).Efficiency = GenPower/(Out.FuelFlow*(block.Fuel.CH4*802952.15+block.Fuel.CO*305200+block.Fuel.H2*240424)); %molar flow rate of fuel inlet
elseif strcmp(string1, 'dY')
    dY = Y*0;
    dY(1) = PowerError*block.IntGain(2);
    dY(2) = FuelError*block.IntGain(3);
    dY(3) = 0;
    Tags.(block.name).dY1 = PowerError;
    Tags.(block.name).dY2 = FuelError;
    Out = dY;
end