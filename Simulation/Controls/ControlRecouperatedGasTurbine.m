function Out = ControlRecouperatedGasTurbine(t,Y, Inlet,block,string1)
% Controls for recoperated micro-turbine
% Three (3) inlets: Turbine exit temperature (TET), RPM, Combustor equivelance
% Three (3) outlets: Generator Power, fuel flow rae and byass vale postion
% Three (3) states: desired RPM, generator power, molar fuel flow
global Tags 
PowerSet = PowerDemandLookup(t);
P_Gain = block.PropGain.*block.Scale;
I_Gain = block.Gain.*block.Scale;

TETerror = (Inlet.TET-block.TETset)/100; 
RPMfromTET =(P_Gain(1)*TETerror + Y(1));
RPMerror = (Inlet.RPM - RPMfromTET)/block.RPMdesign;
GenPower =(P_Gain(2)*RPMerror + Y(2));
FuelError =(PowerSet - GenPower*block.GenEfficiency)/PowerSet;
PowerError = RPMerror;

if strcmp(string1,'Outlet')
    Out.GenPower = GenPower;
    Out.FuelFlow = Y(3) + P_Gain(3)*FuelError;
    Tags.(block.name).TET = Inlet.TET;
    Tags.(block.name).GenPower = Out.GenPower*block.GenEfficiency;
    Tags.(block.name).FuelFlow = Out.FuelFlow;
    Tags.(block.name).SteadyPower = Inlet.WTurbine - Inlet.WCompressor;
    Tags.(block.name).Efficiency = GenPower/(Out.FuelFlow*(block.Fuel.CH4*802952.15+block.Fuel.CO*305200+block.Fuel.H2*240424)); %molar flow rate of fuel inlet
elseif strcmp(string1, 'dY')
    dY = Y*0;
    dY(1) = TETerror*I_Gain(1);
    dY(2) = PowerError*I_Gain(2);
    dY(3) = FuelError*I_Gain(3);
    Tags.(block.name).dY1 = TETerror;
    Tags.(block.name).dY2 = PowerError;
    Tags.(block.name).dY3 = FuelError;
    Out = dY;
end