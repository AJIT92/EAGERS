%This builds a recouperated gas turbine model from the component blocks. 
function Plant = GasTurbine
Plant.NominalPower= 60;

NaturalGas.CH4 = 0.9;
NaturalGas.CO = 0.04;
NaturalGas.CO2 = 0.04;
NaturalGas.H2 = 0;
NaturalGas.H2O = 0;
NaturalGas.N2 = 0.02;
Fuel = NaturalGas;
Air.N2 = .79;
Air.O2 = .21;
    
%% Components
Plant.Components.AirSource.type = 'Source'; %fresh air 
Plant.Components.AirSource.name = 'AirSource';
Plant.Components.AirSource.InitialComposition = Air;
Plant.Components.AirSource.connections = {'';'';'';};

Plant.Components.FuelSource.type = 'Source';
Plant.Components.FuelSource.name = 'FuelSource';
Plant.Components.FuelSource.InitialComposition = Fuel;
Plant.Components.FuelSource.connections = {'';'';'Controller.FuelFlow';};

Plant.Components.Comp.type = 'Compressor';
Plant.Components.Comp.name = 'Comp';
Plant.Components.Comp.Map = 'RadialCompressorT21'; % Loads a saved compressor map
Plant.Components.Comp.connections = {'AirSource.Outlet';'';'HX1.ColdPin';'Shaft.RPM';};
Plant.Components.Comp.Mass = .300;%(kg)

Plant.Components.Comp.PeakEfficiency = 0.7598;
Plant.Components.Comp.Tdesign = 288;%Design temp(K)
Plant.Components.Comp.RPMdesign = 96000;%Design RPM
Plant.Components.Comp.FlowDesign = 0.4894;%design flow rate(kg/Sec)%0.47
Plant.Components.Comp.Pdesign = 4.0889;%design pressure ratio
Plant.Components.Comp.TagInf = {'Flow';'Beta';'NRPM';'Power';'PR';'Nflow';'Temperature';'MassFlow';'Eff'};

Plant.Components.Shaft.type = 'Shaft';
Plant.Components.Shaft.name = 'Shaft';
Plant.Components.Shaft.RPMinit = 96000;
Plant.Components.Shaft.Length = .1;%Shaft Length
Plant.Components.Shaft.Radius = 0.15;
Plant.Components.Shaft.Density = 800;%Shaft Density
Plant.Components.Shaft.connections = {'Turb.PowerTurb';'Comp.Work';'Controller.GenPower'};
Plant.Components.Shaft.TagInf = {'RPM';};

Plant.Components.HX1.type = 'HeatExchanger';
Plant.Components.HX1.name = 'HX1';
Plant.Components.HX1.direction = 'counterflow'; % 'coflow', or 'counterflow' or 'crossflow'
Plant.Components.HX1.columns = 5;
Plant.Components.HX1.rows =1;
Plant.Components.HX1.sizemethod = 'Effectiveness'; %method for sizing HX to initial conditions. Options are: 'fixed' fixed size heat exchanger during intialization, 'ColdT' sizes to match cold exit temp, 'HotT' sizes to match hot ext temp, 'Effectiveness' sizes to match a target effectiveness: % of ideal (infinite area) heat transfer (with no conduction between nodes)
Plant.Components.HX1.Target = 0.7992; %can be numeric or a string of block.property, ex 'Controller.HeaterTarget'. If it can't reach the temperature target it defaults to 98% effectiveness
Plant.Components.HX1.Mass = 5; %mass in kg
Plant.Components.HX1.Vol = .1; % volume in m^3
Plant.Components.HX1.Cold_T_init = 400; %initial guess temperaure of cold inlet
Plant.Components.HX1.Hot_T_init = 900; %initial guess temperaure of hot inlet
Plant.Components.HX1.ColdSpecIn = Air;
Plant.Components.HX1.ColdFlowInit = 0.48/28.84;
Plant.Components.HX1.HotSpecIn = Air;
Plant.Components.HX1.HotFlowInit = 0.48/28.84;
Plant.Components.HX1.connections = {'Comp.Flow';'Turb.Outlet';'Comb.Pin';''};
Plant.Components.HX1.TagInf = {'ColdOut';'HotOut';'Effectiveness';'NetImbalance'};

Plant.Components.Comb.type = 'Combustor';
Plant.Components.Comb.name = 'Comb';
Plant.Components.Comb.connections = {'HX1.ColdOut';'FuelSource.Outlet';'Mix1.Pin';};
Plant.Components.Comb.EquivSet = 0.9; %design equiv ratio
Plant.Components.Comb.MassWall = .150; % Mass of wall between bypass flow and combustor flow.
Plant.Components.Comb.MassCasing = .150; % Mass of casing between ambient and bypass flow.
Plant.Components.Comb.Vol = 0.0628; % volume of combustor in m^3
Plant.Components.Comb.TagInf = {'EquivelanceRatio';'Temperatures';'MassFlow'};

Plant.Components.Mix1.type = 'MixingVolume';
Plant.Components.Mix1.name = 'Mix1';
Plant.Components.Mix1.Vol = 0.1;
Plant.Components.Mix1.inlets = 2;
Plant.Components.Mix1.SpeciesInit.CO2 = .1;
Plant.Components.Mix1.SpeciesInit.CO = 0.01;
Plant.Components.Mix1.SpeciesInit.H2O = .1;
Plant.Components.Mix1.SpeciesInit.N2 = .7;
Plant.Components.Mix1.SpeciesInit.O2 = .1;
Plant.Components.Mix1.Tinit = 1200;
Plant.Components.Mix1.connections = {'Comb.Main';'Comb.Bypass';'Turb.Pin'};
Plant.Components.Mix1.TagInf = {'MassFlow';'Temperature';};

Plant.Components.Turb.type = 'Turbine';
Plant.Components.Turb.name = 'Turb';
Plant.Components.Turb.Map = 'RadialTurbineT23'; % Loads a saved compressor map
Plant.Components.Turb.PeakEfficiency = 0.7704;
Plant.Components.Turb.Tdesign = 1200;%Design temp(K)
Plant.Components.Turb.RPMdesign = 96000;%Design RPM
Plant.Components.Turb.FlowDesign = 0.4894;%design flow rate(kg/sec)%0.47
Plant.Components.Turb.Pdesign = 4.0889;%design pressure ratio
Plant.Components.Turb.Mass = .2;%(kg)
Plant.Components.Turb.connections = {'Mix1.Outlet';'HX1.HotPin';'Shaft.RPM'};
Plant.Components.Turb.TagInf = {'TET';'Power';'PR';'Nflow';'Beta';'NRPM';'Efficiency';'MassFlow'};

%% Controls (note: controls can have specification that depends upon a initialized variable of a component)
Plant.Controls.Controller.type = 'ControlRecouperatedGasTurbine';
Plant.Controls.Controller.name = 'Controller';
Plant.Controls.Controller.NominalPower = Plant.NominalPower;
Plant.Controls.Controller.TETset = 907.8;
Plant.Controls.Controller.RPMdesign = 96000;
Plant.Controls.Controller.GenEfficiency = 0.97;
Plant.Controls.Controller.EstimatedEfficiency = .25;
Plant.Controls.Controller.Fuel = Plant.Components.FuelSource.InitialComposition;
Plant.Controls.Controller.Gain = [4e-4; 4e-3; 5e-4;];%gain for TET control, Power control, fuel control
Plant.Controls.Controller.PropGain = [8e-3; 2e-0; .3;];
Plant.Controls.Controller.TagInf = {'TET';'GenPower';'FuelFlow';'SteadyPower';'Efficiency'};
Plant.Controls.Controller.connections = {'Turb.TET';'Shaft.RPM';'Turb.PowerTurb';'Comp.Work';};

Plant.Scope = {'Controller.FuelFlow';'Shaft.RPM';'Comp.MassFlow';'Turb.TET';}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'Controller.Efficiency';'Controller.GenPower';}];