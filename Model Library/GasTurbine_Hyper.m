%This builds a recouperated gas turbine model from the component blocks. 
global Plant PowerDemand PowerTime RunTime Model_dir Ru
Plant = [];

Ru = 8.314; %Universal gas constant (kJ/kmol*K)
Model_dir=strrep(which('GasTurbine_Hyper.m'),fullfile('Models','GasTurbine_Hyper.m'),'');
AddPaths()
Plant.NetPower= 60;
RunTime =1e6;
PowerDemand    = [1 1 0.9 1 1]*Plant.NetPower;
PowerTime = linspace(0,RunTime,length(PowerDemand));
NaturalGas.CH4 = 0.9;
NaturalGas.CO = 0.04;
NaturalGas.CO2 = 0.04;
NaturalGas.H2 = 0;
NaturalGas.H2O = 0;
NaturalGas.N2 = 0.02;
Fuel = NaturalGas;
Air.N2 = .79;
Air.O2 = .21;

Oxidized.CH4 = 0;
Oxidized.CO = 0;
Oxidized.CO2 = .1;
Oxidized.H2 = 0;
Oxidized.H2O = .1;
Oxidized.N2 = .7;
Oxidized.O2 = .1;

%% Components
Plant.Components.AirSource.type = 'Source'; %fresh air 
Plant.Components.AirSource.name = 'AirSource';
Plant.Components.AirSource.InitialComposition.N2 = 0.79;
Plant.Components.AirSource.InitialComposition.O2 = 0.21;
Plant.Components.AirSource.connections = {'';'';'';};

Plant.Components.FuelSource.type = 'Source';
Plant.Components.FuelSource.name = 'FuelSource';
Plant.Components.FuelSource.InitialComposition = Fuel;
Plant.Components.FuelSource.connections = {'';'';'GTcontrolH.FuelFlow';};

Plant.Components.BleedValve.type = 'LeakageValve';
Plant.Components.BleedValve.name = 'BleedValve';
Plant.Components.BleedValve.connections = {'GTcontrolH.BleedValve','Comp.Flow'};
Plant.Components.BleedValve.leakVal = 0.12;

Plant.Components.ColdBypass.type = 'Valve3Way';
Plant.Components.ColdBypass.name = 'ColdBypass';
Plant.Components.ColdBypass.InitialFlowIn = Air;
Plant.Components.ColdBypass.InitialFlowIn.T = 500;
Plant.Components.ColdBypass.PercOpen = 0; %INITIAL valve position
Plant.Components.ColdBypass.connections = {'GTcontrolH.ColdBypass','BleedValve.Flow1'};

Plant.Components.Comp.type = 'Compressor';
Plant.Components.Comp.name = 'Comp';
Plant.Components.Comp.Map = 'RadialCompressor1'; % Loads a saved compressor map
Plant.Components.Comp.connections = {'AirSource.Outlet';'';'HX1.ColdPin';'Shaft.RPM';};
Plant.Components.Comp.Mass = .300;%(kg)
Plant.Components.Comp.PeakEfficiency = 0.6516;
Plant.Components.Comp.Tdesign = 288;%Design temp(K)
Plant.Components.Comp.RPMdesign = 96000;%Design RPM
Plant.Components.Comp.FlowDesign = 0.4680;%design flow rate(kg/Sec)%0.47
Plant.Components.Comp.Pdesign = 3.9658;%design pressure ratio
Plant.Components.Comp.TagInf = {'Flow';'NRPM';'Power';'PR';'Nflow';'Temperature';'MassFlow';'Eff'};

Plant.Components.Shaft.type = 'Shaft';
Plant.Components.Shaft.name = 'Shaft';
Plant.Components.Shaft.RPMinit = 96000;
Plant.Components.Shaft.Length = .1;%Shaft Length
Plant.Components.Shaft.Radius = 0.15;
Plant.Components.Shaft.Density = 800;%Shaft Density
Plant.Components.Shaft.connections = {'Turb.PowerTurb';'Comp.Work';'GTcontrolH.GenPower'};
Plant.Components.Shaft.TagInf = {'RPM';};

Plant.Components.HX1.type = 'HeatExchanger';
Plant.Components.HX1.name = 'HX1';
Plant.Components.HX1.direction = 2;  % parallel =1, countercurrent =2, crossflow =3
Plant.Components.HX1.columns = 5;
Plant.Components.HX1.rows =1;
Plant.Components.HX1.sizemethod = 'Effectiveness'; %method for sizing HX to initial conditions. Options are: 'fixed' fixed size heat exchanger during intialization, 'ColdT' sizes to match cold exit temp, 'HotT' sizes to match hot ext temp, 'Effectiveness' sizes to match a target effectiveness: % of ideal (infinite area) heat transfer (with no conduction between nodes)
Plant.Components.HX1.Target = 0.95; %can be numeric or a string of block.property, ex 'Controller.HeaterTarget'. If it can't reach the temperature target it defaults to 98% effectiveness
Plant.Components.HX1.Mass = 5; %mass in kg
Plant.Components.HX1.Vol = .1; % volume in m^3
Plant.Components.HX1.Cold_T_init = 400; %initial guess temperaure of cold inlet
Plant.Components.HX1.Hot_T_init = 900; %initial guess temperaure of hot inlet
Plant.Components.HX1.ColdSpecIn = Air;
Plant.Components.HX1.ColdFlowInit = 0.48/28.84;
Plant.Components.HX1.HotSpecIn = Air;
Plant.Components.HX1.HotFlowInit = 0.48/28.84;
Plant.Components.HX1.connections = {'ColdBypass.Out1';'Turb.Outlet';'SOFCVol.Pin';''};
Plant.Components.HX1.TagInf = {'ColdOut';'HotOut';'Effectiveness';'NetImbalance'};

Plant.Components.HotBypass.type = 'Valve3Way';
Plant.Components.HotBypass.name = 'HotBypass';
Plant.Components.bypassValve.InitialFlowIn = Air;
Plant.Components.bypassValve.InitialFlowIn.T = 800;
Plant.Components.bypassValve.PercOpen = 0; %INITIAL valve position
Plant.Components.HotBypass.connections = {'GTcontrolH.HotBypass','HX1.ColdOut'};

Plant.Components.SOFCVol.type = 'PlenumVolume';
Plant.Components.SOFCVol.name = 'SOFCVol';
Plant.Components.SOFCVol.Vol = 2;
Plant.Components.SOFCVol.Tinit = 1200;
Plant.Components.SOFCVol.connections = {'HotBypass.Out1';'Oxidizer.Pin'};

Plant.Components.Oxidizer.type = 'Oxidizer';
Plant.Components.Oxidizer.name = 'Oxidizer';
Plant.Components.Oxidizer.inlets = 5;
Plant.Components.Oxidizer.InitialFlowOut = Oxidized;
Plant.Components.Oxidizer.InitialFlowOut.T = 1200;
Plant.Components.Oxidizer.connections = {'SOFCVol.Outlet','HotBypass.Out2','ColdBypass.Out2','FuelSource.Outlet','BleedValve.Leakage','Turb.Pin'};
Plant.Components.Oxidizer.TagInf = {'EquivelanceRatio';'Temperatures';'MassFlow'};

Plant.Components.Turb.type = 'Turbine';
Plant.Components.Turb.name = 'Turb';
Plant.Components.Turb.Map = 'RadialTurbine1'; % Loads a saved compressor map
Plant.Components.Turb.PeakEfficiency = 0.8293;%0.7615
Plant.Components.Turb.Tdesign = 1200;%Design temp(K)
Plant.Components.Turb.RPMdesign = 96000;%Design RPM
Plant.Components.Turb.FlowDesign = 0.4680;%design flow rate(kg/sec)%0.47
Plant.Components.Turb.Pdesign = 3.9658;%design pressure ratio
Plant.Components.Turb.Mass = .2;%(kg)
Plant.Components.Turb.connections = {'Mix1.Outlet';'HX1.HotPin';'Shaft.RPM'};
Plant.Components.Turb.TagInf = {'TET';'Power';'PR';'Nflow';'NRPM';'Efficiency';'MassFlow'};

%% Controls (note: controls can have specification that depends upon a initialized variable of a component)
Plant.Controls.GTcontrolH.type = 'Control_Hyper';
Plant.Controls.GTcontrolH.name = 'GTcontrolH';
Plant.Controls.GTcontrolH.NominalPower = Plant.NetPower;
Plant.Controls.GTcontrolH.TETset = 907.8;
Plant.Controls.GTcontrolH.RPMdesign = 96000;
Plant.Controls.GTcontrolH.GenEfficiency = 0.97;
Plant.Controls.GTcontrolH.EstimatedEfficiency = .25;
Plant.Controls.GTcontrolH.Fuel = Plant.Components.FuelSource.InitialComposition;
Plant.Controls.GTcontrolH.IntGain = [3e-3; 4e-3; 0; 0;];%Load control, Fuel control, Cold Bypass control, Hot Bypass valve control
Plant.Controls.GTcontrolH.PropGain = [1e-2; 2e-0; 0; 0;];
Plant.Controls.GTcontrolH.TagInf = {'TET';'GenPower';'FuelFlow';'SteadyPower';'Efficiency'};
Plant.Controls.GTcontrolH.connections = {'Turb.TET';'Shaft.RPM';'Turb.PowerTurb';'Comp.Work';};%needs leakage valve interp value

Plant.Plot = {'GTcontrolH.Efficiency';'Turb.MassFlow';'Shaft.RPM';'Turb.TET';'GTcontrolH.GenPower';'GTcontrolH.FuelFlow';};
%% run initialization
global modelParam
tic; BuildModel([]); disp(strcat('Time to build/initialize model:',num2str(toc),' seconds'));

CTmap(modelParam.Comp,modelParam.Turb)
%modify control terms
modelParam.GTcontrol.IntGain = [5e-3; 4e-3; 7e-3; 1e-6];
modelParam.GTcontrol.PropGain = [3e-1; 2e-0; 2e-1; 0];

% Run Model
tic; [T, Y] = RunModel; disp(strcat('Time to run model:',num2str(toc),' seconds'));

CTmap(modelParam.Comp,modelParam.Turb)