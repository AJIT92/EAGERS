%% SOFC with external reformer
function Plant = SOFCsystem
Plant.NominalPower=300;

Steam2Carbon = 2.5;

Air.N2 = .79;
Air.O2 = .21;

NaturalGas.CH4 = 0.9;
NaturalGas.CO = 0.04;
NaturalGas.CO2 = 0.04;
NaturalGas.H2 = 0;
NaturalGas.H2O = 0;
NaturalGas.N2 = 0.02;

%some initial guesses for other streams (don't have to be accurate)
FCexit.CH4 = .0001;
FCexit.CO = .0999;
FCexit.CO2 = .25;
FCexit.H2 = .1;
FCexit.H2O = .54;
FCexit.N2 = .01;

Oxidized.CH4 = 0;
Oxidized.CO = 0.03;
Oxidized.CO2 = .25;
Oxidized.H2 = 0;
Oxidized.H2O = .72;
Oxidized.N2 = .01;

a = Steam2Carbon*(NaturalGas.CH4+.5*NaturalGas.CO)/(FCexit.H2O -Steam2Carbon*(FCexit.CH4 + .5*FCexit.CO));
Humidified.CH4 = (a*FCexit.CH4 + NaturalGas.CH4)/(a+1);
Humidified.CO = (a*FCexit.CO + NaturalGas.CO)/(a+1);
Humidified.CO2 = (a*FCexit.CO2 + NaturalGas.CO2)/(a+1);
Humidified.H2 = (a*FCexit.H2 + NaturalGas.H2)/(a+1);
Humidified.H2O = (a*FCexit.H2O + NaturalGas.H2O)/(a+1);
Humidified.N2 = (a*FCexit.N2 + NaturalGas.N2)/(a+1);


%Fuel Source
Plant.Components.FuelSource.type = 'Source';
Plant.Components.FuelSource.name = 'FuelSource';
Plant.Components.FuelSource.InitialComposition = NaturalGas;
Plant.Components.FuelSource.connections = {300;'';'Controller.FuelFlow';};

%Blower
Plant.Components.Blower.type = 'Blower';
Plant.Components.Blower.name = 'Blower';
Plant.Components.Blower.Map = 'RadialCompressorT6'; % Loads a saved compressor map
Plant.Components.Blower.InitialComposition = Air;
Plant.Components.Blower.PeakEfficiency = 0.65;
Plant.Components.Blower.Tdesign = 300;%Design temp(K)
Plant.Components.Blower.RPMdesign = 30000;%Design RPM
Plant.Components.Blower.FlowDesign = 1.4;%design flow rate(kg/Sec)
Plant.Components.Blower.Pdesign = 1.6;%design pressure ratio
Plant.Components.Blower.connections = {'';'';'';'HX1.ColdPin';'Controller.Blower';}; %no species or temperature connection, stick to IC
Plant.Components.Blower.TagInf = {'Flow';'NRPM';'Power';'PR';'Nflow';'Temperature';'MassFlow';'Eff'};

%recircValve
Plant.Components.recircValve.type = 'Valve3Way';
Plant.Components.recircValve.name = 'recircValve';
Plant.Components.recircValve.InitialFlowIn = FCexit;
Plant.Components.recircValve.InitialFlowIn.T = 1050;
Plant.Components.recircValve.connections = {'FC1.AnodeOut','Controller.AnodeRecirc'};
Plant.Components.recircValve.PercOpen = 0.5; %INITIAL valve position

%Mixing
Plant.Components.Mix1.type = 'MixingVolume';
Plant.Components.Mix1.name = 'Mix1';
Plant.Components.Mix1.Vol = 0.1;
Plant.Components.Mix1.inlets = 2;
Plant.Components.Mix1.SpeciesInit = Humidified;
Plant.Components.Mix1.Tinit = 800;
Plant.Components.Mix1.connections = {'FuelSource.Outlet';'recircValve.Out1';'Reformer.ReformedPin'};
Plant.Components.Mix1.TagInf = {'MassFlow';'Temperature';};

%Oxidizer
Plant.Components.Oxidizer.type = 'Oxidizer';
Plant.Components.Oxidizer.name = 'Oxidizer';
Plant.Components.Oxidizer.inlets = 2;
Plant.Components.Oxidizer.InitialFlowOut = Oxidized;
Plant.Components.Oxidizer.InitialFlowOut.T = 1050;
Plant.Components.Oxidizer.Vol = .1; %volume in m^3
Plant.Components.Oxidizer.connections = {'recircValve.Out2';'FC1.CathodeOut';'Reformer.CooledPin';};
Plant.Components.Oxidizer.TagInf = {'EquivelanceRatio';'Temperature';'MassFlow'};

%Reformer
Plant.Components.Reformer.type = 'Reformer';
Plant.Components.Reformer.name = 'Reformer';
Plant.Components.Reformer.direction = 'counterflow';  % 'co-flow', or 'counterflow'
Plant.Components.Reformer.hotStream = 1; % true if there is a second stream providing heat, false if not
Plant.Components.Reformer.nodes = 5; % # of nodes to simulate
Plant.Components.Reformer.ReformTarget = .75; %    If there is a 2nd stream providing heat, it can use one of the following 4 targets to initialize reformer: Reforming percent, cold outlet temp, hot outlet temp, effectiveness of heat transfer from hot stream. Note that due to the reforming chemistry, a simple heat exchanger effectivess is not applicable.
Plant.Components.Reformer.method = 'RefPerc'; %method for sizing reformer to initial conditions. Options are: 'none' fixed size during intialization, 'RefPerc' sizes for a specific external reforming percent, 'ColdT' sizes to match cold exit temp, 'HotT' sizes to match hot ext temp, 'Effectiveness' sizes to match a target effectiveness: % of ideal (infinite area) heat transfer (with no conduction between nodes)
Plant.Components.Reformer.InletGuess = Humidified;
Plant.Components.Reformer.connections = {'Mix1.Outlet';'FC1.AnodePressureIn';'Oxidizer.Flow';'HX1.HotPin';};
Plant.Components.Reformer.TagInf = {'H2flow';'CH4';};

% Air Bypass
Plant.Components.bypassValve.type = 'Valve3Way';
Plant.Components.bypassValve.name = 'bypassValve';
Plant.Components.bypassValve.InitialFlowIn = Air;
Plant.Components.bypassValve.InitialFlowIn.T = 300;
Plant.Components.bypassValve.connections = {'Blower.Outlet','Controller.HeaterBypass'};
Plant.Components.bypassValve.PercOpen = 0.04; %INITIAL valve position

% Air Heat Exchanger
Plant.Components.HX1.type = 'HeatExchanger';
Plant.Components.HX1.name = 'HX1';
Plant.Components.HX1.direction = 2;  % parallel =1, countercurrent =2, crossflow =3
Plant.Components.HX1.columns = 10;
Plant.Components.HX1.rows =1;
Plant.Components.HX1.sizemethod = 'Effectiveness'; %method for sizing HX to initial conditions. Options are: 'fixed' fixed size heat exchanger during intialization, 'ColdT' sizes to match cold exit temp, 'HotT' sizes to match hot ext temp, 'Effectiveness' sizes to match a target effectiveness: % of ideal (infinite area) heat transfer (with no conduction between nodes)
Plant.Components.HX1.Target = 0.95; %can be numeric or a string of block.property, ex 'Controller.HeaterTarget'. If it can't reach the temperature target it defaults to 98% effectiveness
Plant.Components.HX1.Mass = 25; %mass in kg
Plant.Components.HX1.Vol = 1; % volume in m^3
Plant.Components.HX1.Cold_T_init = 300; %initial guess temperaure of cold inlet
Plant.Components.HX1.Hot_T_init = 1000; %initial guess temperaure of hot inlet
Plant.Components.HX1.ColdSpecIn = Air;
Plant.Components.HX1.HotSpecIn = Oxidized;
Plant.Components.HX1.connections = {'bypassValve.Out2';'Reformer.Cooled';'Mix2.Pin';'';}; 
Plant.Components.HX1.TagInf = {'ColdOut';'HotOut';'Effectiveness';'NetImbalance'};

% Air Bypass mixing
Plant.Components.Mix2.type = 'MixingVolume';
Plant.Components.Mix2.name = 'Mix2';
Plant.Components.Mix2.Vol = 0.1;
Plant.Components.Mix2.inlets = 2;
Plant.Components.Mix2.SpeciesInit= Air;
Plant.Components.Mix2.Tinit = 970;
Plant.Components.Mix2.connections = {'bypassValve.Out1';'HX1.ColdOut';'FC1.CathodePressureIn'};
Plant.Components.Mix2.TagInf = {'MassFlow';'Temperature';};

% fuel cell
Plant.Components.FC1.type = 'FuelCell';
Plant.Components.FC1.name = 'FC1';
Plant.Components.FC1.FCtype = 'SOFC'; %SOFC, or MCFC or oxySOFC
Plant.Components.FC1.Reformer = 'external'; % options are 'internal' for indirect internal reforming, 'direct' for direct internal reforming, 'adiabatic' for external reforming using the heat from the anode (over-rides S2C ratio), 'pox' partial oxidation
Plant.Components.FC1.RefPerc = Plant.Components.Reformer.ReformTarget; % necessary for internal or external reformer. This is the percent towards equilibrium occurin in the reformer plaes
% Plant.Components.FC1.RefSpacing = 1;% necessary for internal reformer. This is the # of active cells between reformer plates
Plant.Components.FC1.direction = 'coflow'; % 'co-flow', or 'counterflow' or 'crossflow'
Plant.Components.FC1.ClosedCathode = 0; %0 means air or some excess flow of O2 in the cathode used as primary means of temerature control (initializations hold to design fuel utilization), 1 means closed end cathode, or simply a fixed oxygen utilization, cooling is done with excess fuel, and the design voltage is met during initialization
Plant.Components.FC1.CoolingStream = 'cathode'; % choices are 'anode' or 'cathode'. Determines which flow is increased to reach desired temperature gradient.
Plant.Components.FC1.PressureRatio = 1.2;
Plant.Components.FC1.columns = 4;
Plant.Components.FC1.rows = 1;
Plant.Components.FC1.RatedStack_kW = 300; %Nominal Stack Power in kW
Plant.Components.FC1.Fuel = NaturalGas; %raw fuel fed to the system (calculates anode recirculation from this to reach S2C)
Plant.Components.FC1.Oxidant = Air; %initial oxidant composition
Plant.Components.FC1.Steam2Carbon = Steam2Carbon;
Plant.Components.FC1.method = 'Achenbach'; %Determines reforming reaction kinetics options: 'Achenbach' , 'Leinfelder' , 'Drescher'   
Plant.Components.FC1.L_Cell= .09;  %Cell length in meters
Plant.Components.FC1.W_Cell = .09;  %Cell Width in meters  
Plant.Components.FC1.DesignTarget = 'power density'; %options are 'power density', 'voltage', or 'current density'
Plant.Components.FC1.DesignTargetValue = 450; % power density specified in mW/cm^2, voltage specified in V/cell, current density specified in A/cm^2
Plant.Components.FC1.Cells = ceil(Plant.Components.FC1.RatedStack_kW*100/(Plant.Components.FC1.L_Cell*Plant.Components.FC1.W_Cell*Plant.Components.FC1.DesignTargetValue)); %# of cells in stack
Plant.Components.FC1.deltaTStack = 70; %temperature difference from cathode inlet to cathode outlet
Plant.Components.FC1.TpenAvg = 1023;% 750 C, average electrolyte operating temperature
Plant.Components.FC1.FuelUtilization =.75; %fuel utilization (net hydrogen consumed/ maximum hydrogen produced with 100% Co and CH4 conversion
Plant.Components.FC1.AnodePdrop = 2; %design anode pressure drop
Plant.Components.FC1.CathodePdrop = 10; %Design cathode pressure drop
Plant.Components.FC1.Map = 'SOFC_map'; %Radiative heat transfer view factors, imported from CAD
Plant.Components.FC1.connections = {'Controller.Current';'Mix2.Outlet';'Reformer.Reformed';'Oxidizer.Pin';'Oxidizer.Pin';};
Plant.Components.FC1.TagInf = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';'TcathOut';'ASR';'O2utilization';}; %Tags to record at each step
Plant.Components.FC1.TagFinal = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';'TcathOut';'ASR';'O2utilization'}; %Tags to record at the final step

% Controller %% order is heater bypass, blower power, anode recirculation
Plant.Controls.Controller.type = 'ControlSOFCsystem';
Plant.Controls.Controller.name = 'Controller';
Plant.Controls.Controller.Target = {'FC1.TpenAvg';'FC1.deltaTStack';'FC1.Steam2Carbon';};
Plant.Controls.Controller.Fuel = NaturalGas;
Plant.Controls.Controller.Cells = 'FC1.Cells';
Plant.Controls.Controller.Utilization = 'FC1.FuelUtilization';
Plant.Controls.Controller.InitialAnodeRecirc = 'FC1.Recirc.Anode';
Plant.Controls.Controller.InitConditions = {'bypassValve.PercOpen';'Blower.NominalPower';'FC1.Current';}; %heater bypass, blower power, FC current  [note: fuel flow and recirculaion are calculated and are not states]
Plant.Controls.Controller.Gain = [1e-3;1e-5;1e-3];
Plant.Controls.Controller.PropGain = [.35;.05;.5];
Plant.Controls.Controller.TagInf = {'Bypass';'Blower';'Recirculation';'FuelFlow';'Current';};
Plant.Controls.Controller.connections = {'FC1.MeasureTcathOut';'Mix2.Temperature';'FC1.MeasureTpen';'FC1.MeasureVoltage';};

Plant.Scope = {'Controller.Blower';'Controller.Bypass';'Controller.Current'}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'FC1.StackdeltaT';'Blower.MassFlow';'Blower.NRPM';'FC1.PENavgT';'FC1.Voltage';'FC1.TcathOut';'Mix2.Temperature';'Oxidizer.Temperature';}];