%% SOFC system with internal or external reformer
function Plant = SOFCsystem
Reformer = 'external';
% options are 'internal' for indirect internal reforming, 
% 'direct' for  direct internal reforming, 
% 'adiabatic' for external reforming using the heat from the anode (over-rides S2C ratio), 
% 'external' for a reformer with heat recovery from the oxidizer, and 
% 'pox' partial oxidation

Plant.NominalPower=300;
Utilization = 0.8; %global fuel utilization
Steam2Carbon = 2.5;

Air.N2 = .79;
Air.O2 = .21;

Fuel.CH4 = 0.9;
Fuel.CO = 0.04;
Fuel.CO2 = 0.04;
Fuel.H2 = 0;
Fuel.H2O = 0;
Fuel.N2 = 0.02;

%some initial guesses for other streams (don't have to be accurate)
FuelOut.CH4 = 0;
FuelOut.CO = 0.4*(Fuel.CH4 + Fuel.CO + Fuel.CO2)/(1 + ((Fuel.CH4 +.5*Fuel.CO)*Steam2Carbon -Fuel.H2O) + 2*Fuel.CH4);
FuelOut.CO2 = 0.6*(Fuel.CH4 + Fuel.CO + Fuel.CO2)/(1 + ((Fuel.CH4 +.5*Fuel.CO)*Steam2Carbon -Fuel.H2O) + 2*Fuel.CH4);
FuelOut.H2 = (1-Utilization)*(4*Fuel.CH4 + Fuel.CO)/(1 + ((Fuel.CH4 +.5*Fuel.CO)*Steam2Carbon -Fuel.H2O) + 2*Fuel.CH4);
FuelOut.N2 = Fuel.N2/(1 + ((Fuel.CH4 +.5*Fuel.CO)*Steam2Carbon -Fuel.H2O) + 2*Fuel.CH4);
FuelOut.H2O = 1 - (FuelOut.CO + FuelOut.CO2 + FuelOut.H2 + FuelOut.N2);

Oxidized.CH4 = 0;
Oxidized.CO = 0;
Oxidized.CO2 = FuelOut.CO2 + FuelOut.CO;
Oxidized.H2 = 0;
Oxidized.H2O = 1 - Oxidized.CO2 - FuelOut.N2;
Oxidized.N2 = FuelOut.N2;

if Fuel.H2O/(Fuel.CH4 +.5*Fuel.CO) <0.9*Steam2Carbon
    a = Steam2Carbon*(Fuel.CH4+.5*Fuel.CO)/(FuelOut.H2O -Steam2Carbon*(FuelOut.CH4 + .5*FuelOut.CO));
    FuelIn.CH4 = (a*FuelOut.CH4 + Fuel.CH4)/(a+1);
    FuelIn.CO = (a*FuelOut.CO + Fuel.CO)/(a+1);
    FuelIn.CO2 = (a*FuelOut.CO2 + Fuel.CO2)/(a+1);
    FuelIn.H2 = (a*FuelOut.H2 + Fuel.H2)/(a+1);
    FuelIn.H2O = (a*FuelOut.H2O + Fuel.H2O)/(a+1);
    FuelIn.N2 = (a*FuelOut.N2 + Fuel.N2)/(a+1);
else
    FuelIn = Fuel;
end

%Fuel Source
Plant.Components.FuelSource.type = 'Source';
Plant.Components.FuelSource.name = 'FuelSource';
Plant.Components.FuelSource.InitialComposition = Fuel;
Plant.Components.FuelSource.connections = {300;'';'Controller.FuelFlow';};

%Blower
Plant.Components.Blower.type = 'Blower';
Plant.Components.Blower.name = 'Blower';
Plant.Components.Blower.Map = 'RadialCompressor1'; % Loads a saved compressor map
Plant.Components.Blower.InitialComposition = Air;
Plant.Components.Blower.PeakEfficiency = 0.65;
Plant.Components.Blower.Tdesign = 300;%Design temp(K)
Plant.Components.Blower.RPMdesign = 30000;%Design RPM
Plant.Components.Blower.FlowDesign = 2;%design flow rate(kg/Sec)
Plant.Components.Blower.Pdesign = 1.25;%design pressure ratio
Plant.Components.Blower.connections = {'';'';'';'HX1.ColdPin';'Controller.Blower';}; %no species or temperature connection, stick to IC
Plant.Components.Blower.TagInf = {'Flow';'NRPM';'Power';'PR';'Nflow';'Temperature';'MassFlow';'Eff'};

%recircValve
Plant.Components.recircValve.type = 'Valve3Way';
Plant.Components.recircValve.name = 'recircValve';
Plant.Components.recircValve.InitialFlowIn = FuelOut;
Plant.Components.recircValve.InitialFlowIn.T = 1050;
Plant.Components.recircValve.connections = {'FC1.Flow1Out','Controller.AnodeRecirc'};
Plant.Components.recircValve.PercOpen = 0.5; %INITIAL valve position

%Mixing
Plant.Components.Mix1.type = 'MixingVolume';
Plant.Components.Mix1.name = 'Mix1';
Plant.Components.Mix1.Vol = 0.1;
Plant.Components.Mix1.inlets = 2;
Plant.Components.Mix1.SpeciesInit = FuelIn;
Plant.Components.Mix1.Tinit = 800;
switch Reformer
    case {'internal';'direct'}
        Plant.Components.Mix1.connections = {'FuelSource.Outlet';'recircValve.Out1';'FC1.Flow1Pin'};
    case {'external';'adiabatic';'pox'}
        Plant.Components.Mix1.connections = {'FuelSource.Outlet';'recircValve.Out1';'Reformer.ReformedPin'};
end
Plant.Components.Mix1.TagInf = {'MassFlow';'Temperature';};

%Oxidizer
Plant.Components.Oxidizer.type = 'Oxidizer';
Plant.Components.Oxidizer.name = 'Oxidizer';
Plant.Components.Oxidizer.inlets = 2;
Plant.Components.Oxidizer.InitialFlowOut = Oxidized;
Plant.Components.Oxidizer.InitialFlowOut.T = 1050;
Plant.Components.Oxidizer.Vol = .1; %volume in m^3
switch Reformer
    case {'internal';'direct'}
        Plant.Components.Oxidizer.connections = {'recircValve.Out2';'FC1.Flow2Out';'HX1.HotPin';};
    case {'external';'adiabatic';'pox'}
        Plant.Components.Oxidizer.connections = {'recircValve.Out2';'FC1.Flow2Out';'Reformer.CooledPin';};
end
Plant.Components.Oxidizer.TagInf = {'EquivelanceRatio';'Temperature';'MassFlow'};

%Reformer
switch Reformer
    case {'internal';'direct'}
        %none
    case 'external'
        Plant.Components.Reformer.type = 'Reformer';
        Plant.Components.Reformer.name = 'Reformer';
        Plant.Components.Reformer.direction = 'counterflow';  % 'co-flow', or 'counterflow', neccessary for external reformer with heat recovery stream
        Plant.Components.Reformer.nodes = 5; % # of nodes to simulate
        Plant.Components.Reformer.ReformTarget = .75; %    If there is a 2nd stream providing heat, it can use one of the following 4 targets to initialize reformer: Reforming percent, cold outlet temp, hot outlet temp, effectiveness of heat transfer from hot stream. Note that due to the reforming chemistry, a simple heat exchanger effectivess is not applicable.
        Plant.Components.Reformer.method = 'RefPerc'; %method for sizing reformer to initial conditions. Options are: 'none' fixed size during intialization, 'RefPerc' sizes for a specific external reforming percent, 'ColdT' sizes to match cold exit temp, 'HotT' sizes to match hot ext temp, 'Effectiveness' sizes to match a target effectiveness: % of ideal (infinite area) heat transfer (with no conduction between nodes)
        Plant.Components.Reformer.InletGuess = FuelIn;
        Plant.Components.Reformer.connections = {'Mix1.Outlet';'FC1.Flow1Pin';'Oxidizer.Flow';'HX1.HotPin';};
        Plant.Components.Reformer.TagInf = {'H2flow';'CH4';};
    case 'adiabatic'
        Plant.Components.Reformer.type = 'Reformer';
        Plant.Components.Reformer.name = 'Reformer';
        Plant.Components.Reformer.nodes = 1; % adiabatic reformer only needs 1 node, it is an equilibrium reactor
        Plant.Components.Reformer.InletGuess = FuelIn;
        Plant.Components.Reformer.connections = {'Mix1.Outlet';'FC1.Flow1Pin';};
        Plant.Components.Reformer.TagInf = {'H2flow';'CH4';};
    case 'pox'
        %have not made reformer handle pox yet
end

% Air Bypass
Plant.Components.bypassValve.type = 'Valve3Way';
Plant.Components.bypassValve.name = 'bypassValve';
Plant.Components.bypassValve.InitialFlowIn = Air;
Plant.Components.bypassValve.InitialFlowIn.T = 300;
Plant.Components.bypassValve.connections = {'Blower.Outlet','Controller.HeaterBypass'};
Plant.Components.bypassValve.PercOpen = 0.08; %INITIAL valve position

% Air Heat Exchanger
Plant.Components.HX1.type = 'HeatExchanger';
Plant.Components.HX1.name = 'HX1';
Plant.Components.HX1.direction = 'counterflow'; % 'coflow', or 'counterflow' or 'crossflow'
Plant.Components.HX1.columns = 10;
Plant.Components.HX1.rows =1;
switch Reformer
    case {'internal';'direct';'adiabatic';'pox'}
        Plant.Components.HX1.sizemethod = 'ColdT'; %method for sizing HX to initial conditions. Options are: 'fixed' fixed size heat exchanger during intialization, 'ColdT' sizes to match cold exit temp, 'HotT' sizes to match hot ext temp, 'Effectiveness' sizes to match a target effectiveness: % of ideal (infinite area) heat transfer (with no conduction between nodes)
        Plant.Components.HX1.Target = 950; %can be numeric or a string of block.property, ex 'Controller.HeaterTarget'. If it can't reach the temperature target it defaults to 98% effectiveness
    case 'external'
        Plant.Components.HX1.sizemethod = 'ColdT'; %method for sizing HX to initial conditions. Options are: 'fixed' fixed size heat exchanger during intialization, 'ColdT' sizes to match cold exit temp, 'HotT' sizes to match hot ext temp, 'Effectiveness' sizes to match a target effectiveness: % of ideal (infinite area) heat transfer (with no conduction between nodes)
        Plant.Components.HX1.Target = 950; %can be numeric or a string of block.property, ex 'Controller.HeaterTarget'. If it can't reach the temperature target it defaults to 98% effectiveness
end
Plant.Components.HX1.Mass = 25; %mass in kg
Plant.Components.HX1.Vol = 1; % volume in m^3
Plant.Components.HX1.Cold_T_init = 300; %initial guess temperaure of cold inlet
Plant.Components.HX1.Hot_T_init = 1000; %initial guess temperaure of hot inlet
Plant.Components.HX1.ColdSpecIn = Air;
Plant.Components.HX1.HotSpecIn = Oxidized;
switch Reformer
    case {'internal';'direct';'adiabatic';'pox'}
        Plant.Components.HX1.connections = {'bypassValve.Out2';'Oxidizer.Flow';'Mix2.Pin';'';}; 
    case 'external'
        Plant.Components.HX1.connections = {'bypassValve.Out2';'Reformer.Cooled';'Mix2.Pin';'';}; 
end
Plant.Components.HX1.TagInf = {'ColdOut';'HotOut';'Effectiveness';'NetImbalance'};

% Air Bypass mixing
Plant.Components.Mix2.type = 'MixingVolume';
Plant.Components.Mix2.name = 'Mix2';
Plant.Components.Mix2.Vol = 0.1;
Plant.Components.Mix2.inlets = 2;
Plant.Components.Mix2.SpeciesInit= Air;
Plant.Components.Mix2.Tinit = 970;
Plant.Components.Mix2.connections = {'bypassValve.Out1';'HX1.ColdOut';'FC1.Flow2Pin'};
Plant.Components.Mix2.TagInf = {'MassFlow';'Temperature';};

% fuel cell
Plant.Components.FC1.type = 'FuelCell';
Plant.Components.FC1.name = 'FC1';
Plant.Components.FC1.FCtype = 'SOFC'; %SOFC, or MCFC or oxySOFC
Plant.Components.FC1.Reformer = Reformer; 
Plant.Components.FC1.direction = 'counterflow'; % 'coflow', or 'counterflow' or 'crossflow'
Plant.Components.FC1.ClosedCathode = 0; %0 means air or some excess flow of O2 in the cathode used as primary means of temerature control (initializations hold to design fuel utilization), 1 means closed end cathode, or simply a fixed oxygen utilization, cooling is done with excess fuel, and the design voltage is met during initialization
Plant.Components.FC1.CoolingStream = 'cathode'; % choices are 'anode' or 'cathode'. Determines which flow is increased to reach desired temperature gradient.
Plant.Components.FC1.PressureRatio = 1.2;
Plant.Components.FC1.columns = 4;
Plant.Components.FC1.rows = 4;
Plant.Components.FC1.RatedStack_kW = 300; %Nominal Stack Power in kW
Plant.Components.FC1.Fuel = Fuel; %raw fuel fed to the system (calculates anode recirculation from this to reach S2C)
Plant.Components.FC1.Flow2Spec = Air; %initial oxidant composition
Plant.Components.FC1.Steam2Carbon = Steam2Carbon;
Plant.Components.FC1.method = 'Achenbach'; %Determines reforming reaction kinetics options: 'Achenbach' , 'Leinfelder' , 'Drescher'   
Plant.Components.FC1.L_Cell= .09;  %Cell length in meters
Plant.Components.FC1.W_Cell = .09;  %Cell Width in meters  
Plant.Components.FC1.Specification = 'power density'; %options are 'cells', 'power density', 'voltage', or 'current density'. Note: careful when specifying cells that it arrives at a feasible power density
Plant.Components.FC1.SpecificationValue = 450; % power density specified in mW/cm^2, voltage specified in V/cell, current density specified in A/cm^2
Plant.Components.FC1.deltaTStack = 70; %temperature difference from cathode inlet to cathode outlet
Plant.Components.FC1.TpenAvg = 1023;% 750 C, average electrolyte operating temperature
Plant.Components.FC1.FuelUtilization = Utilization; %fuel utilization (net hydrogen consumed/ maximum hydrogen produced with 100% Co and CH4 conversion
Plant.Components.FC1.Flow1Pdrop = 2; %design anode pressure drop
Plant.Components.FC1.Flow2Pdrop = 10; %Design cathode pressure drop
Plant.Components.FC1.Map = 'SOFC_map'; %Radiative heat transfer view factors, imported from CAD
switch Reformer
    case 'direct'
        Plant.Components.FC1.connections = {'Controller.Current';'Mix1.Outlet';'Mix2.Outlet';'Oxidizer.Pin';'Oxidizer.Pin';};
    case 'internal'
        Plant.Components.FC1.connections = {'Controller.Current';'Mix1.Outlet';'Mix2.Outlet';'Oxidizer.Pin';'Oxidizer.Pin';};
        Plant.Components.FC1.RefPerc = 0.8;% necessary for internal reformer, proportion of CH4 reforming in the reforming channels
        Plant.Components.FC1.RefSpacing = 1;% necessary for internal reformer. This is the # of active cells between reformer plates
    case 'external';
        Plant.Components.FC1.connections = {'Controller.Current';'Reformer.Reformed';'Mix2.Outlet';'Oxidizer.Pin';'Oxidizer.Pin';};
        Plant.Components.FC1.RefPerc = Plant.Components.Reformer.ReformTarget; % necessary for internal or external reformer. This is the percent towards equilibrium occurin in the reformer plaes
    case 'adiabatic'
        Plant.Components.FC1.connections = {'Controller.Current';'Reformer.Reformed';'Mix2.Outlet';'Oxidizer.Pin';'Oxidizer.Pin';};
    case 'pox'
end
Plant.Components.FC1.TagInf = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';'TcathOut';'ASR';'O2utilization';}; %Tags to record at each step
Plant.Components.FC1.TagFinal = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';'TcathOut';'ASR';'O2utilization'}; %Tags to record at the final step

% Controller %% order is heater bypass, blower power, anode recirculation
Plant.Controls.Controller.type = 'ControlSOFCsystem';
Plant.Controls.Controller.name = 'Controller';
Plant.Controls.Controller.Target = {'FC1.TpenAvg';'FC1.deltaTStack';'FC1.Steam2Carbon';};
Plant.Controls.Controller.Fuel = Fuel;
Plant.Controls.Controller.Cells = 'FC1.Cells';
Plant.Controls.Controller.Utilization = 'FC1.FuelUtilization';
Plant.Controls.Controller.InitialAnodeRecirc = 'FC1.Recirc.Anode';
Plant.Controls.Controller.InitConditions = {'bypassValve.PercOpen';'Blower.NominalPower';'FC1.Current';}; %heater bypass, blower power, FC current  [note: fuel flow and recirculaion are calculated and are not states]
% Plant.Controls.Controller.Gain = [1e-1;1e-3;1e-1];
% Plant.Controls.Controller.PropGain = [.75;4;1];
Plant.Controls.Controller.Gain = [1e-1;1e-3;];
Plant.Controls.Controller.PropGain = [.75;4;];
Plant.Controls.Controller.TagInf = {'Bypass';'Blower';'Recirculation';'FuelFlow';'Current';'Utilization'};
Plant.Controls.Controller.connections = {'FC1.MeasureTflow2';'Mix2.Temperature';'FC1.MeasureVoltage';};

Plant.Scope = {'FC1.PENavgT';'FC1.TcathOut';'Mix2.Temperature';'Controller.Bypass';'Controller.Utilization';'FC1.StackdeltaT';'Controller.Blower';'Blower.MassFlow';'Blower.NRPM';}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'FC1.Voltage';'Controller.Current';'Oxidizer.Temperature';}];