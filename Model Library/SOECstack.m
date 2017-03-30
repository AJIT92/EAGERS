%% SOEC stack with basic control of H2O and air flow
%if design voltage is <1.3, then heat must be supplied via the anode side
%if the design voltage is >1.3, then heat must be removed via the anode flow
function Plant = SOECstack
Plant.NominalPower=300;

Air.N2 = .79;
Air.O2 = .21;

Steam.H2O = 1;

%% Components
Plant.Components.AirSource.type = 'Source'; %fresh air 
Plant.Components.AirSource.name = 'AirSource';
Plant.Components.AirSource.InitialComposition = Air;
Plant.Components.AirSource.connections = {'Controller.OxidantTemp';'';'Controller.OxidantFlow';};

Plant.Components.SteamSource.type = 'Source';
Plant.Components.SteamSource.name = 'SteamSource';
Plant.Components.SteamSource.InitialComposition = Steam;
Plant.Components.SteamSource.connections = {'Controller.SteamTemp';'';'Controller.SteamFlow';};

Plant.Components.EC1.type = 'Electrolyzer';
Plant.Components.EC1.name = 'EC1';
Plant.Components.EC1.FCtype = 'SOEC'; %SOEC, or MCEC 
Plant.Components.EC1.Reformer = 'none'; % options are 'none' or 'methanator' for a seperate set of plates wher CO2 is injected to cause methanation
Plant.Components.EC1.direction = 'coflow'; % 'coflow', or 'counterflow' or 'crossflow'
Plant.Components.EC1.ClosedCathode = 0; %0 means air or some excess flow of O2 in the cathode used as primary means of temerature control (initializations hold to design fuel utilization), 1 means closed end cathode, or simply a fixed oxygen utilization, cooling is done with excess fuel, and the design voltage is met during initialization
Plant.Components.EC1.CoolingStream = 'none'; % choices are 'none' or 'cathode'. Determines which flow is increased to reach desired temperature gradient.
Plant.Components.EC1.PressureRatio = 1.2;
Plant.Components.EC1.columns = 5;
Plant.Components.EC1.rows = 5;
Plant.Components.EC1.RatedStack_kW = 300; %Nominal Stack Power in kW
Plant.Components.EC1.Flow1Spec = Steam; %initial fuel composition at inlet
Plant.Components.EC1.Flow2Spec = Air; %initial oxidant composition if there is a dilution on the anode side added to the O2 production (temperature regualtion)
Plant.Components.EC1.L_Cell= .1;  %Cell length in meters
Plant.Components.EC1.W_Cell = .1;  %Cell Width in meters  
if ~Plant.Components.EC1.ClosedCathode
    Plant.Components.EC1.Specification = 'power density';%options are 'cells', 'power density', 'voltage', or 'current density'. Note: careful when specifying cells that it arrives at a feasible power density
    Plant.Components.EC1.SpecificationValue = 2000; % power density specified in mW/cm^2, voltage specified in V/cell, current density specified in A/cm^2
end
Plant.Components.EC1.deltaTStack = 50; %temperature difference from cathode inlet to cathode outlet
Plant.Components.EC1.TpenAvg = 1023;% 750 C, average electrolyte operating temperature
Plant.Components.EC1.H2O_Utilization =.85; %H2O utilization (net steam consumed/ steam supply)
Plant.Components.EC1.Flow1Pdrop = 2; %design anode pressure drop
Plant.Components.EC1.Flow2Pdrop = 10; %Design cathode pressure drop
Plant.Components.EC1.connections = {'Controller.Current';'SteamSource.Outlet';'AirSource.Outlet';'';'';};
Plant.Components.EC1.TagInf = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2Outilization';'TcathOut';'nCurrent';'nVoltage'}; %Tags to record at each step
Plant.Components.EC1.TagFinal = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2Outilization';'TcathOut';'nCurrent';'nVoltage'}; %Tags to record at the final step

Plant.Controls.Controller.type = 'ControlECstack';
Plant.Controls.Controller.name = 'Controller';
Plant.Controls.Controller.Target = {'EC1.TpenAvg';'EC1.deltaTStack';};
Plant.Controls.Controller.Steam = 'EC1.Flow1Spec';
Plant.Controls.Controller.Cells = 'EC1.Cells';
Plant.Controls.Controller.Utilization = 'EC1.H2O_Utilization';
Plant.Controls.Controller.SteamTemperature = 973;
Plant.Controls.Controller.InitConditions = {'EC1.Flow2.IC';'EC1.Current';}; %oxidant flow rate, net current
Plant.Controls.Controller.Gain = [2e-4;];
Plant.Controls.Controller.PropGain = [3];
Plant.Controls.Controller.TagInf = {'OxidantFlow';'OxidantTemp';'Tsteam';'SteamFlow';'Current';};
Plant.Controls.Controller.connections = {'EC1.MeasureTflow2';'Controller.OxidantTemp';'EC1.MeasureVoltage';'EC1.MeasureTpen';};

Plant.Scope = {'Controller.OxidantFlow';'Controller.Current';'Controller.OxidantTemp';'EC1.Voltage';}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'EC1.StackdeltaT';'EC1.PENavgT';'EC1.TcathOut';}];