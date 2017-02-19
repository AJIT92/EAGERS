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
Plant.Components.EC1.rows = 1;
Plant.Components.EC1.RatedStack_kW = 300; %Nominal Stack Power in kW
Plant.Components.EC1.Steam = Steam; %initial fuel composition at inlet
Plant.Components.EC1.Oxidant = Air; %initial oxidant composition if there is a dilution on the anode side added to the O2 production (temperature regualtion)
Plant.Components.EC1.L_Cell= .1;  %Cell length in meters
Plant.Components.EC1.W_Cell = .1;  %Cell Width in meters  
if Plant.Components.EC1.ClosedCathode
    Plant.Components.EC1.Cells = ceil(Plant.Components.EC1.RatedStack_kW*1000/(1.3*1e4*Plant.Components.EC1.L_Cell*Plant.Components.EC1.W_Cell)); %# of cells in stack (assumes 1 A/cm^2) corrected later
else
    Plant.Components.EC1.DesignTarget = 'power density'; %options are 'power density', 'voltage', or 'current density' (A/cm^2)
    Plant.Components.EC1.DesignTargetValue = 2000; % power density specified in mW/cm^2, voltage specified in V/cell, current density specified in A/cm^2
    if strcmp(Plant.Components.EC1.DesignTarget,'power density')
        Plant.Components.EC1.Cells = ceil(Plant.Components.EC1.RatedStack_kW*100/(Plant.Components.EC1.L_Cell*Plant.Components.EC1.W_Cell*Plant.Components.EC1.DesignTargetValue)); %# of cells in stack
    elseif strcmp(Plant.Components.EC1.DesignTarget,'current density')
        Plant.Components.EC1.Cells = ceil(Plant.Components.EC1.RatedStack_kW*1000/(1.3*1e4*Plant.Components.EC1.L_Cell*Plant.Components.EC1.W_Cell*Plant.Components.EC1.DesignTargetValue)); %# of cells in stack (assumes voltage of 1.3)
    elseif strcmp(Plant.Components.EC1.DesignTarget,'voltage')
        Plant.Components.EC1.Cells = ceil(Plant.Components.EC1.RatedStack_kW*1000/(Plant.Components.EC1.DesignTargetValue*1e4*Plant.Components.EC1.L_Cell*Plant.Components.EC1.W_Cell)); %# of cells in stack (assumes 1 A/cm^2) corrected later
    end 
end
Plant.Components.EC1.deltaTStack = 50; %temperature difference from cathode inlet to cathode outlet
Plant.Components.EC1.TpenAvg = 1023;% 750 C, average electrolyte operating temperature
Plant.Components.EC1.H2O_Utilization =.85; %H2O utilization (net steam consumed/ steam supply)
Plant.Components.EC1.AnodePdrop = 2; %design anode pressure drop
Plant.Components.EC1.CathodePdrop = 10; %Design cathode pressure drop
Plant.Components.EC1.connections = {'Controller.Current';'SteamSource.Outlet';'AirSource.Outlet';'';'';};
Plant.Components.EC1.TagInf = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2Outilization';'TcathOut';'nCurrent';'nVoltage'}; %Tags to record at each step
Plant.Components.EC1.TagFinal = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2Outilization';'TcathOut';'nCurrent';'nVoltage'}; %Tags to record at the final step

Plant.Controls.Controller.type = 'ControlECstack';
Plant.Controls.Controller.name = 'Controller';
Plant.Controls.Controller.Target = {'EC1.TpenAvg';'EC1.deltaTStack';};
Plant.Controls.Controller.Steam = 'EC1.Steam';
Plant.Controls.Controller.Cells = 'EC1.Cells';
Plant.Controls.Controller.Utilization = 'EC1.H2O_Utilization';
Plant.Controls.Controller.SteamTemperature = 973;
Plant.Controls.Controller.InitConditions = {'EC1.AnodeIn.IC';'EC1.Current';}; %oxidant flow rate, net current
Plant.Controls.Controller.Gain = [1e-3;1e-2];
Plant.Controls.Controller.PropGain = [0;1];
Plant.Controls.Controller.TagInf = {'OxidantFlow';'OxidantTemp';'Tsteam';'SteamFlow';'Current';};
Plant.Controls.Controller.connections = {'EC1.MeasureTanodeOut';'Controller.OxidantTemp';'EC1.MeasureTpen';'EC1.MeasureVoltage';};

Plant.Scope = {'Controller.OxidantFlow';'Controller.Current';'Controller.OxidantTemp';'EC1.Voltage';}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'EC1.StackdeltaT';'EC1.PENavgT';'EC1.TcathOut';}];