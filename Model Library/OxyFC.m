%% SOFC stack with basic control of fuel and air flow
function Plant = OxyFC
Plant.NominalPower=300;%power in kW

Steam2Carbon = 2.0;
Oxidant.O2 = 1;
% Oxidant.O2 = .333; %MCFC
% Oxidant.CO2 = .667;

Fuel.CH4 = 0.9;
Fuel.CO = 0.04;
Fuel.CO2 = 0.04;
Fuel.H2 = 0;
Fuel.H2O = 0;
Fuel.N2 = 0.02;

%% Components
Plant.Components.O2Source.type = 'Source'; 
Plant.Components.O2Source.name = 'O2Source';
Plant.Components.O2Source.InitialComposition = Oxidant;
Plant.Components.O2Source.connections = {930;'';'Controller.OxidantFlow';};

Plant.Components.FuelSource.type = 'Source';
Plant.Components.FuelSource.name = 'FuelSource';
Plant.Components.FuelSource.InitialComposition = Fuel;
Plant.Components.FuelSource.connections = {300;'';'Controller.FuelFlow';};

S2C = Fuel.H2O/(Fuel.CH4+.5*Fuel.CO);
if S2C<Steam2Carbon %add anode recirculation
    FCexit.CH4 = .0001;
    FCexit.CO = .0999;
    FCexit.CO2 = .25;
    FCexit.H2 = .1;
    FCexit.H2O = .54;
    FCexit.N2 = .01;

    PartiallyRef.CH4 = .07;
    PartiallyRef.CO = .07;
    PartiallyRef.CO2 = .07;
    PartiallyRef.H2 = .5;
    PartiallyRef.H2O = .28;
    PartiallyRef.N2 = .01;
    
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
    Plant.Components.Mix1.SpeciesInit = PartiallyRef;
    Plant.Components.Mix1.Tinit = 800;
    Plant.Components.Mix1.connections = {'FuelSource.Outlet';'recircValve.Out1';'FC1.AnodePressureIn'};
    Plant.Components.Mix1.TagInf = {'MassFlow';'Temperature';};
end

Plant.Components.FC1.type = 'FuelCell';
Plant.Components.FC1.name = 'FC1';
Plant.Components.FC1.FCtype = 'SOFC'; %SOFC, or MCFC 
Plant.Components.FC1.Reformer = 'internal'; % options are 'internal' for indirect internal reforming, 'direct' for direct internal reforming, 'adiabatic' for external reforming using the heat from the anode (over-rides S2C ratio), 'pox' partial oxidation
Plant.Components.FC1.RefPerc = .75; % necessary for internal reformer. This is the percent towards equilibrium occurin in the reformer plates
Plant.Components.FC1.RefSpacing = 1;% necessary for internal reformer. This is the # of active cells between reformer plates
Plant.Components.FC1.direction = 'coflow'; % 'coflow', or 'counterflow' or 'crossflow'
Plant.Components.FC1.ClosedCathode = 1; %0 means air or some excess flow of O2 in the cathode used as primary means of temerature control (initializations hold to design fuel utilization), 1 means closed end cathode, or simply a fixed oxygen utilization, cooling is done with excess fuel, and the design voltage is met during initialization
Plant.Components.FC1.CoolingStream = 'anode'; % choices are 'anode' or 'cathode'. Determines which flow is increased to reach desired temperature gradient.
Plant.Components.FC1.PressureRatio = 1.2;
Plant.Components.FC1.columns = 4;
Plant.Components.FC1.rows = 1;
Plant.Components.FC1.RatedStack_kW = 300; %Nominal Stack Power in kW
Plant.Components.FC1.Fuel = Fuel; %initial fuel composition at inlet
Plant.Components.FC1.Oxidant = Oxidant;
Plant.Components.FC1.Steam2Carbon = Steam2Carbon; %steam to carbon ratio that fuel or recirculaton is controlled to
Plant.Components.FC1.method = 'Achenbach'; %Determines reforming reaction kinetics options: 'Achenbach' , 'Leinfelder' , 'Drescher'   
Plant.Components.FC1.L_Cell= .09;  %Cell length in meters
Plant.Components.FC1.W_Cell = .09;  %Cell Width in meters  
Plant.Components.FC1.DesignTarget = 'voltage'; %options are 'power density', 'voltage', or 'current density'
Plant.Components.FC1.DesignTargetValue = 0.86; % power density specified in mW/cm^2, voltage specified in V/cell, current density specified in A/cm^2
Plant.Components.FC1.Cells = ceil(Plant.Components.FC1.RatedStack_kW*1000/(Plant.Components.FC1.L_Cell*Plant.Components.FC1.W_Cell*Plant.Components.FC1.DesignTargetValue*5000)); %# of cells in stack (assuming 0.5 A/cm^2 --> 5,000 A/m^2)
Plant.Components.FC1.deltaTStack = 50; %temperature difference from cathode inlet to cathode outlet
Plant.Components.FC1.TpenAvg = 1023;% 750 C, average electrolyte operating temperature
Plant.Components.FC1.FuelUtilization = Plant.Components.FC1.DesignTargetValue*1.6 - .728;% ; %fuel utilization (net hydrogen consumed/ maximum hydrogen produced with 100% Co and CH4 conversion (initial guess, will be iterated)
Plant.Components.FC1.AnodePdrop = 2; %design anode pressure drop
Plant.Components.FC1.CathodePdrop = 10; %Design cathode pressure drop
if S2C<Steam2Carbon %add anode recirculation
    Plant.Components.FC1.connections = {'Controller.Current';'O2Source.Outlet';'Mix1.Outlet';'';'';};
else
    Plant.Components.FC1.connections = {'Controller.Current';'O2Source.Outlet';'FuelSource.Outlet';'';'';};
end
Plant.Components.FC1.TagInf = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';'TcathOut';'LocalNernst';};
Plant.Components.FC1.TagFinal = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';};

Plant.Controls.Controller.type = 'ControlFCstack';
Plant.Controls.Controller.name = 'Controller';
Plant.Controls.Controller.Target = {'FC1.TpenAvg';'FC1.deltaTStack';'FC1.Steam2Carbon'};
Plant.Controls.Controller.OxyFC = 'yes';
Plant.Controls.Controller.OxidantUtilization = 'FC1.OxidantUtilization';
Plant.Controls.Controller.Fuel = 'FC1.Fuel';
Plant.Controls.Controller.Oxidant = 'FC1.Oxidant';
Plant.Controls.Controller.Cells = 'FC1.Cells';
Plant.Controls.Controller.InitConditions = {'FC1.Recirc.Anode';'FC1.AnodeIn.IC';'FC1.Current';}; %Recirculation, fuel flow rate, net current
Plant.Controls.Controller.Gain = [0;1e-1;1e0];
Plant.Controls.Controller.PropGain = [0;0;1];
Plant.Controls.Controller.TagInf = {'OxidantFlow';'OxidantTemp';'FuelFlow';'Current';'Recirculation';'Power';};
Plant.Controls.Controller.connections = {'FC1.MeasureTanodeOut';'Controller.OxidantTemp';'FC1.MeasureTpen';'FC1.MeasureVoltage';};

Plant.Scope = {'Controller.FuelFlow';'Controller.Current';'Controller.Recirculation';'Controller.Power';}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'FC1.StackdeltaT';'FC1.PENavgT';'FC1.Voltage';'FC1.TcathOut';'FC1.LocalNernst';}];