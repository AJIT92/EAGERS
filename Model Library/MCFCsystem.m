%% MCFC with indirect internal reformer
function Plant = MCFCsystem
Plant.NominalPower=300;
Steam2Carbon = 2.0;
Humidify = 1; % 0 = pre-humidify fuel, 1 = humidify with anode reciculation
Air.N2 = .79;
Air.O2 = .21;

NaturalGas.CH4 = 0.9;
NaturalGas.CO = 0.04;
NaturalGas.CO2 = 0.04;
NaturalGas.H2 = 0;
NaturalGas.H2O = 0;
NaturalGas.N2 = 0.02;
if Humidify==1
    Fuel = NaturalGas;
else
    N = 1 + Steam2Carbon*(NaturalGas.CH4+.5*NaturalGas.CO);
    speciesNames = fieldnames(NaturalGas);
    for i = 1:1:length(speciesNames)
        Fuel.(speciesNames(i)) = NaturalGas.(speciesNames(i))/N;
    end
    Fuel.H2O = (N-1)/N;
end

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
Plant.Components.Blower.FlowDesign = 1.5;%design flow rate(kg/Sec)
Plant.Components.Blower.Pdesign = 1.3;%design pressure ratio
Plant.Components.Blower.connections = {'';'';'';'HX1.ColdPin';'Controller.Blower';}; %no species or temperature connection, stick to IC
Plant.Components.Blower.TagInf = {'Flow';'NRPM';'Power';'PR';'Nflow';'Temperature';'MassFlow';'Eff'};

Plant.Components.FC1.type = 'FuelCell';
Plant.Components.FC1.name = 'FC1';
Plant.Components.FC1.FCtype = 'MCFC'; %SOFC, or MCFC or oxySOFC
Plant.Components.FC1.Reformer = 'internal'; % options are 'internal' for indirect internal reforming, 'direct' for direct internal reforming, 'adiabatic' for external reforming using the heat from the anode (over-rides S2C ratio), 'pox' partial oxidation
% Plant.Components.FC1.RefPerc = .75; % necessary for internal reformer. This is the percent towards equilibrium occurin in the reformer plaes
% Plant.Components.FC1.RefSpacing = 1;% necessary for internal reformer. This is the # of active cells between reformer plates
Plant.Components.FC1.direction = 'coflow'; % 'co-flow', or 'counterflow' or 'crossflow'
Plant.Components.FC1.ClosedCathode = 0; %0 means air or some excess flow of O2 in the cathode used as primary means of temerature control (initializations hold to design fuel utilization), 1 means closed end cathode, or simply a fixed oxygen utilization, cooling is done with excess fuel, and the design voltage is met during initialization
Plant.Components.FC1.CoolingStream = 'cathode'; % choices are 'anode' or 'cathode'. Determines which flow is increased to reach desired temperature gradient.
Plant.Components.FC1.PressureRatio = 1.2;
Plant.Components.FC1.columns = 5;
Plant.Components.FC1.rows = 1;
Plant.Components.FC1.RatedStack_kW = 300; %Nominal Stack Power in kW
Plant.Components.FC1.Fuel = Fuel; %initial fuel composition at inlet
Plant.Components.FC1.Oxidant = Air; %initial oxidant composition (with exhaust for CO2)
Plant.Components.FC1.Humidify = Humidify; % 0 = pre-humidify fuel, 1 = humidify with anode reciculation
Plant.Components.FC1.Steam2Carbon = Steam2Carbon;
Plant.Components.FC1.method = 'Achenbach'; %Determines reforming reaction kinetics options: 'Achenbach' , 'Leinfelder' , 'Drescher'   
Plant.Components.FC1.L_Cell= 1.1;  %Cell length in meters
Plant.Components.FC1.W_Cell = .8;  %Cell Width in meters  

Plant.Components.FC1.DesignTarget = 'power density'; %options are 'power density', 'voltage', or 'current density'
Plant.Components.FC1.DesignTargetValue = 120; % power density specified in mW/cm^2, voltage specified in V/cell, current density specified in A/cm^2
Plant.Components.FC1.Cells = ceil(Plant.Components.FC1.RatedStack_kW*100/(Plant.Components.FC1.L_Cell*Plant.Components.FC1.W_Cell*Plant.Components.FC1.DesignTargetValue)); %# of cells in stack
Plant.Components.FC1.deltaTStack = 50; %temperature difference from cathode inlet to cathode outlet
Plant.Components.FC1.TpenAvg = 923;% 650 C, average electrolyte operating temperature
Plant.Components.FC1.FuelUtilization =.75; %fuel utilization (net hydrogen consumed/ maximum hydrogen produced with 100% Co and CH4 conversion
Plant.Components.FC1.RefSpacing = 10;% # of active cells between reformer plates
Plant.Components.FC1.AnodePdrop = 2; %design anode pressure drop
Plant.Components.FC1.CathodePdrop = 10; %Design cathode pressure drop
Plant.Components.FC1.Map = 'SOFC_map'; %Radiative heat transfer view factors, imported from CAD
Plant.Components.FC1.connections = {'Controller.Current';'Mix2.Outlet';'Reformer.Reformed';'Oxidizer.Pin';'Oxidizer.Pin';};
Plant.Components.FC1.TagInf = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';'TcathOut';'ASR';'O2utilization';}; %Tags to record at each step
Plant.Components.FC1.TagFinal = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';'TcathOut';'ASR';'O2utilization'}; %Tags to record at the final step

Plant.Components.Valve1.type = 'Valve3Way'; %Split stream to anode oxidizer
Plant.Components.Valve1.name = 'Valve1';
Plant.Components.Valve1.valveIC = 0.92;
Plant.Components.Valve1.connections = {'AirSource.TXN';'Controller.ValvePerc';};

Plant.Components.Oxidizer.type = 'Oxidizer';
Plant.Components.Oxidizer.name = 'Oxidizer';
Plant.Components.Oxidizer.InitialComposition = Fuel;
Plant.Components.Oxidizer.connections = {'Valve1.TXN2';'FC1.AnodeOut';'FC1.CathodePressureIn';};
Plant.Components.Oxidizer.Vol = 1; %volume of oxidizer
Plant.Components.Oxidizer.Pdrop = 1; %kPa drop across oxidizer

Plant.Components.AirInletMixing.type = 'MixingVolume'; %combine streams
Plant.Components.AirInletMixing.name = 'AirInletMixing';
Plant.Components.AirInletMixing.Vol = 0.1;
Plant.Components.AirInletMixing.P = 120;
Plant.Components.AirInletMixing.inlets = 2;
Plant.Components.AirInletMixing.connections = {'Valve1.TXN1';'Oxidizer.TXN';};

Plant.Controls.Controller.type = 'ControlMCFCstack';
Plant.Controls.Controller.name = 'Controller';
Plant.Controls.Controller.Target = {'FC1.TpenAvg';'FC1.deltaTStack';0.5};
Plant.Controls.Controller.Fuel = 'FC1.Fuel';
Plant.Controls.Controller.Cells = 'FC1.Cells';
Plant.Controls.Controller.Utilization = 'FC1.FuelUtilization';
Plant.Controls.Controller.DesignTemp = 700;
Plant.Controls.Controller.DesignFlow = 'FC1.AirFlow';
Plant.Controls.Controller.valveIC = 0.8;
if Humidify==1
    Plant.Controls.Controller.Gain = [1e-3;1e-4;1e-4;];
    Plant.Controls.Controller.PropGain = [.5;.1;0];
else
    Plant.Controls.Controller.Gain = [1e-3;1e-4;0];
    Plant.Controls.Controller.PropGain = [5;.1;0];
end
Plant.Controls.Controller.TagInf = {'AirFlow';'Toxidant';'Tfuel';'FuelFlow';'ValvePerc';};
Plant.Controls.Controller.connections = {'FC1.MeasureTcathOut';'AirInletMixing.Temperature';'FC1.MeasureTpen';'FC1.MeasureCurrent';'Oxidizer.EquivelanceRatio';};

Plant.Scope = {'Controller.Blower';'Controller.Bypass';'Controller.Current'}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'FC1.StackdeltaT';'Controller.AirFlow';'FC1.PENavgT';'Controller.Toxidant';'FC1.Voltage';'FC1.TcathOut'}];
