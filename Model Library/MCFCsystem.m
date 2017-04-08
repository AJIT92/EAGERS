%% MCFC with indirect internal reformer (is there anode recirculation?)
function Plant = MCFCsystem
Plant.NominalPower=300;
Steam2Carbon = 2.0;

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
Plant.Components.Blower.FlowDesign = 1.5;%design flow rate(kg/Sec)
Plant.Components.Blower.Pdesign = 1.3;%design pressure ratio
Plant.Components.Blower.connections = {'';'';'';'HX1.ColdPin';'Controller.Blower';}; %no species or temperature connection, stick to IC
Plant.Components.Blower.TagInf = {'Flow';'NRPM';'Power';'PR';'Nflow';'Temperature';'MassFlow';'Eff'};

Plant.Components.Valve1.type = 'Valve3Way'; %Split stream to anode oxidizer
Plant.Components.Valve1.name = 'Valve1';
Plant.Components.Valve1.valveIC = 0.92;
Plant.Components.Valve1.connections = {'Blower.Outlet';'Controller.HeaterBypass'};

Plant.Components.Oxidizer.type = 'Oxidizer';
Plant.Components.Oxidizer.name = 'Oxidizer';
Plant.Components.Oxidizer.InitialComposition = Oxidized;
Plant.Components.Oxidizer.connections = {'Valve1.Out2';'FC1.AnodeOut';'FC1.CathodePressureIn';};
Plant.Components.Oxidizer.Vol = 1; %volume of oxidizer
Plant.Components.Oxidizer.Pdrop = 1; %kPa drop across oxidizer

Plant.Components.Mix2.type = 'MixingVolume'; %combine streams
Plant.Components.Mix2.name = 'Mix2';
Plant.Components.Mix2.Vol = 0.1;
Plant.Components.Mix2.P = 120;
Plant.Components.Mix2.inlets = 2;
Plant.Components.Mix2.connections = {'Valve1.Out1';'Oxidizer.Flow';};

Plant.Components.FC1.type = 'FuelCell';
Plant.Components.FC1.name = 'FC1';
Plant.Components.FC1.FCtype = 'MCFC'; %SOFC, or MCFC or oxySOFC
Plant.Components.FC1.Reformer = 'internal'; % options are 'internal' for indirect internal reforming, 'direct' for direct internal reforming, 'adiabatic' for external reforming using the heat from the anode (over-rides S2C ratio), 'pox' partial oxidation
Plant.Components.FC1.RefPerc = .75; % necessary for internal reformer. This is the percent towards equilibrium occurin in the reformer plaes
Plant.Components.FC1.RefSpacing = 10;% necessary for internal reformer. This is the # of active cells between reformer plates
Plant.Components.FC1.direction = 'coflow'; % 'co-flow', or 'counterflow' or 'crossflow'
Plant.Components.FC1.ClosedCathode = 0; %0 means air or some excess flow of O2 in the cathode used as primary means of temerature control (initializations hold to design fuel utilization), 1 means closed end cathode, or simply a fixed oxygen utilization, cooling is done with excess fuel, and the design voltage is met during initialization
Plant.Components.FC1.CoolingStream = 'cathode'; % choices are 'anode' or 'cathode'. Determines which flow is increased to reach desired temperature gradient.
Plant.Components.FC1.PressureRatio = 1.2;
Plant.Components.FC1.columns = 5;
Plant.Components.FC1.rows = 1;
Plant.Components.FC1.RatedStack_kW = 300; %Nominal Stack Power in kW
Plant.Components.FC1.Fuel = Fuel; %initial fuel composition at inlet
Plant.Components.FC1.Flow2Spec = Air; %initial oxidant composition (with exhaust for CO2)
Plant.Components.FC1.Steam2Carbon = Steam2Carbon;
Plant.Components.FC1.method = 'Achenbach'; %Determines reforming reaction kinetics options: 'Achenbach' , 'Leinfelder' , 'Drescher'   
Plant.Components.FC1.L_Cell= 1.1;  %Cell length in meters
Plant.Components.FC1.W_Cell = .8;  %Cell Width in meters  
Plant.Components.FC1.Specification = 'power density'; %options are 'power density', 'voltage', or 'current density'
Plant.Components.FC1.SpecificationValue = 120; % power density specified in mW/cm^2, voltage specified in V/cell, current density specified in A/cm^2
Plant.Components.FC1.deltaTStack = 50; %temperature difference from cathode inlet to cathode outlet
Plant.Components.FC1.TpenAvg = 923;% 650 C, average electrolyte operating temperature
Plant.Components.FC1.FuelUtilization =.75; %fuel utilization (net hydrogen consumed/ maximum hydrogen produced with 100% Co and CH4 conversion
Plant.Components.FC1.Flow1Pdrop = 2; %design anode pressure drop
Plant.Components.FC1.Flow2Pdrop = 10; %Design cathode pressure drop
Plant.Components.FC1.Map = 'MCFC_map'; %Radiative heat transfer view factors, imported from CAD
Plant.Components.FC1.connections = {'Controller.Current';'Mix2.Outlet';'Reformer.Reformed';'Oxidizer.Pin';'Oxidizer.Pin';};
Plant.Components.FC1.TagInf = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';'TcathOut';'ASR';'O2utilization';}; %Tags to record at each step
Plant.Components.FC1.TagFinal = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';'TcathOut';'ASR';'O2utilization'}; %Tags to record at the final step


Plant.Controls.Controller.type = 'ControlSOFCsystem';
Plant.Controls.Controller.name = 'Controller';
Plant.Controls.Controller.Target = {'FC1.TpenAvg';'FC1.deltaTStack';'FC1.Steam2Carbon';};
Plant.Controls.Controller.Fuel = Fuel;
Plant.Controls.Controller.Cells = 'FC1.Cells';
Plant.Controls.Controller.Utilization = 'FC1.FuelUtilization';
Plant.Controls.Controller.InitialAnodeRecirc = 'FC1.Recirc.Anode';
Plant.Controls.Controller.InitConditions = {'bypassValve.PercOpen';'Blower.NominalPower';'FC1.Current';}; %heater bypass, blower power, FC current  [note: fuel flow and recirculaion are calculated and are not states]
Plant.Controls.Controller.Gain = [1e-1;1e-3;];
Plant.Controls.Controller.PropGain = [.75;4;];
Plant.Controls.Controller.TagInf = {'Bypass';'Blower';'Recirculation';'FuelFlow';'Current';'Utilization'};
Plant.Controls.Controller.TagInf = {'AirFlow';'Toxidant';'Tfuel';'FuelFlow';'ValvePerc';};
Plant.Controls.Controller.connections = {'FC1.MeasureTflow2';'Mix2.Temperature';'FC1.MeasureVoltage'};

Plant.Scope = {'Controller.Blower';'Controller.Bypass';'Controller.Current'}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'FC1.StackdeltaT';'Controller.AirFlow';'FC1.PENavgT';'Controller.Toxidant';'FC1.Voltage';'FC1.TcathOut'}];