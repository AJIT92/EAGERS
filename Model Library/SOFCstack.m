%% SOFC stack with basic control of fuel and air flow
function Plant = SOFCstack
Reformer = 'external'; % options are 'internal' for indirect internal reforming, 'direct' for direct internal reforming, 'adiabatic' for external reforming using the heat from the anode (over-rides S2C ratio),'external' for an external reformer with heat captured from oxidixed anode exhaust, 'pox' partial oxidation reformer
Plant.NominalPower=300;
Utilization = 0.8; %global fuel utilization
Steam2Carbon = 3.0;
Air.N2 = .79;
Air.O2 = .21;

Fuel.CH4 = 0.9;
Fuel.CO = 0.04;
Fuel.CO2 = 0.04;
Fuel.H2 = 0;
Fuel.H2O = 0;
Fuel.N2 = 0.02;
S2C = Fuel.H2O/(Fuel.CH4+.5*Fuel.CO);

%some initial guesses for other streams (don't have to be accurate)
%these would need to be adjusted for a non CH4 fuel
if Fuel.H2O/(Fuel.CH4 +.5*Fuel.CO) <0.9*Steam2Carbon
    FuelOut.CH4 = 0;
    FuelOut.CO = 0.4*(Fuel.CH4 + Fuel.CO + Fuel.CO2)/(1 + ((Fuel.CH4 +.5*Fuel.CO)*Steam2Carbon -Fuel.H2O) + 2*Fuel.CH4);
    FuelOut.CO2 = 0.6*(Fuel.CH4 + Fuel.CO + Fuel.CO2)/(1 + ((Fuel.CH4 +.5*Fuel.CO)*Steam2Carbon -Fuel.H2O) + 2*Fuel.CH4);
    FuelOut.H2 = (1-Utilization)*(4*Fuel.CH4 + Fuel.CO)/(1 + ((Fuel.CH4 +.5*Fuel.CO)*Steam2Carbon -Fuel.H2O) + 2*Fuel.CH4);
    FuelOut.N2 = Fuel.N2/(1 + ((Fuel.CH4 +.5*Fuel.CO)*Steam2Carbon -Fuel.H2O) + 2*Fuel.CH4);
    FuelOut.H2O = 1 - (FuelOut.CO + FuelOut.CO2 + FuelOut.H2 + FuelOut.N2);

    a = Steam2Carbon*(Fuel.CH4+.5*Fuel.CO)/(FuelOut.H2O -Steam2Carbon*(FuelOut.CH4 + .5*FuelOut.CO));
    FuelIn.CH4 = (a*FuelOut.CH4 + Fuel.CH4)/(a+1);
    FuelIn.CO = (a*FuelOut.CO + Fuel.CO)/(a+1);
    FuelIn.CO2 = (a*FuelOut.CO2 + Fuel.CO2)/(a+1);
    FuelIn.H2 = (a*FuelOut.H2 + Fuel.H2)/(a+1);
    FuelIn.H2O = (a*FuelOut.H2O + Fuel.H2O)/(a+1);
    FuelIn.N2 = (a*FuelOut.N2 + Fuel.N2)/(a+1);
else
    FuelIn = Fuel;
    FuelOut.CH4 = 0;
    FuelOut.CO = 0.4*(Fuel.CH4 + Fuel.CO + Fuel.CO2)/(1 + 2*Fuel.CH4);
    FuelOut.CO2 = 0.6*(Fuel.CH4 + Fuel.CO + Fuel.CO2)/(1 + 2*Fuel.CH4);
    FuelOut.H2 = (1-Utilization)*(4*Fuel.CH4 + Fuel.CO)/(1 + 2*Fuel.CH4);
    FuelOut.N2 = Fuel.N2/(1 + 2*Fuel.CH4);
    FuelOut.H2O = 1 - (FuelOut.CO + FuelOut.CO2 + FuelOut.H2 + FuelOut.N2);
end
Oxidized.CH4 = 0;
Oxidized.CO = 0;
Oxidized.CO2 = FuelOut.CO2 + FuelOut.CO;
Oxidized.H2 = 0;
Oxidized.H2O = 1 - Oxidized.CO2 - FuelOut.N2;
Oxidized.N2 = FuelOut.N2;
%% Components
Plant.Components.AirSource.type = 'Source'; %fresh air 
Plant.Components.AirSource.name = 'AirSource';
Plant.Components.AirSource.InitialComposition = Air;
Plant.Components.AirSource.connections = {'Controller.OxidantTemp';'';'Controller.OxidantFlow';};

Plant.Components.FuelSource.type = 'Source';
Plant.Components.FuelSource.name = 'FuelSource';
Plant.Components.FuelSource.InitialComposition = Fuel;
if S2C<Steam2Carbon %there will be anode recirculation to heat incoming fuel
    Plant.Components.FuelSource.connections = {300;'';'Controller.FuelFlow';};
else Plant.Components.FuelSource.connections = {900;'';'Controller.FuelFlow';};
end

if S2C<Steam2Carbon || strcmp(Reformer,'adiabatic')%add anode recirculation
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
    Plant.Components.Mix1.connections = {'FuelSource.Outlet';'recircValve.Out1';'FC1.Flow1Pin'};
    Plant.Components.Mix1.TagInf = {'MassFlow';'Temperature';};
end
switch Reformer
    case 'external'
        %Oxidizer needed to have correct temperature into hot side
        Plant.Components.Oxidizer.type = 'Oxidizer';
        Plant.Components.Oxidizer.name = 'Oxidizer';
        Plant.Components.Oxidizer.inlets = 2;
        Plant.Components.Oxidizer.InitialFlowOut = Oxidized;
        Plant.Components.Oxidizer.InitialFlowOut.T = 1050;
        Plant.Components.Oxidizer.Vol = .1; %volume in m^3
        Plant.Components.Oxidizer.TagInf = {'EquivelanceRatio';'Temperature';'MassFlow'};
        %%reformer
        Plant.Components.Reformer.type = 'Reformer';
        Plant.Components.Reformer.name = 'Reformer';
        Plant.Components.Reformer.direction = 'counterflow';  % 'co-flow', or 'counterflow', neccessary for external reformer with heat recovery stream
        Plant.Components.Reformer.nodes = 5; % # of nodes to simulate
        Plant.Components.Reformer.ReformTarget = .75; %    If there is a 2nd stream providing heat, it can use one of the following 4 targets to initialize reformer: Reforming percent, cold outlet temp, hot outlet temp, effectiveness of heat transfer from hot stream. Note that due to the reforming chemistry, a simple heat exchanger effectivess is not applicable.
        Plant.Components.Reformer.method = 'RefPerc'; %method for sizing reformer to initial conditions. Options are: 'none' fixed size during intialization, 'RefPerc' sizes for a specific external reforming percent, 'ColdT' sizes to match cold exit temp, 'HotT' sizes to match hot ext temp, 'Effectiveness' sizes to match a target effectiveness: % of ideal (infinite area) heat transfer (with no conduction between nodes)
        Plant.Components.Reformer.InletGuess = FuelIn;
        Plant.Components.Reformer.TagInf = {'H2flow';'CH4';};
        
        if S2C<Steam2Carbon %add anode recirculation
            Plant.Components.Reformer.connections = {'Mix1.Outlet';'FC1.Flow1Pin';'Oxidizer.Flow';'';};
            Plant.Components.Oxidizer.connections = {'recircValve.Out2';'FC1.Flow2Out';'Reformer.CooledPin';};
        else
            Plant.Components.Reformer.connections = {'FuelSource.Outlet';'FC1.Flow1Pin';'Oxidizer.Flow';'';};
            Plant.Components.Oxidizer.connections = {'FC1.Flow1Out';'FC1.Flow2Out';'Reformer.CooledPin';};
        end
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

Plant.Components.FC1.type = 'FuelCell';
Plant.Components.FC1.name = 'FC1';
Plant.Components.FC1.FCtype = 'SOFC'; %SOFC, or MCFC 
Plant.Components.FC1.Reformer = Reformer; 
Plant.Components.FC1.direction = 'crossflow'; % 'coflow', or 'counterflow' or 'crossflow'
Plant.Components.FC1.ClosedCathode = 0; %0 means air or some excess flow of O2 in the cathode used as primary means of temerature control (initializations hold to design fuel utilization), 1 means closed end cathode, or simply a fixed oxygen utilization, cooling is done with excess fuel, and the design voltage is met during initialization
Plant.Components.FC1.CoolingStream = 'cathode'; % choices are 'anode' or 'cathode'. Determines which flow is increased to reach desired temperature gradient.
Plant.Components.FC1.PressureRatio = 1.2;
Plant.Components.FC1.columns = 5;
Plant.Components.FC1.rows = 5;
Plant.Components.FC1.RatedStack_kW = 300; %Nominal Stack Power in kW
Plant.Components.FC1.Fuel = Fuel; %initial fuel composition at inlet
Plant.Components.FC1.Flow2Spec = Plant.Components.AirSource.InitialComposition; %initial oxidant composition
Plant.Components.FC1.Steam2Carbon = Steam2Carbon;
Plant.Components.FC1.method = 'Achenbach'; %Determines reforming reaction kinetics options: 'Achenbach' , 'Leinfelder' , 'Drescher'   
Plant.Components.FC1.L_Cell = .09;  %Cell length in meters
Plant.Components.FC1.W_Cell = .09;  %Cell Width in meters  
Plant.Components.FC1.Specification = 'power density'; %options are 'cells', 'power density', 'voltage', or 'current density'. Note: careful when specifying cells that it arrives at a feasible power density
Plant.Components.FC1.SpecificationValue = 450; % power density specified in mW/cm^2, voltage specified in V/cell, current density specified in A/cm^2
Plant.Components.FC1.deltaTStack = 50; %temperature difference from cathode inlet to cathode outlet
Plant.Components.FC1.TpenAvg = 1023;% 750 C, average electrolyte operating temperature
Plant.Components.FC1.FuelUtilization = Utilization; %fuel utilization (net hydrogen consumed/ maximum hydrogen produced with 100% Co and CH4 conversion
Plant.Components.FC1.Flow1Pdrop = 2; %design anode pressure drop
Plant.Components.FC1.Flow2Pdrop = 10; %Design cathode pressure drop
Plant.Components.FC1.Map = 'SOFC_map'; %Radiative heat transfer view factors, imported from CAD
switch Reformer
    case {'direct';'internal'}
        if S2C<Steam2Carbon %add anode recirculation
            Plant.Components.FC1.connections = {'Controller.Current';'Mix1.Outlet';'AirSource.Outlet';'';'';};
        else
            Plant.Components.FC1.connections = {'Controller.Current';'FuelSource.Outlet';'AirSource.Outlet';'';'';};
        end
        if strcmp(Reformer,'internal')
            Plant.Components.FC1.RefPerc = 0.8;% necessary for internal reformer, proportion of CH4 reforming in the reforming channels
            Plant.Components.FC1.RefSpacing = 1;% necessary for internal reformer. This is the # of active cells between reformer plates
        end
    case 'external';
        Plant.Components.FC1.connections = {'Controller.Current';'Reformer.Reformed';'AirSource.Outlet';'Oxidizer.Pin';'Oxidizer.Pin';};
        Plant.Components.FC1.RefPerc = Plant.Components.Reformer.ReformTarget; % necessary for internal or external reformer. This is the percent towards equilibrium occurin in the reformer plaes
    case 'adiabatic'
        Plant.Components.FC1.connections = {'Controller.Current';'Reformer.Reformed';'AirSource.Outlet';'';'';};
    case 'pox'
end
Plant.Components.FC1.TagInf = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';'TcathOut';'ASR';'O2utilization';'nCurrent';'nVoltage';}; %Tags to record at each step
Plant.Components.FC1.TagFinal = {'Power';'Current';'Voltage';'PENavgT';'StackdeltaT';'H2utilization';'O2utilization';'TcathOut';'ASR';'O2utilization'}; %Tags to record at the final step

Plant.Controls.Controller.type = 'ControlFCstack';
Plant.Controls.Controller.name = 'Controller';
Plant.Controls.Controller.Target = {'FC1.TpenAvg';'FC1.deltaTStack';'FC1.Steam2Carbon'};
Plant.Controls.Controller.Fuel = Fuel;
Plant.Controls.Controller.Cells = 'FC1.Cells';
Plant.Controls.Controller.Utilization = 'FC1.FuelUtilization';
Plant.Controls.Controller.InitialAnodeRecirc = 'FC1.Recirc.Anode';
Plant.Controls.Controller.InitConditions = {'FC1.StackCathTin';'FC1.Flow2.IC';'FC1.Current';}; %OxidentTemp, oxidant flow rate, net current
Plant.Controls.Controller.Gain = [0;1e-3];%[1e-2;1e-3;1e-1];
Plant.Controls.Controller.PropGain = [0;0]; %[1;0;1];
Plant.Controls.Controller.TagInf = {'OxidantTemp';'OxidantFlow';'FuelFlow';'Current';'Recirculation'};
Plant.Controls.Controller.connections = {'FC1.MeasureTflow2';'Controller.OxidantTemp';'FC1.MeasureVoltage';};

Plant.Scope = {'Controller.OxidantFlow';'Controller.Current';'Controller.OxidantTemp';}; %must be in TagInf of the corresponding block to work here
Plant.Plot = [Plant.Scope;{'FC1.StackdeltaT';'FC1.PENavgT';'FC1.Voltage';'FC1.TcathOut';}];