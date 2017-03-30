function EC = InitializeElectrolyzer(varargin)
%electrolyzer model with many states: Temperatures (oxidizer plate, cathode, electrolyte, anode, fuel plate [methanator]), Cathode species, anode species, [methanator species] Current, cathode pressure, anode pressure
% Five (5) inlets: {'NetCurrent','Flow1','Flow2','Flow1Pout','Flow2Pout'}
% for an elextrolyzer, Flow 1 is the cathode (steam), and Flow2 is the anode (cooling/heating air)
% Current is negative
global F Ru
F=96485.339; % %Faraday's constant in Coulomb/mole
Ru = 8.314472; % Universal gas constant in kJ/K*kmol
EC = varargin{1};
if length(varargin)==1 % first initialization
    EC.nodes = EC.rows*EC.columns;
    %% Secondary Design Variables
    %%--Geometric Varibles  %
    EC.t_plate1 = 0.003;                    % [m] thickness of plate1
    EC.t_plate1_wall =0.005;                % [m] Thickness of the channel wall of the plate
    EC.H_plate1channels = 0.002;            % [m] height of cathode channel
    EC.W_plate1channels = 0.005;            % [m] width of channel                       
    EC.H_plate2channels = 0.002;            % [m] height of anode channel
    EC.W_plate2channels = 0.005;            % [m] width of  channel
    EC.t_plate2 = 0.003;                    % [m] Thickness ofPlate
    EC.t_plate2_wall=0.005;                 % [m] Thickness of the channel wall 
    EC.Nu_flow1 = 4;                        %Nusselt # for square channel aspect ratio =3 & uniform temp
    EC.Nu_flow2 = 4;                        %Nusselt # for square channel aspect ratio =3 & uniform temp

    %%Electrochemical parameters %%%
    % H2 +1/2O2 --> H2O (Nernst Eo)
    %SOFC Conductivity - Implemented Equation  sigma = A*e^(-deltaG/RT)/T
    switch EC.FCtype
        case {'SOFC';'SOEC'}
            EC.ElecConst = 2e3; %(K/ohm*m) Electrolyte Constant  %default SOFC  = 9e7
            EC.deltaG = 8.0e3; %(kJ/kmol)
            EC.t_Membrane = 18e-6;                     % [m] thickness of membrane
            EC.t_Cath = 800e-6;                        % [m] thickness of cathode structure
            EC.t_An = 50e-6;                           % [m] thickness of Anode structure
            EC.t_Elec = EC.t_Membrane+EC.t_Cath+EC.t_An;        % [m] thickness of complete electrolyte
    end

    %%Electrolyte
    EC.k_Elec =6.19;                                % [W/m K]  Conductivity of the Electrolyte
    EC.Density_Elec = 375;                          % [kg/m3] Density of Electrolyte
    EC.C_Elec = .800;                                  % [kJ/(kg K)] specific heat of electrolyte 

    %%--- Plate 1 -------%
    EC.Density_plate1 = 2000;                                % [kg/m3]     density 
    EC.C_plate1 = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    EC.k_plate1 = 5;   %25                                    % [W/(m K)]   conductivity of Fuel Seperator Plate
    %%-----Anode Gas Stream-----------%
    EC.k_flow1 = 67E-3;                                              % (W/m*K) Thermal Conductivity of air at 1000K

    %%-------Cathode Gas Stream---------%
    EC.k_flow2 = 259E-3;                                 % (W/m*K) Thermal Conductivity of 50%H2 & 50%H2O at 1000K

    %%---- Plate 2 -------%   
    EC.Density_plate2 = 2000;                                % [kg/m3]     density 
    EC.C_plate2 = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    EC.k_plate2 = 5;%25;                                % [W/(m K)]   conductivity of Fuel Seperator Plate


%%%---%%% end of user defined variables    
    %Dimensions
    EC.A_Cell = EC.L_Cell*EC.W_Cell; %Cell Area
    EC.A_Node = EC.A_Cell/EC.nodes; %node Area
    EC.L_node = EC.L_Cell/EC.columns; %Node length in meters
    EC.W_node = EC.W_Cell/EC.rows;  %Node Width in meters
    EC.Dh_Flow1 = 4*(EC.H_plate1channels*EC.W_plate1channels)/(2*(EC.H_plate1channels+EC.W_plate1channels)); %(m) hydraulic diameter of channel
    EC.Dh_Flow2 = 4*(EC.H_plate2channels*EC.W_plate2channels)/(2*(EC.H_plate2channels+EC.W_plate2channels)); %(m) hydraulic diameter of channel
    EC.CH_Flow1 = EC.W_node/(EC.W_plate1channels+EC.t_plate1_wall); %Number of channels in each node of the anode
    EC.CH_Flow2 = EC.W_node/(EC.W_plate2channels+EC.t_plate2_wall); %Number of channels in each node of the cathode 
    %% --- Plate 1 -------%
    EC.A_plate1_elecCond = EC.t_plate1_wall*EC.L_node*EC.CH_Flow1;   % [m2] Conduction area between the fuel seperator plate and the electrolyte
    EC.A_plate1_heatCond = (EC.H_plate1channels*EC.t_plate1_wall + (EC.W_plate1channels+EC.t_plate1_wall)*EC.t_plate1)*EC.CH_Flow1; %[m^2] conduction area between nodes
    EC.L_plate1_heatCond = EC.H_plate1channels;                                     % [m] Lenght of conduction between the fuel seperator plate and electrolyte
    EC.Mass_plate1 = (EC.H_plate1channels*EC.t_plate1_wall + (EC.W_plate1channels+EC.t_plate1_wall)*EC.t_plate1)*EC.CH_Flow1*EC.L_node*EC.Density_plate1;
    %% ----- Cathode Gas Stream-----------%
    EC.flow1_crossArea = EC.H_plate1channels*EC.W_plate1channels*EC.CH_Flow1;              % [m2] Crossectional Area of Anode entrance
    EC.h_flow1          = EC.Nu_flow1*EC.k_flow1/EC.Dh_Flow1;                     % [W/m2/K]  Convection coefficient between the anode gas and the Fuel Seperator plate
    EC.A_flow1_plate1      = (2*EC.H_plate1channels + EC.W_plate1channels)*EC.L_node*EC.CH_Flow1;       % [m2]  Area in common between Anode stream and Sep Plate for convection
    EC.A_flow1_elec     = (EC.W_plate1channels)*EC.L_node*EC.CH_Flow1;                  % [m2]  Area in common between Anode stream and Electrolyte for convection
    EC.Vol_flow1        = EC.H_plate1channels*EC.W_plate1channels*EC.L_node*EC.CH_Flow1;               % [m3]  control volume
    %% --------Electrolyte-------------%
    EC.A_Elec_Cond =  EC.W_node*EC.t_Elec;                   % [m2] Conduction surface area of electrolyte
    EC.A_Elec_Heat_Cond = EC.W_node*EC.t_Elec;                    % [m2] Conduction surface area of electrolyte
    EC.Vol_Elec = EC.t_Elec*EC.L_node*EC.W_node;              % [m3] volume of electrolyte    
    %% ------- Anode Gas Stream---------%
    EC.flow2_crossArea= EC.H_plate2channels*EC.W_plate2channels*EC.CH_Flow2;       % [m2] Crossectional Area of Cathode entrance
    EC.h_flow2= EC.Nu_flow2*EC.k_flow2/EC.Dh_Flow2;                 % [W/m2/K]  Convection coefficient between the Anode gas and the Fuel Seperator plate
    EC.A_flow2_plate2 = (2*EC.H_plate2channels + EC.W_plate2channels)*EC.L_node*EC.CH_Flow2;    % [m2]  Area in common between Cathode stream and Sep Plate for convection
    EC.A_flow2_elec = EC.W_plate2channels*EC.CH_Flow2*EC.L_node;                 % [m2]  Area in common between Cathode stream and Electrolyte for convection
    EC.Vol_flow2 = EC.H_plate2channels*EC.W_plate2channels*EC.CH_Flow2*EC.L_node;            % [m3]  control volume Cathode
    %% ---- Plate 2 -------%   
    EC.A_plate2_elecCond = EC.t_plate2_wall*EC.L_node*EC.CH_Flow2;                % [m2] Conduction area between the fuel seperator plate and the electrolyte
    EC.A_plate2_heatCond = (EC.H_plate2channels*EC.t_plate2_wall + (EC.W_plate2channels+EC.t_plate2_wall)*EC.t_plate2)*EC.CH_Flow2; %[m^2] conduction area between nodes
    EC.L_plate2_heatCond=EC.H_plate2channels;                                    % [m] Length of conduction between the fuel seperator plate and electrolyte	
    EC.Mass_plate2 = (EC.H_plate2channels*EC.t_plate2_wall + (EC.W_plate2channels+EC.t_plate2_wall)*EC.t_plate2)*EC.L_node*EC.CH_Flow2*EC.Density_plate2;
    %% Pressure
    EC.Flow1_Pout = EC.PressureRatio*101;
    EC.Flow1_Pinit = EC.Flow1_Pout + EC.Flow1Pdrop;
    EC.Flow2_Pout =  EC.PressureRatio*101;
    EC.Flow2_Pinit = EC.Flow2_Pout + EC.Flow2Pdrop;
    
    switch EC.Reformer %% Load mask parameters (flow direction)
        case 'methanator'
            EC = FlowDir(EC,3);
        case {'none'}
            EC = FlowDir(EC,2);
    end
    %estimate # of cells
    if EC.ClosedCathode
        EC.Cells = ceil(EC.RatedStack_kW*1000/(1.3*1e4*EC.L_Cell*EC.W_Cell)); %# of cells in stack (assumes 1 A/cm^2) corrected later
    else
        if strcmp(EC.Specification,'cells')
            EC.Cells = EC.SpecificationValue;
            EC.Specification = 'power density';
            EC.SpecificationValue = EC.RatedStack_kW*100/(EC.L_Cell*EC.W_Cell*EC.Cells);
        elseif strcmp(EC.Specification,'power density')
            EC.Cells = ceil(EC.RatedStack_kW*100/(EC.L_Cell*EC.W_Cell*EC.SpecificationValue)); %# of cells in stack
        elseif strcmp(EC.Specification,'current density')
            EC.Cells = ceil(EC.RatedStack_kW*1000/(1.3*1e4*EC.L_Cell*EC.W_Cell*EC.SpecificationValue)); %# of cells in stack (assumes voltage of 1.3)
        elseif strcmp(EC.Specification,'voltage')
            EC.Cells = ceil(EC.RatedStack_kW*1000/(EC.SpecificationValue*1e4*EC.L_Cell*EC.W_Cell)); %# of cells in stack (assumes 1 A/cm^2) corrected later
        end 
    end
     %% Estimate heat generated
    [h,~] = enthalpy(EC.TpenAvg,{'H2','H2O','O2'});
    h_rxn3 = h.H2+.5*h.O2-h.H2O;
    Vbalance = 1./(2*F).*h_rxn3; %voltage that balances heat
    if EC.ClosedCathode
        EC.SpecificationValue = Vbalance;
        EC.Voltage = EC.SpecificationValue;
        EC.Specification = 'voltage';
    else
        if strcmp(EC.Specification,'power density')
            EC.Voltage = 1.3;
        elseif strcmp(EC.Specification,'current density')
            EC.Voltage = EC.RatedStack_kW/EC.Cells*1000/(EC.A_Cell*(100^2))/EC.SpecificationValue; %convert kW to W/cm^2, then divide by A/cm^2 to get V
        elseif strcmp(EC.Specification,'voltage')
            EC.Voltage = EC.SpecificationValue;
        end
    end
    
    %% %% 1st guess at Initial Condition
    EC.Current = zeros(EC.nodes,1);
    if strcmp(EC.Specification,'power density')
        i_avg = -EC.SpecificationValue/EC.Voltage/1000; %convert mW/cm^2 to A/cm^2, assume an initial guess voltage of 0.85
    elseif strcmp(EC.Specification,'current density')
        i_avg = -EC.SpecificationValue;
    elseif strcmp(EC.Specification,'voltage')
        i_avg = -EC.RatedStack_kW/EC.Cells*1000/(EC.A_Cell*(100^2))/EC.Voltage; %convert kW to W/cm^2, then divide by V to get A/cm^2
    end
    for j = 1:1:EC.rows
        EC.Current(1+EC.columns*(j-1):EC.columns*j) =linspace(2,1,EC.columns)/sum(linspace(2,1,EC.columns))*i_avg*(100^2)*EC.A_Cell/EC.rows; %make the initial current guess low to not overutilize H2 in 1st iteration of solution
    end

    EC.StackTempIn = EC.TpenAvg-.9*EC.deltaTStack;
    EC.T.Cath = zeros(EC.nodes,1) + EC.TpenAvg;
    EC.T.Elec = zeros(EC.nodes,1) + EC.TpenAvg;
    EC.T.Anode = zeros(EC.nodes,1) + EC.TpenAvg;
   

    AnIn = EC.Flow2Spec;
    AnIn.T = EC.TpenAvg;
    Cp = SpecHeat(AnIn);
    Q = (EC.Voltage-Vbalance)*sum(abs(EC.Current))*EC.Cells/1000; %estimated kW of heat generated in stack
    EC.AirFlow = abs(Q)/(Cp*EC.deltaTStack); %estimate of cooling air flow

    EC.F2Spec = unique([{'O2';};fieldnames(EC.Flow2Spec)]);
    
    EC.SteamFlow  = sum(abs(EC.Current))/(2*F*1000)/(EC.H2O_Utilization*EC.Flow1Spec.H2O)*EC.Cells; % H2O flow rate,  current/(2*F*1000) = kmol H2
    SteamSpec = fieldnames(EC.Flow1Spec);
    EC.F1Spec = unique([{'H2';'H2O';};SteamSpec]);
    
    Inlet = InletFlow(EC);
    if strcmp(EC.Specification,'current density')
        Inlet.NetCurrent = i_avg*(EC.A_Cell*(100^2));
    end
    Inlet.Flow1.T = EC.TpenAvg -.75*EC.deltaTStack;
    Inlet.Flow2.T = EC.TpenAvg -.75*sign(Q)*EC.deltaTStack;
    
    %% Run Initial Condition
    [Cathode, Anode,EC,Inlet] = solveInitCond(Inlet,EC,1);
    
    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    EC.PortNames = {'NetCurrent','Flow1','Flow2','Flow1Pout','Flow2Pout','Flow1Out','Flow2Out','Flow1Pin','Flow2Pin','MeasureVoltage','MeasurePower','MeasureTpen','MeasureTflow1','MeasureTflow2'};
    EC.NetCurrent.type = 'in';
    EC.NetCurrent.IC = sum(EC.Current);
    
    EC.Flow1.type = 'in';
    EC.Flow1.IC = Inlet.Flow1;

    EC.Flow2.type = 'in';
    EC.Flow2.IC =  Inlet.Flow2;

    EC.Flow1Pout.type = 'in';
    EC.Flow1Pout.IC = EC.Flow1_Pout;
    EC.Flow1Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    EC.Flow2Pout.type = 'in';
    EC.Flow2Pout.IC = EC.Flow2_Pout;
    EC.Flow2Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    EC.Flow1Out.type = 'out';
    EC.Flow1Out.IC  = MergeLastColumn(Cathode.Outlet,EC.Flow1Dir,EC.Cells);

    EC.Flow2Out.type = 'out';
    EC.Flow2Out.IC =  MergeLastColumn(Anode.Outlet,EC.Flow2Dir,EC.Cells);

    EC.Flow1Pin.type = 'out';
    EC.Flow1Pin.IC = EC.Flow1_Pinit;
    EC.Flow1Pin.Pstate = length(EC.Scale)-1; %identifies the state # of the pressure state if this block has one

    EC.Flow2Pin.type = 'out';
    EC.Flow2Pin.IC = EC.Flow2_Pinit;
    EC.Flow2Pin.Pstate = length(EC.Scale); %identifies the state # of the pressure state if this block has one

    EC.MeasureVoltage.type = 'out';
    EC.MeasureVoltage.IC = EC.Voltage;

    EC.MeasurePower.type = 'out';
    EC.MeasurePower.IC = sum(abs(EC.Current)*EC.Voltage*EC.Cells)/1000;%power in kW

    EC.MeasureTpen.type = 'out';
    EC.MeasureTpen.IC = EC.T.Elec;

    EC.MeasureTflow1.type = 'out';
    EC.MeasureTflow1.IC = EC.T.Cath(EC.Flow1Dir(:,end));

    EC.MeasureTflow2.type = 'out';
    EC.MeasureTflow2.IC = EC.T.Anode(EC.Flow2Dir(:,end));

    EC.P_Difference = {'Flow1Pin','Flow1Pout'; 'Flow2Pin', 'Flow2Pout';};
    
    for i = 1:1:length(EC.PortNames)
        if length(EC.connections)<i || isempty(EC.connections{i})
            EC.(EC.PortNames{i}).connected={};
        else
            if ischar(EC.connections{i})
                EC.(EC.PortNames{i}).connected = EC.connections(i);
            else
                EC.(EC.PortNames{i}).IC = EC.connections{i};
                EC.(EC.PortNames{i}).connected={};
            end
        end
    end
end
if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    if strcmp(EC.Specification,'power density')
        EC.Specification = 'current density'; %convergence of power done by controller initialization
    end
    
    Flow1SpecNew = fieldnames(Inlet.Flow1);
    Flow1SpecAll = unique([EC.F1Spec;Flow1SpecNew]);
    Flow1SpecAll = Flow1SpecAll(~strcmp('T',Flow1SpecAll));
    for i = 1:1:length(Flow1SpecAll)
        if ~ismember(Flow1SpecAll{i},Flow1SpecNew)
            Inlet.Flow1.(Flow1SpecAll{i})=0;
        end
    end
    EC.F1Spec = Flow1SpecAll;
    
    Flow2SpecNew = fieldnames(Inlet.Flow2);
    Flow2SpecAll = unique([EC.F2Spec;Flow2SpecNew]);
    Flow2SpecAll = Flow2SpecAll(~strcmp('T',Flow2SpecAll));
    for i = 1:1:length(Flow2SpecAll)
        if ~ismember(Flow2SpecAll{i},Flow2SpecNew)
            Inlet.Flow2.(Flow2SpecAll{i})=0;
        end
    end
    EC.F2Spec = Flow2SpecAll;

    EC.Flow1_Pinit = Inlet.Flow1Pout + EC.Flow1Pdrop;
    EC.Flow2_Pinit = Inlet.Flow2Pout + EC.Flow2Pdrop;
    EC.Flow1Pout.IC = Inlet.Flow1Pout;
    EC.Flow2Pout.IC = Inlet.Flow2Pout;
    EC.StackPower = abs(Inlet.NetCurrent)*EC.Cells/1000*EC.Voltage;
    %%--%%
    [Cathode, Anode,EC,~] = solveInitCond(Inlet,EC,2);
    %%%
    EC.Flow2Pin.Pstate = length(EC.Scale); %identifies the state # of the pressure state if this block has one
    EC.Flow1Pin.Pstate = length(EC.Scale)-1; %identifies the state # of the pressure state if this block has one
    EC.Flow1Out.IC  = MergeLastColumn(Cathode.Outlet,EC.Flow1Dir,EC.Cells);
    EC.Flow2Out.IC = MergeLastColumn(Anode.Outlet,EC.Flow2Dir,EC.Cells);
    EC.Flow1Pin.IC = EC.Flow1_Pinit;
    EC.Flow2Pin.IC = EC.Flow2_Pinit;
    EC.MeasureVoltage.IC = EC.Voltage;
    EC.MeasurePower.IC = sum(abs(EC.Current)*EC.Voltage*EC.Cells)/1000;%power in kW
    EC.MeasureTpen.IC = EC.T.Elec;
    EC.MeasureTflow1.IC = EC.T.Cath(EC.Flow1Dir(:,end));
    EC.MeasureTflow2.IC = EC.T.Anode(EC.Flow2Dir(:,end));
end

function [Cathode, Anode,EC,Inlet] = solveInitCond(Inlet,EC,firstSolve)
global Tags F
OldCurrent = EC.Current;
Cp.cath = SpecHeat(Inlet.Flow1);
Cp.anode = SpecHeat(Inlet.Flow2);
%% loop to converge voltage 
% SOEC & MCEC, adjust current and H2 flow to achieve power density and H2O utilization
error = 1;
Tol = 1e-3;
count = 1;
Methanator = []; %Currently not set up for methane production internally
EC.R_CH4 = zeros(EC.nodes,1);
EC.R_WGS = zeros(EC.nodes,1);
while abs(error)>Tol %iterate to find the initial current (holding power or voltage fixed)
    Cathode = FCin2Out(EC.T.Cath,Inlet.Flow1,EC.Flow1Dir, EC.FCtype,EC.Cells,EC.Current,[],'cathode');
    Anode = FCin2Out(EC.T.Cath,Inlet.Flow2,EC.Flow2Dir, EC.FCtype,EC.Cells,EC.Current,[],'anode');
%     Utilization = (sum(abs(EC.Current))*EC.Cells/(2000*F))/(Inlet.Flow1.H2O);
    if count==1 && firstSolve==1
        [EC.Tstates,EC.HTcond,EC.HTconv]= SteadyTemps(EC,Inlet);
    else
        [~, Y] = ode15s(@(t,y) DynamicECtemps(t,y,EC,Anode,Cathode,Methanator), [0, 1e5], EC.Tstates);
        EC.Tstates = Y(end,:)';
    end
    %organize temperatures
    
    EC.T.Elec =  EC.Tstates(2*EC.nodes+1:3*EC.nodes);
    Tcorrection = EC.TpenAvg - mean(EC.T.Elec);
    EC.T.Cath =  EC.Tstates(1*EC.nodes+1:2*EC.nodes)+Tcorrection;
    EC.T.Anode = EC.Tstates(3*EC.nodes+1:4*EC.nodes)+Tcorrection;

    %% Nernst & Losses
    n_an_in = NetFlow(Anode.Inlet);
    n_an_out = NetFlow(Anode.Outlet);
    n_cath_in = NetFlow(Cathode.Inlet);
    n_cath_out = NetFlow(Cathode.Outlet);
    AvgX.O2 = (Anode.Outlet.O2+Anode.Inlet.O2)./(n_an_in+n_an_out);
    AvgX.H2 = (Cathode.Outlet.H2+Cathode.Inlet.H2)./(n_cath_in+n_cath_out);
    AvgX.H2O = (Cathode.Outlet.H2O+Cathode.Inlet.H2O)./(n_cath_in+n_cath_out);

    %% Calculate local voltages
    normTemp = EC.T.Elec+Tcorrection; %assume you will get to the desired temperature (this avoids oscilations in voltage and helps convergence)
    FuelCellNernst(EC.Current,normTemp,EC.Flow2_Pinit,AvgX,EC);
    
    OldVoltage = EC.Voltage;
    EC.Voltage = sum(Tags.(EC.name).nVoltage)/EC.nodes;
    localR = Tags.(EC.name).LocalOhmic./EC.Current;
    
    if strcmp(EC.Specification,'power density')
        error = (EC.RatedStack_kW - Tags.(EC.name).Power)/EC.RatedStack_kW;
        if count>1
            if firstSolve==1 && error>1e-3
                dP_di = max(.7*(Tags.(EC.name).Power-sum(abs(OldCurrent))*OldVoltage*EC.Cells/1000)/(TotCurrent-sum(abs(OldCurrent))),.7*Tags.(EC.name).Power/(TotCurrent));%change in power with change in current
            else
                dP_di = max(1.15*((Tags.(EC.name).Power-sum(abs(OldCurrent))*OldVoltage*EC.Cells/1000)/(TotCurrent-sum(abs(OldCurrent)))),.7*Tags.(EC.name).Power/(TotCurrent));
            end
            scale = 1+ error*EC.RatedStack_kW/dP_di/TotCurrent;
        else % first time through
            scale = (EC.RatedStack_kW*1000/EC.Cells/EC.Voltage)/sum(-EC.Current); %total current it should have at this new voltage/ total current specified right now
        end
    elseif strcmp(EC.Specification,'voltage')
        Tol = 1e-3;
        error = (EC.SpecificationValue - EC.Voltage)/EC.Voltage;
        if count>1
            dV_di = -Tags.(EC.name).ASR/EC.A_Node/100^2;
            scale = 1 + (EC.Voltage - EC.SpecificationValue)/dV_di/TotCurrent;
        else % first time through
            scale =1+sum((Tags.(EC.name).nVoltage - EC.SpecificationValue)./localR)/sum(EC.Current);
        end
    elseif strcmp(EC.Specification,'current density')
        scale = Inlet.NetCurrent/sum(EC.Current);
        error = (sum(EC.Current) - Inlet.NetCurrent)/sum(EC.Current);
    end
    
    OldCurrent = EC.Current;
    EC.Current = redistributeCurrent(EC.Current,scale,Tags.(EC.name).nVoltage,localR,EC.Voltage); %% start by maintaining same current in each row, then allow row voltages to balance (prevents fuel starvation in a row during initialization)
    TotCurrent = sum(abs(EC.Current));
    if firstSolve ==1
        if strcmp(EC.Specification,'voltage') || strcmp(EC.Specification,'current density')
            EC.Cells = ceil((EC.RatedStack_kW*1000/EC.Voltage)/(sum(-EC.Current))); %re-calculate the # of cells
        end
        EC.SteamFlow  = TotCurrent/(2*F*1000)/(EC.H2O_Utilization*EC.Flow1Spec.H2O)*EC.Cells; % Fresh fuel flow rate,  current/(2*F*1000) = kmol H2

        if EC.ClosedCathode
            EC.AirFlow = 0;
            Inlet = InletFlow(EC);
            Inlet.Flow2.T = EC.TpenAvg; % no anode flow in
        else
            Q_cathode = EC.Cells*sum(NetFlow(Cathode.Outlet).*SpecHeat(Cathode.Outlet).*(Cathode.Outlet.T - Cathode.Inlet.T));
            [h,~] = enthalpy(EC.TpenAvg,{'H2','H2O','O2'});
            h_rxn3 = h.H2+.5*h.O2-h.H2O;
            Vbalance = 1/(2*F)*h_rxn3; %voltage that balances heat
            EC.AirFlow = abs((EC.Cells*(EC.Voltage - Vbalance)*TotCurrent/1000) - Q_cathode)/(33*50); %air flow is extra heat / Cp* deltaT
            Inlet = InletFlow(EC);
            if ((EC.Cells*(EC.Voltage - Vbalance)*TotCurrent/1000) - Q_cathode)>0
                Inlet.Flow2.T = EC.TpenAvg-100; %cooling stack
            else
                Inlet.Flow2.T= EC.TpenAvg+100;%heating stack
            end
        end
        Inlet.Flow1.T = EC.TpenAvg -.75*EC.deltaTStack;
        
    end
    count= count+1;
end
%% Finish initialization by organizing state variables
EC.PfactorAnode = EC.AirFlow/EC.Flow2Pdrop;
EC.PfactorCath = EC.SteamFlow/EC.Flow1Pdrop;
EC = Set_IC(EC,Anode,Cathode,Methanator);


function Inlet = InletFlow(EC) %only used 1st time through initialization (before we know what is connected to inlet
% Anode (oxygen production)
if EC.AirFlow>0
    for i = 1:1:length(EC.F2Spec)
        if isfield(EC.Flow2Spec,EC.F2Spec{i})
            Inlet.Flow2.(EC.F2Spec{i}) = EC.Flow2Spec.(EC.F2Spec{i})*EC.AirFlow;%flow rate of every species entering the anode (or reformer if there is one)
        else Inlet.Flow2.(EC.F2Spec{i}) = 0;
        end
    end
else % no dilution air
    for i = 1:1:length(EC.F2Spec)
        Inlet.Flow2.(EC.F2Spec{i}) = 0;
    end
end

%Cathode H2 production
switch EC.FCtype
    case {'MCEC';'SOEC'}
        for i = 1:1:length(EC.F1Spec)
            if isfield(EC.Flow1Spec,EC.F1Spec{i})
                Inlet.Flow1.(EC.F1Spec{i}) = EC.Flow1Spec.(EC.F1Spec{i})*EC.SteamFlow;
            else Inlet.Flow1.(EC.F1Spec{i}) = 0;
            end
        end
end

function EC = Set_IC(EC,Anode,Cathode,Methanator)
global Ru
Cp.cath = SpecHeat(Cathode.Inlet);
Cp.an = SpecHeat(Anode.Outlet);
if ~isempty(Methanator)
    Cp.meth = SpecHeat(Methanator.Outlet);
end
switch EC.Reformer
    case 'none'
        NumOfStates = (5 + length(EC.F2Spec) + length(EC.F1Spec) + 1)*EC.nodes + 2; % 5 temperatures, anode species & cathode species & current at each node and 2 states for anode/cathode pressure 
    case 'methanator'
end

EC.Scale = ones(NumOfStates,1); %2 states for anode/cathode pressure
EC.IC = EC.Scale;
EC.tC = EC.IC; % time constant for derivative dY
EC.Scale(1:5*EC.nodes) = EC.Tstates(1:5*EC.nodes);%temperature (K)

EC.tC(1:EC.nodes) = (EC.Mass_plate2*EC.C_plate2);
EC.tC(1+EC.nodes:2*EC.nodes) = (EC.Vol_flow2*Cp.cath*EC.Flow1_Pinit./(Ru*EC.T.Cath));
EC.tC(2*EC.nodes+1:3*EC.nodes) = (EC.Vol_Elec*EC.Density_Elec*EC.C_Elec);
EC.tC(3*EC.nodes+1:4*EC.nodes) = (EC.Vol_flow1*Cp.an*EC.Flow2_Pinit./(Ru*EC.T.Anode));
EC.tC(4*EC.nodes+1:5*EC.nodes) = (EC.Mass_plate1*EC.C_plate1);

n = 5*EC.nodes;
for i = 1:1:length(EC.F1Spec)
    EC.tC(n+1:n+EC.nodes) = (EC.Vol_flow2*EC.Flow1_Pinit)./(EC.T.Cath*Ru);  % cathode 
    if any(Cathode.Outlet.(EC.F1Spec{i})==0)
        Flow = NetFlow(Cathode.Outlet);
        EC.IC(n+1:n+EC.nodes) = Cathode.Outlet.(EC.F1Spec{i})./Flow;
        EC.Scale(n+1:n+EC.nodes) = Flow; n = n+EC.nodes; %cathode flows
    else
        EC.Scale(n+1:n+EC.nodes) = Cathode.Outlet.(EC.F1Spec{i}); n = n+EC.nodes; %cathode flows
    end
end

Flow = NetFlow(Anode.Outlet);
for i = 1:1:length(EC.F2Spec)
    EC.tC(n+1:n+EC.nodes) = (EC.Vol_flow1*EC.Flow2_Pinit)./(EC.T.Anode*Ru); %anode
    if any(Anode.Outlet.(EC.F2Spec{i})==0) %concentration less than 1%
        EC.IC(n+1:n+EC.nodes) = Anode.Outlet.(EC.F2Spec{i})./Flow;%concentration
        EC.Scale(n+1:n+EC.nodes) = Flow; %anode flow
    else
        EC.Scale(n+1:n+EC.nodes) = Anode.Outlet.(EC.F2Spec{i}); %individual species flow
    end   
    n = n+EC.nodes;
end

switch EC.Reformer
    case 'methanator'
        for i = 1:1:length(EC.F1Spec)
            EC.tC(n+1:n+EC.nodes) = (EC.Vol_Methanator*EC.FuelPinit)./(EC.T.Methanator*Ru);
            if any(Methanator.Outlet.(EC.F1Spec{i})==0)
                Flow = NetFlow(Methanator.Outlet);
                EC.IC(n+1:n+EC.nodes) = Methanator.Outlet.(EC.F1Spec{i})./Flow;
                EC.Scale(n+1:n+EC.nodes) = Flow; 
            else
                EC.Scale(n+1:n+EC.nodes) = Methanator.Outlet.(EC.F1Spec{i}); 
            end
            n = n+EC.nodes;
        end
end
% note sign convention of current is negative for electrolyzers !!
EC.IC(n+1:n+EC.nodes) = -1;  %current
EC.tC(n+1:n+EC.nodes) = 1;  %current
EC.Scale(n+1:n+EC.nodes) = abs(EC.Current);  n = n+EC.nodes; %current

EC.tC(n+1) = (EC.Vol_flow2*EC.nodes*EC.Cells);  %pressure
EC.tC(n+2) = (EC.Vol_flow1*EC.nodes*EC.Cells); %pressure
EC.Scale(n+1) = EC.Flow1_Pinit;%pressure
EC.Scale(n+2) = EC.Flow2_Pinit;%pressure


function dY = DynamicECtemps(t,Y,EC,Anode,Cathode,Methanator)
global F Ru
dY = 0*Y;
nodes = EC.nodes;

if isfield(EC,'tC')
    tC = EC.tC(1:5*nodes);
else
    Cp.cath = 42;% kJ/kmol*K
    Cp.an = 33;% kJ/kmol*K
    tC(1:nodes,1) = (EC.Mass_plate2*EC.C_plate2);
    tC(1+nodes:2*nodes,1) = (EC.Vol_flow2*Cp.cath*EC.Flow1_Pinit./(Ru*EC.T.Cath));
    tC(2*nodes+1:3*nodes,1) = (EC.Vol_Elec*EC.Density_Elec*EC.C_Elec);
    tC(3*nodes+1:4*nodes,1) = (EC.Vol_flow1*Cp.an*EC.Flow2_Pinit./(Ru*EC.T.Anode));
    tC(4*nodes+1:5*nodes,1) = (EC.Mass_plate1*EC.C_plate1);
end

[h,~] = enthalpy(Y(1+2*nodes:3*nodes),{'H2','H2O','O2','CO2'});
Power = EC.Voltage*EC.Current/1000; %cell power in kW
Qgen = EC.Current/(2000*F).*(h.H2+.5*h.O2-h.H2O) - Power;%kW of heat generated by electrochemistry (per node & per cell)
switch EC.FCtype%ion transport across membrane (total enthalpy)
    case {'SOFC';'SOEC'}
        Qion = EC.Current/(4000*F).*h.O2; %O2 ion crossing over (kW)
    case {'MCFC';'MCEC'}
        Qion = EC.Current/(4000*F).*h.O2 + EC.Current/(2000*F).*h.CO2;% O2 & CO2 ion crossing over
end

QT = EC.HTcond*Y + EC.HTconv*Y;

Cathode.Outlet.T = Y(nodes+1:2*nodes);
for j = 1:1:length(EC.Flow1Dir(1,:));%1:columns
    k = EC.Flow1Dir(:,j);
    if j~=1
        Cathode.Inlet.T(k,1) = Cathode.Outlet.T(kprev);
    end
    kprev = k;
end

Anode.Outlet.T = Y(3*nodes+1:4*nodes);
for j = 1:1:length(EC.Flow2Dir(1,:));%1:columns
    k = EC.Flow2Dir(:,j);
    if j~=1
        Anode.Inlet.T(k,1) = Anode.Outlet.T(kprev);
    end
    kprev = k;
end

%energy flows & sepcific heats
HoutCath = enthalpy(Cathode.Outlet);
HinCath = enthalpy(Cathode.Inlet);
HoutAnode = enthalpy(Anode.Outlet);
HinAnode = enthalpy(Anode.Inlet);

switch EC.Reformer
    case 'methanator'
        Methanator.Outlet.T = Y(5*nodes+1:6*nodes);
        for j = 1:1:length(EC.Flow3Dir(1,:));%1:columns
            k = EC.Flow3Dir(:,j);
            if j~=1
                Methanator.Inlet.T(k,1) = Methanator.Outlet.T(kprev);
            end
            kprev = k;
        end
        HoutMeth = enthalpy(Methanator.Outlet);
        HinMeth = enthalpy(Methanator.Inlet);
end

if EC.ClosedCathode %%energy balance
    Qimbalance = sum((HinCath(EC.Flow1Dir(:,1))) - sum(HoutCath(EC.Flow1Dir(:,end)))) + sum(HinAnode(EC.Flow2Dir(:,1)))  - sum(HoutAnode(EC.Flow2Dir(:,end))) - sum(Power);
    Power = Power + Qimbalance*Power./sum(Power);
    Qgen = EC.Current/(2000*F).*(h.H2+.5*h.O2-h.H2O) - Power;%kW of heat generated by electrochemistry (per node & per cell)
end

dY(1:nodes)= QT(1:nodes)./tC(1:nodes);  %Ox Sep Plate
dY(1+nodes:2*nodes)= (QT(1+nodes:2*nodes) + HinCath - HoutCath - Power - Qgen + Qion)./tC(1+nodes:2*nodes); %Cathode
dY(1+2*nodes:3*nodes)= (QT(1+2*nodes:3*nodes) + Qgen)./tC(2*nodes+1:3*nodes); %Electrolyte Plate
dY(1+3*nodes:4*nodes)= (QT(1+3*nodes:4*nodes) + HinAnode - HoutAnode - Qion )./tC(1+3*nodes:4*nodes);  %Anode
dY(1+4*nodes:5*nodes)= QT(1+4*nodes:5*nodes)./tC(4*nodes+1:5*nodes);  %Fuel Sep Plate
switch EC.Reformer
    case 'methanator'
        dY(1+5*nodes:6*nodes)= (EC.RefSpacing*QT(1+5*nodes:6*nodes) + HinMeth - HoutMeth)./tC(1+5*nodes:6*nodes);  %Fuel Reformer Channels
end