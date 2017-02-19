function EC = InitializeElectrolyzer(varargin)
%electrolyzer model with many states: Temperatures (oxidizer plate, cathode, electrolyte, anode, fuel plate [methanator]), Cathode species, anode species, [methanator species] Current, cathode pressure, anode pressure
% Five (5) inlets: {'PowerSupply','CathodeIn','AnodeIn','CathodePressureOut','AnodePressureOut'}
global F Ru
F=96485.339; % %Faraday's constant in Coulomb/mole
Ru = 8.314472; % Universal gas constant in kJ/K*kmol
EC = varargin{1};
if length(varargin)==1 % first initialization

    %% Load mask parameters (flow direction)
    EC.nodes = EC.rows*EC.columns;
    EC = FlowDir(EC);
    Inlet = [];

    EC.SpecCathIn = EC.Steam;

    %% Secondary Design Variables
    %%--Geometric Varibles  %
    EC.t_An_plate = 0.003;                   % [m] thickness of seperator plate
    EC.t_An_plate_wall =0.005;                  % [m] Thickness of the channel wall of the Fuel Seperator Plate
    EC.H_Anode = 0.002;                       % [m] height of anode channel
    EC.W_Anode = 0.005;                       % [m] width of channel                       
    EC.H_Cathode = 0.002;                      % [m] height of Cathode channel
    EC.W_Cathode = 0.005;                      % [m] width of Cathode channel
    EC.t_Cath_plate = 0.003;                     % [m] Thickness of Oxidant Seperator Plate
    EC.t_Cath_plate_wall=0.005;                     % [m] Thickness of the channel wall of the Oxidant Seperator Plate
    EC.Nu_Anode = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp
    EC.Nu_Cathode = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp

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

    %%---Anode half of interconnect plate-------%
    EC.Density_An_plate = 2000;                                % [kg/m3]     density 
    EC.C_An_plate = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    EC.k_An_plate = 5;   %25                                    % [W/(m K)]   conductivity of Fuel Seperator Plate
    %%-----Anode Gas Stream-----------%
    EC.k_Anode = 67E-3;                                              % (W/m*K) Thermal Conductivity of air at 1000K

    %%-------Cathode Gas Stream---------%
    EC.k_Cathode = 259E-3;                                 % (W/m*K) Thermal Conductivity of 50%H2 & 50%H2O at 1000K

    %%----Cathode half of interconnect plate-------%   
    EC.Density_Cath_plate = 2000;                                % [kg/m3]     density 
    EC.C_Cath_plate = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    EC.k_Cath_plate = 5;%25;                                % [W/(m K)]   conductivity of Fuel Seperator Plate


%%%---%%% end of user defined variables    
    %Dimensions
    EC.A_Cell = EC.L_Cell*EC.W_Cell; %Cell Area
    EC.A_Node = EC.A_Cell/EC.nodes; %node Area
    EC.L_node = EC.L_Cell/EC.columns; %Node length in meters
    EC.W_node = EC.W_Cell/EC.rows;  %Node Width in meters
    EC.Dh_Anode = 4*(EC.H_Anode*EC.W_Anode)/(2*(EC.H_Anode+EC.W_Anode)); %(m) hydraulic diameter of channel
    EC.Dh_Cathode = 4*(EC.H_Cathode*EC.W_Cathode)/(2*(EC.H_Cathode+EC.W_Cathode)); %(m) hydraulic diameter of channel
    EC.CH_Anode = EC.W_node/(EC.W_Anode+EC.t_An_plate_wall); %Number of channels in each node of the anode
    EC.CH_Cathode = EC.W_node/(EC.W_Cathode+EC.t_Cath_plate_wall); %Number of channels in each node of the cathode 
    %% ---Fuel Seperator plate-------%
    EC.A_An_plate_Elec_Cond = EC.t_An_plate_wall*EC.L_node*EC.CH_Anode;   % [m2] Conduction area between the fuel seperator plate and the electrolyte
    EC.A_An_plate_Heat_Cond = (EC.H_Anode*EC.t_An_plate_wall + (EC.W_Anode+EC.t_An_plate_wall)*EC.t_An_plate)*EC.CH_Anode; %[m^2] conduction area between nodes
    EC.L_An_plate_Elec = EC.H_Anode;                                     % [m] Lenght of conduction between the fuel seperator plate and electrolyte
    EC.Mass_An_plate = (EC.H_Anode*EC.t_An_plate_wall + (EC.W_Anode+EC.t_An_plate_wall)*EC.t_An_plate)*EC.CH_Anode*EC.L_node*EC.Density_An_plate;
    %% -----Anode Gas Stream-----------%
    EC.A_Anode_Cross= EC.H_Anode*EC.W_Anode*EC.CH_Anode;              % [m2] Crossectional Area of Anode entrance
    EC.h_Anode_Sep=EC.Nu_Anode*EC.k_Anode/EC.Dh_Anode;                     % [W/m2/K]  Convection coefficient between the anode gas and the Fuel Seperator plate
    EC.A_Anode_Sep = (2*EC.H_Anode + EC.W_Anode)*EC.L_node*EC.CH_Anode;       % [m2]  Area in common between Anode stream and Sep Plate for convection
    EC.h_Anode_Elec=EC.Nu_Anode*EC.k_Anode/EC.Dh_Anode;                    % [W/m2/K]  Convection coefficient between the anode gas and the Electrolyte
    EC.A_Anode_Elec = (EC.W_Anode)*EC.L_node*EC.CH_Anode;                  % [m2]  Area in common between Anode stream and Electrolyte for convection
    EC.Vol_Anode = EC.H_Anode*EC.W_Anode*EC.L_node*EC.CH_Anode;               % [m3]  control volume
    EC.A_Node_Surf= EC.A_Node;                                    % [m^2] Surface Area for internal refoming

    %% --------Electrolyte-------------%
    EC.A_Elec_Cond =  EC.W_node*EC.t_Elec;                   % [m2] Conduction surface area of electrolyte
    EC.A_Elec_Heat_Cond = EC.W_node*EC.t_Elec;                    % [m2] Conduction surface area of electrolyte
    EC.Vol_Elec = EC.t_Elec*EC.L_node*EC.W_node;              % [m3] volume of electrolyte    
    %% -------Cathode Gas Stream---------%
    EC.A_Cathode_Cross= EC.H_Cathode*EC.W_Cathode*EC.CH_Cathode;       % [m2] Crossectional Area of Cathode entrance
    EC.h_Cathode_Sep= EC.Nu_Cathode*EC.k_Cathode/EC.Dh_Cathode;                 % [W/m2/K]  Convection coefficient between the FCdesTempode gas and the Fuel Seperator plate
    EC.A_Cathode_Sep = (2*EC.H_Cathode + EC.W_Cathode)*EC.L_node*EC.CH_Cathode;    % [m2]  Area in common between Cathode stream and Sep Plate for convection
    EC.h_Cathode_Elec= EC.Nu_Cathode*EC.k_Cathode/EC.Dh_Cathode;                % [W/m2/K]  Convection coefficient between the Cathode gas and the Electrolyte
    EC.A_Cathode_Elec = EC.W_Cathode*EC.CH_Cathode*EC.L_node;                 % [m2]  Area in common between Cathode stream and Electrolyte for convection
    EC.Vol_Cathode = EC.H_Cathode*EC.W_Cathode*EC.CH_Cathode*EC.L_node;            % [m3]  control volume Cathode
    %% ----Oxidant Seperator plate-------%   
    EC.A_Cath_plate_Elec_Cond = EC.t_Cath_plate_wall*EC.L_node*EC.CH_Cathode;                % [m2] Conduction area between the fuel seperator plate and the electrolyte
    EC.A_Cath_plate_Heat_Cond = (EC.H_Cathode*EC.t_Cath_plate_wall + (EC.W_Cathode+EC.t_Cath_plate_wall)*EC.t_Cath_plate)*EC.CH_Cathode; %[m^2] conduction area between nodes
    EC.L_Cath_plate_Cond=EC.H_Cathode;                                    % [m] Length of conduction between the fuel seperator plate and electrolyte	
    EC.Mass_Cath_plate = (EC.H_Cathode*EC.t_Cath_plate_wall + (EC.W_Cathode+EC.t_Cath_plate_wall)*EC.t_Cath_plate)*EC.L_node*EC.CH_Cathode*EC.Density_Cath_plate;
    %% Pressure
    EC.CathPout = EC.PressureRatio*101;
    EC.CathPinit = EC.CathPout + EC.CathodePdrop;
    EC.AnPout =  EC.PressureRatio*101;
    EC.AnPinit = EC.AnPout + EC.AnodePdrop;
    
     %% Estimate heat generated
    [h,~] = enthalpy(EC.TpenAvg,{'H2','H2O','O2'});
    h_rxn3 = h.H2+.5*h.O2-h.H2O;
    Vbalance = 1./(2*F).*h_rxn3; %voltage that balances heat
    if EC.ClosedCathode
        EC.DesignTargetValue = Vbalance;
        EC.Voltage = EC.DesignTargetValue;
        EC.DesignTarget = 'voltage';
    else
        if strcmp(EC.DesignTarget,'power density')
            EC.Voltage = 1.3;
        elseif strcmp(EC.DesignTarget,'current density')
            EC.Voltage = EC.RatedStack_kW/EC.Cells*1000/(EC.A_Cell*(100^2))/EC.DesignTargetValue; %convert kW to W/cm^2, then divide by A/cm^2 to get V
        elseif strcmp(EC.DesignTarget,'voltage')
            EC.Voltage = EC.DesignTargetValue;
        end
    end
    
    %% %% 1st guess at Initial Condition
    EC.Current = zeros(EC.nodes,1);
    if strcmp(EC.DesignTarget,'power density')
        i_avg = -EC.DesignTargetValue/EC.Voltage/1000; %convert mW/cm^2 to A/cm^2, assume an initial guess voltage of 0.85
    elseif strcmp(EC.DesignTarget,'current density')
        i_avg = -EC.DesignTargetValue;
    elseif strcmp(EC.DesignTarget,'voltage')
        i_avg = -EC.RatedStack_kW/EC.Cells*1000/(EC.A_Cell*(100^2))/EC.Voltage; %convert kW to W/cm^2, then divide by V to get A/cm^2
    end
    for j = 1:1:EC.rows
        EC.Current(1+EC.columns*(j-1):EC.columns*j) =linspace(2,1,EC.columns)/sum(linspace(2,1,EC.columns))*i_avg*(100^2)*EC.A_Cell/EC.rows; %make the initial current guess low to not overutilize H2 in 1st iteration of solution
    end

    EC.StackTempIn = EC.TpenAvg-.9*EC.deltaTStack;
    EC.T.Cath = zeros(EC.nodes,1) + EC.TpenAvg;
    EC.T.Elec = zeros(EC.nodes,1) + EC.TpenAvg;
    EC.T.Anode = zeros(EC.nodes,1) + EC.TpenAvg;
   

    AnIn = EC.Oxidant;
    AnIn.T = EC.TpenAvg;
    Cp = SpecHeat(AnIn);
    Q = (EC.Voltage-Vbalance)*sum(abs(EC.Current))*EC.Cells/1000; %estimated kW of heat generated in stack
    EC.AirFlow = abs(Q)/(Cp*EC.deltaTStack); %estimate of cooling air flow
    
    OxSpec = fieldnames(EC.Oxidant);
    EC.AnSpec = unique([{'O2';};OxSpec]);
    
    EC.SteamFlow  = sum(abs(EC.Current))/(2*F*1000)/(EC.H2O_Utilization*EC.Steam.H2O)*EC.Cells; % H2O flow rate,  current/(2*F*1000) = kmol H2
    SteamSpec = fieldnames(EC.Steam);
    EC.CathSpec = unique([{'H2';'H2O';};SteamSpec]);
    
    Inlet = InletFlow(EC);
    if strcmp(EC.DesignTarget,'current density')
        Inlet.NetCurrent = i_avg*(EC.A_Cell*(100^2));
    end
    Inlet.CathodeIn.T = EC.TpenAvg -.75*EC.deltaTStack;
    Inlet.AnodeIn.T = EC.TpenAvg -.75*sign(Q)*EC.deltaTStack;
    
    %% Run Initial Condition
    [Cathode, Anode,EC,Inlet] = solveInitCond(Inlet,EC,1);
    
    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    EC.PortNames = {'NetCurrent','CathodeIn','AnodeIn','CathodePressureOut','AnodePressureOut','CathodeOut','AnodeOut','CathodePressureIn','AnodePressureIn','MeasureVoltage','MeasurePower','MeasureTpen','MeasureTcathOut','MeasureTanodeOut'};
    EC.NetCurrent.type = 'in';
    EC.NetCurrent.IC = sum(EC.Current);
    
    EC.CathodeIn.type = 'in';
    EC.CathodeIn.IC = Inlet.CathodeIn;

    EC.AnodeIn.type = 'in';
    EC.AnodeIn.IC =  Inlet.AnodeIn;

    EC.CathodePressureOut.type = 'in';
    EC.CathodePressureOut.IC = EC.CathPout;
    EC.CathodePressureOut.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    EC.AnodePressureOut.type = 'in';
    EC.AnodePressureOut.IC = EC.AnPout;
    EC.AnodePressureOut.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    EC.CathodeOut.type = 'out';
    EC.CathodeOut.IC  = MergeLastColumn(Cathode,EC.CathodeFlowDir,EC.Cells);

    EC.AnodeOut.type = 'out';
    EC.AnodeOut.IC =  MergeLastColumn(Anode,EC.AnodeFlowDir,EC.Cells);

    EC.CathodePressureIn.type = 'out';
    EC.CathodePressureIn.IC = EC.CathPinit;
    EC.CathodePressureIn.Pstate = length(EC.Scale)-1; %identifies the state # of the pressure state if this block has one

    EC.AnodePressureIn.type = 'out';
    EC.AnodePressureIn.IC = EC.AnPinit;
    EC.AnodePressureIn.Pstate = length(EC.Scale); %identifies the state # of the pressure state if this block has one

    EC.MeasureVoltage.type = 'out';
    EC.MeasureVoltage.IC = EC.Voltage;

    EC.MeasurePower.type = 'out';
    EC.MeasurePower.IC = sum(abs(EC.Current)*EC.Voltage*EC.Cells)/1000;%power in kW

    EC.MeasureTpen.type = 'out';
    EC.MeasureTpen.IC = EC.T.Elec;

    EC.MeasureTcathOut.type = 'out';
    EC.MeasureTcathOut.IC = EC.T.Cath(EC.CathodeFlowDir(:,end));

    EC.MeasureTanodeOut.type = 'out';
    EC.MeasureTanodeOut.IC = EC.T.Anode(EC.AnodeFlowDir(:,end));

    EC.P_Difference = {'CathodePressureIn','CathodePressureOut'; 'AnodePressureIn', 'AnodePressureOut';};
    
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
    if strcmp(EC.DesignTarget,'power density')
        EC.DesignTarget = 'current density'; %convergence of power done by controller initialization
    end
    AnSpecNew = fieldnames(Inlet.AnodeIn);
    AnSpecAll = unique([EC.AnSpec;AnSpecNew]);
    AnSpecAll = AnSpecAll(~strcmp('T',AnSpecAll));
    for i = 1:1:length(AnSpecAll)
        if ~ismember(AnSpecAll{i},AnSpecNew)
            Inlet.AnodeIn.(AnSpecAll{i})=0;
        end
    end
    EC.AnSpec = AnSpecAll;

    CathSpecNew = fieldnames(Inlet.CathodeIn);
    CathSpecAll = unique([EC.CathSpec;CathSpecNew]);
    CathSpecAll = CathSpecAll(~strcmp('T',CathSpecAll));
    for i = 1:1:length(CathSpecAll)
        if ~ismember(CathSpecAll{i},CathSpecNew)
            Inlet.CathodeIn.(CathSpecAll{i})=0;
        end
    end
    EC.CathSpec = CathSpecAll;
    
    EC.AnPinit = Inlet.AnodePressureOut + EC.AnodePdrop;
    EC.CathPinit = Inlet.CathodePressureOut + EC.CathodePdrop;
    EC.CathodePressureOut.IC = Inlet.CathodePressureOut;
    EC.AnodePressureOut.IC = Inlet.AnodePressureOut;
    EC.StackPower = abs(Inlet.NetCurrent)*EC.Cells/1000*EC.Voltage;
    %%--%%
    [Cathode, Anode,EC,~] = solveInitCond(Inlet,EC,2);
    %%%
    EC.AnodePressureIn.Pstate = length(EC.Scale); %identifies the state # of the pressure state if this block has one
    EC.CathodePressureIn.Pstate = length(EC.Scale)-1; %identifies the state # of the pressure state if this block has one
    EC.CathodeOut.IC  = MergeLastColumn(Cathode,EC.CathodeFlowDir,EC.Cells);
    EC.AnodeOut.IC = MergeLastColumn(Anode,EC.AnodeFlowDir,EC.Cells);
    EC.CathodePressureIn.IC = EC.CathPinit;
    EC.AnodePressureIn.IC = EC.AnPinit;
    EC.MeasureVoltage.IC = EC.Voltage;
    EC.MeasurePower.IC = sum(abs(EC.Current)*EC.Voltage*EC.Cells)/1000;%power in kW
    EC.MeasureTpen.IC = EC.T.Elec;
    EC.MeasureTcathOut.IC = EC.T.Cath(EC.CathodeFlowDir(:,end));
    EC.MeasureTanodeOut.IC = EC.T.Anode(EC.AnodeFlowDir(:,end));
end

function [Cathode, Anode,EC,Inlet] = solveInitCond(Inlet,EC,firstSolve)
global Tags F
OldCurrent = EC.Current;
Cp.cath = SpecHeat(Inlet.CathodeIn);
Cp.anode = SpecHeat(Inlet.AnodeIn);
%% loop to converge voltage 
% SOEC & MCEC, adjust current and H2 flow to achieve power density and H2O utilization
error = 1;
Tol = 1e-3;
count = 1;
Methanator = []; %Currently not set up for methane production internally
while abs(error)>Tol %iterate to find the initial current (holding power or voltage fixed)
    Cathode = In2Out(EC.T.Cath,Inlet.CathodeIn,EC.CathodeFlowDir, EC.FCtype,EC.Cells,-EC.Current,'cathode');
    Anode = In2Out(EC.T.Cath,Inlet.AnodeIn,EC.AnodeFlowDir, EC.FCtype,EC.Cells,-EC.Current,'anode');
%     Utilization = (sum(abs(EC.Current))*EC.Cells/(2000*F))/(Inlet.CathodeIn.H2O);
    if count==1 && firstSolve==1
        [EC.Tstates,EC.HTmatrix]= SteadyECtemps(EC,Anode,Cathode,Methanator,Inlet);
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
    FuelCellNernst(EC.Current,normTemp,EC.AnPinit,AvgX,EC);
    
    OldVoltage = EC.Voltage;
    EC.Voltage = sum(Tags.(EC.name).nVoltage)/EC.nodes;
    localR = Tags.(EC.name).LocalOhmic./EC.Current;
    
    if strcmp(EC.DesignTarget,'power density')
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
    elseif strcmp(EC.DesignTarget,'voltage')
        Tol = 1e-3;
        error = (EC.DesignTargetValue - EC.Voltage)/EC.Voltage;
        if count>1
            dV_di = -Tags.(EC.name).ASR/EC.A_Node/100^2;
            scale = 1 + (EC.Voltage - EC.DesignTargetValue)/dV_di/TotCurrent;
        else % first time through
%             localR = 3*localR;
            scale =1+sum((Tags.(EC.name).nVoltage - EC.DesignTargetValue)./localR)/sum(EC.Current);
        end
    elseif strcmp(EC.DesignTarget,'current density')
        scale = Inlet.NetCurrent/sum(EC.Current);
        error = (sum(EC.Current) - Inlet.NetCurrent)/sum(EC.Current);
    end
    
    OldCurrent = EC.Current;
    EC.Current = redistributeCurrent(EC.Current,scale,Tags.(EC.name).nVoltage,localR,EC.Voltage); %% start by maintaining same current in each row, then allow row voltages to balance (prevents fuel starvation in a row during initialization)
    TotCurrent = sum(abs(EC.Current));
    if firstSolve ==1
        if strcmp(EC.DesignTarget,'voltage') || strcmp(EC.DesignTarget,'current density')
            EC.Cells = ceil((EC.RatedStack_kW*1000/EC.Voltage)/(sum(-EC.Current))); %re-calculate the # of cells
        end
        EC.SteamFlow  = TotCurrent/(2*F*1000)/(EC.H2O_Utilization*EC.Steam.H2O)*EC.Cells; % Fresh fuel flow rate,  current/(2*F*1000) = kmol H2

        if EC.ClosedCathode
            EC.AirFlow = 0;
            Inlet = InletFlow(EC);
            Inlet.AnodeIn.T = EC.TpenAvg; % no anode flow in
        else
            Q_cathode = EC.Cells*sum(NetFlow(Cathode.Outlet).*SpecHeat(Cathode.Outlet).*(Cathode.Outlet.T - Cathode.Inlet.T));
            [h,~] = enthalpy(EC.TpenAvg,{'H2','H2O','O2'});
            h_rxn3 = h.H2+.5*h.O2-h.H2O;
            Vbalance = 1/(2*F)*h_rxn3; %voltage that balances heat
%             if EC.Voltage > Vbalance
%                 TavgError = (mean(EC.T.Elec) - EC.TpenAvg)/EC.deltaTStack; %cooling stack % too hot = increase flow
%             else
%                 TavgError = (EC.TpenAvg - mean(EC.T.Elec))/EC.deltaTStack; %heating stack  %too hot = reduce flow
%             end
%             if count == 1
                EC.AirFlow = abs((EC.Cells*(EC.Voltage - Vbalance)*TotCurrent/1000) - Q_cathode)/(33*50); %air flow is extra heat / Cp* deltaT
%             else
%                 EC.AirFlow = EC.AirFlow*(1 + .8*TavgError);
%             end
            Inlet = InletFlow(EC);
            if ((EC.Cells*(EC.Voltage - Vbalance)*TotCurrent/1000) - Q_cathode)>0
                Inlet.AnodeIn.T = EC.TpenAvg-100; %cooling stack
            else
                Inlet.AnodeIn.T= EC.TpenAvg+100;%heating stack
            end
        end
        Inlet.CathodeIn.T = EC.TpenAvg -.75*EC.deltaTStack;
        
    end
    count= count+1;
end
%% Finish initialization by organizing state variables
EC.PfactorAnode = EC.AirFlow/EC.AnodePdrop;
EC.PfactorCath = EC.SteamFlow/EC.CathodePdrop;
EC = Set_IC(EC,Anode,Cathode,Methanator);


function Inlet = InletFlow(EC) %only used 1st time through initialization (before we know what is connected to inlet
% Anode (oxygen production)
if EC.AirFlow>0
    for i = 1:1:length(EC.AnSpec)
        if isfield(EC.Oxidant,EC.AnSpec{i})
            Inlet.AnodeIn.(EC.AnSpec{i}) = EC.Oxidant.(EC.AnSpec{i})*EC.AirFlow;%flow rate of every species entering the anode (or reformer if there is one)
        else Inlet.AnodeIn.(EC.AnSpec{i}) = 0;
        end
    end
else % no dilution air
    for i = 1:1:length(EC.AnSpec)
        Inlet.AnodeIn.(EC.AnSpec{i}) = 0;
    end
end

%Cathode H2 production
switch EC.FCtype
    case {'MCEC';'SOEC'}
        for i = 1:1:length(EC.CathSpec)
            if isfield(EC.Steam,EC.CathSpec{i})
                Inlet.CathodeIn.(EC.CathSpec{i}) = EC.Steam.(EC.CathSpec{i})*EC.SteamFlow;
            else Inlet.CathodeIn.(EC.CathSpec{i}) = 0;
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
        NumOfStates = (5 + length(EC.AnSpec) + length(EC.CathSpec) + 1)*EC.nodes + 2; % 5 temperatures, anode species & cathode species & current at each node and 2 states for anode/cathode pressure 
    case 'methanator'
end

EC.Scale = ones(NumOfStates,1); %2 states for anode/cathode pressure
EC.IC = EC.Scale;
EC.tC = EC.IC; % time constant for derivative dY
EC.Scale(1:5*EC.nodes) = EC.Tstates(1:5*EC.nodes);%temperature (K)

EC.tC(1:EC.nodes) = (EC.Mass_Cath_plate*EC.C_Cath_plate);
EC.tC(1+EC.nodes:2*EC.nodes) = (EC.Vol_Cathode*Cp.cath*EC.CathPinit./(Ru*EC.T.Cath));
EC.tC(2*EC.nodes+1:3*EC.nodes) = (EC.Vol_Elec*EC.Density_Elec*EC.C_Elec);
EC.tC(3*EC.nodes+1:4*EC.nodes) = (EC.Vol_Anode*Cp.an*EC.AnPinit./(Ru*EC.T.Anode));
EC.tC(4*EC.nodes+1:5*EC.nodes) = (EC.Mass_An_plate*EC.C_An_plate);

n = 5*EC.nodes;
for i = 1:1:length(EC.CathSpec)
    EC.tC(n+1:n+EC.nodes) = (EC.Vol_Cathode*EC.CathPinit)./(EC.T.Cath*Ru);  % cathode 
    if any(Cathode.Outlet.(EC.CathSpec{i})==0)
        Flow = NetFlow(Cathode.Outlet);
        EC.IC(n+1:n+EC.nodes) = Cathode.Outlet.(EC.CathSpec{i})./Flow;
        EC.Scale(n+1:n+EC.nodes) = Flow; n = n+EC.nodes; %cathode flows
    else
        EC.Scale(n+1:n+EC.nodes) = Cathode.Outlet.(EC.CathSpec{i}); n = n+EC.nodes; %cathode flows
    end
end

Flow = NetFlow(Anode.Outlet);
for i = 1:1:length(EC.AnSpec)
    X = Anode.Outlet.(EC.AnSpec{i})./Flow;%concentration
    EC.tC(n+1:n+EC.nodes) = (EC.Vol_Anode*EC.AnPinit)./(EC.T.Anode*Ru); %anode
    if any(X<.01) %concentration less than 1%
        EC.IC(n+1:n+EC.nodes) = X;
        EC.Scale(n+1:n+EC.nodes) = Flow; %anode flow
    else
        EC.Scale(n+1:n+EC.nodes) = Anode.Outlet.(EC.AnSpec{i}); %individual species flow
    end   
    n = n+EC.nodes;
end

switch EC.Reformer
    case 'methanator'
        for i = 1:1:length(EC.CathSpec)
            EC.tC(n+1:n+EC.nodes) = (EC.Vol_Methanator*EC.FuelPinit)./(EC.T.Methanator*Ru);
            if any(Methanator.Outlet.(EC.CathSpec{i})==0)
                Flow = NetFlow(Methanator.Outlet);
                EC.IC(n+1:n+EC.nodes) = Methanator.Outlet.(EC.CathSpec{i})./Flow;
                EC.Scale(n+1:n+EC.nodes) = Flow; 
            else
                EC.Scale(n+1:n+EC.nodes) = Methanator.Outlet.(EC.AnSpec{i}); 
            end
            n = n+EC.nodes;
        end
end
% note sign convention of current is negative for electrolyzers !!
EC.IC(n+1:n+EC.nodes) = -1;  %current
EC.tC(n+1:n+EC.nodes) = 1;  %current
EC.Scale(n+1:n+EC.nodes) = abs(EC.Current);  n = n+EC.nodes; %current

EC.tC(n+1) = (EC.Vol_Cathode*EC.nodes*EC.Cells);  %pressure
EC.tC(n+2) = (EC.Vol_Anode*EC.nodes*EC.Cells); %pressure
EC.Scale(n+1) = EC.CathPinit;%pressure
EC.Scale(n+2) = EC.AnPinit;%pressure

function Current = redistributeCurrent(Current,scale,Voltage,localR,SetVoltage)
currentPercOfAvg = Current/sum(Current)*length(Current);
dCurrent = (Voltage-SetVoltage)./localR.*currentPercOfAvg; %change in current to balance voltage
if min(abs(Current./dCurrent))<1
    a = .5*min(abs(Current./dCurrent));%ensure dCurrent is never more than half of a step towards zero
else a = .5;
end
dCurrent = a*dCurrent;
scale2  = (scale*sum(Current))/sum(Current+dCurrent); %change in current to get to new power
Current = (Current+dCurrent)*scale2;%re-distribute current to achieve average voltage, then scale to new total current

function Flow = In2Out(T,Inlet,Dir,Type,Cells,Current,c_or_a)
global F
Flow.Outlet.T = T;
k = Dir(:,1);
r = length(k);
spec = fieldnames(Inlet);
for i = 1:1:length(spec)
    if strcmp(spec{i},'T')
        Flow.Inlet.T(k,1) = Inlet.T;
    else
        Flow.Inlet.(spec{i})(k,1) = Inlet.(spec{i})/Cells/r;
    end
end
for j= 1:1:length(Dir(1,:))
    if j>1
        k2 = Dir(:,j);
        for i = 1:1:length(spec)
            if strcmp(spec{i},'T')
                Flow.Inlet.T(k2,1) = Flow.Outlet.T(k);
            else
                Flow.Inlet.(spec{i})(k2,1) = Flow.Outlet.(spec{i})(k);
            end
        end
        k = k2;
    end
    if strcmp(c_or_a,'anode') %anode side (O2 production)
        for i = 1:1:length(spec)
            if strcmp(spec{i},'CO2')
                switch Type
                    case 'SOEC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k); %CO2 species concentration
                    case 'MCEC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k) - 2*Current(k)/(4*F*1000); %CO2 species concentration
                end
            elseif strcmp(spec{i},'O2')
                Flow.Outlet.O2(k,1) = Flow.Inlet.O2(k) + Current(k)/(4*F*1000); %O2 species concentration
            elseif ~strcmp(spec{i},'T')
                Flow.Outlet.(spec{i})(k,1) = Flow.Inlet.(spec{i})(k);
            end
        end
    elseif strcmp(c_or_a,'cathode') % cathode side (H2 production, and H2O use)
        for i = 1:1:length(spec)
            if strcmp(spec{i},'CO2')
                switch Type
                    case 'SOEC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k); %CO2 flow
                    case 'MCEC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k) + 2*Current(k)/(4*F*1000); %CO2 flow
                end
            elseif strcmp(spec{i},'H2O')
                Flow.Outlet.H2O(k,1) = Flow.Inlet.H2O(k) - Current(k)/(2*F*1000); %H2O flow
            elseif strcmp(spec{i},'H2')
                Flow.Outlet.H2(k,1) = Flow.Inlet.H2(k) + Current(k)/(2*F*1000); %H2O flow
            elseif ~strcmp(spec{i},'T')
                Flow.Outlet.(spec{i})(k,1) = Flow.Inlet.(spec{i})(k);
            end
        end
    end
end

function Out = MergeLastColumn(Flow, Dir, Cells)
spec = fieldnames(Flow.Outlet);
for i = 1:1:length(spec)
    if strcmp(spec{i},'T')
        Out.T  = mean(Flow.Outlet.T(Dir(:,end))); %temperature 
    else
        Out.(spec{i}) = max(0,sum(Flow.Outlet.(spec{i})(Dir(:,end)))*Cells);%avoid sending negative outflows to other blocks
    end
end


function [T,HTmatrix]= SteadyECtemps(EC,Anode,Cathode,Methanator,Inlet)
%Solve problem of form xdot = Ax-b for xdot =0.
%final solution is x = A\b;
%states represent  heat transfer into each node/layer, the temperatures of each node/layers, and inlet cathode and anode temperatures, Qerror term associated with the small error in air flow rate so that the deltaT and Tavg constraints can both be satisfied
%like the heat exchanger this averages the inlet and exit temperature for the gas streams, and assumes the solid temperature states correspond to the average for the node
global F
a = .5; %weighting of previous node temperature on convection HT calculation
[h,h_s] = enthalpy(EC.T.Elec,{'H2','H2O','O2','CO','CO2','CH4'});
h_rxn1 = h.CO+3*h.H2-h.CH4-h.H2O;
h_rxn2 = h.CO2+h.H2-h.CO-h.H2O;
h_rxn3 = h.H2+.5*h.O2-h.H2O;
Qgen = (EC.Current/(2*F).*h_rxn3-EC.Current*EC.Voltage)/1000;%kW of heat generated by electrochemistry (per node & per cell)
% Qdirect = h_rxn1.*EC.R_CH4+h_rxn2.*EC.R_WGS; %kW of cooling per cell;
switch EC.Reformer
    case 'methanator'
        Qindirect = h_rxn1.*EC.R_CH4ref+h_rxn2.*EC.R_WGSref; %kW of cooling per cell
end
%ion transport across membrane (sensible enthalpy
switch EC.FCtype
    case {'SOFC';'SOEC'}
        Qion = EC.Current/(4*F*1000).*h_s.O2; %O2 ion crossing over (kW)
    case {'MCFC';'MCEC'}
        Qion = EC.Current/(4*F*1000).*h_s.O2 + EC.Current/(2*F*1000).*h_s.CO2;% O2 & CO2 ion crossing over
end

nodes = EC.nodes;
Aflow = mean(NetFlow(Anode.Outlet));
Cflow = mean(NetFlow(Cathode.Outlet));
switch EC.Reformer
    case 'methanator'
        s=6;
        Cells2Reformer=EC.RefSpacing/2; % set equal to 1 to model the cell next to the reformer
        C8 = (EC.h_Ref_Channel*EC.A_Node)/Cells2Reformer/1000;
        Cp_meth = mean(SpecHeat(Methanator.Outlet));
        Rflow = mean(NetFlow(Methanator.Outlet));
    case 'none'
        s=5;
end
states = 2*s*nodes+2;

A = zeros(states+2,states);
b = zeros(states+2,1);

Cp_cath = mean(SpecHeat(Cathode.Outlet));
Cp_anode = mean(SpecHeat(Anode.Outlet));

%all heat transfer coefficients converted to kW/K: thus Q = C*(T1-T2) is in kW
C1 = (EC.k_Cath_plate*EC.L_node*EC.W_node/EC.t_Cath_plate)/1000;
C2 = (EC.h_Cathode_Sep*EC.A_Cathode_Sep)/1000;
C3 = (EC.h_Cathode_Elec*EC.A_Cathode_Elec)/1000;
C4 = (EC.k_Cath_plate*EC.A_Cath_plate_Elec_Cond/EC.L_Cath_plate_Cond)/1000;
C5 = (EC.k_An_plate*EC.A_An_plate_Elec_Cond/EC.L_An_plate_Elec)/1000;
C6 = (EC.h_Anode_Elec*EC.A_Anode_Elec)/1000;
C7 = (EC.h_Anode_Sep*EC.A_Anode_Sep)/1000;

%horizontal
H1_pn = (EC.k_Cath_plate*EC.A_Cath_plate_Heat_Cond/(EC.L_node/2))/1000; %heat transfer coefficient between previous and next node of oxidant plate
H1_lr  = (EC.k_Cath_plate*EC.A_Cath_plate_Heat_Cond/(EC.W_node/2))/1000; %heat transfer coefficient between left and right adjacent nodes of oxidant plate
H2_pn  = (EC.k_Elec*EC.A_Elec_Heat_Cond/(EC.L_node/2))/1000; %heat transfer coefficient between previous and next node of fuel cell electrolyte assembly
H2_lr  = (EC.k_Elec*EC.A_Elec_Heat_Cond/(EC.W_node/2))/1000; %heat transfer coefficient between left and right adjacent nodes  of fuel cell electrolyte assembly
H3_pn  = (EC.k_An_plate*EC.A_An_plate_Heat_Cond/(EC.L_node/2))/1000; %heat transfer coefficient between previous and next node of fuel plate
H3_lr  = (EC.k_An_plate*EC.A_An_plate_Heat_Cond/(EC.W_node/2))/1000; %heat transfer coefficient between left and right adjacent nodes of fuel plate

for k = 1:1:nodes
    %QT1 : heat transfer into oxidizer plate
    A(k,k) = -1;
    if s==5
        A(k,k+s*nodes) = -(C1+C2+C4);
        A(k,k+(s+4)*nodes) = C1;
    elseif s==6
        A(k,k+s*nodes) = -(C8+C2+C4);
        A(k,k+(s+5)*nodes) = (1-a)*C8;
    end
    A(k,k+(s+1)*nodes) = (1-a)*C2; %averaged with inlet temp or previous node temp
    A(k,k+(s+2)*nodes) = C4;

    %QT2 : heat transfer into cathode
    A(k+nodes,k+nodes) = -1;
    A(k+nodes,k+(s+1)*nodes) = -(1-a)*(C2+C3);
    A(k+nodes,k+s*nodes) = C2;
    A(k+nodes,k+(s+2)*nodes) = C3;

    [i,j] = find(EC.CathodeFlowDir==k);
    if j==1 %first column averaged with inlet temperature
        A(k,2*s*nodes+1) = a*C2;
        A(k+nodes,2*s*nodes+1) = -a*(C2+C3);
        A(k+2*nodes,2*s*nodes+1) = a*C3;
    else % other columns averaged with previous one
        k2 = EC.CathodeFlowDir(i,j-1);
        A(k,k2+(s+1)*nodes) = a*C2;
        A(k+nodes,k2+(s+1)*nodes) = -a*(C2+C3);
        A(k+2*nodes,k2+(s+1)*nodes) = a*C3;
    end

    %QT3 : heat transfer into electrolyte
    A(k+2*nodes,k+2*nodes) = -1;
    A(k+2*nodes,k+(s+2)*nodes) = -(C3+C4+C5+C6);
    A(k+2*nodes,k+(s+1)*nodes) = (1-a)*C3; %averaged with inlet temp or previous node temp
    A(k+2*nodes,k+s*nodes) = C4;
    A(k+2*nodes,k+(s+4)*nodes) = C5;
    A(k+2*nodes,k+(s+3)*nodes) = (1-a)*C6;

    %QT4 : heat transfer into anode
    A(k+3*nodes,k+3*nodes) = -1;
    A(k+3*nodes,k+(s+3)*nodes) = -(1-a)*(C6+C7);
    A(k+3*nodes,k+(s+2)*nodes) = C6;
    A(k+3*nodes,k+(s+4)*nodes) = C7;

    [i,j] = find(EC.AnodeFlowDir==k);
    if j==1 %first column averaged with inlet temperature
        if s==6 %pull last temperature from reformer
            k2 = EC.MethanatorFlowDir(i,end);
            A(k+2*nodes,k2+(s+5)*nodes) = a*C6;
            A(k+3*nodes,k2+(s+5)*nodes) = -a*(C6+C7);
            A(k+4*nodes,k2+(s+5)*nodes) = a*C7;
        else
            A(k+2*nodes,2*s*nodes+2) = a*C6;
            A(k+3*nodes,2*s*nodes+2) = -a*(C6+C7);
            A(k+4*nodes,2*s*nodes+2) = a*C7;
        end
    else % other columns averaged with previous one
        k2 = EC.AnodeFlowDir(i,j-1);
        A(k+2*nodes,k2+(s+3)*nodes) = a*C6;
        A(k+3*nodes,k2+(s+3)*nodes) = -a*(C6+C7);
        A(k+4*nodes,k2+(s+3)*nodes) = a*C7;
    end

    %QT5 : heat transfer into anode plate
    A(k+4*nodes,k+4*nodes) = -1;
    if s==5
        A(k+4*nodes,k+(s+4)*nodes) = -(C1+C5+C7);
        A(k+4*nodes,k+s*nodes) = C1;
    elseif s==6
        A(k+4*nodes,k+(s+4)*nodes) = -(C8+C5+C7);
        A(k+4*nodes,k+(s+5)*nodes) = A(k+4*nodes,k+(s+5)*nodes) + (1-a)*C8;
    end
    A(k+4*nodes,k+(s+2)*nodes) = C5;
    A(k+4*nodes,k+(s+3)*nodes) = (1-a)*C7;

    if s==6 %QT6 : heat transfer into methanator gas
       A(k+5*nodes,k+5*nodes) = -1;
       A(k+5*nodes,k+(s+5)*nodes) = -C8; 
       A(k+5*nodes,k+s*nodes) = C8; 
       A(k+5*nodes,k+(s+4)*nodes) = C8; 

       [i,j] = find(EC.MethanatorFlowDir==k);
        if j==1 %first column averaged with inlet temperature
            A(k,2*s*nodes+2) = a*C8;
            A(k+4*nodes,2*s*nodes+2) = a*C8;
            A(k+5*nodes,2*s*nodes+2) = -C8;
        else % other columns averaged with previous one
            k2 = EC.MethanatorFlowDir(i,j-1);
            A(k,k2+(s+5)*nodes) = a*C8;
            A(k+4*nodes,k2+(s+5)*nodes) = a*C8;
            A(k+5*nodes,k2+(s+5)*nodes) = -C8;
        end
    end

    %Tox : Temeraure of oxidizer plate (net heat transfer into plate = 0)
    A(k+s*nodes,k) = 1;

    %Tcath: Temperature of cathode
    A(k+(s+1)*nodes,k+nodes) = 1;
    A(k+(s+1)*nodes,k+(s+1)*nodes) = -Cp_cath*Cflow;
    [i,j] = find(EC.CathodeFlowDir==k);
    if j==1 %first column receives fresh air
        A(k+(s+1)*nodes,2*s*nodes+1) = Cp_cath*Cflow;
    else
        index = EC.CathodeFlowDir(i,j-1);
        A(k+(s+1)*nodes,index +(s+1)*nodes) = Cp_cath*Cflow;
    end
    b(k+(s+1)*nodes) = Qion(k);

    %Telec: Temperature of electrolyte
    A(k+(s+2)*nodes,k+2*nodes) = 1;
    b(k+(s+2)*nodes) = (-Qgen(k));

    %Tan: Temperature of anode
    A(k+(s+3)*nodes,k+3*nodes) = 1;
    A(k+(s+3)*nodes,k+(s+3)*nodes) = -Cp_anode*Aflow; 
    [i,j] = find(EC.AnodeFlowDir==k);
    if j==1 %first column receives fresh fuel
        if s==5
            A(k+(s+3)*nodes,2*s*nodes+2) = Cp_anode*Aflow; %fresh inlet
        elseif s==6
            index = EC.MethanatorFlowDir(i,end);
            A(k+(s+3)*nodes,index+(s+5)*nodes) = Cp_anode*Aflow; %reformer out
        end
    else
        index = EC.AnodeFlowDir(i,j-1);
        A(k+(s+3)*nodes,index+(s+3)*nodes) = Cp_anode*Aflow;
    end
    b(k+(s+3)*nodes) = -Qion(k);

    %Tanode plate: Temperature of anode plate (net heat transfer into plate = 0)
    A(k+(s+4)*nodes,k+4*nodes) = 1;

    if s==6 %Tmethanator : Temperature of methanator gas
       A(k+(s+5)*nodes,k+5*nodes) = 1;
       A(k+(s+5)*nodes,k+(s+5)*nodes) = -Cp_meth*Rflow;
       [i,j] = find(EC.MethanatorFlowDir==k);
        if j==1 %first column receives fresh fuel
            A(k+(s+5)*nodes,2*s*nodes+2) = Cp_meth*Rflow;%fresh inlet
        else
            index = EC.MethanatorFlowDir(i,j-1);
            A(k+(s+5)*nodes,index+(s+5)*nodes) = Cp_meth*Rflow;
        end 
       b(k+(s+5)*nodes) = Qindirect(k);
    end
end

%% left and right, prev and next
prev = EC.HTadjacent(:,1);
next = EC.HTadjacent(:,2);
left = EC.HTadjacent(:,3);
right = EC.HTadjacent(:,4);
for k = 1:1:nodes
    if prev(k) ~=k
        A(k,k+s*nodes) = A(k,k+s*nodes)-H1_pn;
        A(k,prev(k)+s*nodes) = A(k,prev(k)+s*nodes)+H1_pn;
        A(k+2*nodes,k+(s+2)*nodes) = A(k+2*nodes,k+(s+2)*nodes)-H2_pn;
        A(k+2*nodes,prev(k)+(s+2)*nodes) = A(k+2*nodes,prev(k)+(s+2)*nodes)+H2_pn;
        A(k+4*nodes,k+(s+4)*nodes) = A(k+4*nodes,k+(s+4)*nodes)-H3_pn;
        A(k+4*nodes,prev(k)+(s+4)*nodes) = A(k+4*nodes,prev(k)+(s+4)*nodes)+H3_pn;
    end
    if next(k) ~=k
        A(k,k+s*nodes) = A(k,k+s*nodes)-H1_pn;
        A(k,next(k)+s*nodes) = A(k,next(k)+s*nodes)+H1_pn;
        A(k+2*nodes,k+(s+2)*nodes) = A(k+2*nodes,k+(s+2)*nodes)-H2_pn;
        A(k+2*nodes,next(k)+(s+2)*nodes) = A(k+2*nodes,next(k)+(s+2)*nodes)+H2_pn;
        A(k+4*nodes,k+(s+4)*nodes) = A(k+4*nodes,k+(s+4)*nodes)-H3_pn;
        A(k+4*nodes,next(k)+(s+4)*nodes) = A(k+4*nodes,next(k)+(s+4)*nodes)+H3_pn;
    end
    if left(k) ~=k
        A(k,k+s*nodes) = A(k,k+s*nodes)-H1_lr;
        A(k,left(k)+s*nodes) = A(k,left(k)+s*nodes)+H1_lr;
        A(k+2*nodes,k+(s+2)*nodes) = A(k+2*nodes,k+(s+2)*nodes)-H2_lr;
        A(k+2*nodes,left(k)+(s+2)*nodes) = A(k+2*nodes,left(k)+(s+2)*nodes)+H2_lr;
        A(k+4*nodes,k+(s+4)*nodes) = A(k+4*nodes,k+(s+4)*nodes)-H3_lr;
        A(k+4*nodes,left(k)+(s+4)*nodes) = A(k+4*nodes,left(k)+(s+4)*nodes)+H3_lr;
    end
    if right(k) ~=k
        A(k,k+s*nodes) = A(k,k+s*nodes)-H1_lr;
        A(k,right(k)+s*nodes) = A(k,right(k)+s*nodes)+H1_lr;
        A(k+2*nodes,k+(s+2)*nodes) = A(k+2*nodes,k+(s+2)*nodes)-H2_lr;
        A(k+2*nodes,right(k)+(s+2)*nodes) = A(k+2*nodes,right(k)+(s+2)*nodes)+H2_lr;
        A(k+4*nodes,k+(s+4)*nodes) = A(k+4*nodes,k+(s+4)*nodes)-H3_lr;
        A(k+4*nodes,right(k)+(s+4)*nodes) = A(k+4*nodes,right(k)+(s+4)*nodes)+H3_lr;
    end   
end

%remove temperature averaging at cathode inlet node
k1 = EC.CathodeFlowDir(:,1); %first column
for n =1:1:length(k1)
    k = k1(n);
    A(k,k+(s+1)*nodes) = A(k,k+(s+1)*nodes) + a*C2; %HT to ox plate from cathode
    A(k,2*s*nodes+1) = A(k,2*s*nodes+1) - a*C2; %HT to ox plate from cathode

    A(k+nodes,k+(s+1)*nodes) = A(k+nodes,k+(s+1)*nodes) - a*(C2+C3); %HT from ox plate and electrolyte to cathode
    A(k+nodes,2*s*nodes+1) = A(k+nodes,2*s*nodes+1) + a*(C2+C3);   %HT from ox plate and electrolyte to cathode

    A(k+2*nodes,k+(s+1)*nodes) = A(k+2*nodes,k+(s+1)*nodes) + a*C3;%HT from the cathode to electrolyte
    A(k+2*nodes,2*s*nodes+1) = A(k+2*nodes,2*s*nodes+1) - a*C3; %HT from the cathode to electrolyte
end
if s==6
    k1 = EC.MethanatorFlowDir(:,1); %first column
    for n =1:1:length(k1)
        k = k1(n);
        A(k+5*nodes,k+(s+5)*nodes) = A(k+5*nodes,k+(s+5)*nodes) - C8; %HT to reform gas from fuel and ox plates
        A(k+5*nodes,2*s*nodes+2) = A(k+5*nodes,2*s*nodes+2) + C8; %HT to reform gas from fuel and ox plates

        A(k,k+(s+5)*nodes) = A(k,k+(s+5)*nodes) + a*C8; %HT from reformer to ox plate
        A(k,2*s*nodes+2) = A(k,2*s*nodes+2) - a*C8;   %HT from reformer to ox plate

        A(k+4*nodes,k+(s+5)*nodes) = A(k+4*nodes,k+(s+5)*nodes) + a*C8;%HT into the fuel plate from reformer
        A(k+4*nodes,2*s*nodes+2) = A(k+4*nodes,2*s*nodes+2) - a*C8; %HT into the fuel plate from reformer
    end

else
    k1 = EC.AnodeFlowDir(:,1); %first column
    for n =1:1:length(k1)
        k = k1(n);
        A(k+3*nodes,k+(s+3)*nodes) = A(k+3*nodes,k+(s+3)*nodes) -a*(C6+C7); %HT to anode from fuel plate & electrolyte
        A(k+3*nodes,2*s*nodes+2) = A(k+3*nodes,2*s*nodes+2) + a*(C6+C7);   %HT to anode from fuel plate & electrolyte

        A(k+4*nodes,k+(s+3)*nodes) = A(k+4*nodes,k+(s+3)*nodes) + a*C7; %HT from anode to fuel plate
        A(k+4*nodes,2*s*nodes+2) = A(k+4*nodes,2*s*nodes+2) - a*C7;   %HT from anode to fuel plate

        A(k+2*nodes,k+(s+3)*nodes) = A(k+2*nodes,k+(s+3)*nodes) + a*C6; %HT into the electrolyte from anode
        A(k+2*nodes,2*s*nodes+2) = A(k+2*nodes,2*s*nodes+2) - a*C6;%HT into the electrolyte from anode
    end
end
HTmatrix = A(1:s*nodes,s*nodes+1:2*s*nodes);%matrix of coefficients to multiply by vector of temperature and get the heat transfer by conduction & convection between layers and nodes
% Qscale = b(2*s*nodes+3)/max(abs(b(s*nodes+1:2*s*nodes)))/10;%scaling so more emphasis is put on balancing temperature than on balancing Q
% b(s*nodes+1:2*s*nodes) = Qscale*b(s*nodes+1:2*s*nodes);
% A(1:s*nodes,s*nodes+1:2*s*nodes+2) = Qscale*A(1:s*nodes,s*nodes+1:2*s*nodes+2);

%Cathode Temperature
A(2*s*nodes+1,2*s*nodes+1) = 1; 
b(2*s*nodes+1) = Inlet.CathodeIn.T;
%Anode Temperature
A(2*s*nodes+2,2*s*nodes+2) = 1;
b(2*s*nodes+2) = Inlet.AnodeIn.T;

%Average electrolyte temp  constraint 
for j = 1:1:nodes %Average all the electrolyte temps to equal Tpenavg
    A(2*s*nodes+3,j+(s+2)*nodes) = 1/nodes;
end 
b(2*s*nodes+3) = EC.TpenAvg;

%cathode dT constraint: T3 - T2 = T2 - T1
c = length(EC.CathodeFlowDir(1,:));% # of columns
r = length(EC.CathodeFlowDir(:,1)); % # of rows
for i =1:1:c-1
    if i ==1
        A(2*s*nodes+3+i,2*s*nodes+1) = 1; %cathode inlet (T1)
    else A(2*s*nodes+3+i,EC.CathodeFlowDir(:,i-1)+(s+1)*nodes) = 1/r;%first column (T1)
    end
    A(2*s*nodes+3+i,EC.CathodeFlowDir(:,i)+(s+1)*nodes) = -2/r; %middle column (T2)
    A(2*s*nodes+3+i,EC.CathodeFlowDir(:,i+1)+(s+1)*nodes) = 1/r; %last column (T3)
    b(2*s*nodes+3+i) = 0;
end
x= A\b;
T = x(s*nodes+1:2*s*nodes);

function dY = DynamicECtemps(t,Y,EC,Anode,Cathode,Methanator)
global F Ru
dY = 0*Y;
nodes = EC.nodes;

if isfield(EC,'tC')
    tC = EC.tC(1:5*nodes);
else
    Cp.cath = 42;% kJ/kmol*K
    Cp.an = 33;% kJ/kmol*K
    tC(1:nodes,1) = (EC.Mass_Cath_plate*EC.C_Cath_plate);
    tC(1+nodes:2*nodes,1) = (EC.Vol_Cathode*Cp.cath*EC.CathPinit./(Ru*EC.T.Cath));
    tC(2*nodes+1:3*nodes,1) = (EC.Vol_Elec*EC.Density_Elec*EC.C_Elec);
    tC(3*nodes+1:4*nodes,1) = (EC.Vol_Anode*Cp.an*EC.AnPinit./(Ru*EC.T.Anode));
    tC(4*nodes+1:5*nodes,1) = (EC.Mass_An_plate*EC.C_An_plate);
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

QT = EC.HTmatrix*Y;

Cathode.Outlet.T = Y(nodes+1:2*nodes);
for j = 1:1:length(EC.CathodeFlowDir(1,:));%1:columns
    k = EC.CathodeFlowDir(:,j);
    if j~=1
        Cathode.Inlet.T(k,1) = Cathode.Outlet.T(kprev);
    end
    kprev = k;
end

Anode.Outlet.T = Y(3*nodes+1:4*nodes);
for j = 1:1:length(EC.AnodeFlowDir(1,:));%1:columns
    k = EC.AnodeFlowDir(:,j);
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
        for j = 1:1:length(EC.MethanatorFlowDir(1,:));%1:columns
            k = EC.MethanatorFlowDir(:,j);
            if j~=1
                Methanator.Inlet.T(k,1) = Methanator.Outlet.T(kprev);
            end
            kprev = k;
        end
        HoutMeth = enthalpy(Methanator.Outlet);
        HinMeth = enthalpy(Methanator.Inlet);
end

if EC.ClosedCathode %%energy balance
    Qimbalance = sum((HinCath(EC.CathodeFlowDir(:,1))) - sum(HoutCath(EC.CathodeFlowDir(:,end)))) + sum(HinAnode(EC.AnodeFlowDir(:,1)))  - sum(HoutAnode(EC.AnodeFlowDir(:,end))) - sum(Power);
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

function EC = FlowDir(EC)
% Script which orients the nodes relative to the flow directions
% direction = 1: co-flow, direction = 2: counter-flow, direction = 3: cross-flow
nodes = EC.nodes;
columns = EC.columns;
rows = EC.rows;
switch EC.direction
    case 'coflow'
        for j = 1:1:columns
            EC.AnodeFlowDir(:,j) = (j:columns:nodes)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
        end
    case 'counterflow'
        for j = 1:1:columns
            EC.AnodeFlowDir(:,j) = (columns-j+1:columns:nodes)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
        end
    case 'crossflow'
        for j = 1:1:rows
            EC.AnodeFlowDir(:,j) = (1+columns*(j-1):j*columns)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
        end
end
for j = 1:1:columns
    EC.CathodeFlowDir(:,j) = (j:columns:nodes)';
end
EC.HTadjacent = zeros(nodes,4);
for i = 1:1:nodes
    EC.HTadjacent(i,1) = i-1;%previous node
    EC.HTadjacent(i,2) = i+1;%next node
    EC.HTadjacent(i,3) = i-columns;%node to left
    EC.HTadjacent(i,4) = i+columns;%node to right
end
EC.HTadjacent(1:columns:end,1)=linspace(1,nodes-columns+1,rows);%first node in each row has nothing before it
EC.HTadjacent(columns:columns:end,2)=linspace(columns,nodes,rows)';%last node in each row has nothing after it
EC.HTadjacent(1:columns,3)=linspace(1,columns,columns);%first row has nothing to left
EC.HTadjacent(end-columns+1:end,4)=linspace(nodes-columns+1,nodes,columns);%last row has nothing to right