function FC = InitializeFuelCell(varargin)
%FC model with many states: Temperatures (oxidizer plate, cathode, electrolyte, anode, fuel plate [possibly a reformer]), Cathode species, anode species (can include anode recirculation for humidification, [reformer: indirect, adiabatic or none], Current, cathode pressure, anode pressure
% Five (5) inlets: {'Net Current','CathodeIn','AnodeIn','CathodePressureOut','AnodePressureOut'}
global F Ru
F=96485.339; % %Faraday's constant in Coulomb/mole
Ru = 8.314472; % Universal gas constant in kJ/K*kmol
FC = varargin{1};
if length(varargin)==1 % first initialization
    FC.nodes = FC.rows*FC.columns;
    FC = FlowDir(FC); %% Load mask parameters (flow direction)
    FC.AnPercEquilib = 1; %CH4 reforming reaches equilibrium at anode exit.
    Inlet = [];
%%%---%%% User defined variables
    %%--Geometric Varibles  %
    FC.t_An_plate = 0.003;                   % [m] thickness of seperator plate
    FC.t_An_plate_wall =0.005;                  % [m] Thickness of the channel wall of the Fuel Seperator Plate
    FC.H_Anode = 0.002;                       % [m] height of anode channel
    FC.W_Anode = 0.005;                       % [m] width of channel                       
    FC.H_Cathode = 0.002;                      % [m] height of Cathode channel
    FC.W_Cathode = 0.005;                      % [m] width of Cathode channel
    FC.t_Cath_plate = 0.003;                     % [m] Thickness of Oxidant Seperator Plate
    FC.t_Cath_plate_wall=0.005;                     % [m] Thickness of the channel wall of the Oxidant Seperator Plate
    FC.Nu_Anode = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp
    FC.Nu_Cathode = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp

    %%Electrochemical parameters %%%
    % H2 +1/2O2 --> H2O (Nernst Eo)
    %SOFC Conductivity - Implemented Equation  sigma = A*e^(-deltaG/RT)/T
    switch FC.FCtype
        case 'SOFC'
            FC.ElecConst = 2e3; %(K/ohm*m) Electrolyte Constant  %default SOFC  = 9e7
            FC.deltaG = 8.0e3; %(kJ/kmol)
            %%unused:
%             FC.Io = 5000;          % [A/m2] activation current density %default SOFC = 1000
%             FC.DeffO2 = 4.0e-5; %(m^2/s)
%             FC.alpha=.7;
            FC.t_Membrane = 18e-6;                     % [m] thickness of membrane
            FC.t_Cath = 800e-6;                        % [m] thickness of cathode structure
            FC.t_An = 50e-6;                           % [m] thickness of Anode structure
            FC.t_Elec = FC.t_Membrane+FC.t_Cath+FC.t_An;        % [m] thickness of complete electrolyte
        case 'MCFC'
            FC.Io = 500;            % [A/m2] activation current density
            FC.alpha=.4;      
            FC.J_L = 6000;            % [A/m2] Limiting  current density   
            FC.Cr0 = 4.7833e-4;%4.25e-5;%
            FC.Cr1 = -6.6667e-7;%-5e-8;%
            FC.t_Elec = 0.003;        % [m] thickness of complete electrolyte
    end

    %%Electrolyte
    FC.k_Elec =6.19;                                % [W/m K]  Conductivity of the Electrolyte
    FC.Density_Elec = 375;                          % [kg/m3] Density of Electrolyte
    FC.C_Elec = .800;                                  % [kJ/(kg K)] specific heat of electrolyte 

    %%---Anode half of interconnect plate-------%
    FC.Density_An_plate = 2000;                                % [kg/m3]     density 
    FC.C_An_plate = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    FC.k_An_plate = 5;   %25                                    % [W/(m K)]   conductivity of Fuel Seperator Plate
    %%-----Anode Gas Stream-----------%
    FC.k_Anode = 259E-3;                                 % (W/m*K) Thermal Conductivity of 50%H2 & 50%H2O at 1000K

    %%-------Cathode Gas Stream---------%
    FC.k_Cathode = 67E-3;                                              % (W/m*K) Thermal Conductivity of air at 1000K

    %%----Cathode half of interconnect plate-------%   
    FC.Density_Cath_plate = 2000;                                % [kg/m3]     density 
    FC.C_Cath_plate = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    FC.k_Cath_plate = 5;%25;                                % [W/(m K)]   conductivity of Fuel Seperator Plate

    %%----- / ReformerChannel---
    FC.H_Reform = 0.005;                   % (m) height of Reformer channel
    FC.W_Reform = 0.003;                  % (m) width of Reformer channel 
    FC.Nu_Reform = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp
    FC.t_Ref_CH_Wall = .001;               % (m)Thickness of reformer channel wall
    FC.k_Ref_Channel = 112.52E-3;                                          % (W/m*K) Thermal Conductivity of 5%CO, 35%CO2 15%H2, 45%H20 at 900K
%%%---%%% end of user defined variables    
    %Dimensions
    FC.A_Cell = FC.L_Cell*FC.W_Cell; %Cell Area
    FC.A_Node = FC.A_Cell/FC.nodes; %node Area
    FC.L_node = FC.L_Cell/FC.columns; %Node length in meters
    FC.W_node = FC.W_Cell/FC.rows;  %Node Width in meters
    FC.Dh_Anode = 4*(FC.H_Anode*FC.W_Anode)/(2*(FC.H_Anode+FC.W_Anode)); %(m) hydraulic diameter of channel
    FC.Dh_Cathode = 4*(FC.H_Cathode*FC.W_Cathode)/(2*(FC.H_Cathode+FC.W_Cathode)); %(m) hydraulic diameter of channel
    FC.CH_Anode = FC.W_node/(FC.W_Anode+FC.t_An_plate_wall); %Number of channels in each node of the anode
    FC.CH_Cathode = FC.W_node/(FC.W_Cathode+FC.t_Cath_plate_wall); %Number of channels in each node of the cathode 
    %% ---Fuel Seperator plate-------%
    FC.A_An_plate_Elec_Cond = FC.t_An_plate_wall*FC.L_node*FC.CH_Anode;   % [m2] Conduction area between the fuel seperator plate and the electrolyte
    FC.A_An_plate_Heat_Cond=(FC.H_Anode*FC.t_An_plate_wall + (FC.W_Anode+FC.t_An_plate_wall)*FC.t_An_plate)*FC.CH_Anode; %[m^2] conduction area between nodes
    FC.L_An_plate_Elec=FC.H_Anode;                                     % [m] Lenght of conduction between the fuel seperator plate and electrolyte
    FC.Mass_An_plate = (FC.H_Anode*FC.t_An_plate_wall + (FC.W_Anode+FC.t_An_plate_wall)*FC.t_An_plate)*FC.CH_Anode*FC.L_node*FC.Density_An_plate;
    %% -----Anode Gas Stream-----------%
    FC.A_Anode_Cross= FC.H_Anode*FC.W_Anode*FC.CH_Anode;              % [m2] Crossectional Area of Anode entrance
    FC.h_Anode_Sep=FC.Nu_Anode*FC.k_Anode/FC.Dh_Anode;                     % [W/m2/K]  Convection coefficient between the anode gas and the Fuel Seperator plate
    FC.A_Anode_Sep = (2*FC.H_Anode + FC.W_Anode)*FC.L_node*FC.CH_Anode;       % [m2]  Area in common between Anode stream and Sep Plate for convection
    FC.h_Anode_Elec=FC.Nu_Anode*FC.k_Anode/FC.Dh_Anode;                    % [W/m2/K]  Convection coefficient between the anode gas and the Electrolyte
    FC.A_Anode_Elec = (FC.W_Anode)*FC.L_node*FC.CH_Anode;                  % [m2]  Area in common between Anode stream and Electrolyte for convection
    FC.Vol_Anode = FC.H_Anode*FC.W_Anode*FC.L_node*FC.CH_Anode;               % [m3]  control volume
    FC.A_Node_Surf= FC.A_Node;                                    % [m^2] Surface Area for internal refoming

    %% --------Electrolyte-------------%
    FC.A_Elec_Cond =  FC.W_node*FC.t_Elec;                   % [m2] Conduction surface area of electrolyte
    FC.A_Elec_Heat_Cond = FC.W_node*FC.t_Elec;                    % [m2] Conduction surface area of electrolyte
    FC.Vol_Elec = FC.t_Elec*FC.L_node*FC.W_node;              % [m3] volume of electrolyte    
    %% -------Cathode Gas Stream---------%
    FC.A_Cathode_Cross= FC.H_Cathode*FC.W_Cathode*FC.CH_Cathode;       % [m2] Crossectional Area of Cathode entrance
    FC.h_Cathode_Sep= FC.Nu_Cathode*FC.k_Cathode/FC.Dh_Cathode;                 % [W/m2/K]  Convection coefficient between the FCdesTempode gas and the Fuel Seperator plate
    FC.A_Cathode_Sep = (2*FC.H_Cathode + FC.W_Cathode)*FC.L_node*FC.CH_Cathode;    % [m2]  Area in common between Cathode stream and Sep Plate for convection
    FC.h_Cathode_Elec= FC.Nu_Cathode*FC.k_Cathode/FC.Dh_Cathode;                % [W/m2/K]  Convection coefficient between the Cathode gas and the Electrolyte
    FC.A_Cathode_Elec = FC.W_Cathode*FC.CH_Cathode*FC.L_node;                 % [m2]  Area in common between Cathode stream and Electrolyte for convection
    FC.Vol_Cathode = FC.H_Cathode*FC.W_Cathode*FC.CH_Cathode*FC.L_node;            % [m3]  control volume Cathode
    
    %% ----Oxidant Seperator plate-------%   
    FC.A_Cath_plate_Elec_Cond = FC.t_Cath_plate_wall*FC.L_node*FC.CH_Cathode;                % [m2] Conduction area between the fuel seperator plate and the electrolyte
    FC.A_Cath_plate_Heat_Cond = (FC.H_Cathode*FC.t_Cath_plate_wall + (FC.W_Cathode+FC.t_Cath_plate_wall)*FC.t_Cath_plate)*FC.CH_Cathode; %[m^2] conduction area between nodes
    FC.L_Cath_plate_Elec=FC.H_Cathode;                                    % [m] Length of conduction between the fuel seperator plate and electrolyte
    FC.Mass_Cath_plate = (FC.H_Cathode*FC.t_Cath_plate_wall + (FC.W_Cathode+FC.t_Cath_plate_wall)*FC.t_Cath_plate)*FC.L_node*FC.CH_Cathode*FC.Density_Cath_plate;
    switch FC.Reformer
        case 'internal'
            %% ------ Seperate Reformer Channels ---
            FC.Dh_Reform = 4*(FC.H_Reform*FC.W_Reform)/(2*(FC.H_Reform+FC.W_Reform)); %(m) hydraulic diameter of channel
            FC.CH_Reform = FC.W_node/(FC.W_Reform+FC.t_Ref_CH_Wall);    % Number of channels in each node of the reformer
            FC.Vol_Reform=FC.H_Reform*FC.W_Reform*FC.CH_Reform*FC.L_node; % (m^3) Volume of Reformer Channel in cell
            FC.Ref_Channel_Area_Cross=FC.H_Reform*FC.W_Reform*FC.CH_Reform;     % (m^2) Reformer Channel Area per node
            FC.h_Ref_Channel=FC.Nu_Reform*FC.k_Ref_Channel/FC.Dh_Reform;                     % [W/m2/K]  Convection coefficient between the anode gas and the Fuel Seperator plate   
            FC.A_Ref_Node=FC.A_Node;                         % Reformer Channel reformation area
%         case 'adiabatic'
%             FC.Vol_Reform = .01;        % [m^3] volume
%             FC.C_Reformer = .600;       % [kJ/(kg K)] specific heat of fuel seperator plate
%             FC.M_Reformer = 20;         % [kg] mass of reformer
    end

    %% Pressure
    FC.CathPout = FC.PressureRatio*101;
    FC.AirPinit = FC.CathPout + FC.CathodePdrop;
    FC.AnPout =  FC.PressureRatio*101;
    FC.FuelPinit = FC.AnPout + FC.AnodePdrop;

    %% %% 1st guess at Initial Condition
    FC.Current = zeros(FC.nodes,1);
    if strcmp(FC.DesignTarget,'power density')
        FC.Voltage = .85;
        i_avg = FC.DesignTargetValue/FC.Voltage/1000; %convert mW/cm^2 to A/cm^2, assume an initial guess voltage of 0.85
    elseif strcmp(FC.DesignTarget,'current density')
        i_avg = FC.DesignTargetValue;
        FC.Voltage = FC.RatedStack_kW/FC.Cells*1000/(FC.A_Cell*(100^2))/i_avg; %convert kW to W/cm^2, then divide by A/cm^2 to get V
        Inlet.NetCurrent = i_avg*(FC.A_Cell*(100^2));
    elseif strcmp(FC.DesignTarget,'voltage')
        FC.Voltage = FC.DesignTargetValue;
        i_avg = FC.RatedStack_kW/FC.Cells*1000/(FC.A_Cell*(100^2))/FC.Voltage; %convert kW to W/cm^2, then divide by V to get A/cm^2
    end
    for j = 1:1:FC.rows
        FC.Current(1+FC.columns*(j-1):FC.columns*j) =linspace(2,1,FC.columns)/sum(linspace(2,1,FC.columns))*i_avg*(100^2)*FC.A_Cell/FC.rows; %make the initial current guess low to not overutilize H2 in 1st iteration of solution
    end
    FC.StackCathTin  = FC.TpenAvg -.75*FC.deltaTStack;
    FC.T.Cath = zeros(FC.nodes,1) + FC.TpenAvg;
    FC.T.Elec = zeros(FC.nodes,1) + FC.TpenAvg;
    FC.T.Anode = zeros(FC.nodes,1) + FC.TpenAvg;
    FC.T.Reform = zeros(FC.nodes,1) + FC.TpenAvg;
    %% initial guess of reforming cooling
    FC.FuelFlowInit  = sum(FC.Current)/(2*F*1000)/(FC.FuelUtilization*(4*FC.Fuel.CH4+FC.Fuel.CO+FC.Fuel.H2))*FC.Cells; % fuel flow rate,  current/(2*F*1000) = kmol H2
    FC.AnSpec = fieldnames(FC.Fuel);
    FuelTempIn =FC.TpenAvg-FC.deltaTStack;
    switch FC.Reformer
        case 'adiabatic'
            FC.ReformT = 823; %an initial guess temperature for adiabatic reforme, or if uncommenting line 873, this is the setpoint
            FC.Steam2Carbon = 6; %determines anode recirculation, needs to be high to ensure sufficient temperature for some pre-reforming
    end
    R1 = FC.Fuel.CH4*FC.AnPercEquilib*FC.FuelFlowInit;
    switch FC.Reformer
        case 'external'
            CH4ref_ext = FC.RefPerc*R1;
            FC.R_WGSref =  CH4ref_ext*.8;
            FC.R_CH4 = (R1 - CH4ref_ext)/FC.Cells*ones(FC.nodes,1)/FC.nodes;
            FC.R_WGS = FC.R_CH4*.8;
        case 'internal'
            FC.R_CH4ref = FC.RefPerc*R1/FC.Cells*FC.RefSpacing*ones(FC.nodes,1)/FC.nodes;
            FC.R_WGSref =  FC.R_CH4ref*.8;
            FC.R_CH4 = (R1 - sum(FC.R_CH4ref)*FC.Cells/FC.RefSpacing)/FC.Cells*ones(FC.nodes,1)/FC.nodes;
            FC.R_WGS = FC.R_CH4*.8;
        case 'adiabatic'
            FC.R_CH4ref = 0.5*R1;
            FC.R_WGSref =  FC.R_CH4ref*.8;
            FC.R_CH4 = (R1 - FC.R_CH4ref)/FC.Cells*ones(FC.nodes,1)/FC.nodes;
            FC.R_WGS = FC.R_CH4*.8;
        case 'direct'
            FC.R_CH4 = R1/FC.Cells*ones(FC.nodes,1)/FC.nodes;
            FC.R_WGS =  FC.R_CH4*.8;
    end 
%     %% -- get surface areas and radiation view coefficients from file --%%
%     Dir=strrep(which('InitializeFuelCell.m'),fullfile('Components','Initialization','InitializeFuelCell.m'),'FCMaps');
%     load(fullfile(Dir,FC.Map));
%     f = fieldnames(Map);
%     for i = 1:1:length(f)
%         FC.(f{i}) = Map.(f{i});
%     end
%     Sigma = 5.670367e-11;%kW/(m^2 K^4) %all heat transfer coefficients converted to kW/m^2*K^4: thus Q = sigma*Area*(T1^4-T2^4) is in kW
%     FC.RTmatrix = zeros(s*FC.nodes,s*FC.nodes);
%     %% Here is where it needs to agregate a 100x100 view factor map into the rows & columns of this particular FC
%     %% -- %%
%     for j = 1:1:FC.nodes
%         FC.RTmatrix(j,2*FC.nodes+1:3*FC.nodes) = Sigma*FC.ViewFactorCath(j,:)*FC.A_Node; %view factor from cathode plate to electrolyte
%         FC.RTmatrix(j,j) = -Sigma*sum(FC.ViewFactorCath(j,:))*FC.A_Node; % - sum(view factors) for this node
%         
%         FC.RTmatrix(2*FC.nodes+j,1:FC.nodes) = Sigma*FC.ViewFactorCath(j,:)*FC.A_Node; %view factor from electrolyte to cathode plate
%         FC.RTmatrix(2*FC.nodes+j,2*FC.nodes+j) = -Sigma*sum(FC.ViewFactorCath(j,:))*FC.A_Node; % - sum(view factors) for this node
%         
%         FC.RTmatrix(4*FC.nodes+j,2*FC.nodes+1:3*FC.nodes) = Sigma*FC.ViewFactorAn(j,:)*FC.A_Node; %view factor from anode plate to electrolyte
%         FC.RTmatrix(4*FC.nodes+j,4*FC.nodes+j) = -Sigma*sum(FC.ViewFactorAn(j,:))*FC.A_Node; % - sum(view factors) for this node
%         
%         FC.RTmatrix(2*FC.nodes+j,4*FC.nodes+1:5*FC.nodes) = Sigma*FC.ViewFactorAn(j,:)*FC.A_Node; %view factor from electrolyte to anode plate
%         FC.RTmatrix(2*FC.nodes+j,2*FC.nodes+j) = FC.RTmatrix(2*FC.nodes+j,2*FC.nodes+j) -Sigma*sum(FC.ViewFactorAn(j,:))*FC.A_Node; % - sum(view factors) for this node
%         switch FC.Reformer
%             case 'internal'
%                 FC.RTmatrix(4*FC.nodes+j,1:FC.nodes) = Sigma*FC.ViewFactorRef(j,:)*FC.A_Node; %view factor from anode plate to cathode plate, with reformer channels between
%                 FC.RTmatrix(4*FC.nodes+j,4*FC.nodes+j) = FC.RTmatrix(4*FC.nodes+j,4*FC.nodes+j) -Sigma*sum(FC.ViewFactorRef(j,:))*FC.A_Node; % - sum(view factors) for this node
%                 
%                 FC.RTmatrix(j,4*FC.nodes+1:5*FC.nodes) = Sigma*FC.ViewFactorRef(j,:)*FC.A_Node; %view factor from cathode plate to anode plate, with reformer channels between
%                 FC.RTmatrix(j,j) = FC.RTmatrix(j,j) -Sigma*sum(FC.ViewFactorRef(j,:))*FC.A_Node; % - sum(view factors) for this node
%         end
%     end

    
    for i = 1:1:length(FC.AnSpec)
        Inlet.AnodeIn.(FC.AnSpec{i}) = FC.Fuel.(FC.AnSpec{i})*FC.FuelFlowInit;
    end
    switch FC.FCtype
        case {'SOFC';'MCFC'}
            criticalSpecies = {'CH4';'CO';'CO2';'H2';'H2O'};
    end
    
    for i = 1:1:length(criticalSpecies)
        if ~ismember(criticalSpecies{i},FC.AnSpec)
            Inlet.AnodeIn.(criticalSpecies{i}) = 0;
        end
    end
    FC.AnSpec = unique([FC.AnSpec;criticalSpecies]);
    Inlet.AnodeIn.T = FuelTempIn;
    S2C = FC.Fuel.H2O/(FC.Fuel.CH4+.5*FC.Fuel.CO);
    if S2C<FC.Steam2Carbon %add anode recirculation
        Inlet.AnodeIn.T = 300;%mixing provides humidification & pre-heat
        Inlet.FuelMix.CH4 = FC.Fuel.CH4*FC.FuelFlowInit;
        eWGS = .7; %initial guess of effective CO conversion
        r = 0.5; %Initial guess of anode recirculation
        dr = 1e-5;
        error = 1;
        % Inlet = (Inlet + generated - consumed)*r  + New, thus inlet = New/(1-r) + (generated - consumed)*r/(1-r)
        while abs(error)>1e-6
            Inlet.FuelMix.CO = FC.Fuel.CO*FC.FuelFlowInit/(1-r) + (FC.Fuel.CH4 - eWGS*(FC.Fuel.CH4+FC.Fuel.CO))*FC.FuelFlowInit*r/(1-r);
            Inlet.FuelMix.H2O = FC.Fuel.H2O*FC.FuelFlowInit/(1-r) + (FC.Cells*sum(FC.Current)/(2*F*1000) - (FC.Fuel.CH4 + (FC.Fuel.CH4 + FC.Fuel.CO)*eWGS)*FC.FuelFlowInit)*r/(1-r);
            S2C = Inlet.FuelMix.H2O/(Inlet.FuelMix.CH4 + 0.5*Inlet.FuelMix.CO);
            error = FC.Steam2Carbon - S2C;
            r2 = r+dr;
            COin2 = FC.Fuel.CO*FC.FuelFlowInit/(1-r2) + (FC.Fuel.CH4 - eWGS*(FC.Fuel.CH4+FC.Fuel.CO))*FC.FuelFlowInit*r2/(1-r2);
            H2Oin2 = FC.Fuel.H2O*FC.FuelFlowInit/(1-r2) + (FC.Cells*sum(FC.Current)/(2*F*1000) - (FC.Fuel.CH4 + (FC.Fuel.CH4 + FC.Fuel.CO)*eWGS)*FC.FuelFlowInit)*r2/(1-r2);
            S2C2 = H2Oin2/(Inlet.FuelMix.CH4 + 0.5*COin2);
            dSdr = (S2C2 - S2C)/dr;
            r = r + max(-.5*r,min((1-r)/2,error/dSdr));
        end
        FC.Recirc.Anode = r;
        Inlet.FuelMix.CO2 = FC.Fuel.CO2*FC.FuelFlowInit/(1-r) + eWGS*(FC.Fuel.CH4+FC.Fuel.CO)*FC.FuelFlowInit*r/(1-r);
        Inlet.FuelMix.H2 = FC.Fuel.H2*FC.FuelFlowInit/(1-r) + ((3*FC.Fuel.CH4 + (FC.Fuel.CH4 + FC.Fuel.CO)*eWGS)*FC.FuelFlowInit - FC.Cells*sum(FC.Current)/(2*F*1000))*r/(1-r);
        for i = 1:1:length(FC.AnSpec)
            if ~ismember(FC.AnSpec{i},criticalSpecies)
                Inlet.FuelMix.(FC.AnSpec{i}) = FC.Fuel.(FC.AnSpec{i})*FC.FuelFlowInit/(1-r);
            end
            AnOutlet.(FC.AnSpec{i}) = Inlet.FuelMix.(FC.AnSpec{i}) - Inlet.AnodeIn.(FC.AnSpec{i});
        end
        %%find resulting temperature of mixture
        errorT = 1;
        Inlet.FuelMix.T = FuelTempIn;
        AnOutlet.T = FC.TpenAvg + .5*FC.deltaTStack;
        [~,Hin] = enthalpy(Inlet.AnodeIn);
        [~,Hout] = enthalpy(AnOutlet);
        Hnet = Hin + FC.Recirc.Anode*Hout;
        Cp = SpecHeat(AnOutlet);
        NetFlowMix = NetFlow(Inlet.FuelMix);
        while abs(errorT)>1e-3
            [~,Hmix] = enthalpy(Inlet.FuelMix);
            errorT = (Hnet-Hmix)/(Cp*NetFlowMix);
            Inlet.FuelMix.T = Inlet.FuelMix.T + errorT;
        end 
    else
        FC.Recirc.Anode = 0;
        Inlet.FuelMix = Inlet.AnodeIn;
    end
    FC.T.FuelMix = Inlet.FuelMix.T;
    
    
    Inlet.CathodePressureOut = FC.CathPout;
    Inlet.AnodePressureOut = FC.AnPout;
    if ~isfield(FC,'OxidantUtilization')
        if FC.ClosedCathode
            FC.OxidantUtilization = 1;
        elseif strcmp(FC.Reformer,'internal')
            FC.OxidantUtilization = .33;
        else
            FC.OxidantUtilization = .1;
        end
    end
    FC.AirFlow = FC.Cells*sum(FC.Current)/(4*F*FC.Oxidant.O2)/1000/FC.OxidantUtilization;%kmol of oxidant

    FC.CathSpec = fieldnames(FC.Oxidant);
    if FC.ClosedCathode
        FC.CathSpec = {}; %no cathode flow states
        criticalSpecies = {};
    else
        switch FC.FCtype
            case 'SOFC'
                criticalSpecies = {'O2';'N2';};
            case 'MCFC'
                criticalSpecies = {'CO2';'H2O';'O2';'N2';};
        end
    end
    Inlet = InletFlow(FC,Inlet);
    for i = 1:1:length(criticalSpecies)
        if ~ismember(FC.CathSpec,criticalSpecies{i})
            Inlet.CathodeIn.(criticalSpecies{i}) = 0;
        end
    end
    FC.CathSpec = unique([FC.CathSpec,criticalSpecies]);
    %% Run Initial Condition
    [Cathode, Anode,FC,Inlet] = solveInitCond(Inlet,FC,1);
    
    if strcmp(FC.Reformer,'external') || strcmp(FC.Reformer,'adiabatic')
        FC.Reformer = 'direct'; %external and adiabatic reformers handled in seperate block, after 1st initialization
    end
    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    FC.PortNames = {'NetCurrent','CathodeIn','AnodeIn','CathodePressureOut','AnodePressureOut','CathodeOut','AnodeOut','CathodePressureIn','AnodePressureIn','MeasureVoltage','MeasurePower','MeasureTpen','MeasureTcathOut','MeasureTanodeOut'};
    FC.NetCurrent.type = 'in';
    FC.NetCurrent.IC = sum(FC.Current);
    
    FC.CathodeIn.type = 'in';
    FC.CathodeIn.IC = Inlet.CathodeIn;

    FC.AnodeIn.type = 'in';
    FC.AnodeIn.IC = Inlet.AnodeIn; 

    FC.CathodePressureOut.type = 'in';
    FC.CathodePressureOut.IC = Inlet.CathodePressureOut;
    FC.CathodePressureOut.Pstate = []; %identifies the state # of the pressure state if this block has one

    FC.AnodePressureOut.type = 'in';
    FC.AnodePressureOut.IC = Inlet.AnodePressureOut;
    FC.AnodePressureOut.Pstate = []; %identifies the state # of the pressure state if this block has one

    FC.CathodeOut.type = 'out';
    FC.CathodeOut.IC  = MergeFlows(Cathode.Outlet,FC.AirFlowDir(:,end),FC.Cells);

    FC.AnodeOut.type = 'out';
    FC.AnodeOut.IC = MergeFlows(Anode.Outlet,FC.FuelFlowDir(:,end),FC.Cells);

    FC.CathodePressureIn.type = 'out';
    FC.CathodePressureIn.IC = FC.AirPinit;
    FC.CathodePressureIn.Pstate = length(FC.Scale)-1; %identifies the state # of the pressure state if this block has one

    FC.AnodePressureIn.type = 'out';
    FC.AnodePressureIn.IC = FC.FuelPinit;
    FC.AnodePressureIn.Pstate = length(FC.Scale); %identifies the state # of the pressure state if this block has one

    FC.MeasureVoltage.type = 'out';
    FC.MeasureVoltage.IC = FC.Voltage;

    FC.MeasurePower.type = 'out';
    FC.MeasurePower.IC = sum(FC.Current*FC.Voltage*FC.Cells)/1000;%power in kW

    FC.MeasureTpen.type = 'out';
    FC.MeasureTpen.IC = FC.T.Elec;

    FC.MeasureTcathOut.type = 'out';
    FC.MeasureTcathOut.IC = FC.T.Cath(FC.AirFlowDir(:,end));

    FC.MeasureTanodeOut.type = 'out';
    FC.MeasureTanodeOut.IC = FC.T.Anode(FC.FuelFlowDir(:,end));

    FC.P_Difference = {'CathodePressureIn','CathodePressureOut'; 'AnodePressureIn', 'AnodePressureOut';};

    for i = 1:1:length(FC.PortNames)
        if length(FC.connections)<i || isempty(FC.connections{i})
            FC.(FC.PortNames{i}).connected={};
        else
            if ischar(FC.connections{i})
                FC.(FC.PortNames{i}).connected = FC.connections(i);
            else
                FC.(FC.PortNames{i}).IC = FC.connections{i};
                FC.(FC.PortNames{i}).connected={};
            end
        end
    end
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    FC.DesignTarget = 'current density';%converge only to match current density from controller
    FC.Recirc.Anode = 0; %after 1st initialization recirculation is handled in controller, valve & mixing volume
    if Inlet.CathodeIn.T<(FC.TpenAvg - 1.5*FC.deltaTStack)
%         disp('Fuel cell inlet too cold during initialization')
        Inlet.CathodeIn.T  = FC.TpenAvg -.75*FC.deltaTStack;
    end
    Inlet.FuelMix = Inlet.AnodeIn;
    AnSpecNew = fieldnames(Inlet.AnodeIn);
    AnSpecAll = unique([FC.AnSpec;AnSpecNew]);
    AnSpecAll = AnSpecAll(~strcmp('T',AnSpecAll));
    for i = 1:1:length(AnSpecAll)
        if ~ismember(AnSpecAll{i},AnSpecNew)
            Inlet.AnodeIn.(AnSpecAll{i})=0;
        end
    end
    FC.AnSpec = AnSpecAll;

    if FC.ClosedCathode
        FC.CathSpec = {}; %no cathode flow states
    else
        CathSpecNew = fieldnames(Inlet.CathodeIn);
        CathSpecAll = unique([FC.CathSpec;CathSpecNew]);
        CathSpecAll = CathSpecAll(~strcmp('T',CathSpecAll));
        for i = 1:1:length(CathSpecAll)
            if ~ismember(CathSpecAll{i},CathSpecNew)
                Inlet.CathodeIn.(CathSpecAll{i})=0;
            end
        end
        FC.CathSpec = CathSpecAll;
    end
    FC.FuelPinit = Inlet.AnodePressureOut + FC.AnodePdrop;
    FC.AirPinit = Inlet.CathodePressureOut + FC.CathodePdrop;
    FC.CathodePressureOut.IC = Inlet.CathodePressureOut;
    FC.AnodePressureOut.IC = Inlet.AnodePressureOut;
    %%--%%
    [Cathode, Anode,FC,~] = solveInitCond(Inlet,FC,2);
    %%%
    FC.AnodePressureIn.Pstate = length(FC.Scale); %identifies the state # of the pressure state if this block has one
    FC.CathodePressureIn.Pstate = length(FC.Scale)-1; %identifies the state # of the pressure state if this block has one
    FC.CathodeOut.IC  = MergeFlows(Cathode.Outlet,FC.AirFlowDir(:,end),FC.Cells);
    FC.AnodeOut.IC = MergeFlows(Anode.Outlet,FC.FuelFlowDir(:,end),FC.Cells);
    FC.CathodePressureIn.IC = FC.AirPinit;
    FC.AnodePressureIn.IC = FC.FuelPinit;
    FC.MeasureCurrent.IC = sum(FC.Current);
    FC.MeasurePower.IC = sum(FC.Current*FC.Voltage*FC.Cells)/1000;%power in kW
    FC.MeasureTpen.IC = FC.T.Elec;
    FC.MeasureTcathOut.IC = FC.T.Cath(FC.AirFlowDir(:,end));
    FC.MeasureTanodeOut.IC = FC.T.Anode(FC.FuelFlowDir(:,end));
    FC.HumidifiedFuelTemp.IC = FC.T.FuelMix;
end

function [Cathode, Anode,FC,Inlet] = solveInitCond(Inlet,FC,firstSolve)
global Tags F
Cp.cath = SpecHeat(Inlet.CathodeIn);
Cp.anode = SpecHeat(Inlet.AnodeIn);

error = 1;
Tol = 1e-3;
count = 1;
while abs(error)>Tol %iterate to reach target current density, voltage or power density
    Cathode = In2Out(FC.T.Cath,Inlet.CathodeIn,FC.AirFlowDir, FC.FCtype,FC.Cells,FC.Current,[],'cathode');
%     SinglePassUtilization = (sum(FC.Current)*FC.Cells/(2000*F))/(4*Inlet.FuelMix.CH4+Inlet.FuelMix.CO + Inlet.FuelMix.H2);
    [FC,Inlet.FuelMix,Anode,Reformer] = KineticCoef(FC,Inlet,((firstSolve==1) && (count==1)));%% solve for kinetic reaction coefficient which results in this outlet condition (match R_CH4 & R_WGS)
    Offset = 0;
    if count==1 && firstSolve==1
        [FC.Tstates,FC.HTmatrix]= SteadyFCtemps(FC,Anode,Cathode,Reformer,Inlet);
    else
        [~, Y] = ode15s(@(t,y) DynamicFCtemps(t,y,FC,Anode,Cathode,Reformer,Inlet), [0, 1e4], FC.Tstates);
        FC.Tstates = Y(end,:)';
        if firstSolve==1 && abs(mean(FC.Tstates(2*FC.nodes+1:3*FC.nodes))-FC.TpenAvg)>10 %temperature in solve dynamic is diverging to much (1300K), and messing up reforming solution (this is a temporary fix
            Offset = (mean(FC.Tstates(2*FC.nodes+1:3*FC.nodes))-FC.TpenAvg);
        end
    end
    %organize temperatures
    FC.T.Cath =  FC.Tstates(1*FC.nodes+1:2*FC.nodes) - Offset;
    FC.T.Elec =  FC.Tstates(2*FC.nodes+1:3*FC.nodes) - Offset;
    FC.T.Anode =  FC.Tstates(3*FC.nodes+1:4*FC.nodes) - Offset;
    FC.T.FuelMix = Inlet.FuelMix.T;
    switch FC.Reformer
        case 'internal'
            FC.T.Reform = FC.Tstates(5*FC.nodes+1:6*FC.nodes) - Offset;
        case 'adiabatic'
            FC.T.Reform =  Reformer.Outlet.T; %1 temperature state of reformer before anode outlet
    end

    NernstCalc(Anode,Cathode,FC);
    %% calculate the change in current to converge to the desired power density, voltage, or current
    OldVoltage = FC.Voltage;
    FC.Voltage = sum(Tags.(FC.name).nVoltage)/FC.nodes;
    localR = Tags.(FC.name).LocalOhmic./FC.Current;
    if strcmp(FC.DesignTarget,'power density')
        error = (FC.RatedStack_kW - Tags.(FC.name).Power)/FC.RatedStack_kW;
        if count>1
            if firstSolve==1 && error>1e-3
                dP_di = max(.7*(Tags.(FC.name).Power-sum(OldCurrent)*OldVoltage*FC.Cells/1000)/(TotCurrent-sum(OldCurrent)),.7*Tags.(FC.name).Power/(TotCurrent));%change in power with change in current
%                 dP_di_guess = .7*(Tags.(FC.name).Power-sum(OldCurrent)*OldVoltage*FC.Cells/1000)/(TotCurrent-sum(OldCurrent))
%                 dP_di_constraint = .7*Tags.(FC.name).Power/(TotCurrent)
            else
                dP_di = max(1.15*((Tags.(FC.name).Power-sum(OldCurrent)*OldVoltage*FC.Cells/1000)/(TotCurrent-sum(OldCurrent))),.7*Tags.(FC.name).Power/(TotCurrent));
            end
            scale = 1+ error*FC.RatedStack_kW/dP_di/TotCurrent;
        else % first time through
            scale = (FC.RatedStack_kW*1000/FC.Cells/FC.Voltage)/sum(FC.Current); %total current it should have at this new voltage/ total current specified right now
        end
    elseif strcmp(FC.DesignTarget,'voltage')
        Tol = 1e-3;
        error = (FC.DesignTargetValue - FC.Voltage)/FC.Voltage;
        if count>1
            dV_di = -Tags.(FC.name).ASR/FC.A_Node/100^2;
            scale = 1 + (FC.DesignTargetValue - FC.Voltage)/dV_di/TotCurrent;
        else % first time through
            localR = 3*localR;
            scale =1+sum((Tags.(FC.name).nVoltage-FC.DesignTargetValue)./localR)/sum(FC.Current);
        end
        FC.Cells = ceil((FC.RatedStack_kW*1000/FC.DesignTargetValue)/(scale*sum(FC.Current))); %re-calculate the # of cells
    elseif strcmp(FC.DesignTarget,'current density')
        scale = Inlet.NetCurrent/sum(FC.Current);
        error = (sum(FC.Current) - Inlet.NetCurrent)/sum(FC.Current);
    end
    
    OldCurrent = FC.Current;
    FC.Current = redistributeCurrent(FC.Current,scale,Tags.(FC.name).nVoltage,localR,FC.Voltage); %% start by maintaining same current in each row, then allow row voltages to balance (prevents fuel starvation in a row during initialization)
    TotCurrent = sum(FC.Current);
    if firstSolve ==1 %solving FC block to convergence without other blocks or controller
        FC.R_CH4 = scale*FC.R_CH4;
        FC.R_WGS = scale*FC.R_WGS;
        switch FC.Reformer
            case {'internal';'adiabatic'}
                FC.R_CH4ref = scale*FC.R_CH4ref;% fuel flow scales with current so assume reforming will
                FC.R_WGSref = scale*FC.R_WGSref;
        end
        if strcmp(FC.CoolingStream,'cathode')%air flow balances heat generation to maintain deltaT, heat transfer to anode and any fuel reforming is accounted for
            FC.FuelFlowInit  = TotCurrent/(2*F*1000)/(FC.FuelUtilization*(4*FC.Fuel.CH4+FC.Fuel.CO+FC.Fuel.H2))*FC.Cells; % Fresh fuel flow rate,  current/(2*F*1000) = kmol H2
            k = FC.AirFlowDir(:,end);
            dTerror = (mean((FC.Tstates(k+FC.nodes))- Inlet.CathodeIn.T)/FC.deltaTStack-1);
            FC.AirFlow = FC.AirFlow*(1 + dTerror)*scale^2;
            TavgError = (FC.TpenAvg-mean(FC.Tstates(2*FC.nodes+1:3*FC.nodes)))/FC.deltaTStack;
            FC.StackCathTin = FC.StackCathTin + (TavgError + 0.75*dTerror)*FC.deltaTStack;
        elseif strcmp(FC.CoolingStream,'anode') %oxidant flow rate determined by current, fuel flow rate is now determined by thermal balancing
            FC.AirFlow = FC.Cells*TotCurrent/(4000*F*FC.Oxidant.O2)/FC.OxidantUtilization;%kmol of oxidant
            if FC.ClosedCathode %%energy balance
                HinAnode = enthalpy(Anode.Inlet);
                HoutAnode = enthalpy(Anode.Outlet);
                HinCath = enthalpy(Cathode.Inlet);
                HoutCath = enthalpy(Cathode.Outlet);
                HinReform = enthalpy(Reformer.Inlet);
                HoutReform = enthalpy(Reformer.Outlet);
                Power = FC.Voltage*FC.Current/1000; %cell power in kW
                Qimbalance = sum((HinCath - HoutCath) + (HinAnode - HoutAnode) + (HinReform - HoutReform) - Power);
                h = enthalpy(mean(FC.T.Reform),{'H2','H2O','O2','CO','CO2','CH4'});
                Qreform = (h.CO+3*h.H2-h.CH4-h.H2O) + 0.8*(h.CO2+h.H2-h.CO-h.H2O); %kW of cooling per kmol of fuel
                ExtraFuel = 0.25*Qimbalance*FC.Cells/Qreform/FC.Fuel.CH4;
                error = max(abs(error),abs(Qimbalance/sum(Power)));
            else
                ExtraFuel = 0;
                %need to do something with recirculation
            end
            FC.FuelFlowInit  = FC.Cells*TotCurrent/(2000*F)/(FC.FuelUtilization*(4*FC.Fuel.CH4+FC.Fuel.CO+FC.Fuel.H2)); % re-calculate with revised current
            FC.FuelFlowInit = FC.FuelFlowInit + ExtraFuel;
            FC.FuelUtilization = FC.Cells*sum(FC.Current)/(2000*F)/(FC.FuelFlowInit*(4*FC.Fuel.CH4+FC.Fuel.CO+FC.Fuel.H2));
            % change steam to carbon to affect deltaT?
        end
        Inlet = InletFlow(FC,Inlet);
    end
    count= count+1;
end
Cathode = In2Out(FC.T.Cath,Inlet.CathodeIn,FC.AirFlowDir, FC.FCtype,FC.Cells,FC.Current,[],'cathode');
FC.PfactorAnode = NetFlow(Inlet.AnodeIn)/FC.AnodePdrop;
FC.PfactorCath = NetFlow(Inlet.CathodeIn)/FC.CathodePdrop;
FC = Set_IC(FC,Anode,Cathode,Reformer);

function NernstCalc(Anode,Cathode,FC)
%% Nernst & Losses
n_an_in = NetFlow(Anode.Inlet);
n_an_out = NetFlow(Anode.Outlet);
n_cath_in = NetFlow(Cathode.Inlet);
n_cath_out = NetFlow(Cathode.Outlet);
switch FC.FCtype
    case 'SOFC'
        if FC.ClosedCathode
            AvgX.O2 = ones(FC.nodes,1);
        else
            AvgX.O2 = (Cathode.Outlet.O2+Cathode.Inlet.O2)./(n_cath_in+n_cath_out);
        end
    case 'MCFC'
        if FC.ClosedCathode
            AvgX.O2 = ones(FC.nodes,1)/3;
            AvgX.CO2c = 2*ones(FC.nodes,1)/3;
        else
            AvgX.O2 = (Cathode.Outlet.O2+Cathode.Inlet.O2)./(n_cath_in+n_cath_out);
            AvgX.CO2c = (Cathode.Outlet.CO2+Cathode.Inlet.CO2)./(n_cath_in+n_cath_out);
        end
        AvgX.CO2a = (Anode.Outlet.CO2+Anode.Inlet.CO2)./(n_an_in+n_an_out);
end
AvgX.H2 = (Anode.Outlet.H2+Anode.Inlet.H2)./(n_an_in+n_an_out);
AvgX.H2O = (Anode.Outlet.H2O+Anode.Inlet.H2O)./(n_an_in+n_an_out);

k = FC.FuelFlowDir(:,1);
if min(Anode.Inlet.H2(k))==0 && FC.R_CH4(k)>0 %gives some reformed methane as anode inlet
    AvgX.H2 = (Anode.Outlet.H2 + 0.5*(3*FC.R_CH4(k)+FC.R_WGS(k)))./(n_an_in(k)+n_an_out(k));
    AvgX.H2O = (Anode.Outlet.H2O - 0.5*(FC.R_CH4(k)-FC.R_WGS(k)) + Anode.Inlet.H2O(k))./(n_an_in(k)+n_an_out(k));
end

%% Calculate local voltages
normTemp = FC.TpenAvg + (FC.T.Elec-mean(FC.T.Elec)); %assume you will get to the desired temperature (this avoids oscilations in voltage and helps convergence
FuelCellNernst(FC.Current,normTemp,FC.AirPinit,AvgX,FC);

function Inlet = InletFlow(FC,Inlet) %only used 1st time through initialization (before we know what is connected to inlet
% Anode
for i = 1:1:length(FC.AnSpec)
    if isfield(FC.Fuel,FC.AnSpec{i})
        Inlet.AnodeIn.(FC.AnSpec{i}) = FC.Fuel.(FC.AnSpec{i})*FC.FuelFlowInit;%flow rate of every species entering the anode (or reformer if there is one)
    else Inlet.AnodeIn.(FC.AnSpec{i}) = 0;
    end
end
%Cathode
Inlet.CathodeIn.T = FC.StackCathTin;
switch FC.FCtype
    case 'SOFC'
        if FC.ClosedCathode
            Inlet.CathodeIn.O2 = FC.Oxidant.O2*FC.AirFlow;
        else
            for i = 1:1:length(FC.CathSpec)
                if isfield(FC.Oxidant,FC.CathSpec{i})
                    Inlet.CathodeIn.(FC.CathSpec{i}) = FC.Oxidant.(FC.CathSpec{i})*FC.AirFlow;
                else Inlet.CathodeIn.(FC.CathSpec{i}) = 0;
                end
            end
        end
    case 'MCFC' %recalculate cathode inlet species for MCFC (this is an estimate assuming the 100% of non-recirculated anode gas is oxidized and fed to the cathode)
        Inlet.CathodeIn.CO2 = (Inlet.FuelMix.CH4+Inlet.FuelMix.CO+Inlet.FuelMix.CO2) + sum(FC.Current)/(2*F*1000)*FC.Cells;
        Inlet.CathodeIn.H2O = (4*Inlet.FuelMix.CH4+Inlet.FuelMix.CO+Inlet.FuelMix.H2+Inlet.FuelMix.H2O);
        nonCO2_H2O = (FC.AirFlow - Inlet.CathodeIn.CO2 - Inlet.CathodeIn.H2O);
        for i = 1:1:length(FC.CathSpec)
            if isfield(FC.Oxidant,FC.CathSpec{i})
                if strcmp(FC.CathSpec{i},'CO2')||strcmp(FC.CathSpec{i},'H2O')
                    Inlet.CathodeIn.(FC.CathSpec{i}) = Inlet.CathodeIn.(FC.CathSpec{i}) + FC.Oxidant.(FC.CathSpec{i})*nonCO2_H2O;
                else
                    Inlet.CathodeIn.(FC.CathSpec{i}) = FC.Oxidant.(FC.CathSpec{i})*nonCO2_H2O;
                end
            else Inlet.CathodeIn.(FC.CathSpec{i}) = 0;
            end
        end
end

function FC = Set_IC(FC,Anode,Cathode,Reformer)
global Ru
Cp.cath = SpecHeat(Cathode.Inlet);
Cp.an = SpecHeat(Anode.Outlet);
if ~isempty(Reformer)
    Cp.ref = SpecHeat(Reformer.Outlet);
end
switch FC.Reformer
    case 'internal'
        NumOfStates = (6 + 2*length(FC.AnSpec) + length(FC.CathSpec) + 1)*FC.nodes + 2; % 6 temperatures, anode & cathode & reformer & current at each node and 2 states for anode/cathode pressure
  case {'direct';'external';'adiabatic'}
        NumOfStates = (5 + length(FC.AnSpec) + length(FC.CathSpec) + 1)*FC.nodes +2; % 5 temperatures, anode & cathode & current at each node and 2 states for anode/cathode pressure
end
FC.IC = ones(NumOfStates,1); %
FC.tC = FC.IC; % time constant for derivative dY
FC.Scale = FC.IC;
switch FC.Reformer
    case 'internal'
        n = 6*FC.nodes;
        FC.tC(5*FC.nodes+1:6*FC.nodes) = (FC.Vol_Reform*Cp.ref*FC.FuelPinit./(Ru*FC.T.Reform));
    case {'direct';'external';'adiabatic'}
        n = 5*FC.nodes;
end
FC.Scale = FC.Tstates(1:n);%temperature (K)

FC.tC(1:FC.nodes) = (FC.Mass_Cath_plate*FC.C_Cath_plate);
FC.tC(1+FC.nodes:2*FC.nodes) = (FC.Vol_Cathode*Cp.cath*FC.AirPinit./(Ru*FC.T.Cath));
FC.tC(2*FC.nodes+1:3*FC.nodes) = (FC.Vol_Elec*FC.Density_Elec*FC.C_Elec);
FC.tC(3*FC.nodes+1:4*FC.nodes) = (FC.Vol_Anode*Cp.an*FC.FuelPinit./(Ru*FC.T.Anode));
FC.tC(4*FC.nodes+1:5*FC.nodes) = (FC.Mass_An_plate*FC.C_An_plate);
FC.tC(1:n) = FC.tC(1:n)-diag(FC.HTmatrix); %this accounts for the change in HT as temperature of the control volume changes. The change in HT helps balance the energy equation more than the change in enthalpy leaving.

for i = 1:1:length(FC.CathSpec)
    FC.tC(n+1:n+FC.nodes) = (FC.Vol_Cathode*FC.AirPinit)./(FC.T.Cath*Ru);  % cathode 
    if any(Cathode.Outlet.(FC.CathSpec{i})==0)
        Flow = NetFlow(Cathode.Outlet);
        FC.IC(n+1:n+FC.nodes) = Cathode.Outlet.(FC.CathSpec{i})./Flow;
        FC.Scale(n+1:n+FC.nodes) = Flow; n = n+FC.nodes; %cathode flows
    else
        FC.Scale(n+1:n+FC.nodes) = Cathode.Outlet.(FC.CathSpec{i}); n = n+FC.nodes; %cathode flows
    end
end

Flow = NetFlow(Anode.Outlet);
for i = 1:1:length(FC.AnSpec)
    X = Anode.Outlet.(FC.AnSpec{i})./Flow;%concentration
    FC.tC(n+1:n+FC.nodes) = (FC.Vol_Anode*FC.FuelPinit)./(FC.T.Anode*Ru); %anode
    if any(X<.01) %concentration less than 1%
        FC.IC(n+1:n+FC.nodes) = X;
        FC.Scale(n+1:n+FC.nodes) = Flow; %anode flow
    else
        FC.Scale(n+1:n+FC.nodes) = Anode.Outlet.(FC.AnSpec{i}); %individual species flow
    end   
    n = n+FC.nodes;
end

switch FC.Reformer
    case 'internal'
        for i = 1:1:length(FC.AnSpec)
            FC.tC(n+1:n+FC.nodes) = (FC.Vol_Reform*FC.FuelPinit)./(FC.T.Reform*Ru); % reformer
            if any(Reformer.Outlet.(FC.AnSpec{i})==0)
                Flow = NetFlow(Reformer.Outlet);
                FC.IC(n+1:n+FC.nodes) = Reformer.Outlet.(FC.AnSpec{i})./Flow;
                FC.Scale(n+1:n+FC.nodes) = Flow; n = n+FC.nodes; %anode flows
            else
                FC.Scale(n+1:n+FC.nodes) = Reformer.Outlet.(FC.AnSpec{i}); n = n+FC.nodes; %reformer flows
            end
        end
end
FC.tC(n+1:n+FC.nodes) = FC.nodes/100;%  %current changing for voltage balance
FC.Scale(n+1:n+FC.nodes) = FC.Current;  n = n+FC.nodes; %current

FC.tC(n+1) = (FC.Vol_Cathode*FC.nodes*FC.Cells);  %pressure
FC.tC(n+2) = (FC.Vol_Anode*FC.nodes*FC.Cells); %pressure

FC.Scale(n+1) = FC.AirPinit;%pressure
FC.Scale(n+2) = FC.FuelPinit;%pressure


function Current = redistributeCurrent(Current,scale,Voltage,localR,SetVoltage)
currentPercOfAvg = Current/sum(Current)*length(Current);
dCurrent = (Voltage-SetVoltage)./localR.*currentPercOfAvg; %change in current to balance voltage
if min(abs(Current./dCurrent))<1
    a = .5*min(abs(Current./dCurrent));%ensure dCurrent is never more than half of a step towards zero
else a = .5;
end
dCurrent = a*dCurrent;
scale2  = (scale*sum(Current))/sum(Current+dCurrent); %change in current to get to new power
newCurrent = (Current+dCurrent)*scale2;%re-distribute current to achieve average voltage, then scale to new total current
Current = newCurrent;

function [FC,FuelMix, Anode,Reformer] = KineticCoef(FC,Inlet,first)
%% find the kinetic coefficient which results in the net reforming determined by equilibrium
%% first find equilibrium at outlet
global Ru F
FuelMix = Inlet.FuelMix;
StackAnIn = Inlet.AnodeIn;
Rnet.CH4 = sum(FC.R_CH4);
RCH4old =0;
switch FC.Reformer
    case {'internal','adiabatic'}
        CH4max = min(FuelMix.CH4,FuelMix.H2O)/FC.Cells*FC.RefSpacing;
        CH4min = -min(FuelMix.CO,((FuelMix.H2 + FuelMix.CO)/3)*3/4)/FC.Cells*FC.RefSpacing;
        X0guessRef = (sum(FC.R_CH4ref)/FC.RefPerc - CH4min)/(CH4max -CH4min);
        X0guessRef = max(min(X0guessRef,(1-1e-5)),1e-5);
        X0guess = max(min(FC.AnPercEquilib,(1-1e-5)),1e-5);
    case 'direct'
        CH4max = min(FuelMix.CH4,FuelMix.H2O);
        CH4min = -min(FuelMix.CO,((FuelMix.H2 + FuelMix.CO - sum(FC.Current)/(2*F*1000))/3)*3/4);
        X0guess = (sum(FC.R_CH4)/FC.AnPercEquilib - CH4min)/(CH4max -CH4min);
        X0guess = max(min(X0guess,(1-1e-5)),1e-5);
    case 'external'
        CH4max = min((1-FC.RefPerc)*FuelMix.CH4,FuelMix.H2O-1.8*FC.RefPerc*FuelMix.CH4);
        CH4min = -min(FuelMix.CO+.2*FC.RefPerc*FuelMix.CH4,((FuelMix.H2 + FuelMix.CO + 4*FC.RefPerc*FuelMix.CH4 - sum(FC.Current)/(2*F*1000))/3)*3/4);
        X0guess = (sum(FC.R_CH4)/FC.AnPercEquilib - CH4min)/(CH4max -CH4min);
        X0guess = max(min(X0guess,(1-1e-5)),1e-5);
end
count = 0;
Tol = 1e-3;
while abs((Rnet.CH4-RCH4old)/Rnet.CH4)>Tol% && (FC.Recirc.Anode>0 || count==0)
    RCH4old = Rnet.CH4;
    AnOutlet.T = mean(FC.T.Anode(FC.FuelFlowDir(:,end)));
    switch FC.Reformer
        case 'external'
            AnInlet.T = FuelMix.T;
            for i = 1:1:length(FC.AnSpec)
                AnInlet.(FC.AnSpec{i}) = FuelMix.(FC.AnSpec{i})/FC.Cells;
            end
            AnInlet.CH4 = (1-FC.RefPerc)*FuelMix.CH4/FC.Cells;
            AnInlet.CO = (FuelMix.CO + 0.2*FC.RefPerc*FuelMix.CH4)/FC.Cells;
            AnInlet.CO2 = (FuelMix.CO2 + 0.8*FC.RefPerc*FuelMix.CH4)/FC.Cells;
            AnInlet.H2 = (FuelMix.H2 + 3.8*FC.RefPerc*FuelMix.CH4)/FC.Cells;
            AnInlet.H2O = (FuelMix.H2O - 1.8*FC.RefPerc*FuelMix.CH4)/FC.Cells;
            Reformer =[];
        case 'internal'
            RefInlet.T = FuelMix.T;
            for i = 1:1:length(FC.AnSpec)
                RefInlet.(FC.AnSpec{i}) = FuelMix.(FC.AnSpec{i})/FC.Cells*FC.RefSpacing;%reduce mass flow
            end
            RefOutlet.T = mean(FC.T.Reform(FC.ReformFlowDir(:,end)));
            [RefOutlet,Rref_net] = equilib2D(RefInlet,RefOutlet.T,FC.FuelPinit,0,FC.FCtype,FC.RefPerc,X0guessRef);
            AnInlet.T = RefOutlet.T;
            for i = 1:1:length(FC.AnSpec)
                AnInlet.(FC.AnSpec{i}) = RefOutlet.(FC.AnSpec{i})/FC.RefSpacing;%reduce mass flow
            end
        case 'direct'
            AnInlet.T = FuelMix.T;
            for i = 1:1:length(FC.AnSpec)
                AnInlet.(FC.AnSpec{i}) = FuelMix.(FC.AnSpec{i})/FC.Cells;%reduce mass flow
            end
            Reformer =[];
        case 'adiabatic' %% find recirculation that achieves desired reformer temp
            Tol = 1e-2;
            Reformer.Outlet.T = FC.ReformT;
            %% fixed recirculation
            Reformer.Inlet = FuelMix;
            [Reformer.Outlet,Rref_net,~,RefPerc] = equilibReform(Reformer.Inlet,FC.FuelPinit,0,Reformer.Outlet.T,1,'Q');
            
            K_CH4eq = (Reformer.Outlet.H2.^3.*Reformer.Outlet.CO)./(Reformer.Outlet.CH4.*Reformer.Outlet.H2O).*(FC.FuelPinit./NetFlow(Reformer.Outlet)).^2;
            K_WGSeq = Reformer.Outlet.CO2.*Reformer.Outlet.H2./(Reformer.Outlet.CO.*Reformer.Outlet.H2O);
            a = 4352.2./Reformer.Outlet.T - 3.99;
            K_WGS = exp(a);% Water gas shift equilibrium constant
            K_CH4 = 2459000*exp(-6.187*a);
            FC.scaleK_CH4 = K_CH4eq./K_CH4;
            FC.scaleK_WGS = K_WGSeq./K_WGS;
            
            FC.ReformT = Reformer.Outlet.T; %update to have a better guess next time
            for i = 1:1:length(FC.AnSpec)
                AnInlet.(FC.AnSpec{i}) = Reformer.Outlet.(FC.AnSpec{i})/FC.Cells;
            end   
    end
    [AnOutlet,Rnet] = equilib2D(AnInlet,AnOutlet.T,FC.FuelPinit,sum(FC.Current)/(2*F*1000),FC.FCtype,FC.AnPercEquilib,X0guess);
    if FC.Recirc.Anode >0 %only first time through
        errorR = 1;
        while abs(errorR)>1e-5 %% loop to find anode recirculation that meets steam2carbon design
            for i = 1:1:length(FC.AnSpec)
                FuelMix.(FC.AnSpec{i}) = StackAnIn.(FC.AnSpec{i}) + FC.Recirc.Anode*FC.Cells*AnOutlet.(FC.AnSpec{i});
            end
            S2Ccheck = FuelMix.H2O/(FuelMix.CH4+.5*FuelMix.CO);
            errorR = (FC.Steam2Carbon-S2Ccheck)/FC.Steam2Carbon;
            FC.Recirc.Anode = FC.Recirc.Anode*(1+.9*errorR);
        end
        %%find resulting temperature of mixture
        errorT = 1;
        [~,Hin] = enthalpy(Inlet.AnodeIn);
        [~,Hout] = enthalpy(AnOutlet);
        Hnet = Hin + FC.Recirc.Anode*Hout*FC.Cells;
        Cp = SpecHeat(AnOutlet);
        NetFlowMix = NetFlow(FuelMix);
        while abs(errorT)>1e-3
            [~,Hmix] = enthalpy(FuelMix);
            errorT = (Hnet-Hmix)/(Cp*NetFlowMix);
            FuelMix.T = FuelMix.T + errorT;
        end  
    end
    count = count+1;
end 
% %% From running equilibrium function I can calculate exponential fit to WGS equilibrium: K_WGS = exp(4189.8./T -3.8242) : slightly different than in Paradis paper

%% Indirect Reformer
switch FC.Reformer
    case 'internal'
        if first
            X_CH4in = RefInlet.CH4/NetFlow(RefInlet);
            X_CH4out = (RefInlet.CH4-Rref_net.CH4)/(NetFlow(RefInlet)+2*Rref_net.CH4);
            lambda = log(X_CH4out/X_CH4in)/(-FC.columns); %exponential decay in CH4
            R_cumulative = zeros(length(FC.ReformFlowDir(:,1)),1);
            XCH4 = zeros(FC.nodes,1);
            for i= 1:1:FC.columns
                k = FC.ReformFlowDir(:,i);
                XCH4(k) = X_CH4in*exp(-i*lambda);
                if i == 1 % R = (in flow - outflow) = (Xin*flowin + Xout*(flowin +2*R)) solved for R
                    R.CH4(k,1) = (RefInlet.CH4/length(k) - XCH4(k).*(NetFlow(RefInlet)/length(k) +2*R_cumulative))./(1+2*XCH4(k));
                    R_cumulative = R_cumulative+R.CH4(k);
                else
                    R.CH4(k,1) = (XCH4(kold) - XCH4(k)).*(NetFlow(RefInlet)/length(k)+2*R_cumulative)./(1+2*XCH4(k));
                    R_cumulative = R_cumulative+R.CH4(k);
                end
                kold = k;
            end
            R.WGS = Rref_net.WGS/Rref_net.CH4*R.CH4; %assume same initial distribution for WGS reaction
        else
            for r = 1:1:FC.rows 
                k_r = FC.ReformFlowDir(r,:);
                R.CH4(k_r,1) = FC.R_CH4ref(k_r)*(Rref_net.CH4/FC.rows)/sum(FC.R_CH4ref(k_r));%make sure the total reforming is correct.
                R.WGS(k_r,1) = FC.R_WGSref(k_r)*(Rref_net.WGS/FC.rows)/sum(FC.R_WGSref(k_r)); %assume same initial distribution for WGS reaction
            end
        end
        [R, Reformer, FC.KineticCoeff3] = FindKineticCoef(RefInlet,FC.T.Reform,R,FC.ReformFlowDir,Rref_net.CH4,zeros(length(R.CH4),1),FC.FuelPinit,FC.FCtype,FC.method,1e-5);
        FC.R_CH4ref = R.CH4;
        FC.R_WGSref = R.WGS;  
    case 'adiabatic'
        %% solve for recircualtion that gives desired reformer Toutlet;
%         [Hin,~] = enthalpy(Inlet.AnodeIn); %total energy (not just sensible)
%         [Rref_net,Rnet,FC,Reformer,AnInlet,AnOutlet] = ReformRecirculation(FuelMix,Reformer,Hin,FC);
        %%---%
        R.CH4ref = Rref_net.CH4;
        R.WGSref = Rref_net.WGS;
        nout = NetFlow(Reformer.Outlet);
        X_CH4 = Reformer.Outlet.CH4./nout*FC.FuelPinit*1000; %partial pressures in Pa
        X_H2O = Reformer.Outlet.H2O./nout*FC.FuelPinit*1000; %partial pressures in Pa
        if strcmp(FC.method,'Achenbach')
            FC.KineticCoeff3 = R.CH4ref/(X_CH4.*exp(-8.2e4./(Ru*Reformer.Outlet.T)));
        elseif strcmp(FC.method,'Leinfelder')
            FC.KineticCoeff3 = R.CH4ref/(30.8e10*X_CH4.*X_H2O.*exp(-2.05e5./(Ru*Reformer.Outlet.T)));
        elseif strcmp(FC.method,'Drescher')
            FC.KineticCoeff3 = R.CH4ref/(288.52*X_CH4.*X_H2O.*exp(-1.1e4./(Ru*Reformer.Outlet.T))./(1+16*X_CH4+0.143*X_H2O.*exp(3.9e4./(Ru*Reformer.Outlet.T))));
        end
        FC.R_CH4ref = R.CH4ref;
        FC.R_WGSref = R.WGSref;
end
%% Anode Reforming
if first
    X_CH4in = AnInlet.CH4/NetFlow(AnInlet);
    X_CH4out = (AnInlet.CH4-Rnet.CH4)/(NetFlow(AnInlet)+2*Rnet.CH4);
    lambda = log(X_CH4out/X_CH4in)/(-FC.columns); %exponential decay in CH4
    R_cumulative = zeros(length(FC.FuelFlowDir(:,1)),1);
    XCH4 = zeros(FC.nodes,1);
    for i= 1:1:FC.columns
        k = FC.FuelFlowDir(:,i);
        XCH4(k) = X_CH4in*exp(-i*lambda);
        if i == 1 % R = (in flow - outflow) = (Xin*flowin + Xout*(flowin +2*R)) solved for R
            R.CH4(k,1) = (AnInlet.CH4/length(k) - XCH4(k).*(NetFlow(AnInlet)/length(k) +2*R_cumulative))./(1+2*XCH4(k));
            R_cumulative = R_cumulative+R.CH4(k);
        else
            R.CH4(k,1) = (XCH4(kold) - XCH4(k)).*(NetFlow(AnInlet)/length(k) + 2*R_cumulative)./(1+2*XCH4(k));
            R_cumulative = R_cumulative+R.CH4(k);
        end
        kold = k;
    end
    R.WGS = Rnet.WGS/Rnet.CH4*R.CH4; %assume same initial distribution for WGS reaction
else
    R.CH4 = FC.R_CH4*Rnet.CH4/sum(FC.R_CH4); %keep same disribution as last time, but make the sum equal to the global calculation
    R.WGS = FC.R_WGS*Rnet.WGS/sum(FC.R_WGS); %assume same initial distribution for WGS reaction
end
[R, Anode, FC.KineticCoeff1] = FindKineticCoef(AnInlet,FC.T.Anode,R,FC.FuelFlowDir,Rnet.CH4,FC.Current,FC.FuelPinit,FC.FCtype,FC.method,1e-5);
FC.R_CH4 = R.CH4;
FC.R_WGS = R.WGS;

function [R, Flow, KinCoef] = FindKineticCoef(Inlet,T_Out,R, Dir, referenceR_CH4,Current,Pressure,Type,method,Tol)
global Ru F
specInterest = {'CH4','CO','CO2','H2','H2O'};
h = enthalpy(T_Out,specInterest);
s = entropy(T_Out,specInterest);

Flow = In2Out(T_Out,Inlet,Dir,Type,1,Current,R,'anode');
nout = NetFlow(Flow.Outlet);
X_CH4 = Flow.Outlet.CH4./nout*Pressure*1000; %partial pressures in Pa
X_H2O = Flow.Outlet.H2O./nout*Pressure*1000; %partial pressures in Pa
r = length(Dir(:,1));%rows
if strcmp(method,'Achenbach')
    KinCoef = R.CH4./(X_CH4.*exp(-8.2e4./(Ru*T_Out))); %best guess of KinCoef
elseif strcmp(method,'Leinfelder')
    KinCoef = R.CH4./(X_CH4.*30.8e10*X_H2O.*exp(-2.05e5./(Ru*T_Out)));
elseif strcmp(method,'Drescher')
    KinCoef = R.CH4./(X_CH4.*288.52*X_H2O.*exp(-1.1e4./(Ru*T_Out(Dir(:,1))))/(1+16*X_CH4+0.143*X_H2O.*exp(3.9e4./(Ru*T_Out))));
end
KinCoef = sum(KinCoef.*R.CH4/referenceR_CH4);

spec = fieldnames(Inlet);
spec = spec(~strcmp('T',spec));
XnP = 0;
for i = 1:1:length(spec)
    if ~ismember(spec{i},specInterest)
        XnP = XnP + Inlet.(spec{i})/r;
    end
end
count =0;
error = 1;   
while abs(error)>Tol %iterate to converge on a kinetic coefficients (if less CH4 in exhaust than equilibrium, smaller coefficient)
    R_CH4a = loopConverge(Flow,R,T_Out,Pressure,KinCoef,Dir,method);
    eK = 1e-6*KinCoef;
    R_CH4b = loopConverge(Flow,R,T_Out,Pressure,KinCoef+eK,Dir,method);
    dR_dK = (sum(R_CH4b) - sum(R_CH4a))/eK;
    error = (referenceR_CH4/sum(R_CH4a)-1);
    KinCoef = KinCoef + (referenceR_CH4 - sum(R_CH4a))/dR_dK; %adjust kinetic coefficient
    R.CH4 = loopConverge(Flow,R,T_Out,Pressure,KinCoef,Dir,method);

    R_CH4max = (1-1e-10)*min(Flow.Inlet.CH4,Flow.Inlet.H2O);
    R.CH4 = min(R.CH4,R_CH4max);
    for j= 1:1:length(Dir(1,:))
        k = Dir(:,j);
        for i = 1:1:length(specInterest)
            if j>1
                Flow.Inlet.(specInterest{i})(k) = Flow.Outlet.(specInterest{i})(kOld);
            else
                Flow.Inlet.(specInterest{i})(k) = Inlet.(specInterest{i})/r;
            end
            g0.(specInterest{i}) = h.(specInterest{i})(k)-T_Out(k).*s.(specInterest{i})(k);
        end
        %% inlet to outlet (without WGS)
        X.CH4 = Flow.Inlet.CH4(k) - R.CH4(k);
        X.CO = Flow.Inlet.CO(k) + R.CH4(k);
        X.CO2 = Flow.Inlet.CO2(k);
        X.H2 = Flow.Inlet.H2(k) + 3*R.CH4(k) - Current(k)/(2*F*1000);%hydrogen consumed
        X.H2O = Flow.Inlet.H2O(k) - R.CH4(k) + Current(k)/(2*F*1000);% water produced
        if strcmp(Type,'MCFC')
            X.CO2 = X.CO2 + Current(k)/(2*F*1000); % CO2 brought over
        end
        R_COmin = -min(Flow.Inlet.CO2(k),Flow.Inlet.H2(k) + 3*R.CH4(k) - Current(k)/(2*F*1000));
        R_COmax = min((Flow.Inlet.H2O(k)-R.CH4(k)),(Flow.Inlet.CO(k)+R.CH4(k)));
        y0 = max(0+1e-5,min(1-1e-5,(R.WGS(k)-R_COmin)./(R_COmax-R_COmin)));
        y0 = Newton1D(y0,X,R_COmin,R_COmax,T_Out(k),Pressure,g0,XnP,1e-6);

        R.WGS(k) = R_COmin+y0.*(R_COmax-R_COmin);
        Flow.Outlet.CH4(k) = Flow.Inlet.CH4(k)-R.CH4(k);
        Flow.Outlet.CO(k) = Flow.Inlet.CO(k)+R.CH4(k)-R.WGS(k);
        Flow.Outlet.CO2(k) = Flow.Inlet.CO2(k)+R.WGS(k);
        Flow.Outlet.H2(k) = Flow.Inlet.H2(k)+3*R.CH4(k)+R.WGS(k) - Current(k)/(2*F*1000);
        Flow.Outlet.H2O(k) = Flow.Inlet.H2O(k)-R.CH4(k)-R.WGS(k) + Current(k)/(2*F*1000);
        kOld = k;
    end
    count = count+1;
end
% disp(strcat('FindKineticCoef count is:',num2str(count)))

function R_CH4 = loopConverge(Flow,R,T_Out,Pressure,KinCoef,Dir,method)
global Ru
%% find new reforming reaction rates
k = Dir(:,1);
n_in = NetFlow(Flow.Inlet);
n_in = n_in(k);
H2O_in = Flow.Inlet.H2O(k);
X_CH4in = Flow.Inlet.CH4(k)./n_in;
nout = NetFlow(Flow.Outlet);
X_CH4 = Flow.Outlet.CH4./nout;%initial guess of X_CH4
for j= 1:1:length(Dir(1,:))
    k = Dir(:,j);
    if strcmp(method,'Achenbach')
        C = exp(-8.2e4./(Ru*T_Out(k)));
    elseif strcmp(method,'Leinfelder')
        X_H2O = (H2O_in - R.CH4(k) - R.WGS(k))./(n_in+2*R.CH4(k))*Pressure*1000; %partial pressures in Pa
        C = 30.8e10*X_H2O.*exp(-2.05e5./(Ru*T_Out(k)));
    elseif strcmp(method,'Drescher')
        X_H2O = (H2O_in - R.CH4(k) - R.WGS(k))./(n_in+2*R.CH4(k))*Pressure*1000; %partial pressures in Pa
        C = 288.52*X_H2O.*exp(-1.1e4./(Ru*T_Out))./(1+16*X_CH4(k)*Pressure*1000+0.143*X_H2O.*exp(3.9e4./(Ru*T_Out(k))));
    end

    % use newton method to find X_CH4_out that makes R.CH4 so that R = X_CH4 in *flow in - X_CH4out*flow out
    dX1 = (X_CH4in.*n_in - KinCoef.*X_CH4(k)*Pressure*1000.*C)./(n_in + 2*KinCoef*X_CH4(k)*Pressure*1000.*C) - X_CH4(k);
    eX = 1e-5*X_CH4in;
    dX2 = (X_CH4in.*n_in - KinCoef.*(X_CH4(k)+eX)*Pressure*1000.*C)./(n_in + 2*KinCoef*(X_CH4(k)+eX)*Pressure*1000.*C) - (X_CH4(k)+eX);
    m = (dX2-dX1)./eX;
    b = dX1 - m.*X_CH4(k);
    X_CH4(k) = - b./m; 
    
    R.CH4(k) = KinCoef*X_CH4(k)*Pressure*1000.*C;
    %inlet to the next column
    X_CH4in = (X_CH4in.*n_in - R.CH4(k))./(n_in + 2*R.CH4(k));
    H2O_in = H2O_in - R.CH4(k) - R.WGS(k);
    n_in = n_in + 2*R.CH4(k);
%     errorJ = (X_CH4(k) - X_CH4in)./X_CH4in; %remaining error after 1 newton step
end
R_CH4 = R.CH4;

function y0 = Newton1D(y0,X,R_COmin,R_COmax,T,P,g0,XnP,Tol)
if any((X.H2 + R_COmax)<0)
    y0 = ones(length(T),1);
else
    dY = zeros(length(T),1);
    scale = dY;
    error = dY + max(10*Tol,1e-6);
    count = 0;
    while max(error)>Tol && count<12
        e_y = max(.001*error,1e-6);
        if any(y0+2*e_y>=1)
            e_y(y0+2*e_y>=1) = -.1*abs(e_y(y0+2*e_y>=1));
        end
        WGS = y0.*R_COmax + (1-y0).*R_COmin;
        WGS_plus = (y0+e_y).*R_COmax + (1-(y0+e_y)).*R_COmin;
        WGS_minus = (y0-e_y).*R_COmax + (1-(y0-e_y)).*R_COmin;
        G11 = GibbVal(X,WGS,T,P,g0,XnP);
        a = abs(G11);
        G12 = GibbVal(X,WGS_plus,T,P,g0,XnP)./a;
        G13 = GibbVal(X,WGS_minus,T,P,g0,XnP)./a;
        G11 = G11./a;

        dGdy1 = (G12-G11)./e_y;
        dGdy2 = (G11-G13)./e_y;
        dGdydy = (dGdy1-dGdy2)./e_y;
        dY(dGdydy==0) =0;
        dY(dGdydy~=0) = -dGdy1./dGdydy;
        b = y0 + dY;
        scale(:) = 1;
        if any(b>1) || any(b<0) %take a smaller step if iteration takes it beyond 0 or 1
            scale(b>1) = .75*(1-y0(b>1))./dY(b>1);
            scale(b<0) = .75*y0(b<0)./(-dY(b<0));
        end
        y0 = y0 +dY.*scale; 
        error = abs(dY);
        count = count+1;
    end
end
% disp(strcat('Newton1D count is:',num2str(count)))

function G = GibbVal(X,WGS,T,P,g0,XnP)
global Ru
X.CO = X.CO - WGS;
X.CO2 = X.CO2 + WGS;
X.H2 = X.H2  + WGS;
X.H2O = X.H2O - WGS;
G = 0;
spec = fieldnames(g0);
sumX = 0;
for i = 1:1:length(spec)
    sumX = sumX + X.(spec{i});
end
sumX = sumX + XnP;
for i = 1:1:length(spec)
    G = G+X.(spec{i}).*(g0.(spec{i})./(Ru*T)+log(X.(spec{i})./sumX));
    %     G = G+X.(spec{i})*(g0.(spec{i})/(Ru*T)+log(X.(spec{i})/sumX*P));
    %     G = G+X.(spec{i})*(g0.(spec{i}) + Ru*T*log(P)+Ru*T*log(X.(spec{i})/sumX));
end

function Flow = In2Out(T,Inlet,Dir, Type,Cells,Current,R,c_or_a) 
%flow of species out after steam reformation & water gas shift
%(does not normalize back to 1) this is important for finding G correctly)
global F
k = Dir(:,1);
r = length(k);
Flow.Outlet.T = T;
spec = fieldnames(Inlet);
for i = 1:1:length(spec)
    if strcmp(spec{i},'T')
        Flow.Inlet.T(k,1) = Inlet.T;
    else
        if strcmp(c_or_a,'anode') %already split flow into cells when calculating FuelMix
            Flow.Inlet.(spec{i})(k,1) = Inlet.(spec{i})/r;
        elseif strcmp(c_or_a,'cathode') %converting from CathodeInlet to flow per cell per row
            Flow.Inlet.(spec{i})(k,1) = Inlet.(spec{i})/Cells/r;
        end
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
    
    for i = 1:1:length(spec)
        if strcmp(c_or_a,'anode') %anode side (H2 consumption)
            if strcmp(spec{i},'CO2')
                switch Type
                    case 'SOFC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k)+R.WGS(k); %CO2 flow
                    case 'MCFC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k)+R.WGS(k) + Current(k)/(2*F*1000); %CO2 flow
                end
            elseif strcmp(spec{i},'CH4')
                Flow.Outlet.CH4(k,1) = Flow.Inlet.CH4(k)-R.CH4(k);%CH4 flow
            elseif strcmp(spec{i},'CO')
                Flow.Outlet.CO(k,1) = Flow.Inlet.CO(k)+R.CH4(k)-R.WGS(k); %CO flow
            elseif strcmp(spec{i},'H2')
                Flow.Outlet.H2(k,1) = Flow.Inlet.H2(k)+3*R.CH4(k)+R.WGS(k) - Current(k)/(2*F*1000); %H2 flow
            elseif strcmp(spec{i},'H2O')
                Flow.Outlet.H2O(k,1) = Flow.Inlet.H2O(k)-R.CH4(k)-R.WGS(k) + Current(k)/(2*F*1000); %H2O flow
            elseif ~strcmp(spec{i},'T')
                Flow.Outlet.(spec{i})(k,1) = Flow.Inlet.(spec{i})(k);
            end
        elseif strcmp(c_or_a,'cathode') %cathode side (O2 consumption)
            if strcmp(spec{i},'CO2')
                switch Type
                    case 'SOFC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k); %CO2 flow
                    case 'MCFC'
                        Flow.Outlet.CO2(k,1) = Flow.Inlet.CO2(k) - 2*Current(k)/(4*F*1000); %CO2 flow
                end
            elseif strcmp(spec{i},'O2')
                Flow.Outlet.O2(k,1) = Flow.Inlet.O2(k) - Current(k)/(4*F*1000); %O2 flow
            elseif ~strcmp(spec{i},'T')
                Flow.Outlet.(spec{i})(k,1) = Flow.Inlet.(spec{i})(k);
            end
        end
        if abs(Flow.Outlet.(spec{i})(k,1))<1e-18 %zero
            Flow.Outlet.(spec{i})(k,1) = 0;
        end
    end
end

function Out = MergeFlows(Flow, k, scale)
spec = fieldnames(Flow);
for i = 1:1:length(spec)
    if strcmp(spec{i},'T')
        Out.T  = mean(Flow.T(k)); %temperature 
    else
        Out.(spec{i}) = max(0,sum(Flow.(spec{i})(k))*scale);%avoid sending negative outflows to other blocks
    end
end

function [Rref_net,Rnet,FC,Reformer,AnInlet,AnOutlet] = ReformRecirculation(FuelMix,Reformer,Hin,FC)
% this function works if we want to find the recirculation which achieves a
% desired temperature in an adibatic reformer (requires more than standard iterations)
% It is commented out in line 873
errorR = 100;
errorOld=0;
count2 = 0;
Rnet.CH4 = [];
X0guess = .99;
while abs(errorR)>1
    Reformer.Inlet = FuelMix;
    [Reformer.Outlet,Rref_net,~,RefPerc] = equilibReform(Reformer.Inlet,FC.FuelPinit,0,Reformer.Outlet.T,1,'Q');
    for i = 1:1:length(FC.AnSpec)
        AnInlet.(FC.AnSpec{i}) = Reformer.Outlet.(FC.AnSpec{i})/FC.Cells;
    end
    if~isempty(Rnet.CH4)
        CH4max = min(AnInlet.CH4,AnInlet.H2O);
        CH4min = -min(AnInlet.CO,((AnInlet.H2 + AnInlet.CO - sum(FC.Current)/(2*F*1000))/3)*3/4);
        X0guess = (Rnet.CH4 - CH4min)/(CH4max -CH4min); 
    end
    [AnOutlet,Rnet] = equilib2D(AnInlet,AnOutlet.T,FC.FuelPinit,sum(FC.Current)/(2*F*1000),FC.FCtype,FC.AnPercEquilib,X0guess);
    if FC.ReformT > Reformer.Outlet.T
        errorR = log(1+FC.ReformT - Reformer.Outlet.T);
    else errorR = -real(log(-1+FC.ReformT - Reformer.Outlet.T));
    end
    %note: total energy of the fuel mix needs to be total
    %energy of the reformer outlet (adiabatic)
    Rold = FC.Recirc.Anode;
    if count2>=first
        if abs(errorR)>2
            FC.Recirc.Anode = max(.3,min(.99,FC.Recirc.Anode + errorR/100 - errorOld/300));
        else FC.Recirc.Anode = max(.3,min(.99,FC.Recirc.Anode + errorR/200 - errorOld/600));
        end
    end
    dr = FC.Recirc.Anode - Rold;
    nScale = (1-Rold)/((1-Rold)-dr);
    errorOld = errorR;
    %% find new fuel mix 
    for i = 1:1:length(FC.AnSpec)
        FuelMix.(FC.AnSpec{i}) = (StackAnIn.(FC.AnSpec{i}) + FC.Recirc.Anode*AnOutlet.(FC.AnSpec{i})*FC.Cells*nScale);
    end
    FC.Steam2Carbon = FuelMix.H2O/(FuelMix.CH4+.5*FuelMix.CO);
    if FC.Steam2Carbon<2
        FC.Recirc.Anode = FC.Recirc.Anode*2/FC.Steam2Carbon;
        for i = 1:1:length(FC.AnSpec)
            FuelMix.(FC.AnSpec{i}) = StackAnIn.(FC.AnSpec{i}) + FC.Recirc.Anode*AnOutlet.(FC.AnSpec{i})*FC.Cells*nScale;
        end
        FC.Steam2Carbon = FuelMix.H2O/(FuelMix.CH4+.5*FuelMix.CO);
    end
    %%find resulting temperature of mixture
    errorT = 100;
    Hnet = (Hin + FC.Recirc.Anode*enthalpy(AnOutlet)*FC.Cells*nScale);
    Cp = SpecHeat(AnOutlet);
    NetFlowMix = NetFlow(FuelMix);
    while abs(errorT)>1
        errorT = (Hnet-enthalpy(FuelMix))/(Cp*NetFlowMix);
        FuelMix.T = FuelMix.T + errorT;
    end    
    count2 = count2+1;
end

function [T,HTmatrix]= SteadyFCtemps(FC,Anode,Cathode,Reformer,Inlet)
%Solve problem of form xdot = Ax-b for xdot =0.
%final solution is x = A\b;
%states represent  heat transfer into each node/layer, the temperatures of each node/layers, and inlet cathode and anode temperatures, Qerror term associated with the small error in air flow rate so that the deltaT and Tavg constraints can both be satisfied
%like the heat exchanger this averages the inlet and exit temperature for the gas streams, and assumes the solid temperature states correspond to the average for the node

global F
a = .5; %weighting of previous node temperature on convection HT calculation
[h,h_s] = enthalpy(FC.T.Elec,{'H2','H2O','O2','CO','CO2','CH4'});
h_rxn1 = h.CO+3*h.H2-h.CH4-h.H2O;
h_rxn2 = h.CO2+h.H2-h.CO-h.H2O;
h_rxn3 = h.H2+.5*h.O2-h.H2O;
Qgen = (FC.Current/(2*F).*h_rxn3-FC.Current*FC.Voltage)/1000;%kW of heat generated by electrochemistry (per node & per cell)
Qdirect = h_rxn1.*FC.R_CH4+h_rxn2.*FC.R_WGS; %kW of cooling per cell;
switch FC.Reformer
    case 'internal'
        Qindirect = h_rxn1.*FC.R_CH4ref+h_rxn2.*FC.R_WGSref; %kW of cooling per cell
end
%ion transport across membrane (sensible enthalpy
switch FC.FCtype
    case 'SOFC'
        Qion = FC.Current/(4*F*1000).*h_s.O2; %O2 ion crossing over (kW)
    case 'MCFC'
        Qion = FC.Current/(4*F*1000).*h_s.O2 + FC.Current/(2*F*1000).*h_s.CO2;% O2 & CO2 ion crossing over
end

nodes = FC.nodes;
Aflow = mean(NetFlow(Anode.Outlet));
Cflow = mean(NetFlow(Cathode.Outlet));
switch FC.Reformer
    case 'internal'
        s=6;
        Cells2Reformer=1; % set equal to 1 to model the cell next to the reformer
        C8 = (FC.h_Ref_Channel*FC.A_Node)/Cells2Reformer/1000;
        Cp_reform = mean(SpecHeat(Reformer.Outlet));
        Rflow = mean(NetFlow(Reformer.Outlet))/FC.RefSpacing;
    case {'adiabatic';'direct';'external'}%  no direct heat transfer
        s=5;
end
states = 2*s*nodes+2;

A = zeros(states+1,states);
b = zeros(states+1,1);

Cp_cath = mean(SpecHeat(Cathode.Inlet));
Cp_anode = mean(SpecHeat(Anode.Inlet));

%all heat transfer coefficients converted to kW/K: thus Q = C*(T1-T2) is in kW
C1 = (FC.k_Cath_plate*FC.L_node*FC.W_node/FC.t_Cath_plate)/1000;
C2 = (FC.h_Cathode_Sep*FC.A_Cathode_Sep)/1000;
C3 = (FC.h_Cathode_Elec*FC.A_Cathode_Elec)/1000;
C4 = (FC.k_Cath_plate*FC.A_Cath_plate_Elec_Cond/FC.L_Cath_plate_Elec)/1000;
C5 = (FC.k_An_plate*FC.A_An_plate_Elec_Cond/FC.L_An_plate_Elec)/1000;
C6 = (FC.h_Anode_Elec*FC.A_Anode_Elec)/1000;
C7 = (FC.h_Anode_Sep*FC.A_Anode_Sep)/1000;

%horizontal
H1_pn = (FC.k_Cath_plate*FC.A_Cath_plate_Heat_Cond/(FC.L_node/2))/1000; %heat transfer coefficient between previous and next node of oxidant plate
H1_lr  = (FC.k_Cath_plate*FC.A_Cath_plate_Heat_Cond/(FC.W_node/2))/1000; %heat transfer coefficient between left and right adjacent nodes of oxidant plate
H2_pn  = (FC.k_Elec*FC.A_Elec_Heat_Cond/(FC.L_node/2))/1000; %heat transfer coefficient between previous and next node of fuel cell electrolyte assembly
H2_lr  = (FC.k_Elec*FC.A_Elec_Heat_Cond/(FC.W_node/2))/1000; %heat transfer coefficient between left and right adjacent nodes  of fuel cell electrolyte assembly
H3_pn  = (FC.k_An_plate*FC.A_An_plate_Heat_Cond/(FC.L_node/2))/1000; %heat transfer coefficient between previous and next node of fuel plate
H3_lr  = (FC.k_An_plate*FC.A_An_plate_Heat_Cond/(FC.W_node/2))/1000; %heat transfer coefficient between left and right adjacent nodes of fuel plate

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

    [i,j] = find(FC.AirFlowDir==k);
    if j==1 %first column averaged with inlet temperature
        A(k,2*s*nodes+1) = a*C2;
        A(k+nodes,2*s*nodes+1) = -a*(C2+C3);
        A(k+2*nodes,2*s*nodes+1) = a*C3;
    else % other columns averaged with previous one
        k2 = FC.AirFlowDir(i,j-1);
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

    [i,j] = find(FC.FuelFlowDir==k);
    if j==1 %first column averaged with inlet temperature
        if s==6 %pull last temperature from reformer
            k2 = FC.ReformFlowDir(i,end);
            A(k+2*nodes,k2+(s+5)*nodes) = a*C6;
            A(k+3*nodes,k2+(s+5)*nodes) = -a*(C6+C7);
            A(k+4*nodes,k2+(s+5)*nodes) = a*C7;
        else
            A(k+2*nodes,2*s*nodes+2) = a*C6;
            A(k+3*nodes,2*s*nodes+2) = -a*(C6+C7);
            A(k+4*nodes,2*s*nodes+2) = a*C7;
        end
    else % other columns averaged with previous one
        k2 = FC.FuelFlowDir(i,j-1);
        A(k+2*nodes,k2+(s+3)*nodes) = a*C6;
        A(k+3*nodes,k2+(s+3)*nodes) = -a*(C6+C7);
        A(k+4*nodes,k2+(s+3)*nodes) = a*C7;
    end

    %QT5 : heat transfer into fuel plate
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

    if s==6 %QT6 : heat transfer into reformer gas
       A(k+5*nodes,k+5*nodes) = -1;
       A(k+5*nodes,k+(s+5)*nodes) = -C8; 
       A(k+5*nodes,k+s*nodes) = C8; 
       A(k+5*nodes,k+(s+4)*nodes) = C8; 

       [i,j] = find(FC.ReformFlowDir==k);
        if j==1 %first column averaged with inlet temperature
            A(k,2*s*nodes+2) = a*C8;
            A(k+4*nodes,2*s*nodes+2) = a*C8;
            A(k+5*nodes,2*s*nodes+2) = -C8;
        else % other columns averaged with previous one
            k2 = FC.ReformFlowDir(i,j-1);
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
    [i,j] = find(FC.AirFlowDir==k);
    if j==1 %first column receives fresh air
        A(k+(s+1)*nodes,2*s*nodes+1) = Cp_cath*Cflow;
    else
        index = FC.AirFlowDir(i,j-1);
        A(k+(s+1)*nodes,index +(s+1)*nodes) = Cp_cath*Cflow;
    end
    b(k+(s+1)*nodes) = Qion(k);

    %Telec: Temperature of electrolyte
    A(k+(s+2)*nodes,k+2*nodes) = 1;
    b(k+(s+2)*nodes) = (-Qgen(k)+ Qdirect(k));

    %Tan: Temperature of anode
    A(k+(s+3)*nodes,k+3*nodes) = 1;
    A(k+(s+3)*nodes,k+(s+3)*nodes) = -Cp_anode*Aflow; 
    [i,j] = find(FC.FuelFlowDir==k);
    if j==1 %first column receives fresh fuel
        if s==5
            A(k+(s+3)*nodes,2*s*nodes+2) = Cp_anode*Aflow; %fresh inlet
        elseif s==6
            index = FC.ReformFlowDir(i,end);
            A(k+(s+3)*nodes,index+(s+5)*nodes) = Cp_anode*Aflow; %reformer out
        end
    else
        index = FC.FuelFlowDir(i,j-1);
        A(k+(s+3)*nodes,index+(s+3)*nodes) = Cp_anode*Aflow;
    end
    b(k+(s+3)*nodes) = -Qion(k);

    %Tfuel: Temperature of fuel plate (net heat transfer into plate = 0)
    A(k+(s+4)*nodes,k+4*nodes) = 1;

    if s==6 %Treform : Temperature of reformer gas
       A(k+(s+5)*nodes,k+5*nodes) = 1;
       A(k+(s+5)*nodes,k+(s+5)*nodes) = -Cp_reform*Rflow;
       [i,j] = find(FC.ReformFlowDir==k);
        if j==1 %first column receives fresh fuel
            A(k+(s+5)*nodes,2*s*nodes+2) = Cp_reform*Rflow;%fresh inlet
        else
            index = FC.ReformFlowDir(i,j-1);
            A(k+(s+5)*nodes,index+(s+5)*nodes) = Cp_reform*Rflow;
        end 
       b(k+(s+5)*nodes) = Qindirect(k);
    end
end

%% left and right, prev and next
prev = FC.HTadjacent(:,1);
next = FC.HTadjacent(:,2);
left = FC.HTadjacent(:,3);
right = FC.HTadjacent(:,4);
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
k1 = FC.AirFlowDir(:,1); %first column
for n =1:1:length(k1)
    k = k1(n);
    A(k,k+(s+1)*nodes) = A(k,k+(s+1)*nodes) + a*C2; %HT to ox plate from cathode
    A(k,2*s*nodes+1) = A(k,2*s*nodes+1) - a*C2; %HT to ox plate from cathode

    A(k+nodes,k+(s+1)*nodes) = A(k+nodes,k+(s+1)*nodes) - a*(C2+C3); %HT from ox plate and electrolyte to cathode
    A(k+nodes,2*s*nodes+1) = A(k+nodes,2*s*nodes+1) + a*(C2+C3);   %HT from ox plate and electrolyte to cathode

    A(k+2*nodes,k+(s+1)*nodes) = A(k+2*nodes,k+(s+1)*nodes) + a*C3;%HT from the cathode to electrolyte
    A(k+2*nodes,2*s*nodes+1) = A(k+2*nodes,2*s*nodes+1) - a*C3; %HT from the cathode to electrolyte
end
    
if s==6 %remove temperature averaging at anode inlet node
    k1 = FC.ReformFlowDir(:,1); %first column
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
    k1 = FC.FuelFlowDir(:,1); %first column
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
b(2*s*nodes+2) = Inlet.FuelMix.T;

%Average electrolyte temp  constraint 
for j = 1:1:nodes %Average all the electrolyte temps to equal Tpenavg
    A(2*s*nodes+3,j+(s+2)*nodes) = 1/nodes;
end 
b(2*s*nodes+3) = FC.TpenAvg;
%cathode dT constraint: T3 - T2 = T2 - T1
c = length(FC.AirFlowDir(1,:));% # of columns
r = length(FC.AirFlowDir(:,1)); % # of rows
for i =1:1:c-1
    if i ==1
        A(2*s*nodes+3+i,2*s*nodes+1) = 1; %cathode inlet (T1)
    else A(2*s*nodes+3+i,FC.AirFlowDir(:,i-1)+(s+1)*nodes) = 1/r;%first column (T1)
    end
    A(2*s*nodes+3+i,FC.AirFlowDir(:,i)+(s+1)*nodes) = -2/r; %middle column (T2)
    A(2*s*nodes+3+i,FC.AirFlowDir(:,i+1)+(s+1)*nodes) = 1/r; %last column (T3)
    b(2*s*nodes+3+i) = 0;
end

x= A\b;
T = x(s*nodes+1:2*s*nodes);


% %% Radiative heat transfer %%
% Ar = 0*A;
% for j = 1:s*nodes
%     Ar(1:s*nodes,s*nodes+j) = FC.RTmatrix(1:s*nodes,j).*T(i).^3;
% end
% A = A + Ar;
% 
% error = true;
% a = 1; %convergence gain factor
% count = 0;
% while error % converges to correct temperatures with radiative heat transfer
%     Told = T;
%     x = A\b;
%     T = a*x(s*nodes+1:2*s*nodes) + (1-a)*Told;
%     if max(abs(T-Told))<.1
%         error = false;
%     end
%     count = count+1;
%     if count>10
%         disp('Slow convergence in radiative heat transfer')
%     end
% end
% %     QT = x(1:s*nodes);

function dY = DynamicFCtemps(t,Y,FC,Anode,Cathode,Reformer,Inlet)
global F Ru
dY = 0*Y;
FuelInlet = Inlet.FuelMix;
nodes = FC.nodes;
if isfield(FC,'tC')
    tC = FC.tC(1:6*nodes);
else
    Cp.cath = 33;% kJ/kmol*K
    Cp.an = 42;% kJ/kmol*K
    tC(1:nodes,1) = (FC.Mass_Cath_plate*FC.C_Cath_plate);
    tC(1+nodes:2*nodes,1) = (FC.Vol_Cathode*Cp.cath*FC.AirPinit./(Ru*FC.T.Cath));
    tC(2*nodes+1:3*nodes,1) = (FC.Vol_Elec*FC.Density_Elec*FC.C_Elec);
    tC(3*nodes+1:4*nodes,1) = (FC.Vol_Anode*Cp.an*FC.FuelPinit./(Ru*FC.T.Anode));
    tC(4*nodes+1:5*nodes,1) = (FC.Mass_An_plate*FC.C_An_plate);
    if isfield(FC,'Vol_Reform')
        Cp.ref = 42;% kJ/kmol*K
        tC(5*nodes+1:6*nodes,1) = (FC.Vol_Reform*Cp.ref*FC.FuelPinit./(Ru*FC.T.Reform));
    end
end

[h,~] = enthalpy(Y(1+2*nodes:3*nodes),{'H2','H2O','O2','CO2'});
Power = FC.Voltage*FC.Current/1000; %cell power in kW
Qgen = FC.Current/(2000*F).*(h.H2+.5*h.O2-h.H2O)-Power;%kW of heat generated by electrochemistry (per node & per cell)
switch FC.FCtype%ion transport across membrane (total enthalpy)
    case 'SOFC'
        Qion = FC.Current/(4000*F).*h.O2; %O2 ion crossing over (kW)
    case 'MCFC'
        Qion = FC.Current/(4000*F).*h.O2 + FC.Current/(2000*F).*h.CO2;% O2 & CO2 ion crossing over
end

QT = FC.HTmatrix*Y;
% QT = QT + FC.RTmatrix*(Y.^4);

Cathode.Outlet.T = Y(nodes+1:2*nodes);
for j = 1:1:length(FC.AirFlowDir(1,:));%1:columns
    k = FC.AirFlowDir(:,j);
    if j~=1
        Cathode.Inlet.T(k,1) = Cathode.Outlet.T(kprev);
    end
    kprev = k;
end

Anode.Outlet.T = Y(3*nodes+1:4*nodes);
for j = 1:1:length(FC.FuelFlowDir(1,:));%1:columns
    k = FC.FuelFlowDir(:,j);
    if j~=1
        Anode.Inlet.T(k,1) = Anode.Outlet.T(kprev);
    end
    kprev = k;
end

%energy flows & sepcific heats
HoutCath = enthalpy(Cathode.Outlet);
HinCath = enthalpy(Cathode.Inlet);
HoutAnode = enthalpy(Anode.Outlet);
HfreshFuel = enthalpy(Inlet.AnodeIn);

if FC.Recirc.Anode>0 % Only during the first run with unhumidified fuel, find fuelmix temperature
    error2 = 1;
    Cp = SpecHeat(Anode.Outlet); 
    Cp = Cp(end);
    Hmix = HfreshFuel + sum(HoutAnode(FC.FuelFlowDir(:,end)))*FC.Recirc.Anode*FC.Cells;
    netflow = NetFlow(FuelInlet);
    while abs(error2)>1e-4
        Hinlet = enthalpy(FuelInlet);                       %Guess of enthalpy based on new cold temperature and cold pressure
        error2 = (Hmix - Hinlet)./(Cp*netflow);                             %Adjusting the error in temperature based on known enthalpy and specific heat of the cold side
        FuelInlet.T = FuelInlet.T + .75*error2;                                   %Subtraction of a portion of the T_error from cold outlet temp to get closer to the actual temp
    end
end

switch FC.Reformer
    case 'internal'
        Reformer.Outlet.T = Y(5*nodes+1:6*nodes);
        k = FC.ReformFlowDir(:,1);
        Reformer.Inlet.T(k,1) = FuelInlet.T;
        for j = 1:1:length(FC.ReformFlowDir(1,:));%1:columns
            k = FC.ReformFlowDir(:,j);
            if j~=1
                Reformer.Inlet.T(k,1) = Reformer.Outlet.T(kprev);
            end
            kprev = k;
        end
        HinReform = enthalpy(Reformer.Inlet);
        HoutReform = enthalpy(Reformer.Outlet);
        
        k = FC.FuelFlowDir(:,1);
        k2 = FC.ReformFlowDir(:,end);
        Anode.Inlet.T(k,1) = Reformer.Outlet.T(k2,1);
    case {'adiabatic';'direct';'external';}
        k = FC.FuelFlowDir(:,1);
        Anode.Inlet.T(k,1) = FuelInlet.T;
end
HinAnode = enthalpy(Anode.Inlet);

if FC.ClosedCathode %%energy balance
    Qimbalance = sum((HinCath(FC.AirFlowDir(:,1))) - sum(HoutCath(FC.AirFlowDir(:,end)))) + sum(HinReform(FC.ReformFlowDir(:,1)))  - sum(HoutAnode(FC.FuelFlowDir(:,end))) - sum(Power);
    Power = Power + Qimbalance*Power./sum(Power);
end

dY(1:nodes)= QT(1:nodes)./tC(1:nodes);  %Ox Sep Plate
dY(1+nodes:2*nodes)= (QT(1+nodes:2*nodes) - Qion + HinCath - HoutCath)./tC(1+nodes:2*nodes); %Cathode
dY(1+2*nodes:3*nodes)= (QT(1+2*nodes:3*nodes)+ Qgen)./tC(2*nodes+1:3*nodes); %Electrolyte Plate
dY(1+3*nodes:4*nodes)= (QT(1+3*nodes:4*nodes) + HinAnode - HoutAnode + Qion - Power - Qgen)./tC(1+3*nodes:4*nodes);  %Anode
dY(1+4*nodes:5*nodes)= QT(1+4*nodes:5*nodes)./tC(4*nodes+1:5*nodes);  %Fuel Sep Plate
switch FC.Reformer
    case 'internal'
        dY(1+5*nodes:6*nodes)= (FC.RefSpacing*QT(1+5*nodes:6*nodes) + HinReform - HoutReform)./tC(1+5*nodes:6*nodes);  %Fuel Reformer Channels
end

function FC = FlowDir(FC)
% Script which orients the nodes relative to the flow directions
% direction = 1: co-flow, direction = 2: counter-flow, direction = 3: cross-flow
nodes = FC.nodes;
columns = FC.columns;
rows = FC.rows;
switch FC.direction
    case 'coflow'
        for j = 1:1:columns
            FC.AirFlowDir(:,j) = (j:columns:nodes)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
        end
    case 'counterflow'
        for j = 1:1:columns
            FC.AirFlowDir(:,j) = (columns-j+1:columns:nodes)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
        end
    case 'crossflow'
        for j = 1:1:rows
            FC.AirFlowDir(:,j) = (1+columns*(j-1):j*columns)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
        end
end
for j = 1:1:columns
    FC.FuelFlowDir(:,j) = (j:columns:nodes)';
    FC.ReformFlowDir(:,j) = (columns-j+1:columns:nodes)';
end
FC.HTadjacent = zeros(nodes,4);
for i = 1:1:nodes
    FC.HTadjacent(i,1) = i-1;%previous node
    FC.HTadjacent(i,2) = i+1;%next node
    FC.HTadjacent(i,3) = i-columns;%node to left
    FC.HTadjacent(i,4) = i+columns;%node to right
end
FC.HTadjacent(1:columns:end,1)=linspace(1,nodes-columns+1,rows);%first node in each row has nothing before it
FC.HTadjacent(columns:columns:end,2)=linspace(columns,nodes,rows)';%last node in each row has nothing after it
FC.HTadjacent(1:columns,3)=linspace(1,columns,columns);%first row has nothing to left
FC.HTadjacent(end-columns+1:end,4)=linspace(nodes-columns+1,nodes,columns);%last row has nothing to right