function FC = InitializeFuelCell(varargin)
%FC model with many states: Temperatures (oxidizer plate, cathode, electrolyte, anode, fuel plate [possibly a reformer]), Cathode species, anode species (can include anode recirculation for humidification, [reformer: indirect, adiabatic or none], Current, cathode pressure, anode pressure
% Five (5) inlets: {'NetCurrent','CathodeIn','AnodeIn','Flow1Pout','Flow2Pout'}
% Seven (7) outlets: {'Flow1Out','Flow2Out','Flow1Pin','Flow2Pin','MeasureCurrent','MeasureTpen','MeasureTflow1','MeasureTflow2'}
% for a fuel cell, Flow 1 is the anode (fuel), and Flow2 is the cathode (cooling/heating air)
% Current is positive
global F Ru
F=96485.339; % %Faraday's constant in Coulomb/mole
Ru = 8.314472; % Universal gas constant in kJ/K*kmol
FC = varargin{1};
if length(varargin)==1 % first initialization
    FC.nodes = FC.rows*FC.columns;
    FC.AnPercEquilib = 1; %CH4 reforming reaches equilibrium at anode exit.
%%%---%%% User defined variables
    %%--Geometric Varibles  %
    FC.t_plate1 = 0.003;                   % [m] thickness of seperator plate
    FC.t_plate1_wall =0.005;                  % [m] Thickness of the channel wall of the Fuel Seperator Plate
    FC.H_plate1channels = 0.002;                       % [m] height of anode channel
    FC.W_plate1channels = 0.005;                       % [m] width of channel                       
    FC.H_plate2channels = 0.002;                      % [m] height of Cathode channel
    FC.W_plate2channels = 0.005;                      % [m] width of Cathode channel
    FC.t_plate2 = 0.003;                     % [m] Thickness of Oxidant Seperator Plate
    FC.t_plate2_wall=0.005;                     % [m] Thickness of the channel wall of the Oxidant Seperator Plate
    FC.Nu_flow1 = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp
    FC.Nu_flow2 = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp

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
    FC.Density_plate1 = 2000;                                % [kg/m3]     density 
    FC.C_plate1 = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    FC.k_plate1 = 5;   %25                                    % [W/(m K)]   conductivity of Fuel Seperator Plate
    %%-----Anode Gas Stream-----------%
    FC.k_flow1 = 259E-3;                                 % (W/m*K) Thermal Conductivity of 50%H2 & 50%H2O at 1000K

    %%-------Cathode Gas Stream---------%
    FC.k_flow2 = 67E-3;                                              % (W/m*K) Thermal Conductivity of air at 1000K

    %%----Cathode half of interconnect plate-------%   
    FC.Density_plate2 = 2000;                                % [kg/m3]     density 
    FC.C_plate2 = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    FC.k_plate2 = 5;%25;                                % [W/(m K)]   conductivity of Fuel Seperator Plate

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
    FC.Dh_Flow1 = 4*(FC.H_plate1channels*FC.W_plate1channels)/(2*(FC.H_plate1channels+FC.W_plate1channels)); %(m) hydraulic diameter of channel
    FC.Dh_Flow2 = 4*(FC.H_plate2channels*FC.W_plate2channels)/(2*(FC.H_plate2channels+FC.W_plate2channels)); %(m) hydraulic diameter of channel
    FC.CH_Flow1 = FC.W_node/(FC.W_plate1channels+FC.t_plate1_wall); %Number of channels in each node of the anode
    FC.CH_Flow2 = FC.W_node/(FC.W_plate2channels+FC.t_plate2_wall); %Number of channels in each node of the cathode 
    %% ---Fuel Seperator plate-------%
    FC.A_plate1_elecCond = FC.t_plate1_wall*FC.L_node*FC.CH_Flow1;   % [m2] Conduction area between the fuel seperator plate and the electrolyte
    FC.A_plate1_heatCond=(FC.H_plate1channels*FC.t_plate1_wall + (FC.W_plate1channels+FC.t_plate1_wall)*FC.t_plate1)*FC.CH_Flow1; %[m^2] conduction area between nodes
    FC.L_plate1_heatCond=FC.H_plate1channels;                                     % [m] Lenght of conduction between the fuel seperator plate and electrolyte
    FC.Mass_plate1 = (FC.H_plate1channels*FC.t_plate1_wall + (FC.W_plate1channels+FC.t_plate1_wall)*FC.t_plate1)*FC.CH_Flow1*FC.L_node*FC.Density_plate1;
    %% -----Anode Gas Stream-----------%
    FC.flow1_crossArea= FC.H_plate1channels*FC.W_plate1channels*FC.CH_Flow1;              % [m2] Crossectional Area of Anode entrance
    FC.h_flow1=FC.Nu_flow1*FC.k_flow1/FC.Dh_Flow1;                     % [W/m2/K]  Convection coefficient between the anode gas and the Fuel Seperator plate
    FC.A_flow1_plate1 = (2*FC.H_plate1channels + FC.W_plate1channels)*FC.L_node*FC.CH_Flow1;       % [m2]  Area in common between Anode stream and Sep Plate for convection
    FC.A_flow1_elec = (FC.W_plate1channels)*FC.L_node*FC.CH_Flow1;                  % [m2]  Area in common between Anode stream and Electrolyte for convection
    FC.Vol_flow1 = FC.H_plate1channels*FC.W_plate1channels*FC.L_node*FC.CH_Flow1;               % [m3]  control volume
    FC.A_Node_Surf= FC.A_Node;                                    % [m^2] Surface Area for internal refoming

    %% --------Electrolyte-------------%
    FC.A_Elec_Cond =  FC.W_node*FC.t_Elec;                   % [m2] Conduction surface area of electrolyte
    FC.A_Elec_Heat_Cond = FC.W_node*FC.t_Elec;                    % [m2] Conduction surface area of electrolyte
    FC.Vol_Elec = FC.t_Elec*FC.L_node*FC.W_node;              % [m3] volume of electrolyte    
    %% -------Cathode Gas Stream---------%
    FC.flow2_crossArea= FC.H_plate2channels*FC.W_plate2channels*FC.CH_Flow2;       % [m2] Crossectional Area of Cathode entrance
    FC.h_flow2= FC.Nu_flow2*FC.k_flow2/FC.Dh_Flow2;                 % [W/m2/K]  Convection coefficient between the Cathode gas and the Fuel Seperator plate
    FC.A_flow2_plate2 = (2*FC.H_plate2channels + FC.W_plate2channels)*FC.L_node*FC.CH_Flow2;    % [m2]  Area in common between Cathode stream and Sep Plate for convection
    FC.A_flow2_elec = FC.W_plate2channels*FC.CH_Flow2*FC.L_node;                 % [m2]  Area in common between Cathode stream and Electrolyte for convection
    FC.Vol_flow2 = FC.H_plate2channels*FC.W_plate2channels*FC.CH_Flow2*FC.L_node;            % [m3]  control volume Cathode
    
    %% ----Oxidant Seperator plate-------%   
    FC.A_plate2_elecCond = FC.t_plate2_wall*FC.L_node*FC.CH_Flow2;                % [m2] Conduction area between the fuel seperator plate and the electrolyte
    FC.A_plate2_heatCond = (FC.H_plate2channels*FC.t_plate2_wall + (FC.W_plate2channels+FC.t_plate2_wall)*FC.t_plate2)*FC.CH_Flow2; %[m^2] conduction area between nodes
    FC.L_plate2_heatCond=FC.H_plate2channels;                                    % [m] Length of conduction between the fuel seperator plate and electrolyte
    FC.Mass_plate2 = (FC.H_plate2channels*FC.t_plate2_wall + (FC.W_plate2channels+FC.t_plate2_wall)*FC.t_plate2)*FC.L_node*FC.CH_Flow2*FC.Density_plate2;
    switch FC.Reformer
        case 'internal'
            FC = FlowDir(FC,3); %% Load flow direction
            %% ------ Seperate Reformer Channels ---
            FC.Dh_Flow3 = 4*(FC.H_Reform*FC.W_Reform)/(2*(FC.H_Reform+FC.W_Reform)); %(m) hydraulic diameter of channel
            FC.CH_Flow3 = FC.W_node/(FC.W_Reform+FC.t_Ref_CH_Wall);    % Number of channels in each node of the reformer
            FC.Vol_flow3=FC.H_Reform*FC.W_Reform*FC.CH_Flow3*FC.L_node; % (m^3) Volume of Reformer Channel in cell
            FC.flow3_crossArea=FC.H_Reform*FC.W_Reform*FC.CH_Flow3;     % (m^2) Reformer Channel Area per node
            FC.h_flow3=FC.Nu_Reform*FC.k_Ref_Channel/FC.Dh_Flow3;                     % [W/m2/K]  Convection coefficient between the anode gas and the Fuel Seperator plate   
%         case 'adiabatic'
%             FC.Vol_flow3 = .01;        % [m^3] volume
%             FC.C_Reformer = .600;       % [kJ/(kg K)] specific heat of fuel seperator plate
%             FC.M_Reformer = 20;         % [kg] mass of reformer
        case {'external';'direct';'pox';'adiabatic';}
            FC = FlowDir(FC,2); %% Load flow direction
    end

    %% Pressure
    FC.Flow1_Pout =  FC.PressureRatio*101;
    FC.Flow1_Pinit = FC.Flow1_Pout + FC.Flow1Pdrop;
    FC.Flow2_Pout = FC.PressureRatio*101;
    FC.Flow2_Pinit = FC.Flow2_Pout + FC.Flow2Pdrop;
    
    % number of cells
    if strcmp(FC.Specification,'cells')
        FC.Cells = FC.SpecificationValue;
        FC.Specification = 'power density';
        FC.SpecificationValue = FC.RatedStack_kW*100/(FC.L_Cell*FC.W_Cell*FC.Cells);
    elseif strcmp(FC.Specification,'power density')
        FC.Cells = ceil(FC.RatedStack_kW*100/(FC.L_Cell*FC.W_Cell*FC.SpecificationValue)); %# of cells in stack
    elseif strcmp(FC.Specification,'current density')
        FC.Cells = ceil(FC.RatedStack_kW*1000/(0.8*1e4*FC.L_Cell*FC.W_Cell*FC.SpecificationValue)); %# of cells in stack (assumes voltage of 0.8)
    elseif strcmp(FC.Specification,'voltage')
        FC.Cells = ceil(FC.RatedStack_kW*1000/(FC.SpecificationValue*1e4*FC.L_Cell*FC.W_Cell*0.5)); %# of cells in stack (assumes 0.5 A/cm^2) corrected later
    end 
    %% %% 1st guess at Initial Condition
    FC.Current = zeros(FC.nodes,1);
    if strcmp(FC.Specification,'power density')
        FC.Voltage = .85;
        i_avg = FC.SpecificationValue/FC.Voltage/1000; %convert mW/cm^2 to A/cm^2, assume an initial guess voltage of 0.85
    elseif strcmp(FC.Specification,'current density')
        i_avg = FC.SpecificationValue;
        FC.Voltage = FC.RatedStack_kW/FC.Cells*1000/(FC.A_Cell*(100^2))/i_avg; %convert kW to W/cm^2, then divide by A/cm^2 to get V
        Inlet.NetCurrent = i_avg*(FC.A_Cell*(100^2));
    elseif strcmp(FC.Specification,'voltage')
        FC.Voltage = FC.SpecificationValue;
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
    FC.Spec1 = fieldnames(FC.Fuel);
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

    
    for i = 1:1:length(FC.Spec1)
        Inlet.Flow1.(FC.Spec1{i}) = FC.Fuel.(FC.Spec1{i})*FC.FuelFlowInit;
    end
    switch FC.FCtype
        case {'SOFC';'MCFC'}
            criticalSpecies = {'CH4';'CO';'CO2';'H2';'H2O'};
    end
    
    for i = 1:1:length(criticalSpecies)
        if ~ismember(criticalSpecies{i},FC.Spec1)
            Inlet.Flow1.(criticalSpecies{i}) = 0;
        end
    end
    FC.Spec1 = unique([FC.Spec1;criticalSpecies]);
    Inlet.Flow1.T = FuelTempIn;
    Inlet.AnodeIn = Inlet.Flow1;
    S2C = FC.Fuel.H2O/(FC.Fuel.CH4+.5*FC.Fuel.CO);
    if S2C<FC.Steam2Carbon %add anode recirculation
        Inlet.AnodeIn.T = 300;%mixing provides humidification & pre-heat
        Inlet.Flow1.CH4 = FC.Fuel.CH4*FC.FuelFlowInit;
        eWGS = .7; %initial guess of effective CO conversion
        r = 0.5; %Initial guess of anode recirculation
        dr = 1e-5;
        error = 1;
        % Inlet = (Inlet + generated - consumed)*r  + New, thus inlet = New/(1-r) + (generated - consumed)*r/(1-r)
        while abs(error)>1e-6
            Inlet.Flow1.CO = FC.Fuel.CO*FC.FuelFlowInit/(1-r) + (FC.Fuel.CH4 - eWGS*(FC.Fuel.CH4+FC.Fuel.CO))*FC.FuelFlowInit*r/(1-r);
            Inlet.Flow1.H2O = FC.Fuel.H2O*FC.FuelFlowInit/(1-r) + (FC.Cells*sum(FC.Current)/(2*F*1000) - (FC.Fuel.CH4 + (FC.Fuel.CH4 + FC.Fuel.CO)*eWGS)*FC.FuelFlowInit)*r/(1-r);
            S2C = Inlet.Flow1.H2O/(Inlet.Flow1.CH4 + 0.5*Inlet.Flow1.CO);
            error = FC.Steam2Carbon - S2C;
            r2 = r+dr;
            COin2 = FC.Fuel.CO*FC.FuelFlowInit/(1-r2) + (FC.Fuel.CH4 - eWGS*(FC.Fuel.CH4+FC.Fuel.CO))*FC.FuelFlowInit*r2/(1-r2);
            H2Oin2 = FC.Fuel.H2O*FC.FuelFlowInit/(1-r2) + (FC.Cells*sum(FC.Current)/(2*F*1000) - (FC.Fuel.CH4 + (FC.Fuel.CH4 + FC.Fuel.CO)*eWGS)*FC.FuelFlowInit)*r2/(1-r2);
            S2C2 = H2Oin2/(Inlet.Flow1.CH4 + 0.5*COin2);
            dSdr = (S2C2 - S2C)/dr;
            r = r + max(-.5*r,min((1-r)/2,error/dSdr));
        end
        FC.Recirc.Anode = r;
        Inlet.Flow1.CO2 = FC.Fuel.CO2*FC.FuelFlowInit/(1-r) + eWGS*(FC.Fuel.CH4+FC.Fuel.CO)*FC.FuelFlowInit*r/(1-r);
        Inlet.Flow1.H2 = FC.Fuel.H2*FC.FuelFlowInit/(1-r) + ((3*FC.Fuel.CH4 + (FC.Fuel.CH4 + FC.Fuel.CO)*eWGS)*FC.FuelFlowInit - FC.Cells*sum(FC.Current)/(2*F*1000))*r/(1-r);
        for i = 1:1:length(FC.Spec1)
            if ~ismember(FC.Spec1{i},criticalSpecies)
                Inlet.Flow1.(FC.Spec1{i}) = FC.Fuel.(FC.Spec1{i})*FC.FuelFlowInit/(1-r);
            end
            AnOutlet.(FC.Spec1{i}) = Inlet.Flow1.(FC.Spec1{i}) - Inlet.AnodeIn.(FC.Spec1{i});
        end
        %%find resulting temperature of mixture
        errorT = 1;
        Inlet.Flow1.T = FuelTempIn;
        AnOutlet.T = FC.TpenAvg + .5*FC.deltaTStack;
        Hin = enthalpy(Inlet.AnodeIn);
        Hout = enthalpy(AnOutlet);
        Hnet = Hin + FC.Recirc.Anode*Hout;
        Cp = SpecHeat(AnOutlet);
        NetFlowMix = NetFlow(Inlet.Flow1);
        while abs(errorT)>1e-3
            Hmix = enthalpy(Inlet.Flow1);
            errorT = (Hnet-Hmix)/(Cp*NetFlowMix);
            Inlet.Flow1.T = Inlet.Flow1.T + errorT;
        end 
    else
        FC.Recirc.Anode = 0;
%         Inlet.FuelMix = Inlet.AnodeIn;
    end
    FC.T.FuelMix = Inlet.Flow1.T;
    
    
    Inlet.Flow2Pout = FC.Flow2_Pout;
    Inlet.Flow1Pout = FC.Flow1_Pout;
    if ~isfield(FC,'OxidantUtilization')
        if FC.ClosedCathode
            FC.OxidantUtilization = 1;
        elseif strcmp(FC.Reformer,'internal')
            FC.OxidantUtilization = .33;
        else
            FC.OxidantUtilization = .1;
        end
    end
    FC.AirFlow = FC.Cells*sum(FC.Current)/(4*F*FC.Flow2Spec.O2)/1000/FC.OxidantUtilization;%kmol of oxidant

    FC.Spec2 = fieldnames(FC.Flow2Spec);
    if FC.ClosedCathode
        FC.Spec2 = {}; %no cathode flow states
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
        if ~ismember(FC.Spec2,criticalSpecies{i})
            Inlet.Flow2.(criticalSpecies{i}) = 0;
        end
    end
    FC.Spec2 = unique([FC.Spec2,criticalSpecies]);
    %% Run Initial Condition
    [Cathode, Anode,FC,Inlet] = solveInitCond(Inlet,FC,1);
    
    if strcmp(FC.Reformer,'external') || strcmp(FC.Reformer,'adiabatic')
        FC.Reformer = 'direct'; %external and adiabatic reformers handled in seperate block, after 1st initialization
    end
    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    FC.PortNames = {'NetCurrent','Flow1','Flow2','Flow1Pout','Flow2Pout','Flow1Out','Flow2Out','Flow2Pin','Flow1Pin','MeasureVoltage','MeasurePower','MeasureTpen','MeasureTflow1','MeasureTflow2'};
    FC.NetCurrent.type = 'in';
    FC.NetCurrent.IC = sum(FC.Current);

    FC.Flow1.type = 'in';
    FC.Flow1.IC = Inlet.Flow1; 

    FC.Flow2.type = 'in';
    FC.Flow2.IC = Inlet.Flow2;
    
    FC.Flow1Pout.type = 'in';
    FC.Flow1Pout.IC = Inlet.Flow1Pout;
    FC.Flow1Pout.Pstate = []; %identifies the state # of the pressure state if this block has one

    FC.Flow2Pout.type = 'in';
    FC.Flow2Pout.IC = Inlet.Flow2Pout;
    FC.Flow2Pout.Pstate = []; %identifies the state # of the pressure state if this block has one

    FC.Flow1Out.type = 'out';
    FC.Flow1Out.IC = MergeLastColumn(Anode.Outlet,FC.Flow1Dir,FC.Cells);
    
    FC.Flow2Out.type = 'out';
    FC.Flow2Out.IC  = MergeLastColumn(Cathode.Outlet,FC.Flow2Dir,FC.Cells);

    FC.Flow1Pin.type = 'out';
    FC.Flow1Pin.IC = FC.Flow1_Pinit;
    FC.Flow1Pin.Pstate = length(FC.Scale)-1; %identifies the state # of the pressure state if this block has one

    FC.Flow2Pin.type = 'out';
    FC.Flow2Pin.IC = FC.Flow2_Pinit;
    FC.Flow2Pin.Pstate = length(FC.Scale); %identifies the state # of the pressure state if this block has one

    FC.MeasureVoltage.type = 'out';
    FC.MeasureVoltage.IC = FC.Voltage;

    FC.MeasurePower.type = 'out';
    FC.MeasurePower.IC = sum(FC.Current*FC.Voltage*FC.Cells)/1000;%power in kW

    FC.MeasureTpen.type = 'out';
    FC.MeasureTpen.IC = FC.T.Elec;

    FC.MeasureTflow1.type = 'out';
    FC.MeasureTflow1.IC = FC.T.Anode(FC.Flow1Dir(:,end));
    
    FC.MeasureTflow2.type = 'out';
    FC.MeasureTflow2.IC = FC.T.Cath(FC.Flow2Dir(:,end));

    FC.P_Difference = {'Flow2Pin','Flow2Pout'; 'Flow1Pin', 'Flow1Pout';};

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
    FC.Specification = 'current density';%converge only to match current density from controller
    FC.Recirc.Anode = 0; %after 1st initialization recirculation is handled in controller, valve & mixing volume
%     Inlet.FuelMix = Inlet.AnodeIn;
    Inlet.AnodeIn = Inlet.Flow1;
    Flow1New = fieldnames(Inlet.Flow1);
    Flow1All = unique([FC.Spec1;Flow1New]);
    Flow1All = Flow1All(~strcmp('T',Flow1All));
    for i = 1:1:length(Flow1All)
        if ~ismember(Flow1All{i},Flow1New)
            Inlet.Flow1.(Flow1All{i})=0;
        end
    end
    FC.Spec1 = Flow1All;

    if FC.ClosedCathode
        FC.Spec2 = {}; %no cathode flow states
    else
        Flow2New = fieldnames(Inlet.Flow2);
        Flow2All = unique([FC.Spec2;Flow2New]);
        Flow2All = Flow2All(~strcmp('T',Flow2All));
        for i = 1:1:length(Flow2All)
            if ~ismember(Flow2All{i},Flow2New)
                Inlet.Flow2.(Flow2All{i})=0;
            end
        end
        FC.Spec2 = Flow2All;
    end
    FC.Flow1_Pinit = Inlet.Flow1Pout + FC.Flow1Pdrop;
    FC.Flow2_Pinit = Inlet.Flow2Pout + FC.Flow2Pdrop;
    FC.Flow2Pout.IC = Inlet.Flow2Pout;
    FC.Flow1Pout.IC = Inlet.Flow1Pout;
    %%--%%
    [Cathode, Anode,FC,~] = solveInitCond(Inlet,FC,2);
    %%%
    FC.Flow1Pin.Pstate = length(FC.Scale)-1; %identifies the state # of the pressure state if this block has one
    FC.Flow2Pin.Pstate = length(FC.Scale); %identifies the state # of the pressure state if this block has one
    FC.Flow2Out.IC  = MergeLastColumn(Cathode.Outlet,FC.Flow2Dir,FC.Cells);
    FC.Flow1Out.IC = MergeLastColumn(Anode.Outlet,FC.Flow1Dir,FC.Cells);
    FC.Flow2Pin.IC = FC.Flow2_Pinit;
    FC.Flow1Pin.IC = FC.Flow1_Pinit;
    FC.MeasureCurrent.IC = sum(FC.Current);
    FC.MeasurePower.IC = sum(FC.Current*FC.Voltage*FC.Cells)/1000;%power in kW
    FC.MeasureTpen.IC = FC.T.Elec;
    FC.MeasureTflow2.IC = FC.T.Cath(FC.Flow2Dir(:,end));
    FC.MeasureTflow1.IC = FC.T.Anode(FC.Flow1Dir(:,end));
    FC.HumidifiedFuelTemp.IC = FC.T.FuelMix;
end

function [Cathode, Anode,FC,Inlet] = solveInitCond(Inlet,FC,firstSolve)
global Tags F
Cp.anode = SpecHeat(Inlet.Flow1);
Cp.cath = SpecHeat(Inlet.Flow2);

error = 1;
Tol = 1e-3;
count = 1;
while abs(error)>Tol %iterate to reach target current density, voltage or power density
    Cathode = FCin2Out(FC.T.Cath,Inlet.Flow2,FC.Flow2Dir, FC.FCtype,FC.Cells,FC.Current,[],'cathode');
%     SinglePassUtilization = (sum(FC.Current)*FC.Cells/(2000*F))/(4*Inlet.Flow1.CH4+Inlet.Flow1.CO + Inlet.Flow1.H2);
    [FC,Inlet.Flow1,Anode,Reformer] = KineticCoef(FC,Inlet,((firstSolve==1) && (count==1)));%% solve for kinetic reaction coefficient which results in this outlet condition (match R_CH4 & R_WGS)
    Offset = 0;
    if count==1 && firstSolve==1
        [FC.Tstates,FC.HTcond,FC.HTconv,FC.HTrad]= SteadyTemps(FC,Inlet);
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
    FC.T.FuelMix = Inlet.Flow1.T;
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
    if strcmp(FC.Specification,'power density')
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
    elseif strcmp(FC.Specification,'voltage')
        Tol = 1e-3;
        error = (FC.SpecificationValue - FC.Voltage)/FC.Voltage;
        if count>1
            dV_di = -Tags.(FC.name).ASR/FC.A_Node/100^2;
            scale = 1 + (FC.SpecificationValue - FC.Voltage)/dV_di/TotCurrent;
        else % first time through
            localR = 3*localR;
            scale =1+sum((Tags.(FC.name).nVoltage-FC.SpecificationValue)./localR)/sum(FC.Current);
        end
        FC.Cells = ceil((FC.RatedStack_kW*1000/FC.SpecificationValue)/(scale*sum(FC.Current))); %re-calculate the # of cells
    elseif strcmp(FC.Specification,'current density')
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
            k = FC.Flow2Dir(:,end);
            dTerror = (mean((FC.Tstates(k+FC.nodes))- Inlet.Flow2.T)/FC.deltaTStack-1);
            FC.AirFlow = FC.AirFlow*(1 + dTerror)*scale^2;
            TavgError = (FC.TpenAvg-mean(FC.Tstates(2*FC.nodes+1:3*FC.nodes)))/FC.deltaTStack;
            FC.StackCathTin = FC.StackCathTin + (TavgError + 0.75*dTerror)*FC.deltaTStack;
        elseif strcmp(FC.CoolingStream,'anode') %oxidant flow rate determined by current, fuel flow rate is now determined by thermal balancing
            FC.AirFlow = FC.Cells*TotCurrent/(4000*F*FC.Flow2Spec.O2)/FC.OxidantUtilization;%kmol of oxidant
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
Cathode = FCin2Out(FC.T.Cath,Inlet.Flow2,FC.Flow2Dir, FC.FCtype,FC.Cells,FC.Current,[],'cathode');
FC.PfactorAnode = NetFlow(Inlet.Flow1)/FC.Flow1Pdrop;
FC.PfactorCath = NetFlow(Inlet.Flow2)/FC.Flow2Pdrop;
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

k = FC.Flow1Dir(:,1);
if min(Anode.Inlet.H2(k))==0 && FC.R_CH4(k)>0 %gives some reformed methane as anode inlet
    AvgX.H2 = (Anode.Outlet.H2 + 0.5*(3*FC.R_CH4(k)+FC.R_WGS(k)))./(n_an_in(k)+n_an_out(k));
    AvgX.H2O = (Anode.Outlet.H2O - 0.5*(FC.R_CH4(k)-FC.R_WGS(k)) + Anode.Inlet.H2O(k))./(n_an_in(k)+n_an_out(k));
end

%% Calculate local voltages
normTemp = FC.TpenAvg + (FC.T.Elec-mean(FC.T.Elec)); %assume you will get to the desired temperature (this avoids oscilations in voltage and helps convergence
FuelCellNernst(FC.Current,normTemp,FC.Flow2_Pinit,AvgX,FC);

function Inlet = InletFlow(FC,Inlet) %only used 1st time through initialization (before we know what is connected to inlet
% Anode
for i = 1:1:length(FC.Spec1)
    if isfield(FC.Fuel,FC.Spec1{i})
        Inlet.AnodeIn.(FC.Spec1{i}) = FC.Fuel.(FC.Spec1{i})*FC.FuelFlowInit;%flow rate of every species entering the anode (or reformer if there is one)
    else Inlet.AnodeIn.(FC.Spec1{i}) = 0;
    end
end
%Cathode
Inlet.Flow2.T = FC.StackCathTin;
switch FC.FCtype
    case 'SOFC'
        if FC.ClosedCathode
            Inlet.Flow2.O2 = FC.Flow2Spec.O2*FC.AirFlow;
        else
            for i = 1:1:length(FC.Spec2)
                if isfield(FC.Flow2Spec,FC.Spec2{i})
                    Inlet.Flow2.(FC.Spec2{i}) = FC.Flow2Spec.(FC.Spec2{i})*FC.AirFlow;
                else Inlet.Flow2.(FC.Spec2{i}) = 0;
                end
            end
        end
    case 'MCFC' %recalculate cathode inlet species for MCFC (this is an estimate assuming the 100% of non-recirculated anode gas is oxidized and fed to the cathode)
        Inlet.Flow2.CO2 = (Inlet.Flow1.CH4+Inlet.Flow1.CO+Inlet.Flow1.CO2) + sum(FC.Current)/(2*F*1000)*FC.Cells;
        Inlet.Flow2.H2O = (4*Inlet.Flow1.CH4+Inlet.Flow1.CO+Inlet.Flow1.H2+Inlet.Flow1.H2O);
        nonCO2_H2O = (FC.AirFlow - Inlet.Flow2.CO2 - Inlet.Flow2.H2O);
        for i = 1:1:length(FC.Spec2)
            if isfield(FC.Flow2Spec,FC.Spec2{i})
                if strcmp(FC.Spec2{i},'CO2')||strcmp(FC.Spec2{i},'H2O')
                    Inlet.Flow2.(FC.Spec2{i}) = Inlet.Flow2.(FC.Spec2{i}) + FC.Flow2Spec.(FC.Spec2{i})*nonCO2_H2O;
                else
                    Inlet.Flow2.(FC.Spec2{i}) = FC.Flow2Spec.(FC.Spec2{i})*nonCO2_H2O;
                end
            else Inlet.Flow2.(FC.Spec2{i}) = 0;
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
        NumOfStates = (6 + 2*length(FC.Spec1) + length(FC.Spec2) + 1)*FC.nodes + 2; % 6 temperatures, anode & cathode & reformer & current at each node and 2 states for anode/cathode pressure
  case {'direct';'external';'adiabatic'}
        NumOfStates = (5 + length(FC.Spec1) + length(FC.Spec2) + 1)*FC.nodes +2; % 5 temperatures, anode & cathode & current at each node and 2 states for anode/cathode pressure
end
FC.IC = ones(NumOfStates,1); %
FC.tC = FC.IC; % time constant for derivative dY
FC.Scale = FC.IC;
switch FC.Reformer
    case 'internal'
        n = 6*FC.nodes;
        FC.tC(5*FC.nodes+1:6*FC.nodes) = (FC.Vol_flow3*Cp.ref*FC.Flow1_Pinit./(Ru*FC.T.Reform));
    case {'direct';'external';'adiabatic'}
        n = 5*FC.nodes;
end
FC.Scale = FC.Tstates(1:n);%temperature (K)

FC.tC(1:FC.nodes) = (FC.Mass_plate2*FC.C_plate2);
FC.tC(1+FC.nodes:2*FC.nodes) = (FC.Vol_flow2*Cp.cath*FC.Flow2_Pinit./(Ru*FC.T.Cath));
FC.tC(2*FC.nodes+1:3*FC.nodes) = (FC.Vol_Elec*FC.Density_Elec*FC.C_Elec);
FC.tC(3*FC.nodes+1:4*FC.nodes) = (FC.Vol_flow1*Cp.an*FC.Flow1_Pinit./(Ru*FC.T.Anode));
FC.tC(4*FC.nodes+1:5*FC.nodes) = (FC.Mass_plate1*FC.C_plate1);
FC.tC(1:n) = FC.tC(1:n)-diag(FC.HTconv)-diag(FC.HTcond); %this accounts for the change in HT as temperature of the control volume changes. The change in HT helps balance the energy equation more than the change in enthalpy leaving.

for i = 1:1:length(FC.Spec2)
    FC.tC(n+1:n+FC.nodes) = (FC.Vol_flow2*FC.Flow2_Pinit)./(FC.T.Cath*Ru);  % cathode 
    if any(Cathode.Outlet.(FC.Spec2{i})==0)
        Flow = NetFlow(Cathode.Outlet);
        FC.IC(n+1:n+FC.nodes) = Cathode.Outlet.(FC.Spec2{i})./Flow;
        FC.Scale(n+1:n+FC.nodes) = Flow; n = n+FC.nodes; %cathode flows
    else
        FC.Scale(n+1:n+FC.nodes) = Cathode.Outlet.(FC.Spec2{i}); n = n+FC.nodes; %cathode flows
    end
end

Flow = NetFlow(Anode.Outlet);
for i = 1:1:length(FC.Spec1)
    X = Anode.Outlet.(FC.Spec1{i})./Flow;%concentration
    FC.tC(n+1:n+FC.nodes) = (FC.Vol_flow1*FC.Flow1_Pinit)./(FC.T.Anode*Ru); %anode
    if any(X<.01) %concentration less than 1%
        FC.IC(n+1:n+FC.nodes) = X;
        FC.Scale(n+1:n+FC.nodes) = Flow; %anode flow
    else
        FC.Scale(n+1:n+FC.nodes) = Anode.Outlet.(FC.Spec1{i}); %individual species flow
    end   
    n = n+FC.nodes;
end

switch FC.Reformer
    case 'internal'
        for i = 1:1:length(FC.Spec1)
            FC.tC(n+1:n+FC.nodes) = (FC.Vol_flow3*FC.Flow1_Pinit)./(FC.T.Reform*Ru); % reformer
            if any(Reformer.Outlet.(FC.Spec1{i})==0)
                Flow = NetFlow(Reformer.Outlet);
                FC.IC(n+1:n+FC.nodes) = Reformer.Outlet.(FC.Spec1{i})./Flow;
                FC.Scale(n+1:n+FC.nodes) = Flow; n = n+FC.nodes; %anode flows
            else
                FC.Scale(n+1:n+FC.nodes) = Reformer.Outlet.(FC.Spec1{i}); n = n+FC.nodes; %reformer flows
            end
        end
end
FC.tC(n+1:n+FC.nodes) = FC.nodes/100;%  %current changing for voltage balance
FC.Scale(n+1:n+FC.nodes) = FC.Current;  n = n+FC.nodes; %current

FC.tC(n+1) = (FC.Vol_flow1*FC.nodes*FC.Cells); %pressure
FC.tC(n+2) = (FC.Vol_flow2*FC.nodes*FC.Cells);  %pressure
FC.Scale(n+1) = FC.Flow1_Pinit;%pressure
FC.Scale(n+2) = FC.Flow2_Pinit;%pressure


function dY = DynamicFCtemps(t,Y,FC,Anode,Cathode,Reformer,Inlet)
global F Ru
dY = 0*Y;
FuelInlet = Inlet.Flow1;
nodes = FC.nodes;
if isfield(FC,'tC')
    tC = FC.tC(1:6*nodes);
else
    Cp.cath = 33;% kJ/kmol*K
    Cp.an = 42;% kJ/kmol*K
    tC(1:nodes,1) = (FC.Mass_plate2*FC.C_plate2);
    tC(1+nodes:2*nodes,1) = (FC.Vol_flow2*Cp.cath*FC.Flow2_Pinit./(Ru*FC.T.Cath));
    tC(2*nodes+1:3*nodes,1) = (FC.Vol_Elec*FC.Density_Elec*FC.C_Elec);
    tC(3*nodes+1:4*nodes,1) = (FC.Vol_flow1*Cp.an*FC.Flow1_Pinit./(Ru*FC.T.Anode));
    tC(4*nodes+1:5*nodes,1) = (FC.Mass_plate1*FC.C_plate1);
    if isfield(FC,'Vol_flow3')
        Cp.ref = 42;% kJ/kmol*K
        tC(5*nodes+1:6*nodes,1) = (FC.Vol_flow3*Cp.ref*FC.Flow1_Pinit./(Ru*FC.T.Reform));
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

QT = FC.HTconv*Y + FC.HTcond*Y + FC.HTrad*(Y.^4);

Cathode.Outlet.T = Y(nodes+1:2*nodes);
for j = 1:1:length(FC.Flow2Dir(1,:));%1:columns
    k = FC.Flow2Dir(:,j);
    if j~=1
        Cathode.Inlet.T(k,1) = Cathode.Outlet.T(kprev);
    end
    kprev = k;
end

Anode.Outlet.T = Y(3*nodes+1:4*nodes);
for j = 1:1:length(FC.Flow1Dir(1,:));%1:columns
    k = FC.Flow1Dir(:,j);
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
    Hmix = HfreshFuel + sum(HoutAnode(FC.Flow1Dir(:,end)))*FC.Recirc.Anode*FC.Cells;
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
        k = FC.Flow3Dir(:,1);
        Reformer.Inlet.T(k,1) = FuelInlet.T;
        for j = 1:1:length(FC.Flow3Dir(1,:));%1:columns
            k = FC.Flow3Dir(:,j);
            if j~=1
                Reformer.Inlet.T(k,1) = Reformer.Outlet.T(kprev);
            end
            kprev = k;
        end
        HinReform = enthalpy(Reformer.Inlet);
        HoutReform = enthalpy(Reformer.Outlet);
        
        k = FC.Flow1Dir(:,1);
        k2 = FC.Flow3Dir(:,end);
        Anode.Inlet.T(k,1) = Reformer.Outlet.T(k2,1);
    case {'adiabatic';'direct';'external';}
        k = FC.Flow1Dir(:,1);
        Anode.Inlet.T(k,1) = FuelInlet.T;
end
HinAnode = enthalpy(Anode.Inlet);

if FC.ClosedCathode %%energy balance
    Qimbalance = sum((HinCath(FC.Flow2Dir(:,1))) - sum(HoutCath(FC.Flow2Dir(:,end)))) + sum(HinReform(FC.Flow3Dir(:,1)))  - sum(HoutAnode(FC.Flow1Dir(:,end))) - sum(Power);
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