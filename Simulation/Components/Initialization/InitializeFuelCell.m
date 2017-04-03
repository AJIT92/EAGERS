function block = InitializeFuelCell(varargin)
%FC model with many states: Temperatures (oxidizer plate, cathode, electrolyte, anode, fuel plate [possibly a reformer]), Cathode species, anode species (can include anode recirculation for humidification, [reformer: indirect, adiabatic or none], Current, cathode pressure, anode pressure
% Five (5) inlets: {'NetCurrent','CathodeIn','AnodeIn','Flow1Pout','Flow2Pout'}
% Seven (7) outlets: {'Flow1Out','Flow2Out','Flow1Pin','Flow2Pin','MeasureCurrent','MeasureTpen','MeasureTflow1','MeasureTflow2'}
% for a fuel cell, Flow 1 is the anode (fuel), and Flow2 is the cathode (cooling/heating air)
% Current is positive
global F Ru
F=96485.339; % %Faraday's constant in Coulomb/mole
Ru = 8.314472; % Universal gas constant in kJ/K*kmol
block = varargin{1};
if length(varargin)==1 % first initialization
    block.nodes = block.rows*block.columns;
    block.AnPercEquilib = 1; %CH4 reforming reaches equilibrium at anode exit.
%%%---%%% User defined variables
    %%--Geometric Varibles  %
    block.t_plate1 = 0.003;                   % [m] thickness of seperator plate
    block.t_plate1_wall =0.005;                  % [m] Thickness of the channel wall of the Fuel Seperator Plate
    block.H_plate1channels = 0.002;                       % [m] height of anode channel
    block.W_plate1channels = 0.005;                       % [m] width of channel                       
    block.H_plate2channels = 0.002;                      % [m] height of Cathode channel
    block.W_plate2channels = 0.005;                      % [m] width of Cathode channel
    block.t_plate2 = 0.003;                     % [m] Thickness of Oxidant Seperator Plate
    block.t_plate2_wall=0.005;                     % [m] Thickness of the channel wall of the Oxidant Seperator Plate
    block.Nu_flow1 = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp
    block.Nu_flow2 = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp

    %%Electrochemical parameters %%%
    % H2 +1/2O2 --> H2O (Nernst Eo)
    %SOFC Conductivity - Implemented Equation  sigma = A*e^(-deltaG/RT)/T
    switch block.FCtype
        case 'SOFC'
            block.ElecConst = 2e3; %(K/ohm*m) Electrolyte Constant  %default SOFC  = 9e7
            block.deltaG = 8.0e3; %(kJ/kmol)
            %%unused:
%             block.Io = 5000;          % [A/m2] activation current density %default SOFC = 1000
%             block.DeffO2 = 4.0e-5; %(m^2/s)
%             block.alpha=.7;
            block.t_Membrane = 18e-6;                     % [m] thickness of membrane
            block.t_Cath = 800e-6;                        % [m] thickness of cathode structure
            block.t_An = 50e-6;                           % [m] thickness of Anode structure
            block.t_Elec = block.t_Membrane+block.t_Cath+block.t_An;        % [m] thickness of complete electrolyte
        case 'MCFC'
            block.Io = 500;            % [A/m2] activation current density
            block.alpha=.4;      
            block.J_L = 6000;            % [A/m2] Limiting  current density   
            block.Cr0 = 4.7833e-4;%4.25e-5;%
            block.Cr1 = -6.6667e-7;%-5e-8;%
            block.t_Elec = 0.003;        % [m] thickness of complete electrolyte
    end

    %%Electrolyte
    block.k_Elec =6.19;                                % [W/m K]  Conductivity of the Electrolyte
    block.Density_Elec = 375;                          % [kg/m3] Density of Electrolyte
    block.C_Elec = .800;                                  % [kJ/(kg K)] specific heat of electrolyte 

    %%---Anode half of interconnect plate-------%
    block.Density_plate1 = 2000;                                % [kg/m3]     density 
    block.C_plate1 = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    block.k_plate1 = 5;   %25                                    % [W/(m K)]   conductivity of Fuel Seperator Plate
    %%-----Anode Gas Stream-----------%
    block.k_flow1 = 259E-3;                                 % (W/m*K) Thermal Conductivity of 50%H2 & 50%H2O at 1000K

    %%-------Cathode Gas Stream---------%
    block.k_flow2 = 67E-3;                                              % (W/m*K) Thermal Conductivity of air at 1000K

    %%----Cathode half of interconnect plate-------%   
    block.Density_plate2 = 2000;                                % [kg/m3]     density 
    block.C_plate2 = .600;                                        % [kJ/(kg K)] specific heat of fuel seperator plate
    block.k_plate2 = 5;%25;                                % [W/(m K)]   conductivity of Fuel Seperator Plate

    %%----- / ReformerChannel---
    block.H_Reform = 0.005;                   % (m) height of Reformer channel
    block.W_Reform = 0.003;                  % (m) width of Reformer channel 
    block.Nu_Reform = 4;                           %Nusselt # for square channel aspect ratio =3 & uniform temp
    block.t_Ref_CH_Wall = .001;               % (m)Thickness of reformer channel wall
    block.k_Ref_Channel = 112.52E-3;                                          % (W/m*K) Thermal Conductivity of 5%CO, 35%CO2 15%H2, 45%H20 at 900K
%%%---%%% end of user defined variables    
    %Dimensions
    block.A_Cell = block.L_Cell*block.W_Cell; %Cell Area
    block.A_Node = block.A_Cell/block.nodes; %node Area
    block.L_node = block.L_Cell/block.columns; %Node length in meters
    block.W_node = block.W_Cell/block.rows;  %Node Width in meters
    block.Dh_Flow1 = 4*(block.H_plate1channels*block.W_plate1channels)/(2*(block.H_plate1channels+block.W_plate1channels)); %(m) hydraulic diameter of channel
    block.Dh_Flow2 = 4*(block.H_plate2channels*block.W_plate2channels)/(2*(block.H_plate2channels+block.W_plate2channels)); %(m) hydraulic diameter of channel
    block.CH_Flow1 = block.W_node/(block.W_plate1channels+block.t_plate1_wall); %Number of channels in each node of the anode
    block.CH_Flow2 = block.W_node/(block.W_plate2channels+block.t_plate2_wall); %Number of channels in each node of the cathode 
    %% ---Fuel Seperator plate-------%
    block.A_plate1_elecCond = block.t_plate1_wall*block.L_node*block.CH_Flow1;   % [m2] Conduction area between the fuel seperator plate and the electrolyte
    block.A_plate1_heatCond=(block.H_plate1channels*block.t_plate1_wall + (block.W_plate1channels+block.t_plate1_wall)*block.t_plate1)*block.CH_Flow1; %[m^2] conduction area between nodes
    block.L_plate1_heatCond=block.H_plate1channels;                                     % [m] Lenght of conduction between the fuel seperator plate and electrolyte
    block.Mass_plate1 = (block.H_plate1channels*block.t_plate1_wall + (block.W_plate1channels+block.t_plate1_wall)*block.t_plate1)*block.CH_Flow1*block.L_node*block.Density_plate1;
    %% -----Anode Gas Stream-----------%
    block.flow1_crossArea= block.H_plate1channels*block.W_plate1channels*block.CH_Flow1;              % [m2] Crossectional Area of Anode entrance
    block.h_flow1=block.Nu_flow1*block.k_flow1/block.Dh_Flow1;                     % [W/m2/K]  Convection coefficient between the anode gas and the Fuel Seperator plate
    block.A_flow1_plate1 = (2*block.H_plate1channels + block.W_plate1channels)*block.L_node*block.CH_Flow1;       % [m2]  Area in common between Anode stream and Sep Plate for convection
    block.A_flow1_elec = (block.W_plate1channels)*block.L_node*block.CH_Flow1;                  % [m2]  Area in common between Anode stream and Electrolyte for convection
    block.Vol_flow1 = block.H_plate1channels*block.W_plate1channels*block.L_node*block.CH_Flow1;               % [m3]  control volume
    block.A_Node_Surf= block.A_Node;                                    % [m^2] Surface Area for internal refoming

    %% --------Electrolyte-------------%
    block.A_Elec_Cond =  block.W_node*block.t_Elec;                   % [m2] Conduction surface area of electrolyte
    block.A_Elec_Heat_Cond = block.W_node*block.t_Elec;                    % [m2] Conduction surface area of electrolyte
    block.Vol_Elec = block.t_Elec*block.L_node*block.W_node;              % [m3] volume of electrolyte    
    %% -------Cathode Gas Stream---------%
    block.flow2_crossArea= block.H_plate2channels*block.W_plate2channels*block.CH_Flow2;       % [m2] Crossectional Area of Cathode entrance
    block.h_flow2= block.Nu_flow2*block.k_flow2/block.Dh_Flow2;                 % [W/m2/K]  Convection coefficient between the Cathode gas and the Fuel Seperator plate
    block.A_flow2_plate2 = (2*block.H_plate2channels + block.W_plate2channels)*block.L_node*block.CH_Flow2;    % [m2]  Area in common between Cathode stream and Sep Plate for convection
    block.A_flow2_elec = block.W_plate2channels*block.CH_Flow2*block.L_node;                 % [m2]  Area in common between Cathode stream and Electrolyte for convection
    block.Vol_flow2 = block.H_plate2channels*block.W_plate2channels*block.CH_Flow2*block.L_node;            % [m3]  control volume Cathode
    
    %% ----Oxidant Seperator plate-------%   
    block.A_plate2_elecCond = block.t_plate2_wall*block.L_node*block.CH_Flow2;                % [m2] Conduction area between the fuel seperator plate and the electrolyte
    block.A_plate2_heatCond = (block.H_plate2channels*block.t_plate2_wall + (block.W_plate2channels+block.t_plate2_wall)*block.t_plate2)*block.CH_Flow2; %[m^2] conduction area between nodes
    block.L_plate2_heatCond=block.H_plate2channels;                                    % [m] Length of conduction between the fuel seperator plate and electrolyte
    block.Mass_plate2 = (block.H_plate2channels*block.t_plate2_wall + (block.W_plate2channels+block.t_plate2_wall)*block.t_plate2)*block.L_node*block.CH_Flow2*block.Density_plate2;
    switch block.Reformer
        case 'internal'
            block = FlowDir(block,3); %% Load flow direction
            %% ------ Seperate Reformer Channels ---
            block.Dh_Flow3 = 4*(block.H_Reform*block.W_Reform)/(2*(block.H_Reform+block.W_Reform)); %(m) hydraulic diameter of channel
            block.CH_Flow3 = block.W_node/(block.W_Reform+block.t_Ref_CH_Wall);    % Number of channels in each node of the reformer
            block.Vol_flow3=block.H_Reform*block.W_Reform*block.CH_Flow3*block.L_node; % (m^3) Volume of Reformer Channel in cell
            block.flow3_crossArea=block.H_Reform*block.W_Reform*block.CH_Flow3;     % (m^2) Reformer Channel Area per node
            block.h_flow3=block.Nu_Reform*block.k_Ref_Channel/block.Dh_Flow3;                     % [W/m2/K]  Convection coefficient between the anode gas and the Fuel Seperator plate   
%         case 'adiabatic'
%             block.Vol_flow3 = .01;        % [m^3] volume
%             block.C_Reformer = .600;       % [kJ/(kg K)] specific heat of fuel seperator plate
%             block.M_Reformer = 20;         % [kg] mass of reformer
        case {'external';'direct';'pox';'adiabatic';}
            block = FlowDir(block,2); %% Load flow direction
    end

    %% Pressure
    block.Flow1_Pout =  block.PressureRatio*101;
    block.Flow1_Pinit = block.Flow1_Pout + block.Flow1Pdrop;
    block.Flow2_Pout = block.PressureRatio*101;
    block.Flow2_Pinit = block.Flow2_Pout + block.Flow2Pdrop;
    
    % number of cells
    if strcmp(block.Specification,'cells')
        block.Cells = block.SpecificationValue;
        block.Specification = 'power density';
        block.SpecificationValue = block.RatedStack_kW*100/(block.L_Cell*block.W_Cell*block.Cells);
    elseif strcmp(block.Specification,'power density')
        block.Cells = ceil(block.RatedStack_kW*100/(block.L_Cell*block.W_Cell*block.SpecificationValue)); %# of cells in stack
    elseif strcmp(block.Specification,'current density')
        block.Cells = ceil(block.RatedStack_kW*1000/(0.8*1e4*block.L_Cell*block.W_Cell*block.SpecificationValue)); %# of cells in stack (assumes voltage of 0.8)
    elseif strcmp(block.Specification,'voltage')
        block.Cells = ceil(block.RatedStack_kW*1000/(block.SpecificationValue*1e4*block.L_Cell*block.W_Cell*0.5)); %# of cells in stack (assumes 0.5 A/cm^2) corrected later
    end 
    %% %% 1st guess at Initial Condition
    Current = zeros(block.nodes,1);
    if strcmp(block.Specification,'power density')
        block.Voltage = .85;
        i_avg = block.SpecificationValue/block.Voltage/1000; %convert mW/cm^2 to A/cm^2, assume an initial guess voltage of 0.85
    elseif strcmp(block.Specification,'current density')
        i_avg = block.SpecificationValue;
        block.Voltage = block.RatedStack_kW/block.Cells*1000/(block.A_Cell*(100^2))/i_avg; %convert kW to W/cm^2, then divide by A/cm^2 to get V
        Inlet.NetCurrent = i_avg*(block.A_Cell*(100^2));
    elseif strcmp(block.Specification,'voltage')
        block.Voltage = block.SpecificationValue;
        i_avg = block.RatedStack_kW/block.Cells*1000/(block.A_Cell*(100^2))/block.Voltage; %convert kW to W/cm^2, then divide by V to get A/cm^2
    end
    for j = 1:1:block.rows
        Current(1+block.columns*(j-1):block.columns*j) =linspace(2,1,block.columns)/sum(linspace(2,1,block.columns))*i_avg*(100^2)*block.A_Cell/block.rows; %make the initial current guess low to not overutilize H2 in 1st iteration of solution
    end
    block.Current.CO = 0*Current;
    block.Current.H2 = Current;
    
    block.StackCathTin  = block.TpenAvg -.75*block.deltaTStack;
    block.T.Flow2 = zeros(block.nodes,1) + block.TpenAvg;
    block.T.Elec = zeros(block.nodes,1) + block.TpenAvg;
    block.T.Flow1 = zeros(block.nodes,1) + block.TpenAvg;
    block.T.Flow3 = zeros(block.nodes,1) + block.TpenAvg;
    %% initial guess of reforming cooling
    block.FuelFlowInit  = sum(Current)/(2*F*1000)/(block.FuelUtilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))*block.Cells; % fuel flow rate,  current/(2*F*1000) = kmol H2
    block.Spec1 = fieldnames(block.Fuel);
    FuelTempIn =block.TpenAvg-block.deltaTStack;
    switch block.Reformer
        case 'adiabatic'
            block.ReformT = 823; %an initial guess temperature for adiabatic reforme, or if uncommenting line 873, this is the setpoint
            block.Steam2Carbon = 6; %determines anode recirculation, needs to be high to ensure sufficient temperature for some pre-reforming
    end
    R1 = block.Fuel.CH4*block.AnPercEquilib*block.FuelFlowInit;
    switch block.Reformer
        case 'external'
            CH4ref_ext = block.RefPerc*R1;
            block.R_WGSref =  CH4ref_ext*.8;
            block.R_CH4 = (R1 - CH4ref_ext)/block.Cells*ones(block.nodes,1)/block.nodes;
            block.R_WGS = block.R_CH4*.8;
        case 'internal'
            block.R_CH4ref = block.RefPerc*R1/block.Cells*block.RefSpacing*ones(block.nodes,1)/block.nodes;
            block.R_WGSref =  block.R_CH4ref*.8;
            block.R_CH4 = (R1 - sum(block.R_CH4ref)*block.Cells/block.RefSpacing)/block.Cells*ones(block.nodes,1)/block.nodes;
            block.R_WGS = block.R_CH4*.8;
        case 'adiabatic'
            block.R_CH4ref = 0.5*R1;
            block.R_WGSref =  block.R_CH4ref*.8;
            block.R_CH4 = (R1 - block.R_CH4ref)/block.Cells*ones(block.nodes,1)/block.nodes;
            block.R_WGS = block.R_CH4*.8;
        case 'direct'
            block.R_CH4 = R1/block.Cells*ones(block.nodes,1)/block.nodes;
            block.R_WGS =  block.R_CH4*.8;
    end 
%     %% -- get surface areas and radiation view coefficients from file --%%
%     Dir=strrep(which('InitializeFuelCell.m'),fullfile('Components','Initialization','InitializeFuelCell.m'),'FCMaps');
%     load(fullfile(Dir,block.Map));
%     f = fieldnames(Map);
%     for i = 1:1:length(f)
%         block.(f{i}) = Map.(f{i});
%     end
%     Sigma = 5.670367e-11;%kW/(m^2 K^4) %all heat transfer coefficients converted to kW/m^2*K^4: thus Q = sigma*Area*(T1^4-T2^4) is in kW
%     block.RTmatrix = zeros(s*block.nodes,s*block.nodes);
%     %% Here is where it needs to agregate a 100x100 view factor map into the rows & columns of this particular FC
%     %% -- %%
%     for j = 1:1:block.nodes
%         block.RTmatrix(j,2*block.nodes+1:3*block.nodes) = Sigma*block.ViewFactorCath(j,:)*block.A_Node; %view factor from cathode plate to electrolyte
%         block.RTmatrix(j,j) = -Sigma*sum(block.ViewFactorCath(j,:))*block.A_Node; % - sum(view factors) for this node
%         
%         block.RTmatrix(2*block.nodes+j,1:block.nodes) = Sigma*block.ViewFactorCath(j,:)*block.A_Node; %view factor from electrolyte to cathode plate
%         block.RTmatrix(2*block.nodes+j,2*block.nodes+j) = -Sigma*sum(block.ViewFactorCath(j,:))*block.A_Node; % - sum(view factors) for this node
%         
%         block.RTmatrix(4*block.nodes+j,2*block.nodes+1:3*block.nodes) = Sigma*block.ViewFactorAn(j,:)*block.A_Node; %view factor from anode plate to electrolyte
%         block.RTmatrix(4*block.nodes+j,4*block.nodes+j) = -Sigma*sum(block.ViewFactorAn(j,:))*block.A_Node; % - sum(view factors) for this node
%         
%         block.RTmatrix(2*block.nodes+j,4*block.nodes+1:5*block.nodes) = Sigma*block.ViewFactorAn(j,:)*block.A_Node; %view factor from electrolyte to anode plate
%         block.RTmatrix(2*block.nodes+j,2*block.nodes+j) = block.RTmatrix(2*block.nodes+j,2*block.nodes+j) -Sigma*sum(block.ViewFactorAn(j,:))*block.A_Node; % - sum(view factors) for this node
%         switch block.Reformer
%             case 'internal'
%                 block.RTmatrix(4*block.nodes+j,1:block.nodes) = Sigma*block.ViewFactorRef(j,:)*block.A_Node; %view factor from anode plate to cathode plate, with reformer channels between
%                 block.RTmatrix(4*block.nodes+j,4*block.nodes+j) = block.RTmatrix(4*block.nodes+j,4*block.nodes+j) -Sigma*sum(block.ViewFactorRef(j,:))*block.A_Node; % - sum(view factors) for this node
%                 
%                 block.RTmatrix(j,4*block.nodes+1:5*block.nodes) = Sigma*block.ViewFactorRef(j,:)*block.A_Node; %view factor from cathode plate to anode plate, with reformer channels between
%                 block.RTmatrix(j,j) = block.RTmatrix(j,j) -Sigma*sum(block.ViewFactorRef(j,:))*block.A_Node; % - sum(view factors) for this node
%         end
%     end

    
    for i = 1:1:length(block.Spec1)
        Inlet.Flow1.(block.Spec1{i}) = block.Fuel.(block.Spec1{i})*block.FuelFlowInit;
    end
    switch block.FCtype
        case {'SOFC';'MCFC'}
            criticalSpecies = {'CH4';'CO';'CO2';'H2';'H2O'};
    end
    
    for i = 1:1:length(criticalSpecies)
        if ~ismember(criticalSpecies{i},block.Spec1)
            Inlet.Flow1.(criticalSpecies{i}) = 0;
        end
    end
    block.Spec1 = unique([block.Spec1;criticalSpecies]);
    Inlet.Flow1.T = FuelTempIn;
    Inlet.Mixed = Inlet.Flow1;
    S2C = block.Fuel.H2O/(block.Fuel.CH4+.5*block.Fuel.CO);
    if S2C<block.Steam2Carbon %add anode recirculation
        Inlet.Flow1.T = 300;%mixing provides humidification & pre-heat
        Inlet.Mixed.CH4 = block.Fuel.CH4*block.FuelFlowInit;
        eWGS = .7; %initial guess of effective CO conversion
        r = 0.5; %Initial guess of anode recirculation
        dr = 1e-5;
        error = 1;
        % Inlet = (Inlet + generated - consumed)*r  + New, thus inlet = New/(1-r) + (generated - consumed)*r/(1-r)
        while abs(error)>1e-6
            Inlet.Mixed.CO = block.Fuel.CO*block.FuelFlowInit/(1-r) + (block.Fuel.CH4 - eWGS*(block.Fuel.CH4+block.Fuel.CO))*block.FuelFlowInit*r/(1-r);
            Inlet.Mixed.H2O = block.Fuel.H2O*block.FuelFlowInit/(1-r) + (block.Cells*sum(Current)/(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*eWGS)*block.FuelFlowInit)*r/(1-r);
            S2C = Inlet.Mixed.H2O/(Inlet.Mixed.CH4 + 0.5*Inlet.Mixed.CO);
            error = block.Steam2Carbon - S2C;
            r2 = r+dr;
            COin2 = block.Fuel.CO*block.FuelFlowInit/(1-r2) + (block.Fuel.CH4 - eWGS*(block.Fuel.CH4+block.Fuel.CO))*block.FuelFlowInit*r2/(1-r2);
            H2Oin2 = block.Fuel.H2O*block.FuelFlowInit/(1-r2) + (block.Cells*sum(Current)/(2*F*1000) - (block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*eWGS)*block.FuelFlowInit)*r2/(1-r2);
            S2C2 = H2Oin2/(Inlet.Mixed.CH4 + 0.5*COin2);
            dSdr = (S2C2 - S2C)/dr;
            r = r + max(-.5*r,min((1-r)/2,error/dSdr));
        end
        block.Recirc.Anode = r;
        Inlet.Mixed.CO2 = block.Fuel.CO2*block.FuelFlowInit/(1-r) + eWGS*(block.Fuel.CH4+block.Fuel.CO)*block.FuelFlowInit*r/(1-r);
        Inlet.Mixed.H2 = block.Fuel.H2*block.FuelFlowInit/(1-r) + ((3*block.Fuel.CH4 + (block.Fuel.CH4 + block.Fuel.CO)*eWGS)*block.FuelFlowInit - block.Cells*sum(Current)/(2*F*1000))*r/(1-r);
        for i = 1:1:length(block.Spec1)
            if ~ismember(block.Spec1{i},criticalSpecies)
                Inlet.Mixed.(block.Spec1{i}) = block.Fuel.(block.Spec1{i})*block.FuelFlowInit/(1-r);
            end
            Flow1Outlet.(block.Spec1{i}) = Inlet.Mixed.(block.Spec1{i}) - Inlet.Flow1.(block.Spec1{i});
        end
        %%find resulting temperature of mixture
        errorT = 1;
        Inlet.Mixed.T = FuelTempIn;
        Flow1Outlet.T = block.TpenAvg + .5*block.deltaTStack;
        Hin = enthalpy(Inlet.Flow1);
        Hout = enthalpy(Flow1Outlet);
        Hnet = Hin + block.Recirc.Anode*Hout;
        Cp = SpecHeat(Flow1Outlet);
        NetFlowMix = NetFlow(Inlet.Mixed);
        while abs(errorT)>1e-3
            Hmix = enthalpy(Inlet.Mixed);
            errorT = (Hnet-Hmix)/(Cp*NetFlowMix);
            Inlet.Mixed.T = Inlet.Mixed.T + errorT;
        end 
    else
        block.Recirc.Anode = 0;
    end
    block.T.FuelMix = Inlet.Mixed.T;
    
    Inlet.Flow2Pout = block.Flow2_Pout;
    Inlet.Flow1Pout = block.Flow1_Pout;
    if ~isfield(block,'OxidantUtilization')
        if block.ClosedCathode
            block.OxidantUtilization = 1;
        elseif strcmp(block.Reformer,'internal')
            block.OxidantUtilization = .33;
        else
            block.OxidantUtilization = .1;
        end
    end
    block.AirFlow = block.Cells*sum(Current)/(4*F*block.Flow2Spec.O2)/1000/block.OxidantUtilization;%kmol of oxidant

    block.Spec2 = fieldnames(block.Flow2Spec);
    if block.ClosedCathode
        block.Spec2 = {}; %no cathode flow states
        criticalSpecies = {};
    else
        switch block.FCtype
            case 'SOFC'
                criticalSpecies = {'O2';'N2';};
            case 'MCFC'
                criticalSpecies = {'CO2';'H2O';'O2';'N2';};
        end
    end
    Inlet = InletFlow(block,Inlet);
    for i = 1:1:length(criticalSpecies)
        if ~ismember(block.Spec2,criticalSpecies{i})
            Inlet.Flow2.(criticalSpecies{i}) = 0;
        end
    end
    block.Spec2 = unique([block.Spec2,criticalSpecies]);
    %% Run Initial Condition
    [Flow1,Flow2,block,Inlet] = solveInitCond(Inlet,block,1);
    
    if strcmp(block.Reformer,'external') || strcmp(block.Reformer,'adiabatic')
        block.Reformer = 'direct'; %external and adiabatic reformers handled in seperate block, after 1st initialization
    end
    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    block.PortNames = {'NetCurrent','Flow1','Flow2','Flow1Pout','Flow2Pout','Flow1Out','Flow2Out','Flow2Pin','Flow1Pin','MeasureVoltage','MeasurePower','MeasureTpen','MeasureTflow1','MeasureTflow2'};
    block.NetCurrent.type = 'in';
    block.NetCurrent.IC = sum(block.Current.H2 + block.Current.CO);

    block.Flow1.type = 'in';
    block.Flow1.IC = Inlet.Flow1; 

    block.Flow2.type = 'in';
    block.Flow2.IC = Inlet.Flow2;
    
    block.Flow1Pout.type = 'in';
    block.Flow1Pout.IC = Inlet.Flow1Pout;
    block.Flow1Pout.Pstate = []; %identifies the state # of the pressure state if this block has one

    block.Flow2Pout.type = 'in';
    block.Flow2Pout.IC = Inlet.Flow2Pout;
    block.Flow2Pout.Pstate = []; %identifies the state # of the pressure state if this block has one

    block.Flow1Out.type = 'out';
    block.Flow1Out.IC = MergeLastColumn(Flow1.Outlet,block.Flow1Dir,block.Cells);
    
    block.Flow2Out.type = 'out';
    block.Flow2Out.IC  = MergeLastColumn(Flow2.Outlet,block.Flow2Dir,block.Cells);

    block.Flow1Pin.type = 'out';
    block.Flow1Pin.IC = block.Flow1_Pinit;
    block.Flow1Pin.Pstate = length(block.Scale)-1; %identifies the state # of the pressure state if this block has one

    block.Flow2Pin.type = 'out';
    block.Flow2Pin.IC = block.Flow2_Pinit;
    block.Flow2Pin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one

    block.MeasureVoltage.type = 'out';
    block.MeasureVoltage.IC = block.Voltage;

    block.MeasurePower.type = 'out';
    block.MeasurePower.IC = sum((block.Current.H2 + block.Current.CO)*block.Voltage*block.Cells)/1000;%power in kW

    block.MeasureTpen.type = 'out';
    block.MeasureTpen.IC = block.T.Elec;

    block.MeasureTflow1.type = 'out';
    block.MeasureTflow1.IC = block.T.Flow1(block.Flow1Dir(:,end));
    
    block.MeasureTflow2.type = 'out';
    block.MeasureTflow2.IC = block.T.Flow2(block.Flow2Dir(:,end));

    block.P_Difference = {'Flow2Pin','Flow2Pout'; 'Flow1Pin', 'Flow1Pout';};

    for i = 1:1:length(block.PortNames)
        if length(block.connections)<i || isempty(block.connections{i})
            block.(block.PortNames{i}).connected={};
        else
            if ischar(block.connections{i})
                block.(block.PortNames{i}).connected = block.connections(i);
            else
                block.(block.PortNames{i}).IC = block.connections{i};
                block.(block.PortNames{i}).connected={};
            end
        end
    end
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    block.Specification = 'current density';%converge only to match current density from controller
    block.Recirc.Anode = 0; %after 1st initialization recirculation is handled in controller, valve & mixing volume
    Inlet.Mixed = Inlet.Flow1; %mixing moved to seperate mixing block
    Flow1New = fieldnames(Inlet.Flow1);
    Flow1All = unique([block.Spec1;Flow1New]);
    Flow1All = Flow1All(~strcmp('T',Flow1All));
    for i = 1:1:length(Flow1All)
        if ~ismember(Flow1All{i},Flow1New)
            Inlet.Flow1.(Flow1All{i})=0;
        end
    end
    block.Spec1 = Flow1All;

    if block.ClosedCathode
        block.Spec2 = {}; %no cathode flow states
    else
        Flow2New = fieldnames(Inlet.Flow2);
        Flow2All = unique([block.Spec2;Flow2New]);
        Flow2All = Flow2All(~strcmp('T',Flow2All));
        for i = 1:1:length(Flow2All)
            if ~ismember(Flow2All{i},Flow2New)
                Inlet.Flow2.(Flow2All{i})=0;
            end
        end
        block.Spec2 = Flow2All;
    end
    block.Flow1_Pinit = Inlet.Flow1Pout + block.Flow1Pdrop;
    block.Flow2_Pinit = Inlet.Flow2Pout + block.Flow2Pdrop;
    block.Flow2Pout.IC = Inlet.Flow2Pout;
    block.Flow1Pout.IC = Inlet.Flow1Pout;
    %%--%%
    [Flow1,Flow2,block,~] = solveInitCond(Inlet,block,2);
    %%%
    block.Flow1Pin.Pstate = length(block.Scale)-1; %identifies the state # of the pressure state if this block has one
    block.Flow2Pin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    block.Flow2Out.IC  = MergeLastColumn(Flow2.Outlet,block.Flow2Dir,block.Cells);
    block.Flow1Out.IC = MergeLastColumn(Flow1.Outlet,block.Flow1Dir,block.Cells);
    block.Flow2Pin.IC = block.Flow2_Pinit;
    block.Flow1Pin.IC = block.Flow1_Pinit;
    block.MeasureCurrent.IC = sum(block.Current.H2 + block.Current.CO);
    block.MeasurePower.IC = sum((block.Current.H2 + block.Current.CO)*block.Voltage*block.Cells)/1000;%power in kW
    block.MeasureTpen.IC = block.T.Elec;
    block.MeasureTflow2.IC = block.T.Flow2(block.Flow2Dir(:,end));
    block.MeasureTflow1.IC = block.T.Flow1(block.Flow1Dir(:,end));
    block.HumidifiedFuelTemp.IC = block.T.FuelMix;
end

function [Flow1,Flow2,block,Inlet] = solveInitCond(Inlet,block,firstSolve)
global Tags F
error = 1;
Tol = 1e-3;
count = 1;
while abs(error)>Tol %iterate to reach target current density, voltage or power density
    Flow2 = FCin2Out(block.T.Flow2,Inlet.Flow2,block.Flow2Dir, block.FCtype,block.Cells,block.Current,[],'cathode');
%     SinglePassUtilization = (sum(block.Current)*block.Cells/(2000*F))/(4*Inlet.Mixed.CH4+Inlet.Mixed.CO + Inlet.Mixed.H2);
    [block,Inlet.Mixed,Flow1,Flow3] = KineticCoef(block,Inlet,((firstSolve==1) && (count==1)));%% solve for kinetic reaction coefficient which results in this outlet condition (match R_CH4 & R_WGS)
    Offset = 0;
    if count==1 && firstSolve==1
        [block.Tstates,block.HTcond,block.HTconv,block.HTrad]= SteadyTemps(block,Inlet.Mixed,Inlet.Flow2);
    else
        [~, Y] = ode15s(@(t,y) DynamicTemps(t,y,block,Flow1,Flow2,Flow3,Inlet), [0, 1e4], block.Tstates);
        block.Tstates = Y(end,:)';
        if firstSolve==1 && abs(mean(block.Tstates(2*block.nodes+1:3*block.nodes))-block.TpenAvg)>10 %temperature in solve dynamic is diverging to much (1300K), and messing up reforming solution (this is a temporary fix
            Offset = (mean(block.Tstates(2*block.nodes+1:3*block.nodes))-block.TpenAvg);
        end
    end
    %organize temperatures
    block.T.Flow2 =  block.Tstates(1*block.nodes+1:2*block.nodes) - Offset;
    block.T.Elec =  block.Tstates(2*block.nodes+1:3*block.nodes) - Offset;
    block.T.Flow1 =  block.Tstates(3*block.nodes+1:4*block.nodes) - Offset;
    block.T.FuelMix = Inlet.Mixed.T;
    switch block.Reformer
        case 'internal'
            block.T.Flow3 = block.Tstates(5*block.nodes+1:6*block.nodes) - Offset;
        case 'adiabatic'
            block.T.Flow3 =  Flow3.Outlet.T; %1 temperature state of reformer before anode outlet
    end
    T = block.TpenAvg + (block.T.Elec-mean(block.T.Elec)); %assume you will get to the desired temperature (this avoids oscilations in voltage and helps convergence
    FuelCellNernst(Flow1,Flow2,block.Current,T,block.Flow2_Pinit,block);
    %% calculate the change in current to converge to the desired power density, voltage, or current
    OldVoltage = block.Voltage;
    block.Voltage = sum(Tags.(block.name).nVoltage)/block.nodes;
    block.Current.CO = Tags.(block.name).I_CO;
    block.Current.H2 = Tags.(block.name).I_H2;
    netCurrent = block.Current.H2+block.Current.CO;
    TotCurrent = sum(netCurrent);
    localR = Tags.(block.name).LocalOhmic./(block.Current.H2+block.Current.CO);
    if strcmp(block.Specification,'power density')
        error = (block.RatedStack_kW - Tags.(block.name).Power)/block.RatedStack_kW;
        if count>1
            if firstSolve==1 && error>1e-3
                dP_di = max(.7*(Tags.(block.name).Power - OldCurrent*OldVoltage*block.Cells/1000)/(TotCurrent - OldCurrent),.7*Tags.(block.name).Power/(TotCurrent));%change in power with change in current
            else
                dP_di = max(1.15*((Tags.(block.name).Power - OldCurrent*OldVoltage*block.Cells/1000)/(TotCurrent - OldCurrent)),.7*Tags.(block.name).Power/(TotCurrent));
            end
            scale = 1+ error*block.RatedStack_kW/dP_di/TotCurrent;
        else % first time through
            scale = (block.RatedStack_kW*1000/block.Cells/block.Voltage)/TotCurrent; %total current it should have at this new voltage/ total current specified right now
        end
    elseif strcmp(block.Specification,'voltage')
        Tol = 1e-3;
        error = (block.SpecificationValue - block.Voltage)/block.Voltage;
        if count>1
            dV_di = -Tags.(block.name).ASR/block.A_Node/100^2;
            scale = 1 + (block.SpecificationValue - block.Voltage)/dV_di/TotCurrent;
        else % first time through
            localR = 3*localR;
            scale =1+sum((Tags.(block.name).nVoltage-block.SpecificationValue)./localR)/TotCurrent;
        end
        block.Cells = ceil((block.RatedStack_kW*1000/block.SpecificationValue)/(scale*TotCurrent)); %re-calculate the # of cells
    elseif strcmp(block.Specification,'current density')
        scale = Inlet.NetCurrent/TotCurrent;
        error = (TotCurrent - Inlet.NetCurrent)/TotCurrent;
    end
    
    OldCurrent = sum(netCurrent);
    ratio = block.Current.CO./netCurrent;
    netCurrent = redistributeCurrent(netCurrent,scale,Tags.(block.name).nVoltage,localR,block.Voltage); %% start by maintaining same current in each row, then allow row voltages to balance (prevents fuel starvation in a row during initialization)
    block.Current.CO = ratio.*netCurrent;
    block.Current.H2 = (1-ratio).*netCurrent;
    TotCurrent = sum(netCurrent);
    
    if firstSolve ==1 %solving block to convergence without other blocks or controller
        block.R_CH4 = scale*block.R_CH4;
        block.R_WGS = scale*block.R_WGS;
        switch block.Reformer
            case {'internal';'adiabatic'}
                block.R_CH4ref = scale*block.R_CH4ref;% fuel flow scales with current so assume reforming will
                block.R_WGSref = scale*block.R_WGSref;
        end
        if strcmp(block.CoolingStream,'cathode')%air flow balances heat generation to maintain deltaT, heat transfer to anode and any fuel reforming is accounted for
            block.FuelFlowInit  = TotCurrent/(2*F*1000)/(block.FuelUtilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2))*block.Cells; % Fresh fuel flow rate,  current/(2*F*1000) = kmol H2
            k = block.Flow2Dir(:,end);
            dTerror = (mean((block.Tstates(k+block.nodes))- Inlet.Flow2.T)/block.deltaTStack-1);
            block.AirFlow = block.AirFlow*(1 + dTerror)*scale^2;
            TavgError = (block.TpenAvg-mean(block.Tstates(2*block.nodes+1:3*block.nodes)))/block.deltaTStack;
            block.StackCathTin = block.StackCathTin + (TavgError + 0.75*dTerror)*block.deltaTStack;
        elseif strcmp(block.CoolingStream,'anode') %oxidant flow rate determined by current, fuel flow rate is now determined by thermal balancing
            block.AirFlow = block.Cells*TotCurrent/(4000*F*block.Flow2Spec.O2)/block.OxidantUtilization;%kmol of oxidant
            if block.ClosedCathode %%energy balance
                Hin1 = enthalpy(Flow1.Inlet);
                Hout1 = enthalpy(Flow1.Outlet);
                Hin2 = enthalpy(Flow2.Inlet);
                Hout2 = enthalpy(Flow2.Outlet);
                Hin3 = enthalpy(Flow3.Inlet);
                Hout3 = enthalpy(Flow3.Outlet);
                Power = block.Voltage*netCurrent/1000; %cell power in kW
                Qimbalance = sum((Hin2 - Hout2) + (Hin1 - Hout1) + (Hin3 - Hout3) - Power);
                h = enthalpy(mean(block.T.Flow3),{'H2','H2O','O2','CO','CO2','CH4'});
                Qreform = (h.CO+3*h.H2-h.CH4-h.H2O) + 0.8*(h.CO2+h.H2-h.CO-h.H2O); %kW of cooling per kmol of fuel
                ExtraFuel = 0.25*Qimbalance*block.Cells/Qreform/block.Fuel.CH4;
                error = max(abs(error),abs(Qimbalance/sum(Power)));
            else
                ExtraFuel = 0;
                %need to do something with recirculation
            end
            block.FuelFlowInit  = block.Cells*TotCurrent/(2000*F)/(block.FuelUtilization*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2)); % re-calculate with revised current
            block.FuelFlowInit = block.FuelFlowInit + ExtraFuel;
            block.FuelUtilization = block.Cells*TotCurrent/(2000*F)/(block.FuelFlowInit*(4*block.Fuel.CH4+block.Fuel.CO+block.Fuel.H2));
            % change steam to carbon to affect deltaT?
        end
        Inlet = InletFlow(block,Inlet);
    end
    count= count+1;
end
Flow2 = FCin2Out(block.T.Flow2,Inlet.Flow2,block.Flow2Dir, block.FCtype,block.Cells,block.Current,[],'cathode');
block.Pfactor1 = NetFlow(Inlet.Flow1)/block.Flow1Pdrop;
block.Pfactor2 = NetFlow(Inlet.Flow2)/block.Flow2Pdrop;
block = Set_IC(block,Flow1,Flow2,Flow3);

function Inlet = InletFlow(block,Inlet) %only used 1st time through initialization (before we know what is connected to inlet
% Anode
for i = 1:1:length(block.Spec1)
    if isfield(block.Fuel,block.Spec1{i})
        Inlet.Flow1.(block.Spec1{i}) = block.Fuel.(block.Spec1{i})*block.FuelFlowInit;%flow rate of every species entering the anode (or reformer if there is one)
    else Inlet.Flow1.(block.Spec1{i}) = 0;
    end
end
%Cathode
Inlet.Flow2.T = block.StackCathTin;
switch block.FCtype
    case 'SOFC'
        if block.ClosedCathode
            Inlet.Flow2.O2 = block.Flow2Spec.O2*block.AirFlow;
        else
            for i = 1:1:length(block.Spec2)
                if isfield(block.Flow2Spec,block.Spec2{i})
                    Inlet.Flow2.(block.Spec2{i}) = block.Flow2Spec.(block.Spec2{i})*block.AirFlow;
                else Inlet.Flow2.(block.Spec2{i}) = 0;
                end
            end
        end
    case 'MCFC' %recalculate cathode inlet species for MCFC (this is an estimate assuming the 100% of non-recirculated anode gas is oxidized and fed to the cathode)
        Inlet.Flow2.CO2 = (Inlet.Mixed.CH4+Inlet.Mixed.CO+Inlet.Mixed.CO2) + sum(block.Current.H2+block.Current.CO)/(2*F*1000)*block.Cells;
        Inlet.Flow2.H2O = (4*Inlet.Mixed.CH4+Inlet.Mixed.CO+Inlet.Mixed.H2+Inlet.Mixed.H2O);
        nonCO2_H2O = (block.AirFlow - Inlet.Flow2.CO2 - Inlet.Flow2.H2O);
        for i = 1:1:length(block.Spec2)
            if isfield(block.Flow2Spec,block.Spec2{i})
                if strcmp(block.Spec2{i},'CO2')||strcmp(block.Spec2{i},'H2O')
                    Inlet.Flow2.(block.Spec2{i}) = Inlet.Flow2.(block.Spec2{i}) + block.Flow2Spec.(block.Spec2{i})*nonCO2_H2O;
                else
                    Inlet.Flow2.(block.Spec2{i}) = block.Flow2Spec.(block.Spec2{i})*nonCO2_H2O;
                end
            else Inlet.Flow2.(block.Spec2{i}) = 0;
            end
        end
end

function block = Set_IC(block,Flow1,Flow2,Flow3)
global Ru
Cp.cath = SpecHeat(Flow2.Inlet);
Cp.an = SpecHeat(Flow1.Outlet);
if ~isempty(Flow3)
    Cp.ref = SpecHeat(Flow3.Outlet);
end
switch block.Reformer
    case 'internal'
        NumOfStates = (6 + 2*length(block.Spec1) + length(block.Spec2) + 1)*block.nodes + 2; % 6 temperatures, anode & cathode & reformer & current at each node and 2 states for anode/cathode pressure
  case {'direct';'external';'adiabatic'}
        NumOfStates = (5 + length(block.Spec1) + length(block.Spec2) + 1)*block.nodes +2; % 5 temperatures, anode & cathode & current at each node and 2 states for anode/cathode pressure
end
block.IC = ones(NumOfStates,1); %
block.tC = block.IC; % time constant for derivative dY
block.Scale = block.IC;
switch block.Reformer
    case 'internal'
        n = 6*block.nodes;
        block.tC(5*block.nodes+1:6*block.nodes) = (block.Vol_flow3*Cp.ref*block.Flow1_Pinit./(Ru*block.T.Flow3));
    case {'direct';'external';'adiabatic'}
        n = 5*block.nodes;
end
block.Scale = block.Tstates(1:n);%temperature (K)

block.tC(1:block.nodes) = (block.Mass_plate2*block.C_plate2);
block.tC(1+block.nodes:2*block.nodes) = (block.Vol_flow2*Cp.cath*block.Flow2_Pinit./(Ru*block.T.Flow2));
block.tC(2*block.nodes+1:3*block.nodes) = (block.Vol_Elec*block.Density_Elec*block.C_Elec);
block.tC(3*block.nodes+1:4*block.nodes) = (block.Vol_flow1*Cp.an*block.Flow1_Pinit./(Ru*block.T.Flow1));
block.tC(4*block.nodes+1:5*block.nodes) = (block.Mass_plate1*block.C_plate1);
block.tC(1:n) = block.tC(1:n)-diag(block.HTconv)-diag(block.HTcond); %this accounts for the change in HT as temperature of the control volume changes. The change in HT helps balance the energy equation more than the change in enthalpy leaving.

for i = 1:1:length(block.Spec2)
    block.tC(n+1:n+block.nodes) = (block.Vol_flow2*block.Flow2_Pinit)./(block.T.Flow2*Ru);  % cathode 
    if any(Flow2.Outlet.(block.Spec2{i})==0)
        block.IC(n+1:n+block.nodes) = Flow2.Outlet.(block.Spec2{i})./NetFlow(Flow2.Outlet);
        block.Scale(n+1:n+block.nodes) = NetFlow(Flow2.Outlet); n = n+block.nodes; %cathode flows
    else
        block.Scale(n+1:n+block.nodes) = Flow2.Outlet.(block.Spec2{i}); n = n+block.nodes; %cathode flows
    end
end

for i = 1:1:length(block.Spec1)
    X = Flow1.Outlet.(block.Spec1{i})./NetFlow(Flow1.Outlet);%concentration
    block.tC(n+1:n+block.nodes) = (block.Vol_flow1*block.Flow1_Pinit)./(block.T.Flow1*Ru); %anode
    if any(X<.01) %concentration less than 1%
        block.IC(n+1:n+block.nodes) = X;
        block.Scale(n+1:n+block.nodes) = NetFlow(Flow1.Outlet); %anode flow
    else
        block.Scale(n+1:n+block.nodes) = Flow1.Outlet.(block.Spec1{i}); %individual species flow
    end   
    n = n+block.nodes;
end

switch block.Reformer
    case 'internal'
        for i = 1:1:length(block.Spec1)
            block.tC(n+1:n+block.nodes) = (block.Vol_flow3*block.Flow1_Pinit)./(block.T.Flow3*Ru); % reformer
            if any(Flow3.Outlet.(block.Spec1{i})==0)
                block.IC(n+1:n+block.nodes) = Flow3.Outlet.(block.Spec1{i})./NetFlow(Flow3.Outlet);
                block.Scale(n+1:n+block.nodes) = NetFlow(Flow3.Outlet); n = n+block.nodes; %anode flows
            else
                block.Scale(n+1:n+block.nodes) = Flow3.Outlet.(block.Spec1{i}); n = n+block.nodes; %reformer flows
            end
        end
end
block.tC(n+1:n+block.nodes) = block.nodes/100;%  %current changing for voltage balance
block.Scale(n+1:n+block.nodes) = (block.Current.H2+block.Current.CO);  n = n+block.nodes; %current

block.tC(n+1) = (block.Vol_flow1*block.nodes*block.Cells); %pressure
block.tC(n+2) = (block.Vol_flow2*block.nodes*block.Cells);  %pressure
block.Scale(n+1) = block.Flow1_Pinit;%pressure
block.Scale(n+2) = block.Flow2_Pinit;%pressure


function dY = DynamicTemps(t,Y,block,Flow1,Flow2,Flow3,Inlet)
global F Ru
dY = 0*Y;
nodes = block.nodes;
if isfield(block,'tC')
    tC = block.tC(1:6*nodes);
else
    Cp_2 = 33;% kJ/kmol*K
    Cp_1 = 42;% kJ/kmol*K
    tC(1:nodes,1) = (block.Mass_plate2*block.C_plate2);
    tC(1+nodes:2*nodes,1) = (block.Vol_flow2*Cp_2*block.Flow2_Pinit./(Ru*block.T.Flow2));
    tC(2*nodes+1:3*nodes,1) = (block.Vol_Elec*block.Density_Elec*block.C_Elec);
    tC(3*nodes+1:4*nodes,1) = (block.Vol_flow1*Cp_1*block.Flow1_Pinit./(Ru*block.T.Flow1));
    tC(4*nodes+1:5*nodes,1) = (block.Mass_plate1*block.C_plate1);
    if isfield(block,'Vol_flow3')
        Cp_3 = 42;% kJ/kmol*K
        tC(5*nodes+1:6*nodes,1) = (block.Vol_flow3*Cp_3*block.Flow1_Pinit./(Ru*block.T.Flow3));
    end
end

h = enthalpy(Y(1+2*nodes:3*nodes),{'H2','H2O','O2','CO','CO2'});
Power = block.Voltage*(block.Current.H2+block.Current.CO)/1000; %cell power in kW
Qgen = block.Current.H2/(2000*F).*(h.H2+.5*h.O2-h.H2O) + block.Current.CO/(2000*F).*(h.CO+.5*h.O2-h.CO2)-Power;%kW of heat generated by electrochemistry (per node & per cell)
switch block.FCtype%ion transport across membrane (total enthalpy)
    case 'SOFC'
        Qion = (block.Current.H2+block.Current.CO)/(4000*F).*h.O2; %O2 ion crossing over (kW)
    case 'MCFC'
        Qion = (block.Current.H2+block.Current.CO)*(1/(4000*F).*h.O2 + 1/(2000*F).*h.CO2);% O2 & CO2 ion crossing over
end

QT = block.HTconv*Y + block.HTcond*Y + block.HTrad*(Y.^4);

Flow2.Outlet.T = Y(nodes+1:2*nodes);
for j = 1:1:length(block.Flow2Dir(1,:));%1:columns
    k = block.Flow2Dir(:,j);
    if j~=1
        Flow2.Inlet.T(k,1) = Flow2.Outlet.T(kprev);
    end
    kprev = k;
end

Flow1.Outlet.T = Y(3*nodes+1:4*nodes);
for j = 1:1:length(block.Flow1Dir(1,:));%1:columns
    k = block.Flow1Dir(:,j);
    if j~=1
        Flow1.Inlet.T(k,1) = Flow1.Outlet.T(kprev);
    end
    kprev = k;
end

%energy flows & sepcific heats
Hout2 = enthalpy(Flow2.Outlet);
Hin2 = enthalpy(Flow2.Inlet);
Hout1 = enthalpy(Flow1.Outlet);
HfreshFuel = enthalpy(Inlet.Flow1);

if block.Recirc.Anode>0 % Only during the first run with unhumidified fuel, find fuelmix temperature
    error2 = 1;
    Cp = SpecHeat(Flow1.Outlet); 
    Cp = Cp(end);
    Hmix = HfreshFuel + sum(Hout1(block.Flow1Dir(:,end)))*block.Recirc.Anode*block.Cells;
    netflow = NetFlow(Inlet.Mixed);
    while abs(error2)>1e-4
        error2 = (Hmix - enthalpy(Inlet.Mixed))./(Cp*netflow);                             %Adjusting the error in temperature based on known enthalpy and specific heat of the cold side
        Inlet.Mixed.T = Inlet.Mixed.T + .75*error2;                                   %Subtraction of a portion of the T_error from cold outlet temp to get closer to the actual temp
    end
end

switch block.Reformer
    case 'internal'
        Flow3.Outlet.T = Y(5*nodes+1:6*nodes);
        k = block.Flow3Dir(:,1);
        if block.Recirc.Anode>0
            Flow3.Inlet.T(k,1) = Inlet.Mixed.T;
        else Flow3.Inlet.T(k,1) = Inlet.Flow1.T;
        end
        for j = 1:1:length(block.Flow3Dir(1,:));%1:columns
            k = block.Flow3Dir(:,j);
            if j~=1
                Flow3.Inlet.T(k,1) = Flow3.Outlet.T(kprev);
            end
            kprev = k;
        end
        Hin3 = enthalpy(Flow3.Inlet);
        Hout3 = enthalpy(Flow3.Outlet);
        
        k = block.Flow1Dir(:,1);
        k2 = block.Flow3Dir(:,end);
        Flow1.Inlet.T(k,1) = Flow3.Outlet.T(k2,1);
    case {'adiabatic';'direct';'external';}
        k = block.Flow1Dir(:,1);
        Flow1.Inlet.T(k,1) = Inlet.Flow1.T;
end
Hin1 = enthalpy(Flow1.Inlet);

if block.ClosedCathode %%energy balance
    Qimbalance = sum((Hin2(block.Flow2Dir(:,1))) - sum(Hout2(block.Flow2Dir(:,end)))) + sum(Hin3(block.Flow3Dir(:,1)))  - sum(Hout1(block.Flow1Dir(:,end))) - sum(Power);
    Power = Power + Qimbalance*Power./sum(Power);
end

dY(1:nodes)= QT(1:nodes)./tC(1:nodes);  %Ox Sep Plate
dY(1+nodes:2*nodes)= (QT(1+nodes:2*nodes) + Hin2 - Hout2 - Qion)./tC(1+nodes:2*nodes); %Cathode
dY(1+2*nodes:3*nodes)= (QT(1+2*nodes:3*nodes)+ Qgen)./tC(2*nodes+1:3*nodes); %Electrolyte Plate
dY(1+3*nodes:4*nodes)= (QT(1+3*nodes:4*nodes) + Hin1 - Hout1 + Qion - Power - Qgen)./tC(1+3*nodes:4*nodes);  %Anode
dY(1+4*nodes:5*nodes)= QT(1+4*nodes:5*nodes)./tC(4*nodes+1:5*nodes);  %Fuel Sep Plate
switch block.Reformer
    case 'internal'
        dY(1+5*nodes:6*nodes)= (block.RefSpacing*QT(1+5*nodes:6*nodes) + Hin3 - Hout3)./tC(1+5*nodes:6*nodes);  %Fuel Reformer Channels
end