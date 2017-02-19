function block = InitializeOxidizer(varargin)
global Tags
% an oxidizer with n inlet flows 
% Complete combustion is assumed for all fuels (if enough O2 present)
% two or more (2+) inlets: Flows and outletl pressure
% Two (2) outlets: Flow out and pressure at the inlet
% Many (2+n) states: Outlet Temperature, Species, Inlet Pressure
block = varargin{1};
if length(varargin)==1 % first initialization
    % CH4+1.5O2 --> CO+ 2H2O
    % CO + .5O2 --> CO2
	% H2 + .5O2 --> H2O
    
    block.spec = {'CH4';'CO';'CO2';'H2';'H2O';'O2'};
    block.Scale = block.InitialFlowOut.T;
    SpecNew = fieldnames(block.InitialFlowOut);
    SpecAll = unique([block.spec;SpecNew]);
    block.spec = SpecAll(~strcmp('T',SpecAll));
    for i = 1:1:length(block.spec)
        if ~ismember(block.spec{i},SpecNew)
            block.InitialFlowOut.(block.spec{i}) = 0;
        end
        block.Scale(end+1,1) = block.InitialFlowOut.(block.spec{i});
    end
    
    block.Pdrop = 2; % presure drop across mixing volume (kPa)
    block.Scale(end+1,1) = 101+block.Pdrop;
    block.IC = ones(length(block.Scale),1);
    
    block.PortNames = {};
    for j = 1:1:block.inlets
        name = strcat('Inlet',num2str(j));
        block.PortNames(end+1) = cellstr(name);
        block.(name).type = 'in';
        block.(name).IC.T = 800;
    end
    
    block.PortNames(end+1:end+4) = {'Pout';'Flow';'Pin';'MeasureT'}; %inlet port names must be first
    
    block.Pout.type = 'in';
    block.Pout.IC = 101;
    block.Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.Flow.type = 'out';
    block.Flow.IC = block.InitialFlowOut;
    
    block.Pin.type = 'out';
    block.Pin.IC = block.Pout.IC+block.Pdrop;
    block.Pin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    
    block.MeasureT.type = 'out';
    block.MeasureT.IC = block.InitialFlowOut.T;
    
    block.P_Difference = {'Pin','Pout'};
    %no dMdP or mFlow (fixed pressure drop)

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
    Tags.(block.name).EquivelanceRatio = 0.9;
    Tags.(block.name).Temperature = block.InitialFlowOut.T;
    Tags.(block.name).MassFlow = MassFlow(block.InitialFlowOut);
end
if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    n = block.inlets;
    inlets = fieldnames(Inlet);
    Pin = Inlet.Pout+block.Pdrop;
    %mix flows
    NetIn =[];
    H_in = 0;
    for j = 1:1:n
        H_in = H_in + enthalpy(Inlet.(inlets{j}));
        SpecNew = fieldnames(Inlet.(inlets{j}));
        SpecAll = unique([block.spec;SpecNew]);
        block.spec = SpecAll(~strcmp('T',SpecAll));
        for i = 1:1:length(block.spec)
            if ~isfield(NetIn,block.spec{i})
                NetIn.(block.spec{i}) = 0;
            end
            if ismember(block.spec{i},SpecNew)
                NetIn.(block.spec{i}) = NetIn.(block.spec{i}) + Inlet.(inlets{j}).(block.spec{i});
            end
        end
    end
    
    ReactMix.T = block.Scale(1);
    R.CH4 = NetIn.CH4;
    R.H2 = NetIn.H2;
    if ReactMix.T>1700
        Kp = 0.4589*(ReactMix.T/3000)^14.772;
    else
        Kp = max(0,(ReactMix.T-300)*(6.78e-5/1400));
    end
    X_CO_X_CO2 = Kp/((Pin/101)^.5*(NetIn.O2/NetFlow(NetIn))); %estimate ratio of CO to CO2
    R.CO = ((NetIn.CO+R.CH4)-X_CO_X_CO2*NetIn.CO2)/(1+X_CO_X_CO2);

    sumR = 1.5*R.CH4 + 0.5*R.CO + 0.5*R.H2;
    phi = sumR/NetIn.O2;% O2 necessary/O2 supplied
    r = fieldnames(R);
    lean = min(1,NetIn.O2/sumR);
    for j = 1:1:length(r)
        R.(r{j}) = R.(r{j})*lean; %rich combustion limited by O2
    end
    
    for i = 1:1:length(block.spec)
        if strcmp(block.spec{i},'CH4')
            ReactMix.CH4 = NetIn.CH4 - R.CH4;
        elseif strcmp(block.spec{i},'CO')
            ReactMix.CO = NetIn.CO + R.CH4- R.CO;
        elseif strcmp(block.spec{i},'CO2')
            ReactMix.CO2 = NetIn.CO2 + R.CO;
        elseif strcmp(block.spec{i},'H2')
            ReactMix.H2 = NetIn.H2  - R.H2;
        elseif strcmp(block.spec{i},'H2O')
            ReactMix.H2O = NetIn.H2O + 2*R.CH4 + R.H2;
        elseif strcmp(block.spec{i},'O2')
            ReactMix.O2 = NetIn.O2 - 1.5*R.CH4 - .5*R.CO - .5*R.H2;
        else
            ReactMix.(block.spec{i}) = NetIn.(block.spec{i});
        end
    end
    Flow = NetFlow(ReactMix);
    Terror = 100;
    while abs(Terror)>1e-2
        Cp = SpecHeat(ReactMix);
        Terror = (H_in-enthalpy(ReactMix))/(Cp*Flow);
        ReactMix.T = ReactMix.T+Terror;
    end
    block.Flow.IC = ReactMix;
    block.Pfactor = NetFlow(ReactMix)/block.Pdrop;
    %no dMdP or mFlow to update (fixed pressure drop)
    block.Pout.IC = Inlet.Pout; 
    block.Pin.IC = Inlet.Pout+block.Pdrop;
    block.MeasureT.IC = ReactMix.T;
    block.Scale = ReactMix.T;
    block.IC = 1;
    for i = 1:1:length(block.spec)
        if ReactMix.(block.spec{i})==0
            block.IC(end+1,1) = ReactMix.(block.spec{i})/Flow;
            block.Scale(end+1,1) = Flow;
        else
            block.IC(end+1,1) = 1;
            block.Scale(end+1,1) = ReactMix.(block.spec{i});
        end
    end
    block.Scale(end+1,1) = Inlet.Pout+block.Pdrop;
    block.IC(end+1,1) = 1;
    Tags.(block.name).EquivelanceRatio = phi;
    Tags.(block.name).Temperature = ReactMix.T;
    Tags.(block.name).MassFlow = MassFlow(ReactMix);
end