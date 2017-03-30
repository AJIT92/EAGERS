function block = InitializeMixingVolume(varargin)
% a simple adiabatic mixing model mixing two or more inlet flows 
% Multiple (n+1) inlets: Flows (consisting of Temperature and flow rates of individual species) and P out 
% Two (2) outlets: Flow , and temperature 
% Many (n+2) states: Temperature, any species under consideration, and pressure
block = varargin{1};
if length(varargin)==1 %first initialization
    block.Pdrop = 2; % presure drop across mixing volume (kPa)
    if ischar(block.SpeciesInit)
        block.SpeciesInit = ComponentProperty(block.SpeciesInit);
    end
    spec = fieldnames(block.SpeciesInit);
    spec = spec(~strcmp(spec,'T'));
    block.Scale = zeros(2+length(spec),1);
    block.IC = ones(length(block.Scale),1);
    block.Scale(1) = block.Tinit; 
    for i = 1:1:length(spec)
        if block.SpeciesInit.(spec{i}) ==0
            block.Scale(i+1,1) = NetFlow(block.SpeciesInit);
            block.IC(i+1,1) = 0;
        else
            block.Scale(i+1,1) = block.SpeciesInit.(spec{i});
        end
    end
    block.Scale(end) = 101 + block.Pdrop;
    block.spec = spec;
   
    block.PortNames = {};
    for i = 1:1:block.inlets
        name = strcat('Inlet',num2str(i));
        block.PortNames(end+1) = cellstr(name);
        block.(name).type = 'in';
    end
    
    block.PortNames(end+1:end+4) = {'Pout';'Outlet';'Temperature';'Pin';}; %inlet port names must be first
    
    block.Pout.type = 'in';
    block.Pout.IC = 101;
    block.Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.Outlet.type = 'out';
    block.Outlet.IC = block.SpeciesInit;
    block.Outlet.IC.T = block.Tinit; 
    
    block.Temperature.type = 'out';
    block.Temperature.IC = block.Outlet.IC.T;
    
    block.Pin.type = 'out';
    block.Pin.IC = block.Pout.IC+block.Pdrop;
    block.Pin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    
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
elseif length(varargin)==2%% Have inlets connected, re-initialize
    Inlet =varargin{2};
    n = block.inlets;     
    NetIn = {};
    for j = 1:1:length(block.spec)
        NetIn.(block.spec{j}) = 0;
    end
    inlets = fieldnames(Inlet);
    for i = 1:1:n
        spec2 = fieldnames(Inlet.(inlets{i}));
        spec2 = spec2(~strcmp('T',spec2));
        for j = 1:1:length(spec2)
            if ismember(spec2{j},block.spec)
                NetIn.(spec2{j}) = NetIn.(spec2{j}) + Inlet.(inlets{i}).(spec2{j});
            else
                block.spec(end+1) = spec2(j);
                NetIn.(spec2{j}) = Inlet.(inlets{i}).(spec2{j});
            end
        end
    end
    
    NetFlowIn = NetFlow(NetIn);
    block.Pfactor = NetFlowIn/block.Pdrop;
    
    H = zeros(n,1);
    Tin = zeros(n,1);
    for i = 1:1:n
        [~,H(i)] = enthalpy(Inlet.(inlets{i}));
        Tin(i) = Inlet.(inlets{i}).T;
    end
    
    Hin = sum(H);
    NetOut = NetIn;
    NetOut.T = mean(Tin);
    Terror = 1;
    while abs(Terror)>.01
        [~,Hout] = enthalpy(NetOut);
        Cp = SpecHeat(NetOut);
        Terror = (Hin-Hout)/(Cp*NetFlowIn);
        NetOut.T = NetOut.T+Terror;
    end
    block.Pout.IC = Inlet.Pout; 
    block.Pin.IC = Inlet.Pout+block.Pdrop;
    block.Outlet.IC = NetOut; 
    block.Temperature.IC = NetOut.T; 
    
    block.Scale = zeros(2+length(block.spec),1);
    block.IC = ones(length(block.Scale),1);
    block.Scale(1) = NetOut.T; 
    for i = 1:1:length(block.spec)
        if NetIn.(block.spec{i}) ==0
            block.Scale(i+1,1) = NetFlow(NetOut);
            block.IC(i+1,1) = 0;
        else
            block.Scale(i+1,1) = NetOut.(block.spec{i});
        end
    end
    block.Scale(end) = block.Pin.IC;
end