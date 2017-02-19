function block = InitializeHeatExchanger(varargin)
%Nodal Heat Exchanger model with 3 states per node: Temperature hot, Plate temperature, temperature cold
% Four (4) inlets: {'Cold Flow','Hot Flow','Cold Pout','Hot Pout'}
% Four (4) outlets: {'Cold Flow','Hot Flow','Cold Pin','Hot Pin'}
% Many (3*n+2) states: {'Cold Flow', 'Plate', 'Hot Flow','Cold Pin','Hot Pin'}
block = varargin{1};
tic
if length(varargin)==1 % first initialization
    %% Load mask parameters (flow direction)
    block.nodes = block.rows*block.columns;
    block = HXmask(block);
    
    if ischar(block.ColdSpecIn)
        block.ColdSpecIn = lookupVal(block.ColdSpecIn);
    end
    block.spec1 = fieldnames(block.ColdSpecIn);
    block.spec1 = block.spec1(~strcmp(block.spec1,'T'));
    for i = 1:1:length(block.spec1)
        ColdIn.(block.spec1{i}) = block.ColdSpecIn.(block.spec1{i});
    end
    ColdIn.T = block.Cold_T_init;
    
    if ischar(block.HotSpecIn)
        block.HotSpecIn = lookupVal(block.HotSpecIn);
    end
    block.spec2 = fieldnames(block.HotSpecIn);
    block.spec2 = block.spec2(~strcmp(block.spec2,'T'));
    for i = 1:1:length(block.spec2)
        HotIn.(block.spec2{i}) = block.HotSpecIn.(block.spec2{i});
    end
    HotIn.T = block.Hot_T_init;
    
    %% Heat Exchanger
    block.PdropCold = 1;                                    % (kPa) pressure drop
    block.PdropHot = 1;                                    % (kPa) pressure drop
    block.h_conv= 50;                                       %Convective heat coefficient (W/m^2*K)
    block.t_Plate = .001;                                   %Thickness of Heat Exchanger plate (m)
    block.Length = block.Vol^(1/3);                         %Length of heat Exchanger (m)
    block.Width = block.Vol^(1/3);                                   %Length of heat Exchanger (m)
    block.Solid_SpecHeat=.48;                               %Specific Heat of solid material (kJ/kg*K)
    block.Solid_CondCoef = 0;%15;                              %Conduction coefficient (Steel) (W/m*K)
    block.Solid_Density=8055;                               %Density of solid material (kg/m^3)
    block.Vol_Solid = min(1/3*block.Vol, block.Mass/block.Solid_Density);
    block.Vol_Cold = (block.Vol-block.Vol_Solid)/2;
    block.Vol_Hot = (block.Vol-block.Vol_Solid)/2;
    
    Target = lookupVal(block.Target);
    [block.Area,~,~] = findArea(ColdIn,HotIn,block.h_conv,Target,block.sizemethod);
    block.Convection = block.h_conv*block.Area/block.nodes/1000; %all heat transfer coefficients converted to kW/K: thus Q = C*(T1-T2) is in kW
    %horizontal
    block.L_node = block.Length/block.columns;
    block.W_node =block.Width/block.rows;
    block.AcondLR = block.Vol_Solid/block.Width/block.columns;
    block.AcondPN = block.Vol_Solid/block.Length/block.rows;
    block.ConductionPN = block.Solid_CondCoef*block.AcondPN/(block.L_node/2)/1000; %heat transfer coefficient between previous and next node of plate
    block.ConductionLR  = block.Solid_CondCoef*block.AcondLR/(block.W_node/2)/1000; %heat transfer coefficient between left and right adjacent nodes of oxidant plate

    [block.Scale , block.HTconv, block.HTcond] = SolveTemps(block,ColdIn,HotIn);
    block.Scale(end+1:end+2,1) = [101+block.PdropCold;101+block.PdropHot;];
    block.IC = ones(3*block.nodes+2,1);

    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    block.PortNames = {'ColdIn','HotIn','ColdPout','HotPout','ColdOut','HotOut','ColdPin','HotPin'};
    block.ColdIn.type = 'in';
    block.ColdIn.IC = ColdIn;
    
    block.HotIn.type = 'in';
    block.HotIn.IC = HotIn;
    
    block.ColdPout.type = 'in';
    block.ColdPout.IC = 101; %Atmospheric pressure
    block.ColdPout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.HotPout.type = 'in';
    block.HotPout.IC = 101; %Atmospheric pressure
    block.HotPout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    
    block.ColdOut.type = 'out';
    block.ColdOut.IC = ColdIn;
    block.ColdOut.IC.T  = mean(block.Scale(block.ColdFlowDir(:,end),1));
    
    block.HotOut.type = 'out';
    block.HotOut.IC = HotIn;
    block.HotOut.IC.T  = mean(block.Scale(2*block.nodes+block.HotFlowDir(:,end),1));
    
    block.ColdPin.type = 'out';
    block.ColdPin.IC  = block.ColdPout.IC+block.PdropCold;
    block.ColdPin.Pstate = 3*block.nodes+1; %identifies the state # of the pressure state if this block has one
    
    block.HotPin.type = 'out';
    block.HotPin.IC  = block.HotPout.IC+block.PdropHot;
    block.HotPin.Pstate = 3*block.nodes+2; %identifies the state # of the pressure state if this block has one
    
    block.P_Difference = {'HotPin','HotPout'; 'ColdPin', 'ColdPout';};
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
    
%     %calculate effectiveness
%     ColdMax = ColdIn;
%     ColdMax.T = block.Hot_T_init;
%     HotMin = HotIn;
%     HotMin.T = block.Cold_T_init;
%     QT = enthalpy(block.ColdOut.IC) - enthalpy(ColdIn);
%     maxQT1 = enthalpy(ColdMax) - enthalpy(ColdIn);
%     maxQT2 = enthalpy(HotIn) - enthalpy(HotMin);
%     Effectiveness = QT/min(maxQT1,maxQT2);
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    Target = lookupVal(block.Target);
    block.ColdPout.IC = Inlet.ColdPout;
    block.HotPout.IC  = Inlet.HotPout;
    block.PfactorCold = NetFlow(Inlet.ColdIn)/block.PdropCold;
    block.PfactorHot = NetFlow(Inlet.HotIn)/block.PdropHot;
    block.ColdPin.IC  = Inlet.ColdPout+block.PdropCold;
    block.HotPin.IC  = Inlet.HotPout+block.PdropHot;
    
    specCold = fieldnames(Inlet.ColdIn);
    specHot = fieldnames(Inlet.HotIn);
    %% Cold flow
    r = length(block.ColdFlowDir(:,1));
    for i = 1:1:length(specCold)
        if ~strcmp(specCold{i},'T')
            ColdOutlet.(specCold{i})(1:block.nodes,1) = Inlet.ColdIn.(specCold{i})/r;
        end
    end
    %% Hot flow
    r = length(block.HotFlowDir(:,1));
    for i = 1:1:length(specHot)
        if ~strcmp(specHot{i},'T')
            HotOutlet.(specHot{i})(1:block.nodes,1) = Inlet.HotIn.(specHot{i})/r;
        end
    end
    AreaOld = block.Area;
    %% adjust effective area to achieve desired effectiveness or temperature
    [block.Area,block.Effectiveness,method] = findArea(Inlet.ColdIn,Inlet.HotIn,block.h_conv,Target,block.sizemethod);

    block.HTconv = block.HTconv*block.Area/AreaOld;
    Y = [block.Scale(1:3*block.nodes);1];
    [T, Y] = ode15s(@(t,y) SolveTempsDynamic(t,y,block,ColdOutlet,HotOutlet,Inlet,method), [0, 1e5], Y);
    Y = Y(end,:)';
    block.Area = Y(end)*block.Area;
    block.HTconv = Y(end)*block.HTconv;

    block.Scale(1:3*block.nodes) = Y(1:end-1);
    
    block.Convection = block.h_conv*block.Area/block.nodes/1000; %all heat transfer coefficients converted to kW/K: thus Q = C*(T1-T2) is in kW
    block.ColdOut.IC = Inlet.ColdIn;
    block.ColdOut.IC.T  = mean(block.Scale(block.ColdFlowDir(:,end)));
    block.HotOut.IC = Inlet.HotIn;
    block.HotOut.IC.T  = mean(block.Scale(2*block.nodes+block.HotFlowDir(:,end)));
end

function block = HXmask(block)
% Script which orients the nodes relative to the flow directions
% direction = 1: co-flow, direction = 2: counter-flow, direction = 3: cross-flow
nodes = block.nodes;
columns = block.columns;
rows = block.rows;

for j = 1:1:columns
    if block.direction == 1
        block.ColdFlowDir(:,j) = (j:columns:nodes)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
    elseif block.direction ==2
        block.ColdFlowDir(:,j) = (columns-j+1:columns:nodes)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
    elseif block.direction ==3
        block.ColdFlowDir(:,j) = (1+columns*(j-1):j*columns)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
    end
    block.HotFlowDir(:,j) = (j:columns:nodes)';
end
block.HTadjacent = zeros(nodes,4);
for i = 1:1:nodes
    block.HTadjacent(i,1) = i-1;%previous node
    block.HTadjacent(i,2) = i+1;%next node
    block.HTadjacent(i,3) = i-columns;%node to left
    block.HTadjacent(i,4) = i+columns;%node to right
end
block.HTadjacent(1:columns:end,1)=linspace(1,nodes-columns+1,rows);%first node in each row has nothing before it
block.HTadjacent(columns:columns:end,2)=linspace(columns,nodes,rows)';%last node in each row has nothing after it
block.HTadjacent(1:columns,3)=linspace(1,columns,columns);%first row has nothing to left
block.HTadjacent(end-columns+1:end,4)=linspace(nodes-columns+1,nodes,columns);%last row has nothing to right

function [T , HTconv, HTcond] = SolveTemps(block,Cold,Hot)
%Solve problem of form xdot = Ax-b for xdot =0.
%final solution is x = A\b;
%states represent  heat transfer into each node/layer, the temperatures of each node/layers, and inlet cathode and anode temperatures, Qerror term associated with the small error in air flow rate so that the deltaT and Tavg constraints can both be satisfied

nodes = block.nodes;
states = 6*nodes+2;
%% Convection
hA = block.Convection;

A = zeros(states,states);
b = zeros(states,1);
%Tcold in 
A(6*nodes+1,6*nodes+1) = 1; %cold T
b(6*nodes+1) = Cold.T;
%Thot in
A(6*nodes+2,6*nodes+2) = 1; %hot T
b(6*nodes+2) = Hot.T;

Cold.T = (Cold.T+Hot.T)/2;
Hot.T = Cold.T;

C_cold = SpecHeat(Cold)*NetFlow(Cold);
C_hot = SpecHeat(Hot)*NetFlow(Hot);
for k = 1:1:nodes
    %QT1 : heat transfer into cold flow
    A(k,k) = -1;
    A(k,k+3*nodes) = -.5*hA;
    A(k,k+4*nodes) = hA;

    [i,j] = find(block.ColdFlowDir==k);
    if j==1 %first column averaged with inlet temperature
        A(k,6*nodes+1) = -.5*hA;
        A(k+nodes,6*nodes+1) = .5*hA;
    else % other columns averaged with previous one
        k2 = block.ColdFlowDir(i,j-1);
        A(k,k2+3*nodes) = -.5*hA;
        A(k+nodes,k2+3*nodes) = .5*hA;
    end

    %QT2 : heat transfer into plate
    A(k+nodes,k+nodes) = -1;
    A(k+nodes,k+4*nodes) = -2*hA;
    A(k+nodes,k+3*nodes) = .5*hA;
    A(k+nodes,k+5*nodes) = .5*hA;

    %QT3 : heat transfer into hot flow
    A(k+2*nodes,k+2*nodes) = -1;
    A(k+2*nodes,k+5*nodes) = -.5*hA;
    A(k+2*nodes,k+4*nodes) = hA;

    [i,j] = find(block.HotFlowDir==k);
    if j==1 %first column averaged with inlet temperature
        A(k+nodes,6*nodes+2) = .5*hA;
        A(k+2*nodes,6*nodes+2) = -.5*hA;
    else % other columns averaged with previous one
        k2 = block.HotFlowDir(i,j-1);
        A(k+nodes,k2+5*nodes) = .5*hA;
        A(k+2*nodes,k2+5*nodes) = -.5*hA;
    end

    %Tcold: Temperature of cold flow
    A(k+3*nodes,k) = 1;
    A(k+3*nodes,k+3*nodes) = -C_cold;
    [i,j] = find(block.ColdFlowDir==k);
    if j==1 %first column receives fresh air
        A(k+3*nodes,6*nodes+1) = C_cold;
    else
        index = block.ColdFlowDir(i,j-1);
        A(k+3*nodes,index +3*nodes) = C_cold;
    end

    %Tplate: Temperature of plate
    A(k+4*nodes,k+nodes) = 1;
    b(k+4*nodes) = 0;

    %Thot: Temperature of hot flow
    A(k+5*nodes,k+2*nodes) = 1;
    A(k+5*nodes,k+5*nodes) = -C_hot; 
    [i,j] = find(block.HotFlowDir==k);
    if j==1 %first column receives fresh flow
        A(k+5*nodes,6*nodes+2) = C_hot; %fresh inlet
    else
        index = block.HotFlowDir(i,j-1);
        A(k+5*nodes,index +5*nodes) = C_hot;
    end
end
%remove inlet temperature averaging on first column
k1 = block.ColdFlowDir(:,1); %first column
for n =1:1:length(k1)
    k = k1(n);
    A(k,k+3*nodes) = A(k,k+3*nodes) -.5*hA; %HT to cold flow from plate
    A(k,6*nodes+1) = A(k,6*nodes+1) +.5*hA; %ignore cold flow inlet

    A(k+nodes,k+3*nodes) = A(k+nodes,k+3*nodes) +.5*hA;%HT to plate from cold flow
    A(k+nodes,6*nodes+1) = A(k+nodes,6*nodes+1) -.5*hA;%ignore cold flow inlet
end
k1 = block.HotFlowDir(:,1); %first column
for n =1:1:length(k1)
    k = k1(n);
    A(k+2*nodes,k+5*nodes) = A(k+2*nodes,k+5*nodes) -.5*hA; %HT to hot flow from plate
    A(k+2*nodes,6*nodes+2) = A(k+2*nodes,6*nodes+2) +.5*hA; %ignore hot flow inlet

    A(k+nodes,k+5*nodes) = A(k+nodes,k+5*nodes) +.5*hA;%HT to plate from hot flow
    A(k+nodes,6*nodes+2) = A(k+nodes,6*nodes+2) -.5*hA;%ignore hot flow inlet
end
HTconv = A(1:3*nodes,3*nodes+1:6*nodes);%matrix of coefficients to multiply by vector of temperature and get the heat transfer by convection between layers and nodes

%% Conduction: left and right, prev and next
A2 = zeros(states,states);
prev = block.HTadjacent(:,1);
next = block.HTadjacent(:,2);
left = block.HTadjacent(:,3);
right = block.HTadjacent(:,4);
for k = 1:1:nodes
    A2(k+nodes,k+4*nodes) = A2(k+nodes,k+4*nodes) -2*block.ConductionPN -2*block.ConductionLR;
    A2(k+nodes,prev(k)+4*nodes) = A2(k+nodes,prev(k)+4*nodes)+block.ConductionPN;
    A2(k+nodes,next(k)+4*nodes) = A2(k+nodes,next(k)+4*nodes)+block.ConductionPN;
    A2(k+nodes,left(k)+4*nodes) = A2(k+nodes,left(k)+4*nodes)+block.ConductionLR;
    A2(k+nodes,right(k)+4*nodes) = A2(k+nodes,right(k)+4*nodes)+block.ConductionLR; 
end
HTcond = A2(1:3*nodes,3*nodes+1:6*nodes);%matrix of coefficients to multiply by vector of temperature and get the heat transfer by conduction between layers and nodes
A = A + A2;
x= A\b;
T = x(3*nodes+1:6*nodes);

function dY = SolveTempsDynamic(t,Y,block,Cold,Hot,Inlet,method)
nodes = block.nodes;
block.HTconv = block.HTconv*Y(end);
dY = 0*Y;
QT = block.HTconv*Y(1:3*nodes) + block.HTcond*Y(1:3*nodes);
Cold.T = Y(1:nodes);
ColdInlet = Cold;
ColdInlet.T(block.ColdFlowDir(:,1),1) = Inlet.ColdIn.T;
for j = 1:1:length(block.ColdFlowDir(1,:));%1:columns
    k = block.ColdFlowDir(:,j);
    if j~=1
        ColdInlet.T(k,1) = Cold.T(kprev);
    end
    kprev = k;
end

Hot.T = Y(2*nodes+1:3*nodes);
HotInlet = Hot;
HotInlet.T(block.HotFlowDir(:,1),1) = Inlet.HotIn.T;
for j = 1:1:length(block.HotFlowDir(1,:));%1:columns
    k = block.HotFlowDir(:,j);
    if j~=1
        HotInlet.T(k,1) = Hot.T(kprev);
    end
    kprev = k;
end

%energy flows & sepcific heats
[~,HoutCold] = enthalpy(Cold);
[~,HinCold] = enthalpy(ColdInlet);
[~,HoutHot] = enthalpy(Hot);
[~,HinHot] = enthalpy(HotInlet);
tC = (block.Mass*block.Solid_SpecHeat); %have everything scale with the slower time constant of the plate, until it is initialized

dY(1:nodes)= (QT(1:nodes) + HinCold - HoutCold)./tC; %Cold flow
dY(nodes+1:2*nodes)= QT(nodes+1:2*nodes)/tC;  % Plate
dY(2*nodes+1:3*nodes)= (QT(2*nodes+1:3*nodes) + HinHot - HoutHot)./tC; %Hot flow
if strcmp(method,'ColdT')
    error = (Target - mean(Y(block.ColdFlowDir(:,end),1)))/100;
elseif strcmp(method,'HotT')
    error = (Target - mean(Y(2*block.nodes+block.HotFlowDir(:,end),1)));
elseif strcmp(method,'Effectiveness')
    ColdOut = Inlet.ColdIn;
    ColdOut.T =  mean(Y(block.ColdFlowDir(:,end),1));
    QT = enthalpy(ColdOut) - enthalpy(Inlet.ColdIn);
    ColdMax = Inlet.ColdIn;
    ColdMax.T = Inlet.HotIn.T;
    HotMin = Inlet.HotIn;
    HotMin.T = Inlet.ColdIn.T;
    maxQT1 = enthalpy(ColdMax) - enthalpy(Inlet.ColdIn);
    maxQT2 = enthalpy(Inlet.HotIn) - enthalpy(HotMin);

    error = block.Effectiveness - QT/min(maxQT1,maxQT2);
elseif strcmp(method,'fixed')
    error = 0;
end
dY(end) = error; %slow change in area

function [Area,Effectiveness,method] = findArea(ColdIn,HotIn,h_conv,Target,method)
%%Ideal heat transfer
QinC = enthalpy(ColdIn);
ColdOut = ColdIn;
ColdOut.T = HotIn.T;
QinH = enthalpy(HotIn);
HotOut = HotIn;
HotOut.T = ColdIn.T;
Qcoldhot = enthalpy(ColdOut) - QinC;
Qhotcold = QinH - enthalpy(HotOut);
if strcmp(method,'ColdT')
    ColdOut.T = Target;
    QoutC = enthalpy(ColdOut);
    Effectiveness = (QoutC - QinC)/min(Qcoldhot,Qhotcold);
    if Effectiveness>.98
        Target=.98;
        method = 'Effectiveness';
    else Qnet = QoutC - QinC;
    end
end
if strcmp(method,'HotT')
    HotOut.T = Target;
    QoutH = enthalpy(HotOut);
    Effectiveness = (QinH - QoutH)/min(Qcoldhot,Qhotcold);
    if Effectiveness>.98
        Target=.98;
        method = 'Effectiveness';
    else Qnet = QinH - QoutH;
    end
end
if strcmp(method,'Effectiveness')
    Effectiveness = Target;
    Qnet = Effectiveness*min(Qcoldhot,Qhotcold);
end
if ~strcmp(method,'ColdT')
    QoutC = QinC+Qnet;
    ColdOut.T = ColdIn.T + Effectiveness*(HotIn.T - ColdOut.T);
    C = NetFlow(ColdOut)*SpecHeat(ColdOut); % Cp*flow
    errorT = 1;
    while abs(errorT)>1e-2
        errorT = (QoutC - enthalpy(ColdOut))/C;
        ColdOut.T = ColdOut.T + errorT;
    end
end
if ~strcmp(method,'HotT')
    QoutH = QinH-Qnet;
    HotOut.T = HotIn.T - Effectiveness*(HotIn.T - ColdOut.T);
    C = NetFlow(HotOut)*SpecHeat(HotOut); % Cp*flow
    errorT = 1;
    while abs(errorT)>1e-2
        errorT = (QoutH - enthalpy(HotOut))/C;
        HotOut.T = HotOut.T + errorT;
    end
end
%find the log mean temperature difference to find Qnet then find the effective area for this deltaT
LMTD = ((HotOut.T - ColdIn.T) - (HotIn.T - ColdOut.T))/(log(HotOut.T - ColdIn.T) - log(HotIn.T - ColdOut.T));
Area = 2*1000*Qnet/(LMTD*h_conv); %effective area in m^2

function A = lookupVal(initCond)
global modelParam Tags
if ischar(initCond)
    r = strfind(initCond,'.');
    if ~isempty(r)
        if strcmp(initCond(1:r(1)-1),'Tags')
            A = Tags.(initCond(r(1)+1:r(2)-1)).(initCond(r(2)+1:end));
        else
            r = [r,length(initCond)+1];
            A = modelParam.(initCond(1:r(1)-1));
            for i = 2:1:length(r)
                field = initCond(r(i-1)+1:r(i)-1);
                A = A.(field);
            end
        end
    else
        A = modelParam.(initCond);
    end
else A  = initCond;
end