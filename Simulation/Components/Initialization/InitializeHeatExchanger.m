function block = InitializeHeatExchanger(varargin)
%Nodal Heat Exchanger model with 3 states per node: Temperature hot, Plate temperature, temperature cold
% Four (4) inlets: {'Cold Flow','Hot Flow','Cold Pout','Hot Pout'}
% Four (4) outlets: {'Cold Flow','Hot Flow','Cold Pin','Hot Pin'}
% Many (3*n+2) states: {'Cold Flow', 'Plate', 'Hot Flow','Cold Pin','Hot Pin'}
block = varargin{1};
if length(varargin)==1 % first initialization
    %% Load mask parameters (flow direction)
    block.nodes = block.rows*block.columns;
    block = FlowDir(block,2);
    
    if ischar(block.ColdSpecIn)
        block.ColdSpecIn = ComponentProperty(block.ColdSpecIn);
    end
    block.spec1 = fieldnames(block.ColdSpecIn);
    block.spec1 = block.spec1(~strcmp(block.spec1,'T'));
    for i = 1:1:length(block.spec1)
        Inlet.Flow1.(block.spec1{i}) = block.ColdSpecIn.(block.spec1{i});
    end
    Inlet.Flow1.T = block.Cold_T_init;
    
    if ischar(block.HotSpecIn)
        block.HotSpecIn = ComponentProperty(block.HotSpecIn);
    end
    block.spec2 = fieldnames(block.HotSpecIn);
    block.spec2 = block.spec2(~strcmp(block.spec2,'T'));
    for i = 1:1:length(block.spec2)
        Inlet.Flow2.(block.spec2{i}) = block.HotSpecIn.(block.spec2{i});
    end
    Inlet.Flow2.T = block.Hot_T_init;
    
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
    
    Target = ComponentProperty(block.Target);
    [block.Area,~,~] = findArea(Inlet.Flow1,Inlet.Flow2,block.h_conv,Target,block.sizemethod);
    block.Convection = block.h_conv*block.Area/block.nodes/1000; %all heat transfer coefficients converted to kW/K: thus Q = C*(T1-T2) is in kW
    %horizontal
    block.L_node = block.Length/block.columns;
    block.W_node =block.Width/block.rows;
    block.AcondLR = block.Vol_Solid/block.Width/block.columns;
    block.AcondPN = block.Vol_Solid/block.Length/block.rows;
    block.ConductionPN = block.Solid_CondCoef*block.AcondPN/(block.L_node/2)/1000; %heat transfer coefficient between previous and next node of plate
    block.ConductionLR  = block.Solid_CondCoef*block.AcondLR/(block.W_node/2)/1000; %heat transfer coefficient between left and right adjacent nodes of oxidant plate

    [block.Scale , block.HTcond, block.HTconv] = SteadyTemps(block,Inlet.Flow1,Inlet.Flow2);
    block.Scale(end+1:end+2,1) = [101+block.PdropCold;101+block.PdropHot;];
    block.IC = ones(3*block.nodes+2,1);
    
    ColdOut = Inlet.Flow1;
    ColdOut.T =  mean(block.Scale(block.Flow1Dir(:,end),1));
    HotOut = Inlet.Flow2;
    HotOut.T = mean(block.Scale(2*block.nodes+block.Flow2Dir(:,end),1));
    block.Effectiveness = FindEffectiveness(Inlet.Flow1,Inlet.Flow2,ColdOut,[]);%calculate effectiveness

    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    block.InletPorts = {'Flow1','Flow2','ColdPout','HotPout'};
    block.Flow1.IC = Inlet.Flow1;
    block.Flow2.IC = Inlet.Flow2;
    block.ColdPout.IC = 101; %Atmospheric pressure
    block.ColdPout.Pstate = []; %identifies the state # of the pressure state if this block has one
    block.HotPout.IC = 101; %Atmospheric pressure
    block.HotPout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.OutletPorts = {'ColdOut','HotOut','ColdPin','HotPin'};
    block.ColdOut.IC = ColdOut;
    block.HotOut.IC = HotOut;
    block.ColdPin.IC  = block.ColdPout.IC+block.PdropCold;
    block.ColdPin.Pstate = 3*block.nodes+1; %identifies the state # of the pressure state if this block has one
    block.HotPin.IC  = block.HotPout.IC+block.PdropHot;
    block.HotPin.Pstate = 3*block.nodes+2; %identifies the state # of the pressure state if this block has one
    
    block.P_Difference = {'HotPin','HotPout'; 'ColdPin', 'ColdPout';};
    %no dMdP or mFlow (fixed pressure drop)
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    Target = ComponentProperty(block.Target);
    nodes = block.nodes;
    block.ColdPout.IC = Inlet.ColdPout;
    block.HotPout.IC  = Inlet.HotPout;
    block.PfactorCold = NetFlow(Inlet.Flow1)/block.PdropCold;
    block.PfactorHot = NetFlow(Inlet.Flow2)/block.PdropHot;
    block.ColdPin.IC  = Inlet.ColdPout+block.PdropCold;
    block.HotPin.IC  = Inlet.HotPout+block.PdropHot;
    
    specCold = fieldnames(Inlet.Flow1);
    specHot = fieldnames(Inlet.Flow2);
    %% Cold flow
    r = length(block.Flow1Dir(:,1));
    for i = 1:1:length(specCold)
        if ~strcmp(specCold{i},'T')
            ColdOutlet.(specCold{i})(1:block.nodes,1) = Inlet.Flow1.(specCold{i})/r;
        end
    end
    %% Hot flow
    r = length(block.Flow2Dir(:,1));
    for i = 1:1:length(specHot)
        if ~strcmp(specHot{i},'T')
            HotOutlet.(specHot{i})(1:block.nodes,1) = Inlet.Flow2.(specHot{i})/r;
        end
    end
    AreaOld = block.Area;
    %% adjust effective area to achieve desired effectiveness or temperature
    [block.Area,block.Effectiveness,method] = findArea(Inlet.Flow1,Inlet.Flow2,block.h_conv,Target,block.sizemethod);
    
    %% rescale guessed temperatures to match current inlets and this effectiveness
    ColdOut = Inlet.Flow1;
    ColdOut.T  = Inlet.Flow2.T;
    HotOut = Inlet.Flow2;
    HotOut.T = Inlet.Flow1.T;
    Q1 = enthalpy(ColdOut) - enthalpy(Inlet.Flow1);
    Q2 = enthalpy(Inlet.Flow2) - enthalpy(HotOut);
    if Q1>Q2
        HotOut.T = Inlet.Flow2.T - block.Effectiveness*(Inlet.Flow2.T - HotOut.T);
        Q2 = enthalpy(Inlet.Flow2) - enthalpy(HotOut);
        ColdEffective = Q2/Q1;
        ColdOut.T = Inlet.Flow1.T + ColdEffective*(ColdOut.T - Inlet.Flow1.T);
    else
        ColdOut.T = Inlet.Flow1.T + block.Effectiveness*(ColdOut.T - Inlet.Flow1.T);
        Q1 = enthalpy(ColdOut) - enthalpy(Inlet.Flow1);
        HotEffective = Q1/Q2;
        HotOut.T = Inlet.Flow2.T - HotEffective*(Inlet.Flow2.T - ColdOut.T);
    end
    cDist = (block.Scale(1:nodes) - block.Scale(block.Flow1Dir(1,1)))/(block.Scale(block.Flow1Dir(1,end)) - block.Flow1.IC.T);
    hDist = (block.Scale(2*nodes+1:3*nodes) - block.Scale(2*nodes+block.Flow2Dir(1,1)))/(block.Scale(2*nodes+block.Flow2Dir(1,end)) - block.Flow2.IC.T);
    block.Scale(1:nodes) = Inlet.Flow1.T + cDist*(ColdOut.T - Inlet.Flow1.T);
    block.Scale(2*nodes+1:3*nodes) = Inlet.Flow2.T + hDist*(HotOut.T - Inlet.Flow2.T);
    block.Scale(nodes+1:2*nodes) = (block.Scale(1:nodes) + block.Scale(2*nodes+1:3*nodes))/2;% Average hot and cold side
    %%%
    
    block.HTconv = block.HTconv*block.Area/AreaOld;
    Y = [block.Scale(1:3*nodes);1];
    [T, Y] = ode15s(@(t,y) SolveTempsDynamic(t,y,block,ColdOutlet,HotOutlet,Inlet,method,Target), [0, 1e5], Y);
    Y = Y(end,:)';
    block.Area = Y(end)*block.Area;
    block.HTconv = Y(end)*block.HTconv;

    block.Scale(1:3*block.nodes) = Y(1:end-1);
    
    block.Convection = block.h_conv*block.Area/block.nodes/1000; %all heat transfer coefficients converted to kW/K: thus Q = C*(T1-T2) is in kW
    block.Flow1.IC = Inlet.Flow1;
    block.ColdOut.IC = Inlet.Flow1;
    block.ColdOut.IC.T  = mean(block.Scale(block.Flow1Dir(:,end)));
    block.Flow2.IC = Inlet.Flow2;
    block.HotOut.IC = Inlet.Flow2;
    block.HotOut.IC.T  = mean(block.Scale(2*block.nodes+block.Flow2Dir(:,end)));
    
    [block.Effectiveness,Imbalance] = FindEffectiveness(Inlet.Flow1,Inlet.Flow2,block.ColdOut.IC,block.HotOut.IC);%calculate effectiveness
end


function dY = SolveTempsDynamic(t,Y,block,Flow1,Flow2,Inlet,method,Target)
nodes = block.nodes;
block.HTconv = block.HTconv*Y(end);
dY = 0*Y;
QT = block.HTconv*Y(1:3*nodes) + block.HTcond*Y(1:3*nodes);
Flow1.T = Y(1:nodes);
ColdInlet = Flow1;
ColdInlet.T(block.Flow1Dir(:,1),1) = Inlet.Flow1.T;
for j = 1:1:length(block.Flow1Dir(1,:));%1:columns
    k = block.Flow1Dir(:,j);
    if j~=1
        ColdInlet.T(k,1) = Flow1.T(kprev);
    end
    kprev = k;
end

Flow2.T = Y(2*nodes+1:3*nodes);
HotInlet = Flow2;
HotInlet.T(block.Flow2Dir(:,1),1) = Inlet.Flow2.T;
for j = 1:1:length(block.Flow2Dir(1,:));%1:columns
    k = block.Flow2Dir(:,j);
    if j~=1
        HotInlet.T(k,1) = Flow2.T(kprev);
    end
    kprev = k;
end

%energy flows & sepcific heats
Hout1 = enthalpy(Flow1);
Hin1 = enthalpy(ColdInlet);
Hout2 = enthalpy(Flow2);
Hin2 = enthalpy(HotInlet);
tC = (block.Mass*block.Solid_SpecHeat); %have everything scale with the slower time constant of the plate, until it is initialized

dY(1:nodes)= (QT(1:nodes) + Hin1 - Hout1)./tC; %Cold flow
dY(nodes+1:2*nodes)= QT(nodes+1:2*nodes)/tC;  % Plate
dY(2*nodes+1:3*nodes)= (QT(2*nodes+1:3*nodes) + Hin2 - Hout2)./tC; %Hot flow
if strcmp(method,'ColdT')
    error = (Target - mean(Y(block.Flow1Dir(:,end),1)))/1e5;
elseif strcmp(method,'HotT')
    error = (Target - mean(Y(2*block.nodes+block.Flow2Dir(:,end),1)))/1e5;
elseif strcmp(method,'Effectiveness')
    ColdOut = Inlet.Flow1;
    ColdOut.T =  mean(Y(block.Flow1Dir(:,end),1));
    QT = enthalpy(ColdOut) - enthalpy(Inlet.Flow1);
    ColdMax = Inlet.Flow1;
    ColdMax.T = Inlet.Flow2.T;
    HotMin = Inlet.Flow2;
    HotMin.T = Inlet.Flow1.T;
    maxQT1 = enthalpy(ColdMax) - enthalpy(Inlet.Flow1);
    maxQT2 = enthalpy(Inlet.Flow2) - enthalpy(HotMin);

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