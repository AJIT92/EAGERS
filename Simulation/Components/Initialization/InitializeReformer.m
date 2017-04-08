function block = InitializeReformer(varargin)
%Nodal reformer model with or without heat exchange from second fluid
% Two or Four (2 or 4) inlets: {'Primary Flow', 'Primary Pout', 'Secondary Flow','Secondary Pout'}
global Ru
Ru = 8.314472; % Universal gas constant in kJ/K*kmol
block = varargin{1};
if length(varargin)==1 % first initialization
    %% Load flow direction
    if length(block.connections) ==2
        block.hotStream = 0;
        block = FlowDir(block,1);
    else
        block.hotStream = 1;
        block = FlowDir(block,2);
    end

    block.Pdrop1 = 1;                                    % (kPa) pressure drop
    block.h_conv= 50;                                       %Convective heat coefficient (W/m^2*K)
    block.H_Channel = .01;                                 %height of channels
    block.W_Channel = .01;                                  %Width of Channels
    block.t_solid = .005;                                   %Thickness of reformer solid (m)
    block.Length = 1;                                       %Length of reformer (m)
    block.Width = .75;                                      %Width of reformer (m)
    block.Height = .5;                                      %Height of reformer (m)
    block.Solid_Density=8055;                               %Density of solid material (kg/m^3)
    block.Solid_SpecHeat=40;%480;                               %Specific Heat of solid material (J/kg*K)
    block.Solid_CondCoef = 15;                              %Conduction coefficient (Steel) (W/m*K)
    block.Channels = (block.Width/(block.W_Channel+block.t_solid))*(block.Height/(block.H_Channel+block.t_solid)); %# of Channels assuming an even grid of holes in a LxW block
    block.Vol_1 = block.Length*block.Channels*block.H_Channel*block.W_Channel; %gaseous volume
    block.Vol_Solid = block.Length*block.Width*block.Height - block.Vol_1;
    if block.hotStream ==1 % reformer with heat exchange assume alternating rows of channels with reforming flow and heating flow
        block.Vol_1 = block.Vol_1/2;
        block.Vol_2 = block.Vol_1;
    end

    %all heat transfer coefficients converted to kW/K: thus Q = C*(T1-T2) is in kW
    block.Area = block.Length*block.Channels*(2*block.H_Channel+2*block.W_Channel)/block.nodes; %surface area of each node 
    block.Convection = block.h_conv*block.Area/1000; % h *A, W/m^2*K  x m^2  / 1000 = kW/K
    block.SolidArea = block.Width*block.Height - block.Channels*block.W_Channel*block.H_Channel;
    block.ConductionPN = 0;%block.Solid_CondCoef*block.SolidArea/(block.Length/block.nodes/2)/1000; %heat transfer coefficient between previous and next node
    block.ConductionLR = block.ConductionPN; %left right conduction (only doing a single row)
    
    %for initial flow guess assume 50K between solid and primary flow
    Atotal = block.Length*block.Channels*(2*block.H_Channel+2*block.W_Channel); %surface area of each node 
    Qtrans = 10*block.h_conv*Atotal/1000; %maximum heat transfer in kW
    [h,~] = enthalpy(723,{'H2','H2O','O2','CO','CO2','CH4'});
    h_rxn1 = h.CO+3*h.H2-h.CH4-h.H2O;
    h_rxn2 = h.CO2+h.H2-h.CO-h.H2O;
    R.CH4 = block.InletGuess.CH4;
    R.WGS = R.CH4  + 0.8*block.InletGuess.CO;
    Qref = h_rxn1.*R.CH4+h_rxn2.*R.WGS; %reforming energy per mole inlet
    ScaleFlow = Qtrans/Qref;
    
    block.spec1 = fieldnames(block.InletGuess);
    for i = 1:1:length(block.spec1)
        Inlet.Flow1.(block.spec1{i}) = block.InletGuess.(block.spec1{i})*ScaleFlow;
    end
    criticalSpecies = {'CH4';'CO';'CO2';'H2';'H2O'};
    for i = 1:1:length(criticalSpecies)
        if ~ismember(block.spec1,criticalSpecies{i})
            Inlet.Flow1.(criticalSpecies{i}) = 0;
        end
    end
    block.S2C = Inlet.Flow1.H2O/(Inlet.Flow1.CH4 + .5*Inlet.Flow1.CO);
    block.scaleK_WGS = 1;
    block.scaleK_CH4 = 1;
    Inlet.Flow1.T = 723;
    Inlet.Flow1Pout = 101;
    block.ReformedPin.IC = 101+block.Pdrop1;
    
    if block.hotStream ==1
        block.T.Primary = linspace(700,973,block.nodes)';
        block.Pdrop2 = 1;                                    % (kPa) pressure drop
        block.Vol_2 = block.Vol_1;
        block.spec2 = {'N2'};
        Inlet.Flow2.N2 = 1;
        Inlet.Flow2.T = 1050;
        C = SpecHeat(Inlet.Flow2);
        ScaleFlow2 = Qtrans/(C*50); %assume 50K of temp drop on hot side
        Inlet.Flow2.N2 = ScaleFlow2;
        block.CooledPin.IC = 101+block.Pdrop2;
        Inlet.Flow2Pout = 101;
        
        
    else
        block.T.Primary = linspace(800,500,block.nodes)';
        Inlet.Flow2 = [];
        
        
    end
    %% Run Initial Condition
    [Primary,R] = equilib2D(Inlet.Flow1,block.T.Primary(block.Flow1Dir(end)),Inlet.Flow1Pout,0,0,'Reformer',1,[]);
%     refPerc = R.CH4/Inlet.Flow1.CH4
    
    R.CH4 = ones(length(block.Flow1Dir),1)*R.CH4/length(block.Flow1Dir);
    R.WGS = ones(length(block.Flow1Dir),1)*R.WGS/length(block.Flow1Dir);
    
    [h,~] = enthalpy(block.T.Primary,{'H2','H2O','O2','CO','CO2','CH4'});
    h_rxn1 = h.CO+3*h.H2-h.CH4-h.H2O;
    h_rxn2 = h.CO2+h.H2-h.CO-h.H2O;
    block.Qref = h_rxn1.*R.CH4+h_rxn2.*R.WGS; %heat transfer to reforming reactions;
%     
%     block.Qref = zeros(length(block.Flow1Dir),1);
    
    [block.Tstates,block.HTcond,block.HTconv] = SteadyTemps(block,Inlet.Flow1,Inlet.Flow2);
    block = Set_IC(block,Primary,Inlet.Flow2);

    [Primary, Secondary,block] = solveInitCond(Inlet,block);
    
    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    block.InletPorts = {'Flow1','Flow1Pout'};
    block.Flow1.IC  = Inlet.Flow1;
    block.Flow1Pout.IC = 101; %Atmospheric pressure
    block.Flow1Pout.Pstate = []; %identifies the state # of the pressure state if this block has one

    block.OutletPorts = {'Reformed','ReformedPin','MeasureReformT','MeasureS2C'};
    block.Reformed.IC  = Primary;
    block.ReformedPin.IC  = block.ReformedPin.IC;
    block.ReformedPin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    block.MeasureReformT.IC = Primary.T;
    block.MeasureS2C.IC = block.S2C;
    block.PfactorPrimary = NetFlow(Primary)/block.Pdrop1;

    if block.hotStream==1
        block.InletPorts = {'Flow1','Flow1Pout','Flow2','Flow2Pout'};
        block.Flow2.IC  = Inlet.Flow2;
        block.Flow2Pout.IC = 101; %Atmospheric pressure
        block.Flow2Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
        
        block.OutletPorts = {'Reformed','ReformedPin','Cooled','CooledPin','MeasureReformT','MeasureS2C'};
        block.Cooled.IC  = Secondary;
        block.CooledPin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
        block.ReformedPin.Pstate = length(block.Scale)-1;
        block.CooledPin.IC  = block.CooledPin.IC;
        block.PfactorSecondary = NetFlow(Secondary)/block.Pdrop2;
        block.MeasureCooledT.IC = Secondary.T;
        
        block.P_Difference = {'ReformedPin','Flow1Pout'; 'CooledPin', 'Flow2Pout';};        
    end
elseif length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    N = 2;
    if block.hotStream==1
        N=3;
    end
    block.ReformedPin.IC = Inlet.Flow1Pout+block.Pdrop1;
    SpecNew = fieldnames(Inlet.Flow1);
    SpecNew = SpecNew(~strcmp('T',SpecNew));
    sp = length(block.spec1);
    n = (N+sp)*block.nodes;
    for i = 1:1:length(SpecNew)
        if ~ismember(SpecNew{i},block.spec1)
            block.Scale = [block.Scale(1:n); zeros(block.nodes,1); block.Scale(n+1:end);]; %add zero states for new species
            block.tC = [block.tC(1:n); block.tC(n-block.nodes+1:n); block.tC(n+1:end);]; %add states for new species
            n = n+block.nodes;
            block.spec1{end+1} = SpecNew{i};
        end
    end
    for i = 1:1:length(block.spec1)
        if ~ismember(block.spec1{i},SpecNew)
            Inlet.Flow1.(block.spec1{i})=0;
        end
    end
    block.S2C = Inlet.Flow1.H2O/(Inlet.Flow1.CH4 + .5*Inlet.Flow1.CO);
    if block.hotStream==1
        block.CooledPin.IC = Inlet.Flow2Pout+block.Pdrop2;
        SpecNew = fieldnames(Inlet.Flow2);
        SpecNew = SpecNew(~strcmp('T',SpecNew));
        sp = length(block.spec2);
        n = n + sp*block.nodes;
        for i = 1:1:length(SpecNew)
            if ~ismember(SpecNew{i},block.spec2)
                block.Scale = [block.Scale(1:n); zeros(block.nodes,1); block.Scale(n+1:end);]; %add zero states for new species
                block.tC = [block.tC(1:n); block.tC(n-block.nodes+1:n); block.tC(n+1:end);]; %add states for new species
                n = n+block.nodes;
                block.spec2{end+1} = SpecNew{i};
            end
        end
        for i = 1:1:length(block.spec2)
            if ~ismember(block.spec2{i},SpecNew)
                Inlet.Flow2.(block.spec2{i})=0;
            end
        end
    end
    
    %%-- %%
    [Primary, Secondary,block] = solveInitCond(Inlet,block);
    %%%
    
    block.ReformedPin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    
    block.Primary.IC  = Inlet.Flow1;
    block.PrimaryPout.IC = Inlet.Flow1Pout;
    block.Reformed.IC  = Primary;
    block.ReformedPin.IC  = block.ReformedPin.IC;
    block.MeasureReformT.IC = Primary.T;
    block.MeasureS2C.IC = block.S2C;
    block.PfactorPrimary = NetFlow(Primary)/block.Pdrop1;
    if block.hotStream==1
        block.CooledPin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
        block.ReformedPin.Pstate = length(block.Scale)-1;
        block.Secondary.IC  = Inlet.Flow2;
        block.SecondaryPout.IC = Inlet.Flow2Pout; 
        block.Cooled.IC  = Secondary;
        block.CooledPin.IC  = block.CooledPin.IC;
        block.MeasureCooledT.IC = Secondary.T;        
        block.PfactorSecondary = NetFlow(Secondary)/block.Pdrop2;
    end
end


function [Flow1, Flow2,block] = solveInitCond(Inlet,block)
nodes = block.nodes;
specInterest = {'CH4','CO','CO2','H2','H2O'};
spec = fieldnames(Inlet.Flow1);
spec = spec(~strcmp('T',spec));
for i = 1:1:length(specInterest)
    if ~ismember(specInterest{i},spec)
        Inlet.Flow1.(specInterest{i}) = 0;
    end
end

N = 2;
if block.hotStream ==1
    N=3;
end
n = N*nodes;
sp = length(block.spec1);
scaleFlow = NetFlow(Inlet.Flow1)/sum(block.Scale(n+1:nodes:n+sp*nodes));
for i = 1:1:nodes
    block.Scale(n+i:nodes:n+sp*nodes) = block.Scale(n+i:nodes:n+sp*nodes)*scaleFlow;
end
if block.hotStream ==1
    n = (N+sp)*nodes;
    sp = length(block.spec2);
    scaleFlow = NetFlow(Inlet.Flow2)/sum(block.Scale(n+1:nodes:n+sp*nodes));
    for i = 1:1:nodes
        block.Scale(n+i:nodes:n+sp*nodes) = block.Scale(n+i:nodes:n+sp*nodes)*scaleFlow;
    end
end

Y = [block.Scale;1];
[T, Yt] = ode15s(@(t,y) SolveDynamic(t,y,block,Inlet), [0, 1e2],Y);
Y = Yt(end,:)';
n = 0;
Primary.T = Y(n+1:n+nodes);n = n+2*nodes;
if block.hotStream ==1
    Secondary.T = Y(n+1:n+nodes);n = n+nodes;
else Secondary = [];
end
for i = 1:1:length(block.spec1)
    Primary.(block.spec1{i}) = max(0,Y(n+1:n+nodes));n = n+nodes;
end
if block.hotStream ==1
    for i = 1:1:length(block.spec2)
        Secondary.(block.spec2{i}) = max(0,Y(n+1:n+nodes));n = n+nodes;
    end
end
block.Area = Y(end)*block.Area;
block.HTconv = Y(end)*block.HTconv;
block.Tstates = Y(1:N*nodes);

block = Set_IC(block,Primary,Secondary);

Flow1 = MergeLastColumn(Primary,block.Flow1Dir,1);
if block.hotStream==1
    Flow2 = MergeLastColumn(Secondary,block.Flow2Dir,1);
else Flow2 = [];
end


function block = Set_IC(block,Primary,Secondary)
global Ru
nodes = block.nodes;

N = 2 ;
nS = nodes*(N + length(block.spec1)) + 1;
if block.hotStream ==1
    N = 3;
    nS = nodes*(3 + length(block.spec1) + length(block.spec2)) + 2;
end
block.IC = ones(nS,1); 
block.Scale = block.IC;
block.tC = ones(nS,1); % time constant for derivative dY

block.Scale(1:nodes,1) = block.Tstates(1:nodes);%temperature (K)
Cp = SpecHeat(Primary);
block.tC(1:block.nodes,1) = block.Vol_1*Cp*block.ReformedPin.IC./(Ru*block.Tstates(1:nodes));
n = N*nodes;
Flow = NetFlow(Primary);
for i = 1:1:length(block.spec1)
    if any(Primary.(block.spec1{i})==0) %need to prevent possible divide by zero
        block.IC(n+1:n+nodes,1) = Primary.(block.spec1{i})./Flow;%concentration
        block.Scale(n+1:n+nodes,1) = Flow; 
    else
        block.Scale(n+1:n+nodes,1) = Primary.(block.spec1{i}); 
    end
    block.tC(n+1:n+nodes,1) = (block.Vol_1*block.ReformedPin.IC)./(block.Tstates(1:nodes)*Ru);
    n = n + nodes;
end

block.Scale(1+nodes:2*nodes,1) = block.Tstates(1+nodes:2*nodes);%temperature (K)
block.tC(1+nodes:2*nodes,1) = block.Vol_Solid*block.Solid_Density*block.Solid_SpecHeat; %solid

if block.hotStream ==1
    block.Scale(1+2*nodes:3*nodes,1) = block.Tstates(1+2*nodes:3*nodes);%temperature (K)
    Cp = SpecHeat(Secondary);
    block.tC(1+2*nodes:3*nodes,1) = block.Vol_2*Cp*block.CooledPin.IC./(Ru*block.Tstates(1:nodes));
    Flow = NetFlow(Secondary);
    for i = 1:1:length(block.spec2)
        if any(Secondary.(block.spec2{i})==0) %need to prevent possible divide by zero
            block.IC(n+1:n+nodes,1) = Secondary.(block.spec2{i})./Flow;
            block.Scale(n+1:n+nodes,1) = Flow;
        else
            block.Scale(n+1:n+nodes,1) = Secondary.(block.spec2{i}); 
        end
        block.tC(n+1:n+nodes,1) = (block.Vol_2*block.CooledPin.IC)./(block.Tstates(2*nodes+1:3*nodes)*Ru);
        n = n + nodes; 
    end
end
block.Scale(n+1,1) = block.ReformedPin.IC;%pressure
block.tC(n+1,1) = block.Vol_1;  %pressure
if block.hotStream ==1
    block.Scale(n+2,1) = block.CooledPin.IC;%pressure
    block.tC(n+2,1) = block.Vol_2; %pressure
end

function dY = SolveDynamic(t,Y,block,Inlet)
nodes = block.nodes;
block.HTconv = Y(end)*block.HTconv;
dY = 0*Y;
n = 0;
Primary.Outlet.T = Y(n+1:n+nodes);n = n+2*nodes;
s = 0;
N = 2;
if block.hotStream==1
    SecondaryOut.T = Y(n+1:n+nodes);n = n+nodes;
    s = 1;
    N = 3;
end
%% Cold flow
for i = 1:1:length(block.spec1)
    Primary.Outlet.(block.spec1{i}) = max(0,Y(n+1:n+nodes));n = n+nodes;
end
for j = 1:1:length(block.Flow1Dir(1,:));%1:columns
    k = block.Flow1Dir(:,j);
    r = length(k);
    if j==1
        Primary.Inlet.T(k,1) = Inlet.Flow1.T;
        for i = 1:1:length(block.spec1)
            Primary.Inlet.(block.spec1{i})(k,1) = Inlet.Flow1.(block.spec1{i})/r;
        end
    else
        Primary.Inlet.T(k,1) = Primary.Outlet.T(kprev,1);
        for i = 1:1:length(block.spec1)
            Primary.Inlet.(block.spec1{i})(k,1) = Primary.Outlet.(block.spec1{i})(kprev,1);
        end
    end
    kprev = k;
end
if block.hotStream==1
    
end

a = 4352.2./Primary.Outlet.T - 3.99;
K_WGS = block.scaleK_WGS.*exp(a);% Water gas shift equilibrium constant
K_CH4 = block.scaleK_CH4.*2459000.*exp(-6.187*a);
CH4_eq = (Primary.Outlet.H2.^3.*Primary.Outlet.CO)./(K_CH4.*Primary.Outlet.H2O).*(block.ReformedPin.IC./NetFlow(Primary.Outlet)).^2;
R.CH4 = (Primary.Inlet.CH4)-CH4_eq; %inlet CO + CO from reforming - outlet CO
CO_eq = Primary.Outlet.CO2.*Primary.Outlet.H2./(K_WGS.*Primary.Outlet.H2O);
R.WGS = (Primary.Inlet.CO+R.CH4)-CO_eq; %inlet CO + CO from reforming - outlet CO

%confirm we don't violate anything casuing negative species
R.CH4 = min([R.CH4,Primary.Inlet.CH4,Primary.Inlet.H2O],[],2);
R.CH4 = max([R.CH4,-Primary.Inlet.CO,-Primary.Inlet.H2/3],[],2);
R.WGS = min([R.WGS,Primary.Inlet.CO+R.CH4,Primary.Inlet.H2O-R.CH4],[],2);
R.WGS = max([R.WGS,-Primary.Inlet.CO2,-(Primary.Inlet.H2+3*R.CH4)],[],2);
%identify the target outflow that dY will converge to
RefOut.T = Primary.Outlet.T;
for i = 1:1:length(block.spec1)
    if strcmp(block.spec1{i},'CO2') 
        RefOut.CO2 = (Primary.Inlet.CO2 + R.WGS); 
    elseif strcmp(block.spec1{i},'H2') 
        RefOut.H2 = (Primary.Inlet.H2 + 3*R.CH4 + R.WGS);
    elseif strcmp(block.spec1{i},'H2O') 
        RefOut.H2O = (Primary.Inlet.H2O - R.CH4 - R.WGS); 
    elseif strcmp(block.spec1{i},'CH4') 
        RefOut.CH4 = (Primary.Inlet.CH4 - R.CH4); 
    elseif strcmp(block.spec1{i},'CO') 
        RefOut.CO = (Primary.Inlet.CO + R.CH4 - R.WGS);  
    else
        RefOut.(block.spec1{i}) = Primary.Inlet.(block.spec1{i});  
    end
end 

HoutCold = enthalpy(Primary.Outlet);
scale = NetFlow(Primary.Outlet)./NetFlow(RefOut);
HinCold = enthalpy(Primary.Inlet).*scale;
QT = block.HTconv*Y(1:N*nodes) + block.HTcond*Y(1:N*nodes);

for i=1:1:length(block.Flow1Dir(1,:))
    k = block.Flow1Dir(:,i);
    dY(k)= (QT(k) + HinCold(k) - HoutCold(k))./block.tC(k); %Cold flow
    if i>1
        dY(k) = dY(k)+dY(kprev);
    end
    kprev = k;
end
dY(nodes+1:2*nodes)= QT(nodes+1:2*nodes)./block.tC(nodes+1:2*nodes);  % solid
for j = 1:1:length(block.spec1)
    dY((1+j+s)*nodes+1:(2+j+s)*nodes)= (RefOut.(block.spec1{j})- Primary.Outlet.(block.spec1{j})) ./block.tC((1+j+s)*nodes+1:(2+j+s)*nodes); %all species concentration
end 

if block.hotStream==1
    %% Hot flow
    for i = 1:1:length(block.spec2)
        SecondaryOut.(block.spec2{i}) = max(0,Y(n+1:n+nodes));n = n+nodes;
    end
    for j = 1:1:length(block.Flow2Dir(1,:));%1:columns
        k = block.Flow2Dir(:,j);
        r = length(k);
        if j==1
            SecondaryIn.T(k,1) = Inlet.Flow2.T;
            for i = 1:1:length(block.spec2)
                SecondaryIn.(block.spec2{i})(k,1) = Inlet.Flow2.(block.spec2{i})/r;
            end
        else
            SecondaryIn.T(k,1) = SecondaryOut.T(kprev);
            for i = 1:1:length(block.spec2)
                SecondaryIn.(block.spec2{i})(k,1) = SecondaryOut.(block.spec2{i})(kprev,1);
            end
        end
        kprev = k;
    end
    HoutHot = enthalpy(SecondaryOut);
    scale = NetFlow(SecondaryOut)./NetFlow(SecondaryIn);
    HinHot = enthalpy(SecondaryIn).*scale;
    for i=1:1:length(block.Flow2Dir(1,:))
        k = block.Flow2Dir(:,i);
        dY(2*nodes+k)= (QT(2*nodes+k) + HinHot(k) - HoutHot(k))./block.tC(2*nodes+k); %Hot flow
        if i>1
            dY(2*nodes+k) = dY(2*nodes+k)+dY(2*nodes+kprev);
        end
        kprev = k;
    end
    for j = 1:1:length(block.spec2)
        k = (N+length(block.spec1)+j-1)*nodes+1:(N+length(block.spec1)+j)*nodes;
        dY(k) = (SecondaryIn.(block.spec2{j}) - SecondaryOut.(block.spec2{j}))./block.tC(k);  %all  species concentration
    end
    if strcmp(block.method,'RefPerc')
        RefPerc = (Inlet.Flow1.CH4-Primary.Outlet.CH4(block.Flow1Dir(1,end)))/Inlet.Flow1.CH4;
        error = block.ReformTarget - RefPerc;
    elseif strcmp(block.method,'ColdT')
        error = (block.ReformTarget - Primary.T(block.Flow1Dir(1,end)))/100;
    elseif strcmp(block.method,'HotT')
        error = (Secondary.T(block.Flow2Dir(1,end)) - block.ReformTarget)/100;
    elseif strcmp(block.method,'Effectiveness')
        FlowIdeal = Inlet.Flow2;
        FlowIdeal.T = Primary.Outlet.T(block.Flow2Dir(1,end));
        maxQT = enthalpy(FlowIdeal) - enthalpy(Inlet.Flow2); %assume ideal is at same temperature as 1st node on primary side
        Flow2.T = Secondary.T(block.Flow2Dir(1,end));
        for i = 1:1:length(block.spec2)
            Flow2.(block.spec2{i}) = Secondary.(block.spec2{i})(block.Flow2Dir(1,end));
        end
        QT = enthalpy(Flow2) - enthalpy(Inlet.Flow2);
        error = block.ReformTarget - QT/maxQT;
    elseif strcmp(block.method,'none')
        error = 0;
    end
else
    error = 0; %no surface area adjustment for adiabatic reformer
end
dY(end) = error; %slow change in the area