function block = InitializeReformer(varargin)
%Nodal reformer model with or without heat exchange from second fluid
% Two or Four (2 or 4) inlets: {'Primary Flow', 'Primary Pout', 'Secondary Flow','Secondary Pout'}
global Ru
Ru = 8.314472; % Universal gas constant in kJ/K*kmol
block = varargin{1};
if length(varargin)==1 % first initialization
    %% Load mask parameters (flow direction)
    if length(block.connections) ==2
        block.hotStream = 0;
    else block.hotStream = 1;
    end
    block = FlowDir(block);

    %% Heat Exchanger
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
    %%-- 
    
    %all heat transfer coefficients converted to kW/K: thus Q = C*(T1-T2) is in kW
    block.Area = block.Length*block.Channels*(2*block.H_Channel+2*block.W_Channel)/block.nodes; %surface area of each node 
    block.Convection = block.h_conv*block.Area/1000; % h *A, W/m^2*K  x m^2  / 1000 = kW/K
    block.SolidArea = block.Width*block.Height - block.Channels*block.W_Channel*block.H_Channel;
    block.Conduction = 0;%block.Solid_CondCoef*block.SolidArea/(block.Length/block.nodes/2)/1000; %heat transfer coefficient between previous and next node
    
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
        Inlet.Primary.(block.spec1{i}) = block.InletGuess.(block.spec1{i})*ScaleFlow;
    end
    criticalSpecies = {'CH4';'CO';'CO2';'H2';'H2O'};
    for i = 1:1:length(criticalSpecies)
        if ~ismember(block.spec1,criticalSpecies{i})
            Inlet.Primary.(criticalSpecies{i}) = 0;
        end
    end
    block.S2C = Inlet.Primary.H2O/(Inlet.Primary.CH4 + .5*Inlet.Primary.CO);
    block.scaleK_WGS = 1;
    block.scaleK_CH4 = 1;
    Inlet.Primary.T = 723;
    Inlet.PrimaryPout = 101;
    block.ReformedPin.IC = 101+block.Pdrop1;
    block.T.Primary = linspace(700,973,block.nodes)';
    if block.hotStream ==1
        block.Pdrop2 = 1;                                    % (kPa) pressure drop
        block.Vol_2 = block.Vol_1;
        block.spec2 = {'N2'};
        Inlet.Secondary.N2 = 1;
        Inlet.Secondary.T = 1050;
        C = SpecHeat(Inlet.Secondary);
        ScaleFlow2 = Qtrans/(C*50); %assume 50K of temp drop on hot side
        Inlet.Secondary.N2 = ScaleFlow2;
        block.CooledPin.IC = 101+block.Pdrop2;
        Inlet.SecondaryPout = 101;
    end
    %% Run Initial Condition
    [Primary,R] = equilib2D(Inlet.Primary,block.T.Primary(block.PrimaryDir(end)),Inlet.PrimaryPout,0,'Reformer',1,[]);
%     refPerc = R.CH4/Inlet.Primary.CH4
    
    R.CH4 = ones(length(block.PrimaryDir),1)*R.CH4/length(block.PrimaryDir);
    R.WGS = ones(length(block.PrimaryDir),1)*R.WGS/length(block.PrimaryDir);
    
    [h,~] = enthalpy(block.T.Primary,{'H2','H2O','O2','CO','CO2','CH4'});
    h_rxn1 = h.CO+3*h.H2-h.CH4-h.H2O;
    h_rxn2 = h.CO2+h.H2-h.CO-h.H2O;
    Qref = h_rxn1.*R.CH4+h_rxn2.*R.WGS; %heat transfer to reforming reactions;
    [block.Tstates,block.HTconv,block.HTcond] = SteadyTemps(block,Qref,Inlet);
    block = Set_IC(block,Primary,Inlet.Secondary);

    [Primary, Secondary,block] = solveInitCond(Inlet,block);
    
    %% set up ports : Inlets need to either connected or have initial condition, outlets need an initial condition, and it doesn't matter if they have a connection 
    block.PortNames = {'Primary','PrimaryPout','Reformed','ReformedPin','MeasureReformT','MeasureS2C'};
    block.Primary.type = 'in';
    block.Primary.IC  = Inlet.Primary;
    
    block.PrimaryPout.type = 'in';
    block.PrimaryPout.IC = 101; %Atmospheric pressure
    block.PrimaryPout.Pstate = []; %identifies the state # of the pressure state if this block has one

    block.Reformed.type = 'out';
    block.Reformed.IC  = Primary;
    
    block.ReformedPin.type = 'out';
    block.ReformedPin.IC  = block.ReformedPin.IC;
    block.ReformedPin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    
    block.MeasureReformT.type = 'out';
    block.MeasureReformT.IC = Primary.T;
    
    block.MeasureS2C.type = 'out';
    block.MeasureS2C.IC = block.S2C;
    
    block.P_Difference = {'ReformedPin','PrimaryPout';};
    block.PfactorPrimary = NetFlow(Primary)/block.Pdrop1;

    if block.hotStream==1
        block.PortNames = {'Primary','PrimaryPout','Secondary','SecondaryPout','Reformed','ReformedPin','Cooled','CooledPin','MeasureReformT','MeasureCooledT','MeasureS2C'};
        block.Secondary.type = 'in';
        block.Secondary.IC  = Inlet.Secondary;
        
        block.SecondaryPout.type = 'in';
        block.SecondaryPout.IC = 101; %Atmospheric pressure
        block.SecondaryPout.Pstate = []; %identifies the state # of the pressure state if this block has one
        
        block.Cooled.type = 'out';
        block.Cooled.IC  = Secondary;
        
        block.CooledPin.type = 'out';
        block.CooledPin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
        block.ReformedPin.Pstate = length(block.Scale)-1;
        block.CooledPin.IC  = block.CooledPin.IC;
        block.PfactorSecondary = NetFlow(Secondary)/block.Pdrop2;
        
        block.MeasureCooledT.type = 'out';
        block.MeasureCooledT.IC = Secondary.T;
        
        block.P_Difference = {'ReformedPin','PrimaryPout'; 'CooledPin', 'SecondaryPout';};
        
%         Q = enthalpy(Inlet.Secondary) - enthalpy(Secondary);
%         ideal = Secondary;
%         ideal.T = Inlet.Primary.T;
%         Qideal = enthalpy(Inlet.Secondary) - enthalpy(ideal);
%         effectiveness = Q/Qideal
    end
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
    N = 2;
    if block.hotStream==1
        N=3;
    end
    block.ReformedPin.IC = Inlet.PrimaryPout+block.Pdrop1;
    SpecNew = fieldnames(Inlet.Primary);
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
            Inlet.Primary.(block.spec1{i})=0;
        end
    end
    block.S2C = Inlet.Primary.H2O/(Inlet.Primary.CH4 + .5*Inlet.Primary.CO);
    if block.hotStream==1
        block.CooledPin.IC = Inlet.SecondaryPout+block.Pdrop2;
        SpecNew = fieldnames(Inlet.Secondary);
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
                Inlet.Secondary.(block.spec2{i})=0;
            end
        end
    end
    
    %%-- %%
    [Primary, Secondary,block] = solveInitCond(Inlet,block);
    %%%
    
    block.ReformedPin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
    
    block.Primary.IC  = Inlet.Primary;
    block.PrimaryPout.IC = Inlet.PrimaryPout;
    block.Reformed.IC  = Primary;
    block.ReformedPin.IC  = block.ReformedPin.IC;
    block.MeasureReformT.IC = Primary.T;
    block.MeasureS2C.IC = block.S2C;
    block.PfactorPrimary = NetFlow(Primary)/block.Pdrop1;
    if block.hotStream==1
        block.CooledPin.Pstate = length(block.Scale); %identifies the state # of the pressure state if this block has one
        block.ReformedPin.Pstate = length(block.Scale)-1;
        block.Secondary.IC  = Inlet.Secondary;
        block.SecondaryPout.IC = Inlet.SecondaryPout; 
        block.Cooled.IC  = Secondary;
        block.CooledPin.IC  = block.CooledPin.IC;
        block.MeasureCooledT.IC = Secondary.T;        
        block.PfactorSecondary = NetFlow(Secondary)/block.Pdrop2;
%         Q = enthalpy(Inlet.Secondary) - enthalpy(Secondary);
%         ideal = Secondary;
%         ideal.T = Inlet.Primary.T;
%         Qideal = enthalpy(Inlet.Secondary) - enthalpy(ideal);
%         effectiveness = Q/Qideal
%         disp(Primary.T)
%         disp(Inlet.Secondary.T)
    end
end


function [Flow1, Flow2,block] = solveInitCond(Inlet,block)
nodes = block.nodes;
specInterest = {'CH4','CO','CO2','H2','H2O'};
spec = fieldnames(Inlet.Primary);
spec = spec(~strcmp('T',spec));
for i = 1:1:length(specInterest)
    if ~ismember(specInterest{i},spec)
        Inlet.Primary.(specInterest{i}) = 0;
    end
end

N = 2;
if block.hotStream ==1
    N=3;
end
n = N*nodes;
sp = length(block.spec1);
scaleFlow = NetFlow(Inlet.Primary)/sum(block.Scale(n+1:nodes:n+sp*nodes));
for i = 1:1:nodes
    block.Scale(n+i:nodes:n+sp*nodes) = block.Scale(n+i:nodes:n+sp*nodes)*scaleFlow;
end
if block.hotStream ==1
    n = (N+sp)*nodes;
    sp = length(block.spec2);
    scaleFlow = NetFlow(Inlet.Secondary)/sum(block.Scale(n+1:nodes:n+sp*nodes));
    for i = 1:1:nodes
        block.Scale(n+i:nodes:n+sp*nodes) = block.Scale(n+i:nodes:n+sp*nodes)*scaleFlow;
    end
end

Y = [block.Scale;1];
[T, Y] = ode15s(@(t,y) SolveDynamic(t,y,block,Inlet), [0, 1e2],Y);
Y = Y(end,:)';
n = 0;
Primary.T = Y(n+1:n+nodes);n = n+2*nodes;
if block.hotStream ==1
    Secondary.T = Y(n+1:n+nodes);n = n+nodes;
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
Flow1.T = Primary.T(block.PrimaryDir(1,end));
for i = 1:1:length(block.spec1)
    Flow1.(block.spec1{i}) = Primary.(block.spec1{i})(block.PrimaryDir(1,end));
end
if block.hotStream==1
    Flow2.T = Secondary.T(block.SecondaryDir(1,end));
    for i = 1:1:length(block.spec2)
        Flow2.(block.spec2{i}) = Secondary.(block.spec2{i})(block.SecondaryDir(1,end));
    end
else Flow2 = [];
end


function [T,HTconv,HTcond]= SteadyTemps(block,Qref,Inlet)    
nodes = block.nodes;
N = 2;
s=0;
hA = block.Convection;
Cp1 = mean(SpecHeat(Inlet.Primary));
Flow1 = NetFlow(Inlet.Primary);
if block.hotStream ==1
    N = 3;
    s = 1;
    Cp2 = mean(SpecHeat(Inlet.Secondary));
    Flow2 = NetFlow(Inlet.Secondary);
end
n = 2*N*nodes+2;

A = zeros(n,n);
b = zeros(n,1);
%Tprimary in 
A(n-s,n-s) = 1; %cold T
b(n-s) = Inlet.Primary.T;
if block.hotStream ==1 %T secondary in
    A(n,n) = 1; %hot T
    b(n) = Inlet.Secondary.T;
end

for k = 1:1:nodes
    %% averaging inlet and outlet temp for gaseous nodes
    %QT1 : heat transfer into primary flow
    A(k,k) = -1;
    A(k,k+(N+1)*nodes) = hA;
    A(k,k+N*nodes) = -.5*hA;
    [i,j] = find(block.PrimaryDir==k);
    if j==1 %first column averaged with inlet temperature
        A(k,n-s) = -.5*hA;
    else % other columns averaged with previous one
        k2 = block.PrimaryDir(i,j-1);
        A(k,k2+N*nodes) = -.5*hA;
    end

    %QT2 : heat transfer into solid
    A(k+nodes,k+nodes) = -1;
    if block.hotStream ==0
        A(k+nodes,k+(N+1)*nodes) = -1*hA;
        A(k+nodes,k+N*nodes) = .5*hA;
        if j==1 %first column averaged with inlet temperature
            A(k+nodes,n-s) = .5*hA;
        else % other columns averaged with previous one
            A(k+nodes,k2+N*nodes) = .5*hA;
        end
    elseif block.hotStream ==1
        A(k+nodes,k+(N+1)*nodes) = -2*hA;
        A(k+nodes,k+N*nodes) = .5*hA;
        A(k+nodes,k+(N+2)*nodes) = .5*hA;
        if j==1 %first column averaged with inlet temperature
            A(k+nodes,n-s) = .5*hA;
        else % other columns averaged with previous one
            A(k+nodes,k2+N*nodes) = .5*hA;
        end
        [i,j] = find(block.SecondaryDir==k);
        if j==1 %first column averaged with inlet temperature
            A(k+nodes,n) = .5*hA;
        else
            k2 = block.SecondaryDir(i,j-1);
            A(k+nodes,k2+(N+2)*nodes) = .5*hA;
        end

        %QT3 : heat transfer into hot flow
        A(k+2*nodes,k+2*nodes) = -1;
        A(k+2*nodes,k+(N+1)*nodes) = hA;
        A(k+2*nodes,k+(N+2)*nodes) = -.5*hA;
        if j==1 %first column averaged with inlet temperature
            A(k+2*nodes,n) = -.5*hA;
        else % other columns averaged with previous one
            A(k+2*nodes,k2+(N+2)*nodes) = -.5*hA;
        end
    end
    
    %Tcold: Temperature of cold flow
    A(k+N*nodes,k) = 1;
    A(k+N*nodes,k+N*nodes) = -Cp1*Flow1;
    [i,j] = find(block.PrimaryDir==k);
    if j==1 %first column receives fresh air
        A(k+N*nodes,n-s) = Cp1*Flow1;
    else
        index = block.PrimaryDir(i,j-1);
        A(k+N*nodes,index +N*nodes) = Cp1*Flow1;
    end
    b(k+N*nodes) = Qref(k);

    %Tsolid: Temperature of solid
    A(k+(N+1)*nodes,k+nodes) = 1;
    b(k+(N+1)*nodes) = 0;%Qref(k);
    if block.hotStream ==1
        %Temperature of secondary flow
        A(k+(N+2)*nodes,k+2*nodes) = 1;
        A(k+(N+2)*nodes,k+(N+2)*nodes) = -Cp2*Flow2; 
        [i,j] = find(block.SecondaryDir==k);
        if j==1 %first column receives fresh flow
            A(k+(N+2)*nodes,n) = Cp2*Flow2; %fresh inlet
        else
            index = block.SecondaryDir(i,j-1);
            A(k+(N+2)*nodes,index + (N+2)*nodes) = Cp2*Flow2;
        end
    end
end
% remove gaseous temperature averaging from 1st node
k1 = block.PrimaryDir(:,1); %first column
for j =1:1:length(k1)
    k = k1(j);
    A(k,k+N*nodes) = A(k,k+N*nodes) -.5*hA; %HT to cold flow from solid
    A(k,n-s) = A(k,n-s) +.5*hA; %ignore cold flow inlet

    A(k+nodes,k+N*nodes) = A(k+nodes,k+N*nodes) +.5*hA;%HT to solid from cold flow
    A(k+nodes,n-s) = A(k+nodes,n-s) -.5*hA;%ignore cold flow inlet
end
if block.hotStream ==1
    k1 = block.SecondaryDir(:,1); %first column
    for j =1:1:length(k1)
        k = k1(j);
        A(k+2*nodes,k+(N+2)*nodes) = A(k+2*nodes,k+(N+2)*nodes) -.5*hA; %HT to hot flow from solid
        A(k+2*nodes,n) = A(k+2*nodes,n) +.5*hA; %ignore hot flow inlet

        A(k+nodes,k+(N+2)*nodes) = A(k+nodes,k+(N+2)*nodes) +.5*hA;%HT to solid from hot flow
        A(k+nodes,n) = A(k+nodes,n) -.5*hA;%ignore hot flow inlet
    end
end
HTconv = A(1:N*nodes,N*nodes+1:2*N*nodes);%matrix of coefficients to multiply by vector of temperature and get the heat transfer by convection between layers and nodes
%% Conduction: prev and next
A2 = zeros(n,n);
prev = block.HTadjacent(:,1);
next = block.HTadjacent(:,2);
for k = 1:1:nodes
    A2(k+nodes,k+4*nodes) = A2(k+nodes,k+4*nodes) - 2*block.Conduction;
    A2(k+nodes,prev(k)+4*nodes) = A2(k+nodes,prev(k)+4*nodes)+block.Conduction;
    A2(k+nodes,next(k)+4*nodes) = A2(k+nodes,next(k)+4*nodes)+block.Conduction;
end
HTcond = A2(1:N*nodes,N*nodes+1:2*N*nodes);%matrix of coefficients to multiply by vector of temperature and get the heat transfer by convection between layers and nodes
A = A+A2;
x= A\b;
T = x(N*nodes+1:n-s-1);


function dY = SolveDynamic(t,Y,block,Inlet)
nodes = block.nodes;
block.HTconv = Y(end)*block.HTconv;
dY = 0*Y;
n = 0;
Primary.Outlet.T = Y(n+1:n+nodes);n = n+2*nodes;
s = 0;
N = 2;
if isfield(Inlet,'Secondary')
    SecondaryOut.T = Y(n+1:n+nodes);n = n+nodes;
    s = 1;
    N = 3;
end
%% Cold flow
for i = 1:1:length(block.spec1)
    Primary.Outlet.(block.spec1{i}) = max(0,Y(n+1:n+nodes));n = n+nodes;
end
for j = 1:1:length(block.PrimaryDir(1,:));%1:columns
    k = block.PrimaryDir(:,j);
    r = length(k);
    if j==1
        Primary.Inlet.T(k,1) = Inlet.Primary.T;
        for i = 1:1:length(block.spec1)
            Primary.Inlet.(block.spec1{i})(k,1) = Inlet.Primary.(block.spec1{i})/r;
        end
    else
        Primary.Inlet.T(k,1) = Primary.Outlet.T(kprev,1);
        for i = 1:1:length(block.spec1)
            Primary.Inlet.(block.spec1{i})(k,1) = Primary.Outlet.(block.spec1{i})(kprev,1);
        end
    end
    kprev = k;
end
if isfield(Inlet,'Secondary')
    %% Hot flow
    for i = 1:1:length(block.spec2)
        SecondaryOut.(block.spec2{i}) = max(0,Y(n+1:n+nodes));n = n+nodes;
    end
    for j = 1:1:length(block.SecondaryDir(1,:));%1:columns
        k = block.SecondaryDir(:,j);
        r = length(k);
        if j==1
            SecondaryIn.T(k,1) = Inlet.Secondary.T;
            for i = 1:1:length(block.spec2)
                SecondaryIn.(block.spec2{i})(k,1) = Inlet.Secondary.(block.spec2{i})/r;
            end
        else
            SecondaryIn.T(k,1) = SecondaryOut.T(kprev);
            for i = 1:1:length(block.spec2)
                SecondaryIn.(block.spec2{i})(k,1) = SecondaryOut.(block.spec2{i})(kprev,1);
            end
        end
        kprev = k;
    end
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

if isfield(Inlet,'Secondary')
    HoutHot = enthalpy(SecondaryOut);
    scale = NetFlow(SecondaryOut)./NetFlow(SecondaryIn);
    HinHot = enthalpy(SecondaryIn).*scale;
    QT = block.HTconv*Y(1:N*nodes) + block.HTcond*Y(1:N*nodes);
else
    QT = block.HTconv*Y(1:N*nodes) + block.HTcond*Y(1:N*nodes);
end

for i=1:1:length(block.PrimaryDir(1,:))
    k = block.PrimaryDir(:,i);
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

if isfield(Inlet,'Secondary')
    for i=1:1:length(block.SecondaryDir(1,:))
        k = block.SecondaryDir(:,i);
        dY(2*nodes+k)= (QT(2*nodes+k) + HinHot(k) - HoutHot(k))./block.tC(2*nodes+k); %Hot flow
        if i>1
            dY(2*nodes+k) = dY(2*nodes+k)+dY(2*nodes+kprev);
        end
        kprev = k;
    end
    n = (2+j);
    for j = 1:1:length(block.spec2)
        dY((n+j)*nodes+1:(n+1+j)*nodes) = (SecondaryIn.(block.spec2{j}) - SecondaryOut.(block.spec2{j}))./block.tC((n+j)*nodes+1:(n+1+j)*nodes);  %all  species concentration
    end
end

if strcmp(block.method,'RefPerc')
    RefPerc = (Inlet.Primary.CH4-Primary.Outlet.CH4(block.PrimaryDir(1,end)))/Inlet.Primary.CH4;
    error = block.ReformTarget - RefPerc;
elseif strcmp(block.method,'ColdT')
    error = (block.ReformTarget - Primary.T(block.PrimaryDir(1,end)))/100;
elseif strcmp(block.method,'HotT')
    error = (Secondary.T(block.SecondaryDir(1,end)) - block.ReformTarget)/100;
elseif strcmp(block.method,'Effectiveness')
    FlowIdeal = Inlet.Secondary;
    FlowIdeal.T = Primary.Outlet.T(block.SecondaryDir(1,end));
    maxQT = enthalpy(FlowIdeal) - enthalpy(Inlet.Secondary); %assume ideal is at same temperature as 1st node on primary side
    Flow2.T = Secondary.T(block.SecondaryDir(1,end));
    for i = 1:1:length(block.spec2)
        Flow2.(block.spec2{i}) = Secondary.(block.spec2{i})(block.SecondaryDir(1,end));
    end
    QT = enthalpy(Flow2) - enthalpy(Inlet.Secondary);
    error = block.ReformTarget - QT/maxQT;
elseif strcmp(block.method,'none')
    error = 0;
end
dY(end) = error; %slow change in the area


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
    X = Primary.(block.spec1{i})./Flow;%concentration
    if any(X<.01) %need to prevent possible divide by zero
        block.IC(n+1:n+nodes,1) = X;
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
    for i = 1:1:length(block.spec2)
        if any(Secondary.(block.spec2{i})==0) %need to prevent possible divide by zero
            Flow = NetFlow(Secondary);
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


function block = FlowDir(block)
% Script which orients the nodes relative to the flow directions
% direction = 1: co-flow, direction = 2: counter-flow, direction = 3: cross-flow
nodes = block.nodes;
block.PrimaryDir = (1:nodes);
if block.hotStream ==1
    switch block.direction
        case 'coflow'
            block.SecondaryDir = (1:nodes); % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
        case 'counterflow'
            block.SecondaryDir = (nodes:-1:1); % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
    end
end
block.HTadjacent = zeros(nodes,4);
for i = 1:1:nodes
    block.HTadjacent(i,1) = max(1,i-1);%previous node
    block.HTadjacent(i,2) = min(nodes,i+1);%next node
end