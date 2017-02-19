function block = InitializeBlower(varargin)
% a simple blower model
% Four (4) inlets: air flow,  outlet pressure, inlet pressure, power supply
% One (1) outlet: Flow
% One (1) state: Speed
global Ru Tags
Ru = 8.314472; % Universal gas constant in kJ/K*kmol
block = varargin{1};
if length(varargin)==1 % first initialization
    
    block.Mass = 1; % mass in kg
	block.InnerRadius = 0.05;
    block.OuterRadius = 0.15;
    block.Moment_Inertia = block.Mass/2*(block.InnerRadius^2 + block.InnerRadius^2);%Moment of Inertia for thick walled cylinder
    
    block.Scale = block.RPMdesign*2*pi/60; %speed in RPM converted to rad/s
    block.IC = ones(length(block.Scale),1);
    Dir=strrep(which('InitializeBlower.m'),fullfile('Components','Initialization','InitializeBlower.m'),'CompressorMaps');
    load(fullfile(Dir,block.Map));
    f = fieldnames(map);
    for i = 1:1:length(f)
        block.(f{i}) = map.(f{i});
    end
    Eff = map.Efficiency;
    Eff(Eff<0)=nan;
    block.minEfficiency = min(min(Eff));
    block.Efficiency(block.Efficiency<0) = block.minEfficiency;
    
    block.spec = fieldnames(block.InitialComposition);
    Mmass = MassFlow(block.InitialComposition);%kg/kmol
    Flow.T = block.Tdesign;
    for i = 1:1:length(block.spec)
        Flow.(block.spec{i}) = block.InitialComposition.(block.spec{i})*block.FlowDesign/Mmass;
    end
    Cp = SpecHeat(Flow);
    H1 = enthalpy(Flow);
    T2s = Flow.T*(block.Pdesign)^((1.4 -1)/1.4);
    Flow.T = T2s;
    Cp = (SpecHeat(Flow)+Cp)/2;
    Gamma = Cp/(Cp - Ru);
    T2s = Flow.T*(block.Pdesign)^((Gamma -1)/Gamma);
    Flow.T = T2s;
    H2s = enthalpy(Flow);  
    H2a = H1+(H2s-H1)/block.PeakEfficiency;
    block.NominalPower = H2a - H1;
       
    %% set up dM/dP*Pout - C*Pin = mdot
    i = find(block.RPM>=1,1);
    j = ceil(length(block.NflowGMap(i,:))/2);
    P1 = block.PressRatio(i,j);
    P2 = block.PressRatio(i,j-1);
    M1 = block.NflowGMap(i,j);
    M2 = block.NflowGMap(i,j-1);
    block.mFlow = block.FlowDesign;
    dMdP = block.FlowDesign*(M2 - M1)/(101*(block.Pdesign-1)*(P2 - P1));
    C = (dMdP*101*block.Pdesign-block.mFlow)/101;
    block.dMdP = [dMdP, C]; 

    %%
    block.PortNames = {'Temperature','Species','Pin','Pout','Power','Outlet','Speed'};
    
    block.Temperature.type = 'in';
    block.Temperature.IC = block.Tdesign; 
    
    block.Species.type = 'in';
    block.Species.IC = block.InitialComposition; 
    
    block.Pin.type = 'in';
    block.Pin.IC = 101;%in kPa
    block.Pin.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.Pout.type = 'in';
    block.Pout.IC = block.Pin.IC*block.Pdesign;%in kPa
    block.Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.Power.type = 'in';
    block.Power.IC = block.NominalPower;
    
    block.Outlet.type = 'out';
    block.Outlet.IC = Flow;
      
    block.Speed.type = 'out';
    block.Speed.IC = block.RPMdesign;
    
    block.P_Difference = {'Pout','Pin'};

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
    
    Tags.(block.name).RPM = block.RPMdesign;
    Tags.(block.name).Flow = Flow;
    Tags.(block.name).Power = block.NominalPower;
    Tags.(block.name).MassFlow = MassFlow(Flow);
    Tags.(block.name).Temperature = Flow.T;
    Tags.(block.name).Eff = block.PeakEfficiency;
end
if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    block.spec = fieldnames(Inlet.Species);
    
    %problem is to find speed and flow that matches the PR and power specified by input
    RPM = Tags.(block.name).RPM;
    error = Inlet.Power;
    nRPM = RPM/(block.RPMdesign*(Inlet.Temperature/block.Tdesign)^.5);%normalized RPM
    Tol = (1e-3)*Inlet.Power;
    while abs(error)>Tol
        [Power,Flow,block] = OpPoint(nRPM,Inlet,block);
        if block.RPM(end)<=nRPM
            dRPM = -1e-5*nRPM;
        else
            dRPM = 1e-5*nRPM;
        end
        [Power2,~,~] = OpPoint(nRPM+dRPM,Inlet,block);
        dPdRPM = (Power2-Power)/dRPM;
        error = (Inlet.Power - Power);
        if dPdRPM<0
            dPdRPM = Inlet.Power*2;
        end
        
        if block.RPM(end)<=nRPM && error>0
            error = 0; %at edge of map, can't increase power anymore
        elseif block.RPM(1)>=nRPM && error<0
            error = 0; % at lower edge of map, can't decrease power
        else
            nRPM = max(min(nRPM + error/dPdRPM,block.RPM(end)),block.RPM(1));
        end
        
    end
    RPM = nRPM*(block.RPMdesign*(Inlet.Temperature/block.Tdesign)^.5);
    
    FlowIn = Flow;
    FlowIn.T = Inlet.Temperature;
    H2a = enthalpy(FlowIn)+Inlet.Power;
    Cp = SpecHeat(Flow);
    tC = (Cp*NetFlow(Flow));
    errorT =1;
    while abs(errorT)>.01 %solve for the correct outlet temperature
        Hout = enthalpy(Flow);
        errorT = (H2a - Hout)/tC;
        Flow.T = Flow.T + errorT;
    end
    
    block.IC = (RPM*2*pi/60)/block.Scale; %speed in RPM converted to rad/s
    block.Outlet.IC = Flow;
    block.Speed.IC = RPM;
    block.mFlow = MassFlow(Flow);
    
    %% set up dM/dP*Pout - C*Pin = mdot
    Mscale = block.FlowDesign*(Inlet.Pin/block.P0map)/(MassFlow(Inlet.Species)*(Inlet.Temperature/block.Tdesign)^.5);%scalar maping for mass flow
    CompPR = (Inlet.Pout/Inlet.Pin -1)/(block.Pdesign - 1) + 1;%normalized compressor pressure ratio
    i2 = find(block.RPM >=nRPM,1,'first');
    if isempty(i2)
        i2 = length(block.RPM);
    elseif i2==1
        i2 =2;
    end
    i1 = i2-1;
    r = (nRPM-block.RPM(i1))/(block.RPM(i2)-block.RPM(i1));
    j12 = find(block.PressRatio(i1,:)>=CompPR,1,'first');
    if isempty(j12)
        j12 = length(block.PressRatio(i1,:));
    elseif j12==1
        j12 =2;
    end
    j11 = j12-1;
    j22 = find(block.PressRatio(i2,:)>=CompPR,1,'first');
    if isempty(j22)
        j22 = length(block.PressRatio(i2,:));
    elseif j22==1
        j22 = 2;
    end
    j21 = j22-1;
    P1 = block.PressRatio(i1,j11)*r+ block.PressRatio(i2,j21)*(1-r);
    P2 = block.PressRatio(i1,j12)*r+ block.PressRatio(i2,j22)*(1-r);
    P1 = Inlet.Pin*((P1 -1)*(block.Pdesign -1) + 1);
    P2 = Inlet.Pin*((P2 -1)*(block.Pdesign -1) + 1);
    M1 = block.NflowGMap(i1,j11)*r+ block.NflowGMap(i2,j21)*(1-r);
    M2 = block.NflowGMap(i1,j12)*r+ block.NflowGMap(i2,j22)*(1-r);
    block.dMdP(1,1) = Mscale*(M2 - M1)/(P2 - P1);
    block.dMdP(1,2) =  (block.dMdP(1,1)*Inlet.Pout-block.mFlow)/Inlet.Pin;

    Tags.(block.name).RPM = RPM;
    Tags.(block.name).Flow = Flow;
    Tags.(block.name).Power = Inlet.Power;
    Tags.(block.name).MassFlow = MassFlow(Flow);
    Tags.(block.name).Temperature = Flow.T;
end

function [Power,Flow,block] = OpPoint(nRPM,Inlet,block)
global Ru
CompPR = (Inlet.Pout/Inlet.Pin -1)/(block.Pdesign - 1) + 1;%normalized compressor pressure ratio
Mmass = MassFlow(Inlet.Species);%kg/kmol
Mscale = block.FlowDesign*(Inlet.Pin/block.P0map)/(Mmass*(Inlet.Temperature/block.Tdesign)^.5);%scalar maping for mass flow
%% Find table indecies (i,j) surrounding actual value
% s and r correspond to the fractional coordinate position 
i2 = find(block.RPM >=nRPM,1,'first');
if isempty(i2)
    i2 = length(block.RPM);
elseif i2==1
    i2 =2;
end
i1 = i2-1;
s =(nRPM - block.RPM(i1))/(block.RPM(i2) - block.RPM(i1));
PRVec = block.PressRatio(i1,:)*(1-s) + block.PressRatio(i2,:)*s;
if CompPR>max(PRVec) && nRPM<1
    CompPR = (Inlet.Pout/Inlet.Pin -1)*nRPM/(block.Pdesign - 1) + 1;%renormalized compressor pressure ratio
end
Beta = interp1(PRVec,block.Beta,CompPR,'spline');
Nflow = interp2(block.Beta,block.RPM,block.NflowGMap,Beta,nRPM,'spline');
Eff = interp2(block.Beta,block.RPM,block.Efficiency,Beta,nRPM,'spline')*block.PeakEfficiency;
NetFlowOut = Nflow*Mscale;

Flow.T = Inlet.Temperature;
for i = 1:1:length(block.spec)
    Flow.(block.spec{i}) = Inlet.Species.(block.spec{i})*NetFlowOut;
end
Cp = SpecHeat(Flow);
H1 = enthalpy(Flow);
T2s = Flow.T*(Inlet.Pout/Inlet.Pin)^((1.4 -1)/1.4);
Flow.T = T2s;
Cp = (SpecHeat(Flow)+Cp)/2;
Gamma = Cp/(Cp - Ru);
T2s = Flow.T*(Inlet.Pout/Inlet.Pin)^((Gamma -1)/Gamma);
Flow.T = T2s;
H2s = enthalpy(Flow);  
H2a = H1+(H2s-H1)/Eff;
Power = (H2a - H1);