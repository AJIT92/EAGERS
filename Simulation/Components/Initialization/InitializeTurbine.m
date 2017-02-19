function block = InitializeTurbine(varargin)
% a simple turbine model with 1 inlet flow, using turbine maps to relat mass flow and efficiency to pressure ratio and RPM
% Three (3) inlets: flow,  outlet pressure, and shaft RPM
% Three (3) outlets: Pressure in, Flow , Work output
% Three (3) states: Tgas out, Twall, and pressure
global Ru
Ru = 8.314472; % Universal gas constant in kJ/K*kmol
block = varargin{1};
if length(varargin)==1 % first initialization
    block.AmbConvC = 5;%ambient convection coefficient
    block.ConvCoef = 5;%convection coefficient
    block.Epsilon = .8;%radiation epsilon value
    block.Sigma = 5.67e-8;%Radiation sigma value
    block.SpecHeat = 0.5;%Specific Heat of Turbine kJ/(kg*K)

    block.Length = .03;
    block.Diameter = 0.2;
    block.Volume = pi/4*block.Diameter^2*block.Length;
    block.SurfA = pi*block.Diameter*block.Length;%surface area of compressor

    block.Tamb = 330;%ambient air temp
    T2s = block.Tdesign*(1/block.Pdesign)^((1.4 -1)/1.4);
    T2a = block.Tdesign - block.PeakEfficiency*(block.Tdesign-T2s);
    block.Scale =[T2a; (T2a+block.Tamb)/2; 101*block.Pdesign]; %Tflow, TTurb, Pressure
    block.IC = ones(length(block.Scale),1);
    
    Dir=strrep(which('InitializeTurbine.m'),fullfile('Components','Initialization','InitializeTurbine.m'),'CompressorMaps');
    load(fullfile(Dir,block.Map));
    f = fieldnames(map);
    for i = 1:1:length(f)
        block.(f{i}) = map.(f{i});
    end
    %%
    i = find(block.RPM>=1,1);
    j = find(block.NflowGMap(i,:)>=1,1);
    P1 = block.PressRatio(i,j);
    P2 = block.PressRatio(i,j-1);
    M1 = block.NflowGMap(i,j);
    M2 = block.NflowGMap(i,j-1);
    block.mFlow = block.FlowDesign;
    dMdP = block.FlowDesign*(M2 - M1)/((block.Pdesign-1)*(P2 - P1));
    C = (dMdP*101*block.Pdesign-block.mFlow)/101;
    block.dMdP = [dMdP, C];
    %%
    block.PortNames = {'FlowIn','Pout','RPMin','Pin','Outlet','PowerTurb','TET'};

    block.FlowIn.type = 'in';
    block.FlowIn.IC.T = block.Tdesign;
    block.FlowIn.IC.CO2 = 0.1;
    block.FlowIn.IC.H2O = 0.1;
    block.FlowIn.IC.N2 = 0.7;
    block.FlowIn.IC.O2 = 0.1;
    mMass = MassFlow(block.FlowIn.IC); %molar mass
    specFlow = fieldnames(block.FlowIn.IC);
    for i = 1:1:length(specFlow)%convert to apropriate flow rate
        if ~strcmp(specFlow{i},'T')
            block.FlowIn.IC.(specFlow{i}) = block.FlowIn.IC.(specFlow{i})*block.FlowDesign/mMass;
        end
    end

    block.Pout.type = 'in';
    block.Pout.IC = 101;%in kPa
    block.Pout.Pstate = []; %identifies the state # of the pressure state if this block has one
    
    block.RPMin.type = 'in';
    block.RPMin.IC = block.RPMdesign;

    block.Pin.type = 'out';
    block.Pin.IC = 101*block.Pdesign;%in kPa
    block.Pin.Pstate = 3; %identifies the state # of the pressure state if this block has one

    block.Outlet.type = 'out';
    block.Outlet.IC = block.FlowIn.IC;
    block.Outlet.IC.T = T2a;

    block.PowerTurb.type = 'out';
    [~,H1] = enthalpy(block.FlowIn.IC);
    [~,H2] = enthalpy(block.Outlet.IC);
    block.PowerTurb.IC = H1-H2;

    block.TET.type = 'out';
    block.TET.IC = T2a;

    block.P_Difference = {'Pin','Pout'};
    
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
end
if length(varargin)==2 %% Have inlets connected, re-initialize
    Inlet = varargin{2};
    NetFlowIn = NetFlow(Inlet.FlowIn);
    Mmass = MassFlow(Inlet.FlowIn)/NetFlowIn;
    Pin = block.IC(3)*block.Scale(3);
    nT =Inlet.FlowIn.T/block.Tdesign;% normalized T
    nRPM =Inlet.RPMin/(block.RPMdesign*nT^.5);%normalized RPM
    Mscale = block.FlowDesign*(block.Scale(3)/Inlet.Pout)/block.Pdesign/(Mmass*nT^.5);%scalar maping for mass flow
    Nflow = NetFlow(Inlet.FlowIn)/Mscale;
    TurbPR1 =(Pin/Inlet.Pout-1)/(block.Pdesign - 1) + 1;%Turbine pressure ratio
%% Find table indecies (i,j) surrounding actual value
% s and r correspond to the fractional coordinate position 
    i2 = find(block.RPM>=nRPM,1,'first');
    if isempty(i2)
        i2 = length(block.RPM);
    elseif i2 ==1
        i2=2;
    end
    i1 = i2-1;
    r = (nRPM-block.RPM(i1))/(block.RPM(i2)-block.RPM(i1));
    
    FlowVec = block.NflowGMap(i1,:)*(1-r) + block.NflowGMap(i2,:)*r;
    Nflow = min(FlowVec(end),max(FlowVec(1),Nflow));
    Beta = interp1(FlowVec,block.Beta,Nflow,'spline');
    TurbPR2 = interp2(block.Beta,block.RPM,block.PressRatio,Beta,nRPM,'spline');
    TurbPR = .7*TurbPR1 + .3*TurbPR2;%move in direction of the current inlet mass flow
    Eff = interp2(block.Beta,block.RPM,block.Efficiency,Beta,nRPM,'spline')*block.PeakEfficiency;
    
    j12 = find(block.PressRatio(i1,:)>=TurbPR,1,'first');
    if isempty(j12)
        j12 = length(block.PressRatio(i1,:));
    elseif j12 ==1
        j12=2;
    end
    j11 = j12-1;
    j22 = find(block.PressRatio(i2,:)>=TurbPR,1,'first');
    if isempty(j22)
        j22 = length(block.PressRatio(i2,:));
    elseif j22 ==1
        j22=2;
    end
    j21 = j22-1;
    P1 = block.PressRatio(i1,j11)*r+ block.PressRatio(i2,j21)*(1-r);
    P2 = block.PressRatio(i1,j12)*r+ block.PressRatio(i2,j22)*(1-r);
    P1 = Inlet.Pout*((P1 -1)*(block.Pdesign -1) + 1);
    P2 = Inlet.Pout*((P2 -1)*(block.Pdesign -1) + 1);
    M1 = block.NflowGMap(i1,j11)*r+ block.NflowGMap(i2,j21)*(1-r);
    M2 = block.NflowGMap(i1,j12)*r+ block.NflowGMap(i2,j22)*(1-r);
     %find flow out at this estimated pressure ratio
    PRVec = block.PressRatio(i1,:)*(1-r) + block.PressRatio(i2,:)*r;
    Beta = interp1(PRVec,block.Beta,TurbPR,'spline');
    Nflow = interp2(block.Beta,block.RPM,block.NflowGMap,Beta,nRPM,'spline');
    
    %calculate outlet flow
    NetFlowOut = Nflow*Mscale;
    specname = fieldnames(Inlet.FlowIn);
    for i = 1:1:length(specname)
        if ~strcmp(specname{i},'T')
            FlowOut.(specname{i}) = Inlet.FlowIn.(specname{i})*NetFlowOut/NetFlowIn;
        end
    end
    FlowOut.T = block.IC(1)*block.Scale(1);
    
    Cp = (SpecHeat(Inlet.FlowIn)+SpecHeat(FlowOut))/2;
    Gamma = Cp/(Cp - Ru);
    %update dMdP and mFlow
    block.mFlow = MassFlow(FlowOut);
    dMdP = Mscale*(M2 - M1)/(P2 - P1);
    PinNew = ((TurbPR-1)*(block.Pdesign - 1)+1)*Inlet.Pout;
    C = (dMdP*PinNew-block.mFlow)/Inlet.Pout;
    block.dMdP = [dMdP, C];
    
    %calculate energy balance
    TurbFlow = FlowOut;
    TurbFlow.T = Inlet.FlowIn.T; %same mass flow but at inlet temp
    [~, H1] = enthalpy(TurbFlow);
    T2s = Inlet.FlowIn.T*(1/(PinNew/Inlet.Pout))^((Gamma -1)/Gamma);
    FlowS = TurbFlow;
    FlowS.T = T2s;
    [~,H2s] = enthalpy(FlowS);
  
    errorT =1;
    Wall.T = (block.Tamb + FlowOut.T)/2;
    tC = (Cp*NetFlowIn);
    while abs(errorT)>.01 %solve for the correct outlet temperature
       [~,Hout] = enthalpy(FlowOut);
       %calculate heat transfer
       Q_WallAmbC =(Wall.T - block.Tamb)*block.AmbConvC*block.SurfA/1000;%Convection from wall to ambient
       Q_WallAmbR =(Wall.T^4 - (block.Tamb)^4)*block.Epsilon*block.Sigma*block.SurfA/1000;%Radiation from wall to ambient
       Q_FlowWallC =(FlowOut.T - Wall.T)*block.ConvCoef*block.SurfA/1000;%Convection from flow to wall
       
       H2a = H1-(H1 - H2s)*Eff - Q_FlowWallC;
       Wt = (H1-H2a) - Q_FlowWallC;
        
       errorT = (H2a - Hout)/tC;
       FlowOut.T = FlowOut.T + errorT;
       Wall.T = Wall.T + (Q_FlowWallC - Q_WallAmbC - Q_WallAmbR)/tC;
    end

    %%
    block.Scale = [FlowOut.T; Wall.T; block.Scale(3)]; %Tflow, TTurb, Pressure
    block.Pout.IC = Inlet.Pout;
    block.Pin.IC = Pin; 
    block.Outlet.IC = FlowOut;
    block.FlowIn.IC = Inlet.FlowIn;
    block.PowerTurb.IC = Wt;
    block.TET.IC = FlowOut.T;
end