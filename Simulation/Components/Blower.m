function Out = Blower(t,Y, Inlet,block,string1)
% Four (4) inlets: air flow,  outlet pressure, inlet pressure, power supply
% One (1) outlet: Flow
% One (1) state: Speed
global Ru Tags
Y = Y.*block.Scale;
Mmass = MassFlow(Inlet.Species);%kg/kmol
RPM = Y(1)*60/(2*pi); %convert rad/s back to RPM

%compressor map
nT = Inlet.Temperature/block.Tdesign;%square root of the normalized T
nRPM = RPM/(block.RPMdesign*nT^.5);%normalized RPM
CompPR = (Inlet.Pout/Inlet.Pin -1)/(block.Pdesign - 1) + 1;%normalized compressor pressure ratio
Mscale = block.FlowDesign*(Inlet.Pin/block.P0map)/(Mmass*nT^.5);%scalar maping for mass flow        
i2 = find(block.RPM >=nRPM,1,'first');
if isempty(i2)
    i2 = length(block.RPM);
elseif i2==1
    i2 =2;
end
i1 = i2-1;
s =(nRPM - block.RPM(i1))/(block.RPM(i2) - block.RPM(i1));
PRVec = block.PressRatio(i1,:)*(1-s) + block.PressRatio(i2,:)*s;
if CompPR>PRVec(end)
    disp('Stall Reached')
    Beta = CompPR/PRVec(end);
    Eff = block.minEfficiency;
    flow_RPM = block.NflowGMap(i1,end)*(1-s) + block.NflowGMap(i2,end)*s;
    Nflow = 1/Beta*flow_RPM;
else
    Beta = interp1(PRVec,block.Beta,CompPR,'spline');
    Eff = interp2(block.Beta,block.RPM,block.Efficiency,Beta,nRPM,'spline')*block.PeakEfficiency;
    Nflow = interp2(block.Beta,block.RPM,block.NflowGMap,Beta,nRPM,'spline');
end
%outlet flow
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
NetPower = Inlet.Power - (H2a-H1);

if strcmp(string1,'Outlet')
    Out.Outlet = Flow;
    Out.Speed = RPM;
    Tags.(block.name).RPM = RPM;
    Tags.(block.name).NRPM = nRPM;
    Tags.(block.name).Beta = Beta;
    Tags.(block.name).Flow = Flow;
    Tags.(block.name).Power = Inlet.Power;
    Tags.(block.name).PR = Inlet.Pout/Inlet.Pin;
    Tags.(block.name).Nflow = Nflow;
    Tags.(block.name).MassFlow = MassFlow(Flow);
    Tags.(block.name).Temperature = Flow.T;
    Tags.(block.name).Eff = Eff;
    Tags.(block.name).nMflow = MassFlow(Flow)/block.FlowDesign;
    Tags.(block.name).Pressure = Inlet.Pout;
elseif strcmp(string1,'dY')
    dY = 0*Y;
    dY(1) = NetPower/(block.Moment_Inertia*Y(1)); %dw/dt = P/(w*l) units of rad/s
    dY = dY./block.Scale;
    Out = dY;
end
