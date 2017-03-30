function Out = Turbine(t,Y, Inlet,block,string1)
% a simple turbine model with 1 inlet flow, using turbine maps to relat mass flow and efficiency to pressure ratio and RPM
% Three (3) inlets: air flow,  outlet pressure, and shaft RPM
% Three (3) outlets: Pressure in, Flow, Work output
% Three (3) states: Tgas out, Twall, and pressure
global Ru Tags
NetFlowIn = NetFlow(Inlet.FlowIn);
Mmass = MassFlow(Inlet.FlowIn)/NetFlowIn;

%calculate map position
nT =Inlet.FlowIn.T/block.Tdesign;% normalized T
nRPM =Inlet.RPMin/(block.RPMdesign*nT^.5);%normalized RPM
TurbPR =(Y(3)/Inlet.Pout-1)/(block.Pdesign - 1) + 1;%Turbine pressure ratio
    
i2 = find(block.RPM >=nRPM,1,'first');
if isempty(i2)
    i2 = length(block.RPM);
elseif i2 ==1
    i2=2;
end
i1 = i2-1;
s =(nRPM - block.RPM(i1))/(block.RPM(i2) - block.RPM(i1));

PRVec = block.PressRatio(i1,:)*(1-s) + block.PressRatio(i2,:)*s;
Beta = interp1(PRVec,block.Beta,TurbPR,'spline');
Eff = interp2(block.Beta,block.RPM,block.Efficiency,Beta,nRPM,'spline')*block.PeakEfficiency;
Nflow = interp2(block.Beta,block.RPM,block.NflowGMap,Beta,nRPM,'spline');
Mscale = block.FlowDesign*(Y(3)/Inlet.Pout)/block.Pdesign/(nT^.5*Mmass);

%calculate outlet flow
NetFlowOut = Nflow*Mscale;
specname = fieldnames(Inlet.FlowIn);
for i = 1:1:length(specname)
    if ~strcmp(specname{i},'T')
        FlowOut.(specname{i}) = Inlet.FlowIn.(specname{i})*NetFlowOut/NetFlowIn;
    end
end
FlowOut.T = Y(1);

%calculate heat transfer
Q_WallAmbC =(Y(2) - block.Tamb)*block.AmbConvC*block.SurfA/1000;%Convection from wall to ambient
Q_WallAmbR =((Y(2))^4 - (block.Tamb)^4)*block.Epsilon*block.Sigma*block.SurfA/1000;%Radiation from wall to ambient
Q_FlowWallC =(Y(1) - Y(2))*block.ConvCoef*block.SurfA/1000;%Convection from flow to wall

%% calculate energy balance
Cp = (SpecHeat(Inlet.FlowIn)+SpecHeat(FlowOut))/2;
Gamma = Cp/(Cp - Ru);
TurbFlow = FlowOut;
TurbFlow.T = Inlet.FlowIn.T; %inlet flow at same rate as exit, but inlet temp
H1 = enthalpy(TurbFlow);
T2s = Inlet.FlowIn.T*(1/(Y(3)/Inlet.Pout))^((Gamma -1)/Gamma);
FlowS = TurbFlow;
FlowS.T = T2s;
H2s = enthalpy(FlowS);
H2a = H1-(H1 - H2s)*Eff - Q_FlowWallC;
Wt = (H1-H2a) - Q_FlowWallC;

if strcmp(string1,'Outlet')
    Out.PowerTurb = Wt;
    Out.Outlet = FlowOut;
    Out.TET = FlowOut.T;
    Out.Pin = Y(3);
    Tags.(block.name).TIT = Inlet.FlowIn.T;
    Tags.(block.name).TET = FlowOut.T;
    Tags.(block.name).Power = Wt;
    Tags.(block.name).PR = Y(3)/Inlet.Pout;
    Tags.(block.name).Nflow = Nflow;
    Tags.(block.name).NRPM = nRPM;
    Tags.(block.name).Efficiency = Eff;
    Tags.(block.name).MassFlow = MassFlow(FlowOut);
    Tags.(block.name).Beta = Beta;
    Tags.(block.name).nMflow = MassFlow(FlowOut)/block.FlowDesign;
elseif strcmp(string1,'dY')
    dY = 0*Y;
    Hout = enthalpy(FlowOut);
    dY(1) = (H2a - Hout)*Ru*block.Tdesign/(Inlet.Pout*block.Volume*Cp);
    dY(2) = (Q_FlowWallC - Q_WallAmbC - Q_WallAmbR)/(block.Mass*block.SpecHeat);
    dY(3) = (NetFlowIn - NetFlowOut)*Ru*Y(1)/block.Volume;
    Out = dY;
end