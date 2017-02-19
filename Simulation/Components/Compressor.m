function Out = Compressor(t,Y, Inlet,block,string1)
% Four (4) inlets: air flow,  outlet pressure, inlet pressure, shaft RPM
% Two (2) outlets: Flow, Work input required
% Two (2) states: Tgas out and Twall
global Ru Tags
Y = Y.*block.Scale;
NetFlowIn = NetFlow(Inlet.FlowIn);
Mmass = MassFlow(Inlet.FlowIn)/NetFlowIn;
%compressor map
nT =Inlet.FlowIn.T/block.Tdesign;%square root of the normalized T
nRPM =Inlet.RPMin/(block.RPMdesign*nT^.5);%normalized RPM
CompPR =(Inlet.Pout/Inlet.Pin -1)/(block.Pdesign - 1) + 1;%normalized compressor pressure ratio

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
%     disp('Stall Reached')
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
NetFlowOut = Nflow*Inlet.Pin/block.P0map*block.FlowDesign/Mmass/nT^.5;
specname = fieldnames(Inlet.FlowIn);
for i = 1:1:length(specname)
    if ~strcmp(specname{i},'T')
        FlowOut.(specname{i}) = Inlet.FlowIn.(specname{i})*NetFlowOut/NetFlowIn;
    end
end
FlowOut.T = Y(1);

%heat transfer
Q_WallAmbC =(Y(2) - block.Tamb)* block.AmbConvC*block.SurfA/1000;%Convection from wall to ambient
Q_WallAmbR =((Y(2))^4 - (block.Tamb)^4)*block.Epsilon*block.Sigma*block.SurfA/1000;%Radiation from wall to ambient
Q_FlowWallC =(Y(1) - Y(2)) * block.ConvCoef * block.SurfA/1000;%Convection from flow to wall
    
%energy balance
Cp = (SpecHeat(Inlet.FlowIn)+SpecHeat(FlowOut))/2;
Gamma = Cp/(Cp - Ru);
T2s = Inlet.FlowIn.T*(Inlet.Pout/Inlet.Pin)^((Gamma -1)/Gamma);

FlowComp = FlowOut;
FlowComp.T = Inlet.FlowIn.T;
[~, H1] = enthalpy(FlowComp);
FlowS = FlowComp;
FlowS.T = T2s;
[~,H2s] = enthalpy(FlowS);

H2a = H1+(H2s-H1)/Eff - Q_FlowWallC;
Win = (H2a-H1) + Q_FlowWallC;

if strcmp(string1,'Outlet')
    Out.Flow = FlowOut;
    Out.Work = Win;
    Tags.(block.name).NRPM = nRPM;
    Tags.(block.name).Beta = Beta;
    Tags.(block.name).Flow = FlowOut;
    Tags.(block.name).Power = Win;
    Tags.(block.name).PR = Inlet.Pout/Inlet.Pin;
    Tags.(block.name).Nflow = Nflow;
    Tags.(block.name).MassFlow = MassFlow(FlowOut);
    Tags.(block.name).Temperature = Y(1);
    Tags.(block.name).Eff = Eff;
    Tags.(block.name).nMflow = MassFlow(FlowOut)/block.FlowDesign;
    Tags.(block.name).Pressure = Inlet.Pout;
elseif strcmp(string1,'dY')
    dY = 0*Y;
    [~, Hout] = enthalpy(FlowOut);
    dY(1) = (H2a-Hout)*Ru*Y(1)/(Inlet.Pout*block.Volume*Cp); %dT = dQ / n Cp
    dY(2) = (Q_FlowWallC - Q_WallAmbC - Q_WallAmbR)/(block.Mass*block.SpecHeat);
    Out = dY./block.Scale;
end
