function Out = Combustor(t,Y, Inlet,block,string1)
% a combustor with 2 inlet flows (air and fuel) modeled with 3 reactions
% The completeness of each reaction is adjusted with Rxn1 --3;
% The combustor is divided into a core flow and a bypass flow
% Four (4) inlets: air flow (TXN), fuel flow (TXN), bypass split between core and bypass flow, and outlet pressure
% Three (3) outlets: Bypass flow, core flow and pressure at the inlet
% Six (6) states: Tcombustor, Tbypass, Twall, Tcasing, Inlet Pressure
global Ru Tags
Y = Y.*block.Scale;
Bypass = Y(5);
Pin = Y(6);
NetOut = block.Pfactor*(Pin-Inlet.Pout);%total cold flow out
spec = {'CH4','CO','CO2','H2','H2O'};
for i =1:1:length(spec)
    if ~isfield(Inlet.Fuel,spec{i})
        Inlet.Fuel.(spec{i}) = 0;
    end
end

%merge flows
inlets = {'Air';'Fuel'};
H_in = 0;
block.spec={};
for j = 1:1:2
    H_in = H_in + enthalpy(Inlet.(inlets{j}));
    Spec = fieldnames(Inlet.(inlets{j}));
    Spec = Spec(~strcmp('T',Spec));
    if strcmp(inlets{j},'Air')
        b=Bypass;
    else b= 0;
    end
    for i = 1:1:length(Spec)
        if ismember(Spec{i},block.spec)
            ReactMix.(Spec{i}) = ReactMix.(Spec{i}) + Inlet.(inlets{j}).(Spec{i})*(1-b);
            if strcmp(inlets{j},'Air')
                BypassFlow.(Spec{i}) = Inlet.(inlets{j}).(Spec{i})*b;
            end
        else
            block.spec(end+1) = Spec(i);
            ReactMix.(Spec{i}) = Inlet.(inlets{j}).(Spec{i})*(1-b);
            if strcmp(inlets{j},'Air')
                BypassFlow.(Spec{i}) = Inlet.(inlets{j}).(Spec{i})*b;
            end
        end
    end
end
BypassFlow.T = Inlet.Air.T;

%% 3 reaction:
% CH4 + 1.5O2 --> CO + 2H2O
% CO + .5O2 --> CO2
% H2 + .5O2 --> H2O
R.CH4 = ReactMix.CH4*block.Rxn1;
R.CO = (ReactMix.CO+R.CH4)*block.Rxn2;
R.H2 = ReactMix.H2*block.Rxn3;
sumR = 1.5*R.CH4 + 0.5*R.CO + 0.5*R.H2;
r = fieldnames(R);
phi = sumR/(ReactMix.O2); %rich combustion limited by O2;
EquivError = block.EquivSet - phi;
if phi>1 %rich combustion
    for i = 1:1:length(r)
        R.(r{i}) = R.(r{i})/phi; %rich combustion limited by O2
    end
end

for i = 1:1:length(block.spec)
    if strcmp(block.spec{i},'CH4')
        CombustMix.CH4 = ReactMix.CH4 - R.CH4;
    elseif strcmp(block.spec{i},'CO')
        CombustMix.CO = ReactMix.CO + R.CH4- R.CO;
    elseif strcmp(block.spec{i},'CO2')
        CombustMix.CO2 = ReactMix.CO2 + R.CO;
    elseif strcmp(block.spec{i},'H2')
        CombustMix.H2 = ReactMix.H2  - R.H2;
    elseif strcmp(block.spec{i},'H2O')
        CombustMix.H2O = ReactMix.H2O + 2*R.CH4 + R.H2;
    elseif strcmp(block.spec{i},'O2')
        CombustMix.O2 = ReactMix.O2 - 1.5*R.CH4 - .5*R.CO - .5*R.H2;
    else
        CombustMix.(block.spec{i}) = ReactMix.(block.spec{i});
    end
end
CombustMix.T = Y(1);

FlowFraction = NetOut/(NetFlow(CombustMix)+NetFlow(BypassFlow));%Outflow/Inflow
Out.Main.T = Y(1);
for i = 1:length(block.spec)
    Out.Main.(block.spec{i}) = CombustMix.(block.spec{i})*FlowFraction;%scale outflow by p-factor
end

spec = fieldnames(Inlet.Air);
for i = 1:1:length(spec)
    Out.Bypass.(spec{i}) = BypassFlow.(spec{i})*FlowFraction;
end
Out.Bypass.T = Y(2);

if strcmp(string1,'Outlet')
    Out.Pin = Pin;
    Tags.(block.name).EquivelanceRatio = phi;
    Tags.(block.name).Temperatures = Y(1:4);
    Tags.(block.name).MassFlow = MassFlow(Out.Main)+MassFlow(Out.Bypass);
elseif strcmp(string1,'dY')
    dY = 0*Y;

    Hcombusted = enthalpy(CombustMix);
    Hbypass = enthalpy(BypassFlow);
    
    Cp = SpecHeat(Out.Main);
    Cp2 = SpecHeat(Inlet.Air);
    
    Q_CombWall =(Y(1) - Y(3))*block.AmbConv*block.SurfA/1000;
    Q_BypassCasing = (Y(2) - Y(4))*block.SurfA*block.AmbConv/1000; %convection from bypass to casing
    Q_WallBypassC = (Y(3) - Y(2))*block.SurfA*block.AmbConv/1000;%Convection from Wall to Bypass
    Q_WallBypassR = ((Y(3))^4 - (Y(2))^4)*block.AmbSigma*block.AmbEpsilon*block.SurfA/1000;%Radiation from Wall to Bypass
    Q_CasingAmbR = ((Y(4))^4 - (block.Tamb)^4)*block.AmbSigma*block.AmbEpsilon*block.SurfA/1000; % Heat transfer from casing to ambient due to radiation
    Q_CasingAmbC = (Y(4) - block.Tamb)*block.AmbConv*block.SurfA/1000; % Heat transfer from casing to ambient due to convection
    
    dY(1) = (H_in - Hbypass - Hcombusted - Q_CombWall)*Ru*Y(1)/(Cp*block.Vol*Inlet.Pout); %main air flow through core combustor
    dY(2) = (Hbypass*FlowFraction - enthalpy(Out.Bypass) -Q_BypassCasing + Q_WallBypassR + Q_WallBypassC)*Y(2)*Ru/(Cp2*block.Vol*Inlet.Pout); %total change of bypass temp
    dY(3) =(Q_CombWall - Q_WallBypassR - Q_WallBypassC)/(block.MassWall * block.SpecHeatWall);%total change in wall temp
    dY(4) =(Q_BypassCasing - Q_CasingAmbR - Q_CasingAmbC)/(block.MassCasing * block.SpecHeatCasing); %Total change in Casing temp
    dY(5) = EquivError;
    if (Pin-Inlet.Pout)<0 && NetFlow(ReactMix)<=0
        dY(6) = -(Pin-Inlet.Pout);
    else
        dY(6) = (NetFlow(CombustMix)+NetFlow(BypassFlow)-NetOut)*Ru*Y(1)/(block.Vol);
    end
    Out = dY./block.Scale;
end