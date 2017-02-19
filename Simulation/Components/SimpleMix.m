function Out = SimpleMix(t,Y, Inlet,block,string1) 
% a simple adiabatic mixing model mixing two or more inlet flows 
% Multiple (n+1) inlets: Flows (consisting of Temperature and flow rates of individual species) and P out 
% Three (3) outlets: Pressure, Flow , and Temperature 
% Two (2) states: Temperature and pressure
global Ru Tags
Y = Y.*block.Scale;
  
Pin =Y(end);
NetOut.T = Y(1);
inlets = fieldnames(Inlet);
n = block.inlets;
H = zeros(n,1);
spec ={};
for i = 1:1:n
    spec2 = fieldnames(Inlet.(inlets{i}));
    spec2 = spec2(~strcmp('T',spec2));
    for j = 1:1:length(spec2)
        if ~ismember(spec2{j},spec)
            spec(end+1) = spec2(j);
            NetIn.(spec2{j}) = Inlet.(inlets{i}).(spec2{j});
        else
            NetIn.(spec2{j}) = NetIn.(spec2{j}) + Inlet.(inlets{i}).(spec2{j});
        end
    end
    H(i) = enthalpy(Inlet.(inlets{i}));
end
   
scaleFlow = (block.Pfactor*(Pin-Inlet.Pout))/NetFlow(NetIn);%total cold flow out
for i = 1:1:length(spec)
    NetOut.(spec{i}) = NetIn.(spec{i})*scaleFlow;
end

if strcmp(string1,'Outlet')
    Out.Pin = Pin;
    Out.Outlet = NetOut;
    Out.Temperature = Out.Outlet.T;
elseif strcmp(string1,'dY')
    dY = 0*Y;
    Hin = sum(H);  
    %temperature
    scaleFlow2 = NetFlow(NetIn)/NetFlow(NetOut);
    Hout = enthalpy(NetOut)*scaleFlow2; %scale Hin and Hout to the same flow rate (no reactions so easy)
    Cp = SpecHeat(NetOut);
    dY(1) = (Hin-Hout)/(block.Vol*Cp*Pin).*NetOut.T*Ru;
    %pressure
    dY(end) = (NetFlow(NetIn) - NetFlow(NetOut))*Ru*Y(1)/(block.Vol);
    Out = dY./block.Scale;
    Tags.(block.name).MassFlow = MassFlow(NetOut);
    Tags.(block.name).Temperature = NetOut.T;
end

