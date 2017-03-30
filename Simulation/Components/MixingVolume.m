function Out = MixingVolume(t,Y, Inlet,block,string1) 
% a simple adiabatic mixing model mixing two or more inlet flows 
% Multiple (n+1) inlets: Flows (consisting of Temperature and flow rates of individual species) and P out 
% Three (3) outlets: Pressure, Flow , and Temperature 
% Many (n+2) states: Temperature, any species under consideration, and pressure
global Ru Tags 
Pin =Y(end);
FlowOut = sum(max(0,Y(2:1+length(block.spec))));
scaleFlow = (block.Pfactor*(Pin-Inlet.Pout))/FlowOut;%total flow out    
ActualOut.T = Y(1);
for i = 1:1:length(block.spec)
    ActualOut.(block.spec{i}) = max(0,Y(i+1)*scaleFlow);
end
if strcmp(string1,'Outlet')
    Out.Pin = Pin;
    Out.Outlet = ActualOut;
    Out.Temperature = Out.Outlet.T;
    Tags.(block.name).MassFlow = MassFlow(ActualOut)/scaleFlow;
    Tags.(block.name).Temperature = Y(1);
elseif strcmp(string1,'dY')
    dY = 0*Y;
    n = block.inlets;
    NetIn = {};
    for j = 1:1:length(block.spec)
        NetIn.(block.spec{j}) = 0;
    end
    inlets = fieldnames(Inlet);
    Hin = 0;
    for i = 1:1:n
        spec = fieldnames(Inlet.(inlets{i}));
        spec = spec(~strcmp('T',spec));
        for j = 1:1:length(spec)
            NetIn.(spec{j}) = NetIn.(spec{j}) + Inlet.(inlets{i}).(spec{j});
        end
        Hin = Hin + enthalpy(Inlet.(inlets{i}));
    end
    FlowIn = NetFlow(NetIn);
    FlowError = NetFlow(ActualOut) - FlowOut; 
    %temperature
    scaleFlow2 = NetFlow(NetIn)/FlowOut ;
    Hout = enthalpy(ActualOut)*scaleFlow2/scaleFlow; %scale Hin and Hout to the same flow rate (no reactions so easy)
    Cp = SpecHeat(ActualOut);
    dY(1) = (Hin-Hout)/(block.Vol*Cp*Pin).*Y(1)*Ru;
    %species
    for i = 1:1:length(block.spec)
%         dY(1+i) = (NetIn.(block.spec{i})-Y(i+1)).*Y(1)*Ru/(block.Vol*Pin);%change in stored mass of each species
        dY(1+i) = ((NetIn.(block.spec{i})/FlowIn-Y(i+1)/FlowOut)*FlowOut + FlowError*max(0,Y(i+1))/FlowOut)/(block.Vol*Pin).*Y(1)*Ru;%change in stored mass of each species
    end
    %pressure
    dY(end) = (FlowIn - NetFlow(ActualOut))*Ru*Y(1)/(block.Vol);
    Out = dY;
end

