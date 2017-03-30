function Out = Oxidizer(t,Y, Inlet,block,string1)
% an oxidizer with n inlet flows 
% Complete combustion is assumed for all fuels (if enough O2 present)
% two or more (2+) inlets: Flows and outletl pressure
% Two (2) outlets: Flow out and pressure at the inlet
% Many (2+n) states: Outlet Temperature, Species, Inlet Pressure
global Ru Tags
Pin = Y(end);
%merge flows
Hin = 0;
NetIn =[];
n = block.inlets;
inlets = fieldnames(Inlet);
for j = 1:1:n
    Hin = Hin + enthalpy(Inlet.(inlets{j}));
    Spec = fieldnames(Inlet.(inlets{j}));
    for i = 1:1:length(block.spec)
        if ~isfield(NetIn,block.spec{i})
            NetIn.(block.spec{i}) = 0;
        end
        if ismember(block.spec{i},Spec)
            NetIn.(block.spec{i}) = NetIn.(block.spec{i}) + Inlet.(inlets{j}).(block.spec{i});
        end
    end
end
ReactMix.T = Y(1);
%% 3 reaction:
% CH4 + 1.5O2 --> CO + 2H2O
% CO + .5O2 --> CO2
% H2 + .5O2 --> H2O
R.CH4 = NetIn.CH4;
R.H2 = NetIn.H2;
if ReactMix.T>1700
    Kp = 0.4589*(ReactMix.T/3000)^14.772;
else
    Kp = max(0,(ReactMix.T-300)*(6.78e-5/1400));
end
X_CO_X_CO2 = Kp/((Pin/101)^.5*(NetIn.O2/NetFlow(NetIn))); %estimate ratio of CO to CO2
R.CO = ((NetIn.CO+R.CH4)-X_CO_X_CO2*NetIn.CO2)/(1+X_CO_X_CO2);

sumR = 1.5*R.CH4 + 0.5*R.CO + 0.5*R.H2;
phi = sumR/NetIn.O2;

scaleFlow = (block.Pfactor*(Pin-Inlet.Pout))/sum(max(0,Y(2:1+length(block.spec))));%total flow out
ActualOut.T = Y(1);
for i = 1:length(block.spec)
    ActualOut.(block.spec{i}) = max(0,Y(1+i))*scaleFlow;
end

if strcmp(string1,'Outlet')
    Out.Flow = ActualOut;
    Out.Pin = Pin;
    Out.MeasureT = ActualOut.T;
    Tags.(block.name).EquivelanceRatio = phi;
    Tags.(block.name).Temperature = Y(1);
    Tags.(block.name).MassFlow = MassFlow(Out.Flow);
elseif strcmp(string1,'dY')
    dY = 0*Y;
    
    r = fieldnames(R);
    if phi>1 %rich combustion
        for i = 1:1:length(r)
            R.(r{i}) = R.(r{i})/phi; %rich combustion limited by O2
        end
    end
    for i = 1:1:length(block.spec)
        if strcmp(block.spec{i},'CH4')
            ReactMix.CH4 = NetIn.CH4 - R.CH4;
        elseif strcmp(block.spec{i},'CO')
            ReactMix.CO = NetIn.CO + R.CH4- R.CO;
        elseif strcmp(block.spec{i},'CO2')
            ReactMix.CO2 = NetIn.CO2 + R.CO;
        elseif strcmp(block.spec{i},'H2')
            ReactMix.H2 = NetIn.H2  - R.H2;
        elseif strcmp(block.spec{i},'H2O')
            ReactMix.H2O = NetIn.H2O + 2*R.CH4 + R.H2;
        elseif strcmp(block.spec{i},'O2')
            ReactMix.O2 = NetIn.O2 - 1.5*R.CH4 - .5*R.CO - .5*R.H2;
        else
            ReactMix.(block.spec{i}) = NetIn.(block.spec{i});
        end
    end
    

    ReactFlow = NetFlow(ReactMix);
    FlowOut = sum(Y(2:1+length(block.spec)));
    FlowError = NetFlow(ActualOut) - FlowOut;
    %temperature
    scaleFlow2 = NetFlow(ReactMix)/FlowOut;
    Hout = enthalpy(ActualOut)*scaleFlow2/scaleFlow; %scale Hin and Hout to the same flow rate (no reactions so easy)
    Cp = SpecHeat(ActualOut);
    dY(1) = (Hin-Hout)/(block.Vol*Cp*Pin).*Y(1)*Ru;
    for i = 1:1:length(block.spec)
%         dY(1+i) = (ReactMix.(block.spec{i})-Y(i+1))*Ru*Y(1)/(block.Vol*Pin);%change in stored mass of each species
        dY(1+i) = ((ReactMix.(block.spec{i})/ReactFlow-Y(i+1)/FlowOut)*FlowOut + FlowError*max(0,Y(i+1))/FlowOut)*Ru*Y(1)/(block.Vol*Pin);%change in stored mass of each species
    end
    dY(2+i) = (ReactFlow - NetFlow(ActualOut))*Ru*Y(1)/(block.Vol); %change in Pressure
    Out = dY;
end