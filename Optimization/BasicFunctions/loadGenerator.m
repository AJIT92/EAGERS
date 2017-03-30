function loadGenerator% Loads generators for economic dispatch
%% this function identifies the values that will be used to represent each generator in the quadratic optimizations
global Plant UB LB dX_dt MinThresh 
%UB: the upper limit (capacity) of each generator
%LB: the lower limit (capacity) of each generator when on
%selfDisch: the amount of self dicharging per hour for a storage system
%MinThresh: a minimum buying constraint which may exist for some utilities

nG = length(Plant.Generator);
MinThresh = [];
LB = zeros(1,nG);
UB = zeros(1,nG);
dX_dt = zeros(1,nG);
Plant.optimoptions.Outputs = {};
for i = 1:1:nG
    Plant.Generator(i).OpMatA = {}; %delete this line when you begin using the gui again
    if isempty(Plant.Generator(i).OpMatA)%only load generators that have not been loaded yet. New run, new generator, or edited generator
        typeNoSpace = char(Plant.Generator(i).Type(~isspace(char(Plant.Generator(i).Type))));
        [Plant.Generator(i).OpMatA, Plant.Generator(i).OpMatB, LB(i), UB(i), dX_dt(i),SS] = eval(strcat('load',typeNoSpace,'(Plant.Generator(i))'));
        if ~isempty(SS)
            SSi(i) = SS;
        end
        if strcmp(typeNoSpace,'Utility') && ~isempty(Plant.Generator(i).OpMatB.output) && Plant.Generator(i).OpMatB.X.lb>0
            MinThresh(i) = Plant.Generator(i).OpMatB.X.lb;
        end
    end
    if ~isempty(Plant.Generator(i).OpMatA.output)
        Outs = fieldnames(Plant.Generator(i).OpMatA.output);
        for j = 1:1:length(Outs)
            if nnz(strcmp(Outs{j}, Plant.optimoptions.Outputs))==0 %not in the current list of outputs
                Plant.optimoptions.Outputs(end+1) = Outs(j);
            end
        end
    end
end 
%need to wait until everything else has been loaded to load the buffer,
%because it is reliant on how quickly everything else can charge the system
findbuffer

agregateSSmodel(SSi)


function [OpMatA, OpMatB, LB, UB, dX_dt,SSi] = loadUtility(Gen)
% since utilites are unchanged between the two optimizations OpMatB = OpMatA
util = Gen.VariableStruct;
UB = inf;
dX_dt = inf;
SSi =[];
if strcmp(Gen.Source, 'Electricity')
    OpMatA.states = {'X'};
    OpMatA.output.E = 1;
    OpMatA.X.H = 0;
    OpMatA.X.f = 1;
    OpMatA.X.lb = -inf;
    OpMatA.X.ub = inf;
    OpMatA.xL = 'nS';%length of states
    OpMatA.req = '0';%# of rows taken in Ax = b
    OpMatA.r = '0';%# of rows taken in Ax <= b
    LB = -inf;
    if ~isfield(util,'SellBack') || max(util.SellBack)==0 % no sellback allowed (only 1 state)
        OpMatA.X.lb = util.MinImportThresh;
        LB = util.MinImportThresh;
    elseif util.SellBack == -1 || util.SellBack==1% reversed meter (only 1 state)
        %default (no changes from above)
    else %seperate purchasing and selling states
        OpMatA.states = {'X';'Y'};
        OpMatA.Y.f = -max(util.SellBack,1-1e-6);%ensure less than 1, so no issues with pass through power
        OpMatA.Y.H = 0;
        OpMatA.Y.lb = 0;
        OpMatA.Y.ub = inf;
        OpMatA.xL = '2*nS';
    end
else %if it is a fuel utility
    OpMatA = [];
    OpMatA.output = [];
    OpMatA.states = [];
    OpMatA.xL = '0';
    OpMatA.req = '0';
    OpMatA.r = '0';
    LB = 0;
end
OpMatB = OpMatA;


function [OpMatA, OpMatB, LB, UB, dX_dt,SSi] = loadDistrictHeating(Gen)
%this function loads the parameters for a distric heating supply. It is
%very similar to the electric utility, but there is no sellback capability
%and the lower bound is dependent on the CHP generators and heaters in the
%system (but their individual constraints will boud ho much can be put back
%on the district heating grid.
UB = inf;
LB = 0; %can buy or sell to district heating
dX_dt = inf;
SSi =[];
OpMatA.states = {'X'};
OpMatA.output.H = 1;
OpMatA.X.H = 0; %no quadratic term for the cost here
OpMatA.X.f = 2; %there is a linear term for the cost
OpMatA.X.lb = 0; %LB = -sum(UB(find(GenType==2)).*Hratio(find(GenType==2)))-sum(UB(find(GenType==4)));
OpMatA.X.ub = inf;
OpMatA.xL = 'nS';
OpMatA.req = '0';
OpMatA.r = '0';
OpMatB = OpMatA;

function [OpMatA, OpMatB, LB, UB, dX_dt, SSi] = loadDistrictCooling(Gen)
% this function loads the parameters for a district cooling supply
% it is very similar to the electric utility, but there is no sell back
% capability and the lower bound is dependent on the chillers in the system
%nearly same as district heating/utility with sellback
OpMatA.states = {'X'};
OpMatA.req = '0';%# of rows taken in Ax = b
OpMatA.r = '0';%# of rows taken in Ax <= b
OpMatA.xL = 'nS';%# of states
OpMatA.output.C = 1;
% OpMatA.Z.f = 2;
% OpMatA.Z.H = 0;
% OpMatA.Z.lb = 0;
% OpMatA.Z.ub = inf;
% OpMatA.Y.f = -6;%penalize putting power back on the grid
% OpMatA.Y.H = 0;
% OpMatA.Y.lb = -inf;
% OpMatA.Y.ub = 0;
% OpMatA.link.eq = [1, -1, -1]; %X = Y+Z
% OpMatA.link.beq = 0;
OpMatA.X.f = 2;
OpMatA.X.H = 0;
OpMatA.X.lb = 0;
OpMatA.X.ub = inf;
LB = 0;
UB = inf;
dX_dt = inf;
SSi = [];
OpMatB = OpMatA;


function [OpMatA, OpMatB, LB, UB, dX_dt,SSi] = loadElectricGenerator(Gen)
% this function loads the parameters for an electric generator generators
[OpMatA, OpMatB, LB, UB, dX_dt,SSi] = loadCHPGenerator(Gen);

function [OpMatA, OpMatB, LB, UB, dX_dt,SSi] = loadCHPGenerator(Gen)
% this function loads the parameters for a combined heat and power
% generator, regular electric generator, or chiller
global Plant
UB = Gen.Size;
if isfield(Gen.Output,'Cooling')&& Gen.Output.Cooling(end)>0
    LB = Gen.VariableStruct.Startup.Cooling(1); %chiller
    costTerms = GenCosts(Gen,LB,UB,'Cooling');
    if Plant.optimoptions.sequential == 0
        Cratio = Gen.Output.Cooling(end);
        OpMatA.output.C = 1;
        OpMatB.output.C = 1;
        OpMatA.output.E = -1/Cratio;
        OpMatB.output.E = -1/Cratio;
        costTerms.P = UB;
        costTerms.I = UB;
        costTerms.Intercept(1) =0;
        costTerms.Convex(1) = 0;
    else
        OpMatA.output.C = 1;
        OpMatB.output.C = 1;
    end
else
    LB = Gen.VariableStruct.Startup.Electricity(end); %electric or CHP generator 
    costTerms = GenCosts(Gen,LB,UB,'Electricity');
    OpMatA.output.E = 1;
    OpMatB.output.E = 1;
end
if isfield(Gen.Output,'Heat')&& Gen.Output.Heat(end)>0
    Hratio = Gen.Output.Heat(end)/Gen.Output.Electricity(end);
    OpMatA.output.H = Hratio;
    OpMatB.output.H = Hratio;
else Hratio = [];
end
[dX_dt,SSi] = RampRateCalc(Gen.VariableStruct.StateSpace,LB,UB,Hratio);
if costTerms.P == UB %linear fit use 1 state
    OpMatA.states = {'X'};%x
    OpMatA.xL = '1*nS+1';% 1 states + IC
    OpMatA.req = '1'; %one for the initial condition
    OpMatA.r = '2*nS'; %ramp inequality
    OpMatA.X.H = 0;
    OpMatA.X.f = costTerms.Convex(1);
    OpMatA.X.lb = 0;%relax lower bound until have figured out what will be off or on
    OpMatA.X.ub = UB;
    OpMatA.Ramp.A = [-1, 1; 1, -1;]; %-output1+output2=ramp up %output1-output2=-rampdown
    OpMatA.Ramp.b = [dX_dt;dX_dt];
else %piecewise quadratic fit
    OpMatA.states = {'Y';'Z';};%beta, gamma
    OpMatA.xL = '2*nS+1';% 2 states + IC
    OpMatA.req = '1'; % one equality for the initial condition
    OpMatA.r = '2*nS'; %ramp inequality

    OpMatA.Ramp.A = [-1, 1; 1, -1;]; %-output1+output2=ramp up %output1-output2=-rampdown
    OpMatA.Ramp.b = [dX_dt;dX_dt];

    OpMatA.Y.H = 0;
    OpMatA.Y.f = costTerms.Convex(1);
    OpMatA.Y.lb = 0;%relax lower bound until have figured out what will be off or on
    OpMatA.Y.ub = costTerms.P;

    OpMatA.Z.H = costTerms.Convex(3);
    OpMatA.Z.f = costTerms.Convex(2);
    OpMatA.Z.lb = 0;
    OpMatA.Z.ub = UB-costTerms.P;
end

if costTerms.I == UB %linear fit use 1 state
    OpMatB.constCost = costTerms.Intercept(4);
    OpMatB.states = {'X'};%x
    OpMatB.xL = '1*nS+1';% 1 states + IC
    OpMatB.req = '1'; %one for the initial condition
    OpMatB.r = '2*nS'; %ramp inequality
    OpMatB.X.H = 0;
    OpMatB.X.f = costTerms.Intercept(1);
    OpMatB.X.lb = LB;%relax lower bound until have figured out what will be off or on
    OpMatB.X.ub = UB;
    OpMatB.Ramp.A = [-1, 1; 1, -1;]; %-output1+output2=ramp up %output1-output2=-rampdown
    OpMatB.Ramp.b = [dX_dt;dX_dt];   
else
    OpMatB.constCost = costTerms.Intercept(4);
    OpMatB.states = {'Y';'Z'};%x, beta, gamma
    OpMatB.xL = '2*nS+1';% 2 states + IC
    OpMatB.req = '1'; %link equality, one per timestep, and one for the initial condition
    OpMatB.r = '2*nS'; %ramp inequality

    OpMatB.Ramp.A = [-1, 1; 1, -1;]; %-output1+output2=ramp up %output1-output2=-rampdown
    OpMatB.Ramp.b = [dX_dt;dX_dt];

    OpMatB.Y.H = 0;
    OpMatB.Y.f = costTerms.Intercept(1);
    OpMatB.Y.lb = LB;
    OpMatB.Y.ub = costTerms.I;

    OpMatB.Z.H = costTerms.Intercept(3);
    OpMatB.Z.f = costTerms.Intercept(2);
    OpMatB.Z.lb = 0;
    OpMatB.Z.ub = UB-costTerms.I;  
end

function [OpMatA, OpMatB, LB, UB, dX_dt,SSi] = loadChiller(Gen)
% this function loads the parameters for a chiller generators
% it is very similar to the way Electric Generators are loaded
[OpMatA, OpMatB, LB, UB, dX_dt,SSi] = loadCHPGenerator(Gen);


function [OpMatA, OpMatB, LB, UB, dX_dt,SSi] = loadHeater(Gen)
% this function loads the parameters for a heater generator.
LB = Gen.VariableStruct.Startup.Heat(end);
UB = Gen.Size;
[dX_dt, SSi] = RampRateCalc(Gen.VariableStruct.StateSpace, LB, UB, {});
capacity = Gen.Output.Capacity*UB;
efficiency = Gen.Output.Heat;
costLinear = 1/mean(efficiency(capacity>=LB));%efficiency term  put in cost (will later be scaled by utility cost.

OpMatA.states = {'X'};
OpMatA.output.H = 1;
OpMatA.xL = 'nS+1';% 1 state + IC
OpMatA.req = '1'; 
OpMatA.r = '2*nS'; %ramp inequality

OpMatA.X.H = 0;
OpMatA.X.f = costLinear;
OpMatA.X.lb = 0;
OpMatA.X.ub = UB;
OpMatA.Ramp.A = [-1, 1; 1, -1;];
OpMatA.Ramp.b = [dX_dt; dX_dt];

OpMatB = OpMatA;
OpMatB.X.lb = LB;

function [OpMatA, OpMatB, LB, UB, dX_dt,SSi] = loadElectricStorage(Gen)
%this function loads all the parameters needed for electric storage
[OpMatA, OpMatB, LB, UB, dX_dt,SSi] = Storage(Gen);

function [OpMatA, OpMatB, LB, UB, dX_dt,SSi] = loadThermalStorage(Gen)
%loads all types of thermal storage
[OpMatA, OpMatB, LB, UB, dX_dt,SSi] = Storage(Gen);

function [OpMatA, OpMatB, LB, UB, dX_dt,SSi] = Storage(Gen)
%this function just directs to either hot or cold thermal storage
%if we can get rid of the CS, HS structures then this function can be
%eliminated, because all storage can be handled the same way.
global scaleTime
SSi =[];
Stor.Size = Gen.Size*scaleTime;
Stor.SelfDischarge  = Gen.VariableStruct.SelfDischarge;% SelfDischarge per hour (fraction of total charge)
if isfield(Gen.VariableStruct, 'EnStoreType') %if its thermal storage
    Stor.PeakDisch = (Gen.VariableStruct.DischRatePerc/100*Gen.Size); %Thermal kW out
    Stor.PeakCharge = (Gen.VariableStruct.FillRatePerc/100*Gen.Size); %Thermal kW in
    Stor.ChargeEff = Gen.VariableStruct.ChargeEff;
    Stor.DischEff = Gen.VariableStruct.DischargeEff;
    Stor.UsableSize  = Stor.Size*Stor.DischEff; % usable size 
    if strcmp(Gen.VariableStruct.EnStoreType, 'ColdTES')
        OpMatA.output.C = []; 
    elseif strcmp(Gen.VariableStruct.EnStoreType, 'HotTES')
        OpMatA.output.H = [];
        Stor.ChargeEff = 1; %set to ideal 
        Stor.DischEff = 1;
    elseif strcmp(Gen.VariableStruct.EnStoreType, 'HVAC')
        if strcmp(Gen.VariableStruct.EnStoreType, 'Heat only')
            OpMatA.output.H = [];
        elseif strcmp(Gen.VariableStruct.EnStoreType, 'AC only')
            OpMatA.output.C = [];
        else
            OpMatA.output.C = []; 
            OpMatA.output.H = [];
        end
    end
    
else %electric battery
    OpMatA.output.E = []; 
    Stor.Voltage = Gen.VariableStruct.Voltage;
    DischCurrent = Gen.VariableStruct.PeakDisch.*Gen.Size/Stor.Voltage*1000;
    Stor.DischResistScaled = (100/DischCurrent)*Gen.VariableStruct.DischResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
    DischVoltLoss = DischCurrent*Stor.DischResistScaled; %keep in mind when calculating loss as function of discharge current
    ChargeCurrent = Gen.VariableStruct.PeakCharge*Gen.Size/Stor.Voltage*1000;
    Stor.ChargeResistScaled = (100/ChargeCurrent)*Gen.VariableStruct.ChargeResist*(1/1000); %Scale so the loss of power is equivelant to that specified at 100Amps
    ChargeVoltLoss = ChargeCurrent*Stor.ChargeResistScaled;
    Stor.ChargeEff = Stor.Voltage/(Stor.Voltage+ChargeVoltLoss);
    Stor.DischEff = (Stor.Voltage-DischVoltLoss)/Stor.Voltage;
    Stor.PeakDisch = (DischCurrent*Stor.Voltage*Stor.DischEff/1000); %PeakDischPower kW out
    Stor.PeakCharge = (ChargeCurrent*Stor.Voltage/Stor.ChargeEff/1000); % PeakChargePower kW in
    Stor.UsableSize  = Stor.Size*(Gen.VariableStruct.MaxDOD/100)*Stor.DischEff; % usable size (state x is 0_> usable size, must calculate this initial condition each time from voltage & SOC
    Stor.MaxDODcapacity = Gen.Size*(1-Gen.VariableStruct.MaxDOD/100);
end

OpMatA.Stor = Stor;
Outs = fieldnames(OpMatA.output);
LB = 0;
UB = Stor.UsableSize;
dX_dt = Stor.PeakDisch;% storage discharge constraint in kW
if length(Outs)==2 %HVAC system
    OpMatA.states = {'X'; 'Y'}; %state of charge, charging power, no buffers
    OpMatA.link.ineq = [1 -1/(1-(HS.ChargeEff*HS.DischEff))]; % -SOC(t-1) + SOC(t) -charging <0 ------ Charging is the 1/inefficiency + 1
    OpMatA.link.bineq = 0;
    OpMatA.xL = '2*nS+1';% 2 states + IC
    OpMatA.req = '1';%one line for the IC
    OpMatA.r = '3*nS'; %ramp inequality & link inequality
    
    OpMatA.X.lb = 0;
    OpMatA.X.ub = Stor.UsableSize;
    OpMatA.X.H = 0;
    OpMatA.X.f = 0;
    OpMatA.Ramp.A = [-1,1;1,-1]; % SOC2-SOC1<peakcharge; SOC1-SOC2<peakdischarge
    OpMatA.Ramp.b = [Gen.VariableStruct.FillRatePerc/100*Gen.Size; inf];%since the stored energy is applied directly to the demand site, the max discharge rate is infinite.

    OpMatA.Y.lb = 0;
    OpMatA.Y.ub = inf;%the limit on how much charging power can be delivered is handled by the generators' limits, so put inf here to prevent redundancy
    OpMatA.Y.H = 0;%the cost of the charging power is handled by the generators 
    OpMatA.Y.f = 0;
elseif Stor.ChargeEff*Stor.DischEff>=1 %ideal storage, ignore charging state
    OpMatA.states = {'X';'Z';'W'};%SOC(t+1), charging power, upper buffer, lower buffer
    OpMatA.link.ineq = [-1 0 -1; 1 -1 0];%-SOC(t)-lowerbuffer<-0.2  and %SOC-upperbuffer<0.8
    OpMatA.link.bineq = [0;0;]; %% note: the magnitude of the buffer is set later in findBuffer
    OpMatA.xL = '3*nS+1';% 3 states + IC
    OpMatA.req = '1';%1 for the IC 
    OpMatA.r = '4*nS'; %ramp inequality & link inequalities
    OpMatA.Stor.ChargeEff = 1; %eliminate chargin inefficiencies in CHP plants
    OpMatA.Stor.DischEff = 1;
    
    OpMatA.X.lb = 0;
    OpMatA.X.ub = UB;
    OpMatA.X.H = 0; %no cost directly associated with storage
    OpMatA.X.f = 0;
    OpMatA.Ramp.A = [-1, 1; 1, -1];%SOC2-SOC1<PeakCharge; SOC1-SOC2<PeakDisch
    OpMatA.Ramp.b = [Stor.PeakCharge; Stor.PeakDisch];

    OpMatA.Z.lb = 0;
    OpMatA.Z.ub = 0;%% note: the magnitude of the buffer is set later in findBuffer
    OpMatA.Z.H = 0;
    OpMatA.Z.f = 0;

    OpMatA.W.lb = 0;
    OpMatA.W.ub = 0;%% note: the magnitude of the buffer is set later in findBuffer
    OpMatA.W.H = 0;
    OpMatA.W.f = 0;
else % include charging state
    OpMatA.states = {'X';'Y';'Z';'W'};%SOC(t+1), charging power, upper buffer, lower buffer
    OpMatA.link.ineq = [1, -1/(1-(Stor.ChargeEff*Stor.DischEff)), 0, 0];%-SOC(t-1) + SOC(t) - chargingpower < 0 ------ Charging is the 1/inefficiency + 1
    OpMatA.link.ineq = [OpMatA.link.ineq; -1 0 0 -1];%-SOC(t)-lowerbuffer<-0.2
    OpMatA.link.ineq = [OpMatA.link.ineq; 1 0 -1 0];%SOC-upperbuffer<0.8
    OpMatA.link.bineq = [0;0;0;];
    OpMatA.xL = '4*nS+1';% 4 states + IC
    OpMatA.req = '1';%1 for the IC 
    OpMatA.r = '5*nS'; %ramp inequality & link inequalities
    
    OpMatA.X.lb = 0;
    OpMatA.X.ub = UB;
    OpMatA.X.H = 0; %no cost directly associated with storage
    OpMatA.X.f = 0;
    OpMatA.Ramp.A = [-1, 1; 1, -1];%SOC2-SOC1<PeakCharge; SOC1-SOC2<PeakDisch
    OpMatA.Ramp.b = [Stor.PeakCharge; Stor.PeakDisch];

    OpMatA.Y.lb = 0; %the lower bound of the charging state is 0 (only shows up as load when charging
    OpMatA.Y.ub = Stor.PeakCharge*(1/(Stor.ChargeEff*Stor.DischEff)-1); %the upper bound of the charging power is limited by the generators themselves. put inf here to prevent redundant limits.
    OpMatA.Y.H = 0;
    OpMatA.Y.f = 0;

    OpMatA.Z.lb = 0;
    OpMatA.Z.ub = 0;
    OpMatA.Z.H = 0;
    OpMatA.Z.f = 0;

    OpMatA.W.lb = 0;
    OpMatA.W.ub = 0;
    OpMatA.W.H = 0;
    OpMatA.W.f = 0;
end
OpMatB = OpMatA; %the bounds do not change for the second optimization

function [OpMatA, OpMatB, LB, UB, dX_dt,SSi] = loadSolar(Gen)
OpMatA = [];
OpMatA.output = [];%there are no states or outputs for solar because renewable outputs are handled on the demand side
OpMatA.states = [];
OpMatA.xL = '0';
OpMatA.req = '0';
OpMatA.r = '0';
OpMatB = OpMatA;
LB = 0;
UB = inf;
dX_dt = inf;
SSi = [];