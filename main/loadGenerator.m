function loadGenerator% Loads generators for economic dispatch
%% this function identifies the values that will be used to represent each generator in the quadratic optimizations
global Plant UB LB dischEff MinThresh chargeEff selfDisch SSmpc GenSSindex GendX

%UB: the upper limit (capacity) of each generator
%LB: the lower limit (capacity) of each generator when on
%dischEff: Discharge efficiency of energy storage
%chargeEff: charge efficiency of energy storage
%selfDisch: the amount of self dicharging per hour for a storage system
%MinThresh: a minimum buying constraint which may exist for some utilities
%SSmpc: the agregated state-space model seen by the MPC, the time may be scaled by both Tmpc and scaletime.
%GenSSindex: A matrix used to relate the primary and secondary states of each generator to their index in the states of the agregated state space model SSmpc.
%GendX: The derivative of the generator state, assume steady initial condition and GendX=0. Need to keep track of this for MPC

nG = length(Plant.Generator);
dischEff = [];
chargeEff = [];
MinThresh = [];
LB = zeros(1,nG);
UB = zeros(1,nG);
selfDisch = zeros(1,nG);
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
    if isfield(Plant.Generator(i).OpMatA,'Stor')
        dischEff(i) = Plant.Generator(i).OpMatA.Stor.DischEff;
        chargeEff(i) = Plant.Generator(i).OpMatA.Stor.ChargeEff;
        selfDisch(i) = Plant.Generator(i).OpMatA.Stor.SelfDischarge;
    end
end        
findbuffer
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
Plant.OpMatA = buildMatrices('A',Time); %build quadratic programming matrices for FitA
Plant.OpMatB = buildMatrices('B',Time);%build quadratic programming matrices for FitB
dt = Time - [0, Time(1:end-1)];
n = Plant.optimoptions.thresholdSteps;
Time2 = linspace(1,n,n)*dt(1)/n;
Plant.Threshold = buildMatrices('B',Time2);%build quadratic programming matrices used for threshold calculation with linear spacing from t =0 to t = dt(1).
Plant.OneStep = buildMatrices1Step(dt);%build quadratic programming matrices for 1Step at interval spacing of dt

A.Horizon = Time(1);%Plant.optimoptions.Resolution;%the horizon is not the resolution converted to seconds
A.Resolution = Plant.optimoptions.Topt/3600;%the resolution is the frequency of Topt
A.tspacing = 'constant';
Time3 = buildTimeVector(A);%% set up dt vector of time interval length
Plant.Online = [];
Plant.Online.QP = [];
Plant.Online.Organize = [];
Plant.Online.Timestamp = [];
for t = 1:1:length(Time3)
    Plant.Online(t) = buildMatrices('B',Time3(t:end)); %build the matrix for the onlineOptimLoop using FitB
end

%% Build agregate state space model:
GenSSindex = [];
include = {'CHP Generator', 'Electric Generator','Chiller'};
gen = [];
nCHP = 0;
for i=1:1:nG
    if ismember(cellstr(Plant.Generator(i).Type),include)
        gen(end+1) = i;
        Outs = fieldnames(Plant.Generator(i).OpMatB.output);
        if ismember('H',Outs) && ismember('E',Outs)
            nCHP = nCHP+1;
        end
    end
end
GenSSindex.Primary = zeros(length(gen),2);
GenSSindex.Secondary = zeros(nCHP,2);
A = [];
B = [];
C = [];
r = 0;
w = 0;
chp=0;
for i = 1:1:length(gen)
    sA = size(SSi(gen(i)).A);
    if sA(1) ==2 %SISO 2nd order response (1 output)
        A = [A zeros(w,2); zeros(2,w) SSi(gen(i)).A];
        B = [B zeros(w,1); zeros(2,r) SSi(gen(i)).B];
        C = [C zeros(w/2,2); zeros(1,w) SSi(gen(i)).C];
        GenSSindex.Primary(i,1) = w+1; %1st column gives index in state vector x
        GenSSindex.Primary(i,2) = w/2+1; %2nd column gives index in output (y = Cx)
        w = w+2; %total states of X
    elseif sA(1) ==4 %SIMO 2nd order response (2 outputs)
        chp = chp+1;
        A = [A zeros(w,4); zeros(4,w) SSi(gen(i)).A];
        B = [B zeros(w,1); zeros(4,r) SSi(gen(i)).B];
        C = [C zeros(w/2,4); zeros(2,w) SSi(gen(i)).C];
        GenSSindex.Primary(i,1) = w+1; %1st column gives index in state vector x
        GenSSindex.Primary(i,2) = w/2+1; %2nd column gives index in output (y = Cx)
        GenSSindex.Secondary(chp,1) = w+3;  %1st column gives index in state vector x
        GenSSindex.Secondary(chp,2) = w/2+2; %2nd column gives index in output (y = Cx)
        w = w+4;%total states of X
    end 
    r = r+1; %input #
end
SS.A = A;
SS.B = B;
SS.C = C;
SS.D = zeros(length(C(:,1)),length(B(1,:)));
SS.Dt = 1; %sampling time
SSmpc = changeTimestep(SS,Plant.optimoptions.Tmpc,1);
GendX = zeros(1,nG);


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
    OpMatA.link.ineq = [1 -1/(1-(HS.ChargeEff*HS.DischEff))]; %SOC(t-1) - SOC(t) -charging <0 ------ Charging is the 1/inefficiency + 1
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
    OpMatA.link.bineq = [0;0;];
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
    OpMatA.Z.ub = 0;
    OpMatA.Z.H = 0;
    OpMatA.Z.f = 0;

    OpMatA.W.lb = 0;
    OpMatA.W.ub = 0;
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


function [dX_dt, SS_1] = RampRateCalc(SS,LB,UB,Hratio)
global scaleTime
SS = SS(1,:);
if isfield(SS,'Dt')
    Dt = SS.Dt;
else Dt =1;
    SS.Dt = 1;
end
if Dt~=1 %convert to 1 second sampling time
    convSS = ss(SS.A,SS.B,SS.C,SS.D,Dt);
    newSS = d2d(convSS,1);
    r = length(newSS);
    SS_1.A = newSS(1).A;
    SS_1.B = newSS(1).B;
    SS_1.C = newSS(1).C;
    SS_1.D = newSS(1).D;
    for k = 2:1:r
        SS_1.C(end+1,:) = newSS(k).C;
        SS_1.D(end+1,:) = newSS(k).D;
    end
else SS_1 = SS;
end
x0 = LB;
if ~isempty(Hratio)
    x0(2) = LB*Hratio; %CHP heat produed per unit electricity
end
nS = round(4*3600/Dt)+1; % assume ramping is less than 4 hours (i forget why I made this limit)
t = linspace(0, Dt*(nS-1),nS);
u = UB*linspace(1,1,nS)';
[z,z2] = size(SS.C);
X0 = zeros(z2,1);
for k = 1:1:z
    X0(find(SS.C(k,:),1,'first'))=x0(k);%finds the first non-zero element in SS.C and makes X0 at that index = x0
end
SS = ss(SS.A,SS.B,SS.C,SS.D,Dt);
[y,t] = lsim(SS,u,t,X0);
tRise = t(find(y(:,1)>(.95*u(1)),1,'first'))/3600; %rise time in hours
if isempty(tRise)
    tRise = 4;
end
dX_dt = (UB.*(0.95)./tRise)./scaleTime;

function costTerms = GenCosts(Gen,LB,UB,out)
capacity = Gen.Output.Capacity*UB;
efficiency = Gen.Output.(out);
operationRange = find(capacity>=LB);
x = capacity(operationRange);
options = optimset('Algorithm','active-set','LargeScale','off','MaxIter',50,'Display','none');
c = x./efficiency(operationRange); %cost of generator in terms of input
Index = find(efficiency(operationRange)==max(efficiency(operationRange)),1,'last');

%% FIT A: piecewise convex curve with zero y-intercept (may be same as linear)
costTerms.P = x(Index);
costTerms.Convex(1) = c(Index)/x(Index); %linear segment
if Index < length(c)%if index is equal or greater than c, then there is no convex portion
    alpha = x(Index:end)-x(Index);%this is (x-P). In the paper alpha is gamma
    C=[alpha .5*alpha.^2];
    c_alpha = c(Index:end)-c(Index);
    costTerms.Convex(2:3) = lsqlin(C,c_alpha,[-1 0;0 -1],[-costTerms.Convex(1);0],[C(end,1), C(end,2)],c_alpha(end),[],[],[],options);%%Quadratic fit to cost (made for alpha)
else %if there is no convex portion (Index>= length(c))
    costTerms.Convex(2) = 0;% no quadratic term
    costTerms.Convex(3) = 0; %no quadratic term
end
costTerms.Convex(3) = max(0,costTerms.Convex(3));%ensure a positive H value

%% Fit B: piecewise convex with non-zero y-intercept, first find point I beyond which cost curve is convex
dc_dxi = (c(2:end)-c(1:end-1))./(x(2:end)-x(1:end-1));
convex = and((dc_dxi(2:end)>1.001*(dc_dxi(1:end-1))),(dc_dxi(1:end-1)>0));
k = length(convex);
P_i = find(convex==1,1,'first');
if ~isempty(P_i)%possibly a convex section
    MostlyConvex = 0*convex+1;
    for j = 1:1:k-4
        MostlyConvex(j) = (sum(convex(j:j+4))>3);
    end
    while P_i<k && MostlyConvex(P_i)==0
        P_i = P_i+1;
    end
end
if isempty(P_i) || P_i==k%never convex, make a linear fit
    costTerms.I = x(end);
    costTerms.Intercept = (c(end)-c(1))/(x(end)-x(1)); %beta extends to UB, fit still has a constant term
    if costTerms.Intercept<=0
        costTerms.Intercept = c(end)/x(end);
    end
    costTerms.Intercept(2:3) = 0;
else
    costTerms.I = x(P_i+1);
    alpha = x(P_i+1:end)-x(P_i+1);
    c_alpha = c(P_i+1:end)-c(P_i+1);%cost associated with alpha segment
%     c_alpha = zeros(length(alpha),1);
%     for j = 1:1:length(alpha)
%         c_alpha(j) = sum(c(P_i+1:P_i+j));
%     end
    C =[alpha .5*alpha.^2];
    maxSlopeBeta = c_alpha(end)/alpha(end);
    costTerms.Intercept = min(maxSlopeBeta,(c(P_i+1)-c(1))/(x(P_i+1)-x(1))); %slope of beta segment
    costTerms.Intercept(2:3) = lsqlin(C,c_alpha,[0 -1;-1 0;],[0;-costTerms.Intercept(1)],[C(end,1), C(end,2)],c_alpha(end),[],[],[],options);%%Quadratic fit to cost (made for alpha)
    costTerms.Intercept(3) = max(0,costTerms.Intercept(3));%ensure positive H value
end
costTerms.Intercept(4) = c(1)-x(1)*costTerms.Intercept(1); % constant term (y-intercept)

%% plot for double check
% figure(1)
% plot([0; x],[0;c],'g')
% hold on
% 
% QT = costTerms.Convex;
% beta = [linspace(0,costTerms.P) linspace(costTerms.P,costTerms.P)]';
% alpha = [linspace(0,0) linspace(0,x(end)-costTerms.P)]';
% convexFit =  QT(1)*beta + QT(2)*alpha + .5*QT(3)*alpha.^2;
% plot([0; beta+alpha],[0;convexFit],'r')
% 
% QT = costTerms.Intercept;
% beta = [linspace(0,costTerms.I) linspace(costTerms.I,costTerms.I)]';
% alpha = [linspace(0,0) linspace(0,x(end)-costTerms.I)]';
% betterFit =  QT(1)*beta + QT(2)*alpha + .5*QT(3)*alpha.^2 + QT(4);
% plot([0; beta+alpha],[0;betterFit],'b')


%need to wait until everything else has been loaded to load the buffer,
%because it is reliant on how quickly everything else can charge the system
function findbuffer
%storage needs to be loaded last so that the buffers can be calculated
global Plant UB dischEff
nG = length(Plant.Generator);
BuffPerc = Plant.optimoptionsBuffer; % percentage for buffer on storage
for j = 1:1:nG
    if isfield(Plant.Generator(j).OpMatA,'Stor')
        chargeCapacity = 0;
        if isfield(Plant.Generator(j).OpMatA.output, 'E')
            include = {'Electric Generator','CHP Generator'};
            for i = 1:1:nG
                if ismember(Plant.Generator(i).Type,include)
                    chargeCapacity = chargeCapacity+Plant.Generator(i).Size;
                end
            end
        elseif isfield(Plant.Generator(j).OpMatA.output, 'C')
            include = {'Chiller'};
            for i = 1:1:nG
                if ismember(Plant.Generator(i).Type,include)
                    chargeCapacity = chargeCapacity+Plant.Generator(i).Size;
                end
            end
        elseif isfield(Plant.Generator(j).OpMatA.output, 'H')
            include1 = {'Heater'};
            include2 = {'CHP Generator'};
            for i = 1:1:nG
                if ismember(Plant.Generator(i).Type,include1)
                    chargeCapacity = chargeCapacity+Plant.Generator(i).Size;
                elseif ismember(Plant.Generator(i).Type,include2)
                    chargeCapacity = chargeCapacity+Plant.Generator(i).Size.*Plant.Generator(i).OpMatA.output.H;%OpMatA.output.H is the heat ratio
                end
            end
        end
        Buffer = min((BuffPerc/100)*UB(j), (chargeCapacity))*dischEff(j);
        Plant.Generator(j).OpMatA.link.bineq(end-1) = -Buffer; %lower buffer ineq :  -SOC - W <= -Buffer becomes W>= buffer -SOC
        Plant.Generator(j).OpMatA.link.bineq(end) = UB(j)*dischEff(j)-Buffer; %upper buffer ineq :  SOC - Z <= (UB-Buffer)
        Plant.Generator(j).OpMatA.Z.ub = Buffer;
        Plant.Generator(j).OpMatA.W.ub = Buffer;
        Plant.Generator(j).OpMatB = Plant.Generator(j).OpMatA;
    end
end