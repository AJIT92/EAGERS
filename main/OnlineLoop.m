function OnlineLoop 
%% calculate new set-points for generators & storage systems within the current dispatch time interval
%% use a fixed timestep, but the # of steps reduces each time as it approaches the next dispatch interval
%% can incorporate a revised short-term forecast
%% this function does not re-calculate the dispatch, but may re-allocate the distribution of power between generators & storage
%% this function subjects generators to different upper and lower bounds if they are in the process of ramping up or down according to the dispatch optimization
%% changes in the allocation of storage are penalized via the final state: landing at the same EC as the current dispatch time interval has zero cost, but a marginal cost = margninal cost of generation
global Plant  dischEff LB%variables from loadGenerator
global setGen UBopt LBopt EC % from the dispatch loop
global UBmpc LBmpc ShutdownRamp%editied in this loop & passed down to MPC loop
global CurrentState DateSim OnOff %from the MPC loop
global Threshold

%dischEff: Discharge efficiency of energy storage
%t_mpc: a counter for the MPC loop between sucessive online optimization (index for Operation)
%setGen: the current setting of each generator as determined by the on-line dispatch. The estimated storage output determined by the on-line dispatch converted to power
%UBmpc: upper bound of generators seen by MPC, can be different from the global upper bound due to ramping constraints.
%LBmpc: lower bound of generators seen by MPC, can be different from the global lower bound due to ramping constraints
%ShutdownRamp: if a generator has just been shut down (with OnOff not actualOnOff) this describes the residual power that will be produced during shutdown so that it can be subtracted from the optimization of the remaining generators
%CurrentState: Current state of all generators (kW) & storage (kWh) 
%DateSim: Current time in the simulation.
%OnOff: the desired generator state from the controller (may differ from actual on/off because it leaves generators on as they ramp down to min power, then shuts them off)
for t = 1:1:length(Plant.Online)
    Timestamp(t) = Plant.Online(t).Timestamp;
end
A = (Timestamp>DateSim);
nS = nnz(A); %# of timesteps remaining in Online before next dispatch
index = find(A,1,'first');
nG = length(Plant.Generator);
dX = zeros(1,nG);
dt = Plant.optimoptions.Topt./3600; % interval in hours between succesive Online Loops
TimeDispLoop = buildTimeVector(Plant.optimoptions);
dtDispLoop = TimeDispLoop-[0,TimeDispLoop(1:end-1)];%this is used instead of Plant.optimoptions.Resolution
for i = 1:1:nG
    if isfield(Plant.Generator(i).OpMatB,'Ramp') 
        dX(i) = Plant.Generator(i).OpMatB.Ramp.b(1)*dt;
    end
end
IC = CurrentState.Generators;

stor = [];
UBmpc = UBopt;
downPow = [];
lockOff = zeros(1,nG);
lockon = zeros(1,nG);
LBmpc = LB;
for i = 1:1:nG
    if isfield(Plant.Generator(i).OpMatB,'Stor')
        stor(end+1) = i;
        IC(i)=IC(i).*dischEff(stor(end));%scale storage value by discharge efficiency (this scaling is used in optimizations)
        maxCharge = (UBopt(i)-IC(i))/dtDispLoop(1);%Plant.optimoptions.Resolution;
        maxDischarge = (IC(i)-LBopt(i))/dtDispLoop(1);%Plant.optimoptions.Resolution;
        UBmpc(i) = min(dX(i),maxDischarge);
        LBmpc(i) = max(-dX(i),-maxCharge);
    elseif ~strcmp(Plant.Generator(i).Source, 'Renewable') && ~strcmp(Plant.Generator(i).Type, 'Utility') && ~strcmp(Plant.Generator(i).Type, 'DistrictCooling') && ~strcmp(Plant.Generator(i).Type, 'DistrcitHeating')
        if OnOff(i)>0
            UBmpc(i) = min(UBopt(i),(CurrentState.Generators(i)+dX(i)*nS));
            LBmpc(i) = max(LBopt(i),CurrentState.Generators(i)-dX(i)*nS);
        else
            LBmpc(i) = 0;
            if IC(i)>0 %% if a generator recently shut off, account for that here
                dRamp = ShutdownRamp(i).t2>=Timestamp(index);
                if ~isempty(Plant.Generator(i).Output) && nnz(dRamp)>0
                    downPow = [];
                    Outs = fieldnames(Plant.Generator(i).OpMatA.output);
                    lockOff(i) = i;
                    downPow.(Outs{1}) = zeros(1,length(dRamp));
                    downPow.(Outs{1})(1:nnz(dRamp)) = downPow.(Outs{1})(1:nnz(dRamp))+ShutdownRamp(i).Pow2(dRamp);%forecast power during shutdown
                    %% need to add forecasted heat production for CHP generators
                end
            end
        end
    elseif strcmp(Plant.Generator(i).Type, 'DistrictCooling') || strcmp(Plant.Generator(i).Type, 'DistrcitHeating')
        LBmpc(i) = LBopt(i);
        UBmpc(i) = UBopt(i);
        lockon(i) = i;
    else %if it is renewable, utility, or district heat/cool
        lockon(i) = i;
    end
end
%make sure that if there is a switch, then locked gets adjusted accordingly
Outs = Plant.optimoptions.Outputs;
tswitch = inf;
genswitchon = [];
genswitchoff = [];
for s = 1:1:length(Outs)
    if isfield(Threshold.(Outs{s}),'t')
        time = min(Threshold.(Outs{s}).t);
        if time<tswitch
            genswitchon = Threshold.(Outs{s}).On;
            genswitchoff = Threshold.(Outs{s}).Off;
        end
    else time = inf;
    end
    tswitch = min(time,tswitch);
end
tlocked = nnz(Timestamp(A)<=tswitch); %number that should have the first set locked

lockOff = nonzeros(lockOff);
lockon = nonzeros(lockon);
scaleCost = updateGeneratorCost(DateSim);
marginCost = updateMarginalCost([IC;EC],scaleCost,dt);%give the marginal cost for the power from now until the end of the online loop
QPall = Plant.Online(index).QP;
Organize = Plant.Online(index).Organize;
[QPall,~] = updateMatrices(QPall,Organize,IC,(Timestamp(A)-DateSim)*24,scaleCost,marginCost,EC);
Locked = false(nS+1,nG);
Locked(:,OnOff>0) =true;
Locked(:,lockOff) = false; %if its ramping down lock it off in the QP.
% Locked(2:end,EC>LBopt) = true; %if its ramping up allow it on
Locked(1:end,lockon) = true;
if tlocked<nS
    Locked(tlocked+1:end,genswitchon) = true;
    Locked(tlocked+1:end,genswitchoff) = false;
end
%subtract the ramp down power from the demand
if ~isempty(downPow)
    Outs = fieldnames(downPow);
    for j = 1:1:length(Outs);
        QPall.(Outs{j}).beq(1:length(downPow.(Outs{j}))) = QPall.(Outs{j}).beq(1:length(downPow))-downPow.(Outs{j});
    end
end
[GenDisp,~,~] = DispatchQP(QPall,Organize,Locked);
setGen = zeros(nG,1);
setGen(OnOff>0) = GenDisp(2,(OnOff>0));
setGen(EC>LBopt) = GenDisp(2,(EC>LBopt));
if ~isempty(stor)
    setGen(stor) = (IC(stor)-EC(stor) - GenDisp(2,stor)).*(3600/Plant.optimoptions.Topt); %convert change in energy storage to power
end