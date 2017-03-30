function GenOutput = StepByStepDispatch(Demand,Renewable,scaleCost,dt,IC,limit,FirstProfile)
% Time is the time from the current simulation time (DateSim), positive numbers are forward looking
global Plant UB dX_dt 
nG = length(Plant.Generator);
if isempty(IC)
    nS = 0;
else
    [nS,~] = size(scaleCost);
end
stor = find(Plant.OneStep.Organize.StorageEquality>0);
StorPower = zeros(nS,nG);
GenOutput = zeros(nS+1, nG);%nS should equal 1 in this case (finding IC)
if ~isempty(IC)
    GenOutput(1,:) = IC;
end
TestCombos = zeros(nS,1); %how many QP optimizations needed to be run in eliminate combinations
%note #1: IC of storage is already sclaed by discharge losses
%% note #2:  need to account for self discharging of storage
Outs = fieldnames(Demand);
Outs = Outs(~strcmp(Outs,'T'));
for t = 1:1:max(1,nS) %for every timestep
    for j = 1:1:length(Outs)
        Loads.(Outs{j}) = Demand.(Outs{j})(t,:); %update demand
        netDemand.(Outs{j}) = sum(Loads.(Outs{j}));
    end
    if isempty(FirstProfile) %finding initial conditions
        QP = updateMatrices1Step(Plant.OneStep,Loads,Renewable,scaleCost(t,:),dt(t),[],[],[],[]);
    else
        if strcmp(limit, 'constrained')
            MinPower = max(0,IC-dX_dt*dt(t));
            MaxPower = min(UB,IC+dX_dt*dt(t));
        else
            MinPower = max(0,IC-dX_dt*sum(dt(1:t)));
            MaxPower = min(UB,IC+dX_dt*sum(dt(1:t)));
        end
        StorPower(t,stor) = (FirstProfile(t,stor) - FirstProfile(t+1,stor))/dt(t);
        StorPower(t,stor) = (IC(stor) - FirstProfile(t+1,stor))/dt(t); %expected output of storage in kW to reach the SOC from the 1st dispatch (penalties are always pushing it back on track if it needed more storage than expected somewhere)
        QP = updateMatrices1Step(Plant.OneStep,Loads,Renewable(t,:),scaleCost(t,:),dt(t),FirstProfile(t,:),StorPower(t,:),MinPower,MaxPower);
    end
    [BestDispatch,TestCombos(t)] = eliminateCombinations(QP,netDemand);%% combination elimination loop
    if isempty(BestDispatch)
        disp(['No feasible combination of generators at step' num2str(t)]);
        BestDispatch = IC;
        BestDispatch(stor) = 0;
        TestCombos(t)=0;
    end
    EC = BestDispatch;
    if ~isempty(IC)
        %if charging, subtract charging losses
        dSOC = zeros(1,length(stor));
        for i = 1:1:length(stor)
            if (BestDispatch(stor(i))+StorPower(t,stor(i)))>0 %discharging
                dSOC(i) = -(BestDispatch(stor(i))+StorPower(t,stor(i)))*dt(t);
            else %charging
                eff = Plant.Generator(stor(i)).OpMatA.Stor.DischEff*Plant.Generator(stor(i)).OpMatA.Stor.ChargeEff;
                dSOC(i) = -(BestDispatch(stor(i))+StorPower(t,stor(i)))*dt(t)*eff; 
            end
        end
        EC(stor) = IC(stor) + dSOC; %change in storage for this power output
        if strcmp(limit, 'constrained')%if its constrained but not initially constrained then make the last output the initial condition
            IC = EC;
        else
            IC(stor) = EC(stor); %update the state of storage
        end
    end
    GenOutput(t+(~isempty(IC)),:) = EC;
end
% disp(['Average # of combinations tested is ' num2str(sum(TestCombos)/nS)]);