function [GenDisp,tsim] = DispatchLoop(IC,Forecast,Renewable,PredictDispatch)
%% calculate optimal dispatch over the forecast horizon
global Plant %%  loaded by GUI & load generators
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
dt = Time - [0, Time(1:end-1)];
nG = length(Plant.Generator);
nS = length(Time);
%% Update IC & matrices (including Make forecast & predict renewables, to create forecasted demands)
ren = [];
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Source,'Renewable')
        ren(end+1) = i;
    end
end

scaleCost = updateGeneratorCost(Time); %% All feedstock costs were assumed to be 1 when building matrices 
marginCost = updateMarginalCost(PredictDispatch,scaleCost,Time);%the dispatch is whatever has been dispatched so far, except for the initial condition.
QP = updateMatrices(Plant.OpMatA,IC,Time,scaleCost,marginCost,Forecast,Renewable,[]);
%% Step 1 Determine initial dispatch
Locked = true(nS+1,nG);
for i = 1:1:nG
    if ~Plant.Generator(i).Enabled
        Locked(:,i) = 0;
    end
end
tic
[FirstDisp,~,Feasible] = DispatchQP(QP,Locked);
tsim(1,1) = toc;
for i = 1:1:length(ren)
    FirstDisp(2:end,ren(i))= Renewable(:,ren(i));
end
if ~(Feasible==1)%% hopefully not here
    disp('error: initial dispatch was infeasible. Defaulting to previous dispatch');
    FirstDisp= PredictDispatch;
end

%% Step 2:  dispatch step by step
tic
OptimalState = StepByStepDispatch(Forecast,Renewable,scaleCost,dt,IC,'initially constrained',FirstDisp);
tsim(1,2) = toc;
for i = 1:1:length(ren)
    OptimalState(2:end,ren(i))= Renewable(:,ren(i));%un-scale storage value by discharge efficiency (this scaling is used in optimizations)
end

%% Start with optimal dispatch, and check if feasible
tic
for i = 1:1:nG
    if QP.Organize.Dispatchable(i) ==1
        Locked(OptimalState(:,i)==0,i)=false;
    end
end
marginCost = updateMarginalCost(FirstDisp,scaleCost,Time);
QP = updateMatrices(Plant.OpMatB,IC,Time,scaleCost,marginCost,Forecast,Renewable,[]); %update fit B matrices
OnCost = [zeros(1,length(QP.constCost));ones(nS,1)*QP.constCost];
[GenDisp, Cost, Feasible] = DispatchQP(QP,Locked);%this is the dispatch with fit B

if Feasible ~=1
    [GenDisp, Cost, Feasible] = FindFeasible(QP,Locked);
    if ~(Feasible==1)%% hopefully not here
        disp('error: Cannot Find Feasible Dispatch');
    end
end
Cost = Cost + sum(sum(Locked.*OnCost));
tsim(1,3) = toc;