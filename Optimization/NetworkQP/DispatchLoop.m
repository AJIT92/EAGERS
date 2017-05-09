function [GenDisp,tsim] = DispatchLoop(Date,Time,IC,Forecast,PredictDispatch)
%% calculate optimal dispatch over the forecast horizon
global Plant%%  loaded by GUI & load generators
dt = Time - [0; Time(1:end-1)];
nG = length(Plant.Generator);
nS = length(Time);
%% Update IC & matrices (including Make forecast & predict renewables, to create forecasted demands)
if isempty(PredictDispatch)
    PredictDispatch = ones(length(Time)+1,1)*IC;
end
scaleCost = updateGeneratorCost(Time/24 + Date); %% All feedstock costs were assumed to be 1 when building matrices 
marginCost = updateMarginalCost(PredictDispatch,scaleCost,Time);%the dispatch is whatever has been dispatched so far, except for the initial condition.
QP = updateMatrices(Plant.OpMatA,IC,Date,Time,scaleCost,marginCost,Forecast,[]);
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
if ~(Feasible==1)%% hopefully not here
    disp('error: initial dispatch was infeasible. Defaulting to previous dispatch');
    FirstDisp= PredictDispatch;
end
%% Step 2:  dispatch step by step
if Plant.optimoptions.MixedInteger
    tic
    OptimalState = StepByStepDispatch(Forecast,scaleCost,dt,IC,'initially constrained',FirstDisp);
    tsim(1,2) = toc;
else
    OptimalState = FirstDisp;
end
%% Start with optimal dispatch, and check if feasible
tic
marginCost = updateMarginalCost(OptimalState,scaleCost,Time);
QP = updateMatrices(Plant.OpMatB,IC,Date,Time,scaleCost,marginCost,Forecast,[]); %update fit B matrices
for i = 1:1:nG
    if QP.Organize.Dispatchable(i) ==1
        Locked(OptimalState(:,i)==0,i)=false;
    end
end
[GenDisp, ~, Feasible] = DispatchQP(QP,Locked);%this is the dispatch with fit B
if Feasible ~=1
    [GenDisp, QP, Feasible] = FindFeasible(QP,Locked);
end
if Feasible==1
    if ~Plant.optimoptions.MixedInteger
        GenDisp = FilterGenerators(QP,GenDisp,Locked,[0;Time/24] + Date);
    end
else
    disp('error: Cannot Find Feasible Dispatch');
end
tsim(1,3) = toc;
% Cost = sum(NetCostCalc(GenDisp,[0;Time/24]+Date,'Dispatch'));