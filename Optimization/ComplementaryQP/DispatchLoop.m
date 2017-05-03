function DispatchLoop%% calculate new optimal Dispatch
global Plant RealTime DispatchWaitbar Dispatch Si dischEff scaleTime UB LB timersOn% loaded by GUI & load generators
global CurrentState DateSim LBmpc %changed in MPC loop
global EC targetCharge targetDischarge  %result from dispatch loop used in Online loop
options = Plant.optimoptions;
Time = buildTimeVector(options);%% set up dt vector of time interval length
dt = Time - [0, Time(1:end-1)];
nG = length(Plant.Generator);
nS = length(Time);
%% step 1 Update IC & matrices (including Make forecast & predict renewables, to create forecasted demands)
stor = find(dischEff>0);
if RealTime
    [CurrentState.Generators,~,~,~] = GetCurrentState;
    CurrentState.Generators(stor) = CurrentState.Generators(stor)/(3600)*scaleTime; %convert to kWh & scale storage
    IC = max(0,CurrentState.Generators);
else
    for i = 1:1:nG
        if isfield(Plant.Generator(i).OpMatA.output, 'C') %make sure chillers have positive initial conditions
            if CurrentState.Generators(i)<0
                CurrentState.Generators(i) = LBmpc(i);
            end
        end
    end
    IC = CurrentState.Generators;
end
if ~isempty(dischEff) %must be a storage system
    IC(stor)=max(IC(stor).*(1-(1-dischEff(stor))*dt(1)),LB(stor));%scale storage value by discharge efficiency (this scaling is used in optimizations)
end
IC = min(UB,IC);
scaleCost = updateGeneratorCost(Time/24+DateSim); %% All feedstock costs were assumed to be 1 when building matrices 
if Si<=1
    PredictDispatch = ones(length(Time)+1,1)*IC;
else PredictDispatch = [IC;Dispatch.Predicted(Si).GenDisp(3:end,:);Dispatch.Predicted(Si).GenDisp(end,:)];
end
marginCost = updateMarginalCost(PredictDispatch,scaleCost,Time);%the dispatch is whatever has been dispatched so far, except for the initial condition.
Organize = Plant.OpMatA.Organize;
[QPall,~] = updateMatrices(Plant.OpMatA.QP,Organize,IC,Time,scaleCost,marginCost,[]);
%% Step 2 Record Operation & Summarize the Demand/Dispatch at major intervals (Si)
if Si>1
    recordFromMPCloop
end
Si=Si+1;
%% Step 3 Determine initial dispatch
[FirstDisp,~,Feasible] = DispatchQP(QPall,Organize,[]);
if ~(Feasible==1)%% hopefully not here
    disp('error: initial dispatch was infeasible. Defaulting to previous dispatch');
    FirstDisp= PredictDispatch;
end

%% Step 4 Refine Dispatch & Plot
Out = Plant.optimoptions.Outputs;
for j = 1:1:length(Out);
    GenDemand.(Out{j}) = zeros(nS,1);
    StorPower.(Out{j}) = zeros(nS,1);
    for i = 1:1:nG
        if isfield(Plant.Generator(i).OpMatA.output,Out{j})
            if ~isfield(Plant.Generator(i).OpMatA,'Stor') %don't include the storage use/charging in this new profile;
                GenDemand.(Out{j}) = GenDemand.(Out{j}) + FirstDisp(2:end,i)*Plant.Generator(i).OpMatB.output.(Out{j});
            else  StorPower.(Out{j}) = StorPower.(Out{j}) + (FirstDisp(1:end-1,i)-FirstDisp(2:end,i))./dt'; 
            end
        end
    end
end
tic
OptimalState = StepByStepDispatch(GenDemand,scaleCost,dt,IC,'initially constrained',FirstDisp);
if timersOn && Si==2
    disp(['time for step by step dispatch is ' num2str(toc)])
end
% OptimalState = OptimalMap3Dread(IC,GenDemand); %% pre-calculated tables for generator deployment
% if isempty(OptimalState)
% %% compare to baseline with original demand:
% [marginCost, OptimalState2] = StepByStepDispatch(Forecast,scaleCost,dt,IC,'initially constrained',GenDisp);
% end

%% check if any gen can be completley removed
Locked = true(nS+1,nG);
%% Rule 1 if off-line for entire optimal dispatch, then the generator should be off
for i = 1:1:nG
    if nnz(OptimalState(:,i))==0
        Locked(:,i)=false;
    end
end
tic
[FirstDisp, ~, Feasible] = DispatchQP(QPall,Organize,Locked);%this is the dispatch with fit A
if Si==2 && timersOn
    disp(['time for fit A dispatch is ' num2str(toc)])
end
if Feasible ~=1
    Locked = true(nS+1,nG); %reset to all on default
end
    
Organize = Plant.OpMatB.Organize;
marginCost = updateMarginalCost(FirstDisp,scaleCost,Time);
[QPall,Forecast] = updateMatrices(Plant.OpMatB.QP,Organize,IC,Time,scaleCost,marginCost,[]); %update fit B matrices
tic
[GenDisp, dX] = FilterGenerators(QPall,Organize,IC,Forecast,FirstDisp,OptimalState,scaleCost,Locked);
if Si==2 && timersOn
    disp(['time for fit B dispatch with rules is ' num2str(toc)])
end
LimitUBforShutOff(GenDisp,dX);%Limit UB if generator needs to shut off soon!!!
EC = GenDisp(2,:); %forecasted end condition
targetCharge = max(0,(EC(stor)-IC(stor))/dt(1));
targetDischarge = max(0,(IC(stor)-EC(stor))/dt(1));
plotDispatch(GenDisp,Time)

%% Step 5 identify On/Off threshold for use on-line
ThresholdSet(GenDisp,dt(1));

%% Step 6: Record predicted dispatch && prepare for control over next prediction timestep
% need to make adjustments to ensure we don't hog memory with this (clear
% some of the old predictions during long runs?)
Dispatch.Predicted(Si).Timestamp = Time;
Dispatch.Predicted(Si).Forecast = Forecast;
Dispatch.Predicted(Si).GenDisp = GenDisp;

%% Step 7 Prepare for OnLine Loop
options = Plant.optimoptions;
if Time(1)>=options.Topt/3600 %if the first step is shorter than the amount of time it takes to make one online loop, then make the horizon be one online loop
    options.Horizon = Time(1);%the horizon is not the resolution
    Plant.optimoption.Topt = Time(1).*3600;
else options.Horizon = options.Topt/3600;
end
options.Resolution = options.Topt/3600;%the resolution is the frequency of Topt converted to hours
options.tspacing = 'constant';
Timestamp = buildTimeVector(options)/24+DateSim;%create TimeSim here so that it is constant throughout each dispatch loop. TimeSim timestamps for each run of OptimOptions
n = length(Plant.Online)-length(Timestamp);%if the timestep is short and doesn't allow for normal amount of online loops, then use QP with correct number of steps
for t = 1:1:length(Timestamp)
    Plant.Online(t+n).Timestamp = Timestamp(t);
end
OnlineLoop
waitbar(Si/(Plant.optimoptions.Horizon/dt(1)),DispatchWaitbar,strcat('Running Dispatch'));