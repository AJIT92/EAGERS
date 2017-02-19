function InitialDispatch
global Plant DateSim scaleTime  HistProf LB UB Last24hour RealTimeData CurrentState OnOff Si 
Time = buildTimeVector(Plant.optimoptions);
Outs = Plant.optimoptions.Outputs;
%% Load Historicals
Si = 2;
nS = length(Time);
RealTimeData.Timestamp = [0,Time]./24+DateSim;
RealTimeData.Temperature = ones(1,length(RealTimeData.Timestamp));
HistProf.Timestamp = Time./24+DateSim;
HistProf.Temperature= ones(12,nS);
%create the surface of the demand so that updateMatrices forecasts the
%exact demand that is input
x = zeros(nS*2+2,1);
y = zeros(nS*2+2,1);
x(1:2:nS*2+1,:) = [0,Time];
x(2:2:end) = [0,Time];
y(1:2:end) = 2;%y is the temperature variable and goes from 0 to 2, temp is always 1 in TestDisp
z = zeros(nS*2+2,1);
S = fields(Plant.Data.Demand);
for i = 1:1:length(S)
    ID.(S{i}) = Plant.Data.Demand.(S{i})(1); %Initial Demand
    day1 = interp1((Plant.Data.Timestamp - Plant.Data.Timestamp(1))*24,Plant.Data.Demand.(S{i}),Time);
    dem = [ID.(S{i}), day1]';
    z(1:2:end) = dem;
    z(2:2:end) = dem;
    RealTimeData.Demand.(S{i}) = [ID.(S{i}), Plant.Data.Demand.(S{i})(1:nS)];
    HistProf.(S{i}) = fit([x,y],z,'linearinterp');
end
nG = length(Plant.Generator);

%% Load generators and find initial condition
scaleTime = 1;
loadGenerator 
stor = zeros(1,nG);
for i=1:1:nG
    if isfield(Plant.Generator(i).OpMatA,'Stor')
        stor(i) = i;
    end
end
scaleCost = updateGeneratorCost(Time);
Data_t0 = GetCurrentData(DateSim);
IC = StepByStepDispatch(Data_t0.Demand,scaleCost(1,:),Plant.optimoptions.Resolution, [],'',[]);
IC(stor>0) = .5*UB(stor>0); % IC = halfway charged energy storage
OnOff = logical(IC>=LB);
CurrentState.Generators=IC;
nS_step1 = round(Plant.optimoptions.Horizon/Time(1));
Last24hour =[];%eliminate any old data stored here
%need to have this in terms of the first timestep
Last24hour = updateForecast(DateSim-1,(1:nS_step1).*Time(1));%need to have one extra timestep for prediction
Last24hour.Timestamp = DateSim-1+((0:nS_step1-1).*(Time(1)/24));

%% load fitA, update matrices and run 1st step
Organize = Plant.OpMatA.Organize;
marginCost = updateMarginalCost(UB,scaleCost(1,:),Time);
[QPall,~] = updateMatrices(Plant.OpMatA.QP,Plant.OpMatA.Organize,IC,Time,scaleCost,marginCost,[]); 

[FirstDisp, ~, Feasible] = DispatchQP(QPall,Organize,[]);%this is the dispatch with fit A

%% use modified generator demand profile for Step 2
Out = Plant.optimoptions.Outputs;
%change demand to not include the amount covered by storage
for j = 1:1:length(Out);
    GenDemand.(Out{j}) = zeros(nS,1);
    for i = 1:1:length(Plant.Generator)
        if isfield(Plant.Generator(i).OpMatB.output,Out{j})
            if ~isfield(Plant.Generator(i).OpMatB,'Stor') %don't include the storage use/charging in this new profile;
                GenDemand.(Out{j}) = GenDemand.(Out{j}) + FirstDisp(2:end,i)*Plant.Generator(i).OpMatB.output.(Out{j});
            end
        end
    end
end
tic
OptimalState = StepByStepDispatch(GenDemand,scaleCost,ones(nS,1)*Plant.optimoptions.Resolution,IC, 'initially constrained',FirstDisp);
% disp(['time for step 2 is ' num2str(toc)])
%% check if any gen can be completley removed
Locked = true(nS+1,nG);
%% Rule 1 if off-line for entire optimal dispatch, then the generator should be off
redoStep1 = false;
for i = 1:1:nG
    if nnz(OptimalState(:,i))==0
        Locked(:,i)= false;
        redoStep1 = true;
    end
end
if redoStep1
    [FirstDisp, ~, Feasible] = DispatchQP(QPall,Organize,Locked);%this is the dispatch with fit A
    if Feasible ~=1
        Locked = true(nS+1,nG); %reset to all on default
    end
end
%% load & run fitB for step 3
tic
Organize = Plant.OpMatB.Organize;
marginCost = updateMarginalCost(FirstDisp,scaleCost,Time);
[QPall,Forecast] = updateMatrices(Plant.OpMatB.QP,Organize,IC,Time,scaleCost,marginCost,[]); %update fit B matrices
[GenDisp, dX] = FilterGenerators(QPall,Organize,IC,Forecast,FirstDisp,OptimalState,scaleCost,Locked);
plotDispatch2(GenDisp,Time)