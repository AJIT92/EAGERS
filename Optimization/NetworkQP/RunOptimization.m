global Plant RealTimeData
global Virtual scaleTime DateSim DispatchWaitbar %relates to how optimization proceeds
global CurrentState Si Last24hour OnOff timers%initialized here
global GenAvailTime RestartTime %  Global vairables that need to be reset between each run, but not each loop. Trying to remove these global variables
%Step 1: Step 1 Load generators, & build QP matrices
%Step 2: Initialize Variables & prepare to iterate time steps
%Step 3: Run Through Dispatch Optimization 

%% variable descriptions
%RealTimeData: the high resolution data used for this test

%Virtual: Running a simulation only, set to zero when the end of the test data set is reached
%scaletime: the ratio of emulated time to time in the test-data. For example 24 hours of test data can be run in 3 hours with a scaletime of 8. scaletime enlarges any energy storage and slows down the transient response of any generators
%DateSim: Current time in the simulation.
%NumSteps: the number of dispatch optimiztions that will occur during the entire simulation

%CurrentState: Current state of all generators (kW) & storage (kWh) 
%Dispatch: recorded data at the frequency of the dispatch optimization
%Si: Counter for dispatch loop
%Last24hour: recorded data for the last 24 hours
%OnOff: the desired generator state from the controller (may differ from actual on/off because it leaves generators on as they ramp down to min power, then shuts them off)

%GenAvailTime: The next time a generator is available to turn on
%RestartTime: The amount of time necessary to wait for restarting after a generator has turned off .

%% Step 1 Load generators, & build QP matrices
timers = []; % To record times set to zeros(1,3), to not record set to [];
scaleTime = Plant.optimoptions.scaletime; %should always be 1 for this type of optimization

LoadTestData %% Revise this so you can pull from more than what is loaded in Plant
if ~isfield(Plant.Data,'HistProf')
    calculateHistoricalFit %% calculate fits used in forecasting
end
Xi = nnz(Plant.Data.Timestamp<=DateSim);
Xf = nnz(Plant.Data.Timestamp<=DateSim+Plant.optimoptions.Interval);
RealTimeData = interpolateData(Plant.optimoptions.Tmpc,Xi,Xf,0.00);%create test data at correct frequency

tic
loadGenerator % Loads generators & build optimization matrices, stored in Plant.Generator
build_subNet
Plant.OpMatA = buildMatrices('A'); %build quadratic programming matrices for FitA
Plant.OpMatB = buildMatrices('B');%build quadratic programming matrices for FitB
Plant.OneStep = buildMatrices1Step;%build quadratic programming matrices for 1Step at interval spacing of dt

if strcmp(Plant.optimoptions.method,'Control')
    A.Horizon = Plant.optimoptions.Resolution;%the horizon is the resolution
    A.Resolution = Plant.optimoptions.Topt/3600;%the resolution is the frequency of Topt
    A.tspacing = 'constant';
    OnlineTime = buildTimeVector(A);%% set up dt vector of time interval length
    dt2 = OnlineTime - [0, OnlineTime(1:end-1)];
    for t = 1:1:length(OnlineTime)
        Plant.Online(t) = buildMatrices('B',dt2(t:end)); %build the matrix for the onlineOptimLoop using FitB
    end
end

nG = length(Plant.Generator);
[~,n] = size(Plant.OneStep.organize);
nL = n-nG;
LB = zeros(1,nG);
for i = 1:1:nG
    states = Plant.Generator(i).OpMatB.states;
    for j = 1:1:length(states);
        LB(i) = LB(i) + Plant.Generator(i).OpMatB.(states{j}).lb;
    end
end

%% Step 2 Initialize Variables & prepare to iterate time steps
Data_t0 = GetCurrentData(DateSim);
%only optimize for the demands specified by the user
Demand_t0 = [];
Outs =  Plant.optimoptions.Outputs;
for i = 1:1:length(Outs)
    Demand_t0.(Outs{i}) = Data_t0.Demand.(Outs{i});
end
Data_t0.Demand = Demand_t0;
% K = center_menu('Select Option','Manually specify initial conditions','Automatically determine initial conditions');
K = 2;
if K ==1
    list = {};
    sizes = {};
    Index =[];
    IC = zeros(1,length(Plant.Generator));
    include = {'CHP Generator', 'Electric Generator', 'Chiller','Heater'};
    for i = 1:1:length(Plant.Generator)
        if isfield(Plant.Generator(i).OpMatA,'Stor')
            list(end+1) = {strcat(Plant.Generator(i).Name,' --- % of max capacity')};
            sizes(end+1) = {'50'};
            Index(end+1) = i;
        elseif ismember(cellstr(Plant.Generator(i).Type),include)
            list(end+1) = {strcat(Plant.Generator(i).Name,' --- If greater than 0, IC must be between lower bound(',num2str(LB(i)),') and upper bound(',num2str(Plant.Generator(i).Size),').')};
            sizes(end+1) = {'0'};
            Index(end+1) = i;
        elseif ~isempty(strcmp(Plant.Generator(i).Source,'Renewable'))
            %renewable
        elseif ~isempty(strcmp(Plant.Generator(i).Type,'Utility'))
            %utility
         end
    end
    IC(Index) = str2double(inputdlg(list,'Specify Initial Condition (kW) or State of Charge (%)',1,sizes));
else
    scaleCost = updateGeneratorCost(DateSim);%% All costs were assumed to be 1 when building matrices, update Generator costs for the given time
    Renewable = zeros(1,length(Plant.Generator));
    for i = 1:1:length(Plant.Generator)
        if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
            Renewable(i) = RenewableOutput(Plant.Generator(i).VariableStruct,DateSim,0,'Actual');
        end
    end
    Idispatch = StepByStepDispatch(Data_t0.Demand,Renewable,scaleCost,Plant.optimoptions.Resolution,[],'',[]);
    IC = Idispatch(1:nG);
    for i=1:1:length(IC)
        if isfield(Plant.Generator(i).OpMatA,'Stor') && Plant.Generator(i).Enabled
            IC(i) = 50; % IC = halfway charged energy storage
        elseif isfield(Plant.Generator(i).OpMatA,'Stor')
            IC(i) = 0; %storage that is disabled has 0 IC
        end
    end
end
for i=1:1:length(IC)
    if isfield(Plant.Generator(i).OpMatA,'Stor')
        IC(i) = IC(i)/100*Plant.Generator(i).OpMatA.Stor.UsableSize; % IC = halfway charged energy storage
    end
end
OnOff = logical(IC>=LB);
CurrentState.Generators=IC;
FirstIC = IC;

%% initialize all other variables
NumSteps = Plant.optimoptions.Interval*24/Plant.optimoptions.Resolution+1;
Plant.Dispatch.Temperature = zeros(NumSteps,1);
Plant.Dispatch.Timestamp = zeros(NumSteps,1);
Plant.Dispatch.GeneratorState = zeros(NumSteps,nG+nL);

Plant.RunData.Timestamp = zeros(NumSteps,1);
% Plant.RunData.Timestamp(1)  = Data_t0.Timestamp;
Plant.RunData.GeneratorState = zeros(NumSteps,nG+nL);
% Plant.RunData.GeneratorState(1,:) = IC;

Plant.RunData.GeneratorInput = zeros(NumSteps,nG+nL);

S = fieldnames(Data_t0.Demand);
for i = 1:1:length(S)
    Plant.Dispatch.Demand.(S{i}) = zeros(NumSteps,length(Data_t0.Demand.(S{i})));
    Plant.RunData.Demand.(S{i}) = zeros(NumSteps,length(Data_t0.Demand.(S{i})));
end

nS_step1 = round(24/Plant.optimoptions.Resolution);
Last24hour =[];%eliminate any old data stored here
%need to have this in terms of the first timestep
Xi = nnz(Plant.Data.Timestamp<=DateSim);
Xf = nnz(Plant.Data.Timestamp<DateSim+1);
Last24hour = interpolateData(Plant.optimoptions.Resolution*3600,Xi,Xf,0.00);
Last24hour.Timestamp = DateSim-1+((0:nS_step1-1)'.*(Plant.optimoptions.Resolution/24));
for i = 1:1:nG
    if isfield(Plant.Generator(i).VariableStruct,'RestartTime')
        RestartTime(i) = Plant.Generator(i).VariableStruct.RestartTime/60;%restart time in hours
    else RestartTime(i) = 0;
    end
end
GenAvailTime = ones(1,length(Plant.Generator)).*DateSim;
%% Step 3: Run initial dispatch and start optimization sequence (timers)
PredictDispatch = [];
GenDisp = [];
Si=1; %counter for # of times dispatch loop has run
while Virtual
    RealData = GetCurrentData(DateSim);
    if ~isempty(RealData)
        [Forecast,Renewable] = updateForecast(DateSim,RealData);%% function that creates demand vector with time intervals coresponding to those selected
        %can expand forecast to create high and low estimates to be used in 2nd/3rd call of dispatch loop or step by step
        IC = CurrentState.Generators;
        for i = 1:1:nG
            IC(i) = IC(i)*Plant.Generator(i).Enabled;%remove IC for disabled gens
        end
        if ~isempty(GenDisp)
            PredictDispatch = [IC,zeros(1,nL);Plant.Predicted(Si-1).GenDisp(3:end,:);Plant.Predicted(Si-1).GenDisp(end,:)];
        end
        [GenDisp,Time,tsim] = DispatchLoop(DateSim,IC,Forecast,Renewable,PredictDispatch);
        Plant.Predicted(Si).GenDisp = GenDisp;
        Plant.Predicted(Si).Timestamp = DateSim+[0;Time]/24;
        Plant.Dispatch.GeneratorState(Si,:) = GenDisp(1,:);
        Plant.Dispatch.Timestamp(Si) = DateSim;
        Plant.Dispatch.Temperature(Si) = Forecast.T(1);
        for i = 1:1:length(S)
            Plant.Dispatch.Demand.(S{i})(Si,:) = Forecast.(S{i})(1,:);
        end

        if ~isempty(timers)
            timers(Si,1:3) = tsim;
            disp(strcat('FistDisp:',num2str(tsim(1))));
            disp(strcat('StebByStep:',num2str(tsim(2))));
            disp(strcat('FinalDisp:',num2str(tsim(3))));
        end
        if strcmp(Plant.optimoptions.method,'Control')
            History = Plant.Dispatch.GenDisp(max(1,Si-Plant.optimoptions.Horizon/Plant.optimoptions.Resolution):Si-1,:);
        else
            if Si ==1
                History = [];
            else
                nS = length(GenDisp(:,1))-1;
                backSteps = min(Si-1,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution);
                History = zeros(backSteps+1,nG+nL);
                History(1:2,:) = Plant.Predicted(Si-backSteps).GenDisp(1:2,:);
                for t = 2:1:backSteps
                    History(t+1,:) = Plant.Predicted(Si-backSteps+t-1).GenDisp(2,:);
                end
            end
        end
        plotDispatch(Plant.Predicted(Si),History)

        if strcmp(Plant.optimoptions.method,'Planning')
            CurrentState.Generators = Plant.Predicted(Si).GenDisp(end,1:nG);
            Xi = nnz(Plant.Data.Timestamp<=DateSim);
            Xf = nnz(Plant.Data.Timestamp<DateSim+1);
            Last24hour = interpolateData(Plant.optimoptions.Resolution*3600,Xi,Xf,0.00);
            nS_step1 = round(24/Plant.optimoptions.Resolution);
            Last24hour.Timestamp = DateSim-1+((0:nS_step1-1)'.*(Plant.optimoptions.Resolution/24));
            DateSim = DateSim+Plant.optimoptions.Horizon/24;%%count forward a day
        elseif strcmp(Plant.optimoptions.method,'Dispatch')
            CurrentState.Generators = Plant.Predicted(Si).GenDisp(2,1:nG);
            Outs = fieldnames(Last24hour.Demand);
            Last24hour.Timestamp = [Last24hour.Timestamp(2:end,:);RealData.Timestamp;];
            Last24hour.Temperature = [Last24hour.Temperature(2:end,:);RealData.Temperature;];
            for i = 1:1:length(Outs)
                Last24hour.Demand.(Outs{i}) = [Last24hour.Demand.(Outs{i})(2:end,:);RealData.Demand.(Outs{i})(1,:);];
            end
            DateSim = DateSim+Plant.optimoptions.Resolution/24;%% count forward 1 step
        elseif strcmp(Plant.optimoptions.method,'Control')
            %count forward in time in the control loop
                %% Real-time control 
            if Plant.optimoptions.fastsimulation==0
                if isempty(timerfindall)
                    Timers(Plant.optimoptions)
                else %do nothing;
                end
            else
                %% Virtual Plant
                D = DateSim;
                while DateSim<(D+Plant.optimoptions.Resolution/24)
                    OnlineLoop
                end
            end
        end
        Si=Si+1;
        if ~isempty(DispatchWaitbar)
            waitbar(Si/(Plant.optimoptions.Horizon/Plant.optimoptions.Resolution),DispatchWaitbar,strcat('Running Dispatch'));
        end
    else %final conditions
        Plant.Dispatch.GeneratorState(Si+1,:) = GenDisp(2,:);
        Plant.Dispatch.Timestamp(Si+1) = DateSim;
        Plant.Dispatch.Temperature(Si+1) = Forecast.T(2);
        for i = 1:1:length(S)
            Plant.Dispatch.Demand.(S{i})(Si+1,:) = Forecast.(S{i})(2,:);
        end
    end
end
if strcmp(Plant.optimoptions.method,'Dispatch')
    Plant.NetCost = NetCostCalc(Plant.Dispatch.GeneratorState,Plant.Dispatch.Timestamp,'Dispatch');
elseif strcmp(Plant.optimoptions.method,'Control')
    Plant.NetCost = NetCostCalc(Plant.RunData.GeneratorInput,Plant.RunData.Timestamp,'Input');
end

% Plant.Baseline = RunBaseline(FirstIC); %finish simulation by running baseline