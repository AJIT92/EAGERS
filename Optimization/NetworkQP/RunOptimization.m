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
timers = zeros(1,3); % To record times set to zeros(1,3), to not record set to [];
scaleTime = Plant.optimoptions.scaletime; %should always be 1 for this type of optimization

LoadTestData %% Revise this so you can pull from more than what is loaded in Plant
if ~isfield(Plant.Data,'HistProf')
    calculateHistoricalFit %% calculate fits used in forecasting
end
Xi = nnz(Plant.Data.Timestamp<=DateSim);
Xf = nnz(Plant.Data.Timestamp<=DateSim+Plant.optimoptions.Interval);
if any(strcmp(Plant.optimoptions.method,{'Dispatch';'Planning'}))
    RealTimeData = interpolateData(Plant.optimoptions.Resolution*3600,Xi,Xf,0.00);%create test data at correct frequency
else
    RealTimeData = interpolateData(Plant.optimoptions.Tmpc,Xi,Xf,0.00);%create test data at correct frequency
end

tic
loadGenerator % Loads generators & build optimization matrices, stored in Plant.Generator
build_subNet
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
dt = Time - [0; Time(1:end-1)];
Plant.OpMatA = buildMatrices('OpMatA',dt); %build quadratic programming matrices for FitA
Plant.OpMatB = buildMatrices('OpMatB',dt);%build quadratic programming matrices for FitB
Plant.OneStep = buildMatrices1Step;%build quadratic programming matrices for 1Step at interval spacing of dt

if strcmp(Plant.optimoptions.method,'Control')
    A.Horizon = Plant.optimoptions.Resolution;%the horizon is the resolution
    A.Resolution = Plant.optimoptions.Topt/3600;%the resolution is the frequency of Topt
    A.tspacing = 'constant';
    OnlineTime = buildTimeVector(A);%% set up dt vector of time interval length
    dt2 = OnlineTime - [0, OnlineTime(1:end-1)];
    for t = 1:1:length(OnlineTime)
        Plant.Online(t) = buildMatrices('OpMatB',dt2(t:end)); %build the matrix for the onlineOptimLoop using FitB
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
Outs =  fieldnames(Data_t0.Demand);
for i = 1:1:length(Outs)
    Demand_t0.(Outs{i}) = Data_t0.Demand.(Outs{i});
end
Data_t0.Demand = Demand_t0;
% K = center_menu('Select Option','Manually specify initial conditions','Automatically determine initial conditions');
K = 2;
if K ==1
    IC = manualInitialCondition;
else
    IC = automaticInitialCondition(Data_t0);
end
OnOff = logical(IC(1:nG)>=LB);
CurrentState.Generators=IC(1:nG);
CurrentState.Lines=IC(nG+1:nG+nL);
FirstIC = IC;

%% initialize all other variables
NumSteps = Plant.optimoptions.Interval*24/Plant.optimoptions.Resolution;
Plant.Dispatch.Temperature = zeros(NumSteps,1);
Plant.Dispatch.Timestamp = zeros(NumSteps,1);
Plant.Dispatch.GeneratorState = zeros(NumSteps,nG+nL);

% Plant.RunData.Timestamp = zeros(NumSteps,1);
% % Plant.RunData.Timestamp(1)  = Data_t0.Timestamp;
% Plant.RunData.GeneratorState = zeros(NumSteps,nG+nL);
% % Plant.RunData.GeneratorState(1,:) = IC;
% 
% Plant.RunData.GeneratorInput = zeros(NumSteps,nG+nL);

Plant.Predicted.GenDisp = zeros(round(Plant.optimoptions.Horizon/Plant.optimoptions.Resolution)+1,nG+nL,NumSteps);
Plant.Predicted.Timestamp = zeros(round(Plant.optimoptions.Horizon/Plant.optimoptions.Resolution)+1,NumSteps);
S = fieldnames(Data_t0.Demand);
for i = 1:1:length(S)
    Plant.Dispatch.Demand.(S{i}) = zeros(NumSteps,length(Data_t0.Demand.(S{i})));
    Plant.RunData.Demand.(S{i}) = zeros(NumSteps,length(Data_t0.Demand.(S{i})));
end

nS_step1 = round(24/Plant.optimoptions.Resolution);
%need to have this in terms of the first timestep
Xi = nnz(Plant.Data.Timestamp<=DateSim);
Xf = nnz(Plant.Data.Timestamp<=DateSim+1);
Last24hour = interpolateData(Plant.optimoptions.Resolution*3600,Xi,Xf,0.00);
Last24hour.Timestamp = DateSim-1+((0:nS_step1-1)'.*(Plant.optimoptions.Resolution/24));
for i = 1:1:nG
    if isfield(Plant.Generator(i).VariableStruct,'RestartTime')
        RestartTime(i) = Plant.Generator(i).VariableStruct.RestartTime/60;%restart time in hours
    else
        RestartTime(i) = 0;
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
        Forecast = updateForecast(DateSim,RealData);%% function that creates demand vector with time intervals coresponding to those selected
        %can expand forecast to create high and low estimates to be used in 2nd/3rd call of dispatch loop or step by step
        IC = [CurrentState.Generators, CurrentState.Lines];
        for i = 1:1:nG
            IC(i) = IC(i)*Plant.Generator(i).Enabled;%remove IC for disabled gens
        end
        if ~isempty(GenDisp)
            PredictDispatch = [IC,zeros(1,nL);Plant.Predicted.GenDisp(3:end,:,Si-1);Plant.Predicted.GenDisp(end,:,Si-1)];
        end
        [GenDisp,tsim] = DispatchLoop(DateSim,Time,IC,Forecast,PredictDispatch);
        Plant.Predicted.GenDisp(:,:,Si) = GenDisp;
        Plant.Predicted.Timestamp(:,Si) = DateSim+[0;Time]/24;
        Plant.Dispatch.GeneratorState(Si,:) = GenDisp(1,:);
        Plant.Dispatch.Timestamp(Si) = DateSim;
        Plant.Dispatch.Temperature(Si) = Forecast.T(1);
        for i = 1:1:length(S)
            Plant.Dispatch.Demand.(S{i})(Si,:) = Forecast.(S{i})(1,:);
        end

        if ~isempty(timers)
            timers(Si,1:3) = tsim;
%             disp(strcat('FistDisp:',num2str(tsim(1))));
%             disp(strcat('StebByStep:',num2str(tsim(2))));
%             disp(strcat('FinalDisp:',num2str(tsim(3))));
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
                History(1:2,:) = Plant.Predicted.GenDisp(1:2,:,Si-backSteps);
                for t = 2:1:backSteps
                    History(t+1,:) = Plant.Predicted.GenDisp(2,:,Si-backSteps+t-1);
                end
            end
        end
        updateGUIstatus(Plant.Predicted.GenDisp(1:2,:,Si),History)
        plotDispatch(Plant.Predicted.Timestamp(:,Si),Plant.Predicted.GenDisp(:,:,Si),History)

        if strcmp(Plant.optimoptions.method,'Planning')
            CurrentState.Generators = Plant.Predicted.GenDisp(end,1:nG,Si);
            CurrentState.Lines = Plant.Predicted.GenDisp(end,nG+1:nG+nL,Si);
            Xi = nnz(Plant.Data.Timestamp<=DateSim);
            Xf = nnz(Plant.Data.Timestamp<DateSim+1);
            Last24hour = interpolateData(Plant.optimoptions.Resolution*3600,Xi,Xf,0.00);
            nS_step1 = round(24/Plant.optimoptions.Resolution);
            Last24hour.Timestamp = DateSim-1+((0:nS_step1-1)'.*(Plant.optimoptions.Resolution/24));
            DateSim = round(1e5*(DateSim+Plant.optimoptions.Horizon/24))/1e5;%%count forward by length of the horizon, rounded to nearest second
        elseif strcmp(Plant.optimoptions.method,'Dispatch')
            CurrentState.Generators = Plant.Predicted.GenDisp(2,1:nG,Si);
            CurrentState.Lines = Plant.Predicted.GenDisp(2,nG+1:nG+nL,Si);
            Outs = fieldnames(Last24hour.Demand);
            Last24hour.Timestamp = [Last24hour.Timestamp(2:end);RealData.Timestamp;];
            Last24hour.Temperature = [Last24hour.Temperature(2:end,:);RealData.Temperature;];
            for i = 1:1:length(Outs)
                Last24hour.Demand.(Outs{i}) = [Last24hour.Demand.(Outs{i})(2:end,:);RealData.Demand.(Outs{i})(1,:);];
            end
            DateSim = round(1e5*(DateSim+Plant.optimoptions.Resolution/24))/1e5;%% count forward 1 step, rounded to nearest second
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