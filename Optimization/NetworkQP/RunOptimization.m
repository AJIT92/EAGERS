global Plant RealTimeData LB % Generators Interval Resolution Horizon Tmpc Topt t_mpc spGEN  %passed from GUI
global Virtual scaleTime DateSim%relates to how optimization proceeds
global CurrentState Dispatch Si Last24hour OnOff timers%initialized here
global GenAvailTime RestartTime %  Global vairables that need to be reset between each run, but not each loop. Trying to remove these global variables
global DispatchWaitbar 
%Step 1: Step 1 Load generators, & build QP matrices
%Step 2: Initialize Variables & prepare to iterate time steps
%Step 3: Run Through Dispatch Optimization 

%% variable descriptions
%RealTimeData: the high resolution data used for this test
%UB: the upper limit (capacity) of each generator
%LB: the lower limit (capacity) of each generator when on

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
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
dt = Time - [0, Time(1:end-1)];
Plant.OpMatA = buildMatrices('A',dt); %build quadratic programming matrices for FitA
Plant.OpMatB = buildMatrices('B',dt);%build quadratic programming matrices for FitB
Plant.OneStep = buildMatrices1Step;%build quadratic programming matrices for 1Step at interval spacing of dt

if strcmp(Plant.optimptions.method,'Control')
    A.Horizon = Time(1);%Plant.optimoptions.Resolution;%the horizon is not the resolution converted to seconds
    A.Resolution = Plant.optimoptions.Topt/3600;%the resolution is the frequency of Topt
    A.tspacing = 'constant';
    OnlineTime = buildTimeVector(A);%% set up dt vector of time interval length
    dt2 = OnlineTime - [0, OnlineTime(1:end-1)];
    for t = 1:1:length(OnlineTime)
        Plant.Online(t) = buildMatrices('B',dt2(t:end)); %build the matrix for the onlineOptimLoop using FitB
    end
end

nG = length(Plant.Generator);
nL = length(Plant.subNet.lineNames);
stor = [];
DischEff = [];
ren = [];
for i = 1:1:nG
    if ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
        stor(end+1) = i;
        DischEff(end+1) = Plant.Generator(i).OpMatA.Stor.DischEff;
    elseif strcmp(Plant.Generator(i).Source,'Renewable')
        ren(end+1) = i;
    end
end
%% Step 2 Initialize Variables & prepare to iterate time steps
Data_t0 = GetCurrentData(DateSim);
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
    scaleCost = updateGeneratorCost(0);%% All costs were assumed to be 1 when building matrices, update Generator costs for the given time
    Renewable = zeros(1,length(Plant.Generator));
    for i = 1:1:length(Plant.Generator)
        if strcmp(Plant.Generator(i).Source,'Renewable')
            Renewable(i) = RenewableOutput(Plant.Generator(i).VariableStruct,DateSim,0,'Actual');
        end
    end
    IC = StepByStepDispatch(Data_t0.Demand,Renewable,scaleCost,Plant.optimoptions.Resolution,[],'',[]);
    for i=1:1:length(IC)
        if isfield(Plant.Generator(i).OpMatA,'Stor')
            IC(i) = 50; % IC = halfway charged energy storage
        end
    end
end
for i=1:1:length(IC)
    if isfield(Plant.Generator(i).OpMatA,'Stor')
        IC(i) = IC(i)/100*Plant.Generator(i).OpMatA.Stor.UsableSize/Plant.Generator(i).OpMatA.Stor.DischEff; % IC = halfway charged energy storage
    end
end
OnOff = logical(IC>=LB);
CurrentState.Generators=IC;

%% initialize all other variables
NumSteps = Plant.optimoptions.Interval*24/Plant.optimoptions.Resolution+1;
Dispatch =[];
Dispatch.Dispatch.Temperature = zeros(NumSteps,1);

Dispatch.Dispatch.Timestamp = zeros(NumSteps,1);
Dispatch.Dispatch.Timestamp(1) = Data_t0.Timestamp;

Dispatch.RunData.Timestamp = zeros(NumSteps,1);
Dispatch.RunData.Timestamp(1)  = Data_t0.Timestamp;

Dispatch.Dispatch.GeneratorState = zeros(NumSteps,nG+nL);
Dispatch.Dispatch.GeneratorState(1,:) = IC; %

Dispatch.RunData.GeneratorState = zeros(NumSteps,nG+nL);
Dispatch.RunData.GeneratorState(1,:) = IC;

Dispatch.Dispatch.GeneratorInput = zeros(NumSteps,nG+nL);

S = fieldnames(Data_t0.Demand);
for i = 1:1:length(S)
    Dispatch.Dispatch.Demand.(S{i}) = zeros(NumSteps,length(Data_t0.Demand.(S{i})));
    Dispatch.Dispatch.Demand.(S{i})(1,:) = Data_t0.Demand.(S{i});
    Dispatch.RunData.Demand.(S{i}) = zeros(NumSteps,length(Data_t0.Demand.(S{i})));
    Dispatch.RunData.Demand.(S{i})(1,:) = Data_t0.Demand.(S{i});
end

nS_step1 = round(Plant.optimoptions.Horizon/Time(1));
Last24hour =[];%eliminate any old data stored here
%need to have this in terms of the first timestep
Xi = nnz(Plant.Data.Timestamp<=DateSim);
Xf = nnz(Plant.Data.Timestamp<DateSim+1);
Last24hour = interpolateData(Plant.optimoptions.Resolution*3600,Xi,Xf,0.00);
Last24hour.Timestamp = DateSim-1+((0:nS_step1-1)'.*(Time(1)/24));
for i = 1:1:nG
    if isfield(Plant.Generator(i).VariableStruct,'RestartTime')
        RestartTime(i) = Plant.Generator(i).VariableStruct.RestartTime/60;%restart time in hours
    else RestartTime(i) = 0;
    end
end
GenAvailTime = ones(1,length(Plant.Generator)).*DateSim;

%% Step 3: Run initial dispatch and start optimization sequence (timers)
Si=1; %counter for # of times dispatch loop has run
while Virtual
    [Forecast,Renewable] = updateForecast(DateSim,Time);%% function that creates demand vector with time intervals coresponding to those selected
    %can expand forecast to create high and low estimates to be used in 2nd/3rd call of dispatch loop or step by step
    IC = CurrentState.Generators;
    IC(stor)= IC(stor).*DischEff;%scale storage value by discharge efficiency (this scaling is used in optimizations)
    if Si<=1
        PredictDispatch = ones(length(Time)+1,1)*IC;
    else PredictDispatch = [IC;Dispatch.Predicted(Si-1).GenDisp(3:end,:);Dispatch.Predicted(Si-1).GenDisp(end,:)];
    end
    [GenDisp,tsim] = DispatchLoop(IC,Forecast,Renewable,PredictDispatch);
    for i = 1:1:length(stor)
        GenDisp(:,stor(i))= GenDisp(:,stor(i))./DischEff(i);%un-scale storage value by discharge efficiency (this scaling is used in optimizations)
    end
    for i = 1:1:length(ren)
        GenDisp(2:end,ren(i))= Renewable(:,ren(i));%un-scale storage value by discharge efficiency (this scaling is used in optimizations)
    end
    Dispatch.Predicted(Si).GenDisp = GenDisp;
    if ~isempty(timers)
        timers(Si,3) = tsim;
    end
    if strcmp(Plant.optimptions.method,'Control')
        History = Dispatch.Actual.GenDisp(max(1,Si-Plant.optimoptions.Horizon/Plant.optimoptions.Resolution):Si-1,:);
    else
        if Si ==1
            History = Dispatch.Predicted(Si).GenDisp;
        else
            nS = length(Time);
            backSteps = min(Si-1,Plant.optimoptions.Horizon/Plant.optimoptions.Resolution);
            History = zeros(nS+1+backSteps,nG+nL);
            History(1:2,:) = Dispatch.Predicted(Si-backSteps).GenDisp(1:2,:);
            for t = 2:1:backSteps
                History(t+1,:) = Dispatch.Predicted(Si-backSteps+t-1).GenDisp(2,:);
            end
            History(backSteps+2:end,:) =Dispatch.Predicted(Si).GenDisp(2:end,:);
        end
    end
    plotDispatch(Dispatch.Predicted(Si).GenDisp,History,Time,datevec(DateSim))
    
    PlantD = datevec(DateSim);
    if strcmp(Plant.optimptions.method,'Planning')
        PlantD(4) = PlantD(4)+Plant.optimoptions.Horizon;%%count forward a day
    elseif strcmp(Plant.optimptions.method,'Dispatch')
        PlantD(4) = PlantD(4)+Plant.optimoptions.Resolution;%% count forward 1 step
    end
    DateSim = datenum(PlantD);  
    RealData = GetCurrentData(DateSim);
    if ~isempty(RealData)
        if strcmp(Plant.optimptions.method,'Planning')
            CurrentState.Generators = Dispatch.Predicted(Si).GenDisp(end,1:nG);
            Xi = nnz(Plant.Data.Timestamp<=DateSim);
            Xf = nnz(Plant.Data.Timestamp<DateSim+1);
            Last24hour = interpolateData(Plant.optimoptions.Resolution*3600,Xi,Xf,0.00);
            nS_step1 = round(Plant.optimoptions.Horizon/Time(1));
            Last24hour.Timestamp = DateSim-1+((0:nS_step1-1)'.*(Time(1)/24));
        elseif strcmp(Plant.optimptions.method,'Dispatch')   
            CurrentState.Generators = Dispatch.Predicted(Si).GenDisp(2,1:nG);
            Outs = fieldnames(Last24hour.Demand);
            Last24hour.Timestamp = [Last24hour.Timestamp(2:end,:);RealData.Timestamp;];
            Last24hour.Temperature = [Last24hour.Temperature(2:end,:);RealData.Temperature;];
            for i = 1:1:length(Outs)
                Last24hour.Demand.(Outs{i}) = [Last24hour.Demand.(Outs{i})(2:end,:);RealData.Demand.(Outs{i})(1,:);];
            end
        elseif strcmp(Plant.optimptions.method,'Control')
            %% Real-time control 
            if Plant.optimoptions.fastsimulation==0
                if isempty(timerfindall)
                    Timers(Plant.optimoptions)
                else %do nothing;
                end
            else
                %% Virtual Plant
                D = DateSim;
                while DateSim<(D+Time(1)/24)
                    OnlineLoop
                end
            end
        end
        Si=Si+1;
        if ~isempty(DispatchWaitbar)
            waitbar(Si/(Plant.optimoptions.Horizon/dt(1)),DispatchWaitbar,strcat('Running Dispatch'));
        end
    end
end
Input = Dispatch.Dispatch.GeneratorInput;
Outs = fieldnames(Last24hour);
for i = 1:1:length(IC)
    if strcmp(Plant.Generator(i).Type,'Chiller') %don't include cost, it shows up in generator demand
        Input(:,i) = 0;
    end
end
Time = (Dispatch.Dispatch.Timestamp(2:end)'-Dispatch.Dispatch.Timestamp(1));
Dispatch.NetCost = NetCostCalc(Dispatch.Dispatch.GeneratorState,Time,Input);

% Dispatch.Baseline = RunBaseline; %finish simulation by running baseline
Plant.Dispatch = Dispatch;