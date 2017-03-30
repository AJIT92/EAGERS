global Plant RealTimeData DataLog nLog UB LB% Generators Interval Resolution Horizon Tmpc Topt t_mpc spGEN  %passed from GUI
global RealTime Virtual scaleTime OptimalMap MapWaitbarHandle DateSim%relates to how optimization proceeds
global FanPortWrite %related to communication for Savona E-Hub
global CurrentState Dispatch Operation Si Last24hour OnOff t_mpc timersOn%initialized here
global GenAvailTime RestartTime ShutdownRamp spGEN%  Global vairables that need to be reset between each run, but not each loop. Trying to remove these global variables

%global capacity efficiency SSi%SSi is the state space model. It is pulled from loadGenerator
%Step 1: Step 1 Load generators, & build QP matrices
%Step 2: Initialize Variables & prepare to iterate time steps
%Step 3: Run Initial Dispatch & initial on-line optimization
%Step 4: Start all timers and comence the three loops simultaneously
%Step 5: perform offline plant optimization (provides baseline) & cost Analysis
%% variable descriptions
%RealTimeData: the high resolution data used for this test
%DataLog: storage of the operation and optimization as the simulation runs
%nLog: index in DataLog of current step
%UB: the upper limit (capacity) of each generator
%LB: the lower limit (capacity) of each generator when on

%RealTime: if connected to labview and running hardware
%Virtual: Running a simulation only, set to zero when the end of the test data set is reached
%scaletime: the ratio of emulated time to time in the test-data. For example 24 hours of test data can be run in 3 hours with a scaletime of 8. scaletime enlarges any energy storage and slows down the transient response of any generators
%OptimalMap: Avoids running feed-forwd optimization by using a pre-determined lookup table for optimal dispatch vs. load (doesn't speed things up much with only a few generators)
%MapWaitbarHandle: handle for creating the look-up table
%DateSim: Current time in the simulation.
%NumSteps: the number of dispatch optimiztions that will occur during the entire simulation

%CurrentState: Current state of all generators (kW) & storage (kWh) 
%Dispatch: recorded data at the frequency of the dispatch optimization
%Operation: data between subsequent dispatch optimizations. note that energy storage is specified in terms of output kW not capacity kWh
%Si: Counter for dispatch loop
%Last24hour: recorded data for the last 24 hours
%OnOff: the desired generator state from the controller (may differ from actual on/off because it leaves generators on as they ramp down to min power, then shuts them off)

%GenAvailTime: The next time a generator is available to turn on
%RestartTime: The amount of time necessary to wait for restarting after a generator has turned off .
%ShutdownRamp: if a generator has just been shut down (with OnOff not actualOnOff) this describes the residual power that will be produced during shutdown so that it can be subtracted from the optimization of the remaining generators
%SSmpc: the agregated state-space model seen by the MPC, the time may be scaled by both Tmpc and scaletime.
%GenSSindex: A matrix used to relate the primary and secondary states of each generator to their index in the states of the agregated state space model SSmpc. 

%% Step 1 Load generators, & build QP matrices
options = Plant.optimoptions;
scaleTime = options.scaletime;
Time = buildTimeVector(options);
if options.Topt/3600>Time(1) %prevent RunHMPC from miscounting time if Topt>the first timestep
    options.Topt = Time(1)*3600;
    Plant.optimoptions.Topt = Time(1)*3600;
end

LoadTestData %% Revise this so you can pull from more than what is loaded in Plant
%need to imterpolate this data as it is running to handle the variable timestep
Xi = nnz(Plant.Data.Timestamp<=DateSim);
Xf = nnz(Plant.Data.Timestamp<=DateSim+options.Interval);
RealTimeData = interpolateData(options.Tmpc,Xi,Xf,0.00);%create test data at correct frequency

timersOn = false;
tic
loadGenerator % Loads generators & build optimization matrices, stored in Plant.Generator
if timersOn
    disp(['time for loading components and building matrices is ' num2str(toc)])
end

Data_t0 = GetCurrentData(DateSim);

dt1 = Time(1);%the first timestep dt may varry from the resolution if you chose a linear or manual option
OptimalMap =[];
MapWaitbarHandle =[];

NumSteps = length(Time)+1;
nMPC = round(3600/options.Tmpc);%nMPC is the number of times MPC loop is run per hour
%% load/open comunication ports
if RealTime
    openPorts
    %initilaize DataLog
    nLog = 0;
    DataLog.Timestamp = zeros(NumSteps*nMPC,1);
    DataLog.FanTherm = zeros(NumSteps*nMPC,1);
    DataLog.AmbTemp = zeros(NumSteps*nMPC,1);
    DataLog.mGTpow = zeros(NumSteps*nMPC,1);
    DataLog.mGTheat = zeros(NumSteps*nMPC,1);
    DataLog.mGTstate = zeros(NumSteps*nMPC,1);
    DataLog.mGTfuel = zeros(NumSteps*nMPC,1);
    DataLog.ICEpow = zeros(NumSteps*nMPC,1);
    DataLog.ICEheat = zeros(NumSteps*nMPC,1);
    DataLog.ICEstate = zeros(NumSteps*nMPC,1);
    DataLog.ICEfuel = zeros(NumSteps*nMPC,1);
    DataLog.TES_SOC = zeros(NumSteps*nMPC,1); %kJ
    %DataLog.HVAC_SOC = zeros(Numsteps*nMPC,1); %would this be able to be
    %lumped in with the Thermal Energy Storage datalog?
    K =1;
else K = center_menu('Select Option','Manually specify initial conditions','Automatically determine initial conditions');
end
%% Step 2 Initialize Variables & prepare to iterate time steps
nG = length(Plant.Generator);
IC = zeros(1,nG);
include = {'CHP Generator', 'Electric Generator', 'Chiller','Heater'};
stor = zeros(1,nG);
Gen = zeros(1,nG);
for i=1:1:nG
    if isfield(Plant.Generator(i).OpMatA,'Stor')
        stor(i) = i;
    elseif ismember(cellstr(Plant.Generator(i).Type),include)
        Gen(i) = i;
    end
end

if K ==1
    list = {};
    sizes = {};
    Index =[];
    for i = 1:1:length(Plant.Generator)
        if stor(i)>0
            list(end+1) = cellstr(strcat(Plant.Generator(i).Name,' --- Max IC is:',num2str(UB(i))));
            sizes(end+1) = cellstr(num2str(UB(i)/2));
            Index(end+1) = i;
        elseif Gen(i)>0
            list(end+1) = cellstr(strcat(Plant.Generator(i).Name,' --- If greater than 0, IC must be between lower bound(',num2str(LB(i)),') and upper bound(',num2str(UB(i)),').'));
            sizes(end+1) = cellstr('0');
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
    IC = StepByStepDispatch(Data_t0.Demand,scaleCost,options.Resolution,[],'',[]);
    IC(stor>0) = .5*UB(stor>0); % IC = halfway charged energy storage
end

OnOff = logical(IC>=LB);
CurrentState.Generators=IC;
if RealTime %send initial conditions
    SetGeneratorOnOff(OnOff>0,OnOff==0);
    SetGeneratorOutput(IC);
    %specific to running Genoa e-Hub
    fwrite(FanPortWrite,num2str(Data_t0.Demand.H),'char');%%SET fans
end

%% initialize all other variables
Dispatch =[];
Dispatch.Dispatch.GeneratorState = zeros(NumSteps,length(UB));
Dispatch.RunData.Timestamp = zeros(options.Interval*24*nMPC+1,1);
Dispatch.RunData.GeneratorState = zeros(options.Interval*24*nMPC+1,length(UB));
S = fieldnames(Data_t0.Demand);
for i = 1:1:length(S)
    Dispatch.Dispatch.Demand.(S{i}) = zeros(NumSteps,1);
    Dispatch.Dispatch.Demand.(S{i})(1) = Data_t0.Demand.(S{i});
    Dispatch.RunData.Demand.(S{i}) = zeros(options.Interval*24*nMPC+1,1);
    Dispatch.RunData.Demand.(S{i})(1) = Data_t0.Demand.(S{i});
end
%make the heating demand smoothed by adding the demands for the heat together here.
if options.nsSmooth>0 && isstruct(Dispatch.Dispatch.Demand.H)
    for i = 2:1:length(Dispatch.Dispatch.Demand.H)-options.nsSmooth
        Dispatch.Dispatch.Demand.H(i) = sum(Dispatch.Dispatch.Demand.H(i:i+options.nsSmooth));
    end
    for i = length(Dispatch.Dispatch.Demand.H)-options.nsSmooth+1:length(Dispatch.Dispatch.Demand.H)
        Dispatch.Dispatch.Demand.H(i)= sum(Dispatch.Dispatch.Demand.H(i:end));
    end
end

Dispatch.Dispatch.GeneratorInput = zeros(NumSteps,length(UB));
Dispatch.Dispatch.Timestamp = zeros(NumSteps,1);
Dispatch.Dispatch.Temperature = zeros(NumSteps,1);
Dispatch.Dispatch.Timestamp(1) = Data_t0.Timestamp;
Dispatch.RunData.Timestamp(1)  = Data_t0.Timestamp;
Dispatch.Dispatch.GeneratorState(1,:) = IC; %
Dispatch.RunData.GeneratorState(1,:) = IC;
Dispatch.Predicted = Dispatch.Dispatch;
Dispatch.Predicted.Timestamp = linspace(0,options.Horizon,options.Horizon/dt1+1)./24+DateSim;
Dispatch.Predicted.GenDisp = zeros(length(Dispatch.Predicted.Timestamp),length(IC));
Dispatch.Predicted.GenDisp(1,:) = IC;

%operation is always going to be the current timestep, and the current
%timestep may not have nMPC MPC loops, so calculate the nMPC_step1
nMPC_step1 = floor(dt1*3600/options.Tmpc); %this is the same as nMPC for constant one hour timesteps
Operation.Timestamp = zeros(nMPC_step1,1);
Operation.Temperature = zeros(nMPC_step1,1);
Operation.Cogen = zeros(nMPC_step1,length(Plant.Generator));
Operation.GeneratorState = zeros(nMPC_step1,length(Plant.Generator));
Operation.Demand = [];
for i = 1:1:length(S)
    Operation.Demand.(S{i}) = zeros(1,nMPC_step1);
end

nS_step1 = round(options.Horizon/Time(1));
Last24hour =[];%eliminate any old data stored here
%need to have this in terms of the first timestep
Xi = nnz(Plant.Data.Timestamp<=DateSim);
Xf = nnz(Plant.Data.Timestamp<DateSim+1);
Last24hour = interpolateData(options.Resolution*3600,Xi,Xf,0.00);
Last24hour.Timestamp = DateSim-1+((0:nS_step1-1).*(Time(1)/24));
for i = 1:1:length(LB)
    if isfield(Plant.Generator(i).VariableStruct,'RestartTime')
        RestartTime(i) = Plant.Generator(i).VariableStruct.RestartTime/60;%restart time in hours
    else RestartTime(i) = 0;
    end
end
GenAvailTime = ones(1,length(Plant.Generator)).*DateSim;

spGEN = [];%this is the setpoint for generators in MPCloop
ShutdownRamp = [];


%% Step 3: Run initial dispatch and optimization
Si=1; %counter for # of times dispatch loop has run
DispatchLoop
t_mpc=1;
%% Step 4: set up timers for dispatch loop, optimization loop and MPC loop
%(also write thermal load timer)
% Create and start a timer
if RealTime || Plant.optimoptions.fastsimulation==0 %run slowly with timers so that MPC loop runs simultaneously as dispatch is recalculated
    mpcTimer     = timer('TimerFcn'  ,@MPCloop, ...
                       'StartDelay'   , 0            , ...
                       'Period'       , options.Tmpc/options.scaleTime, ...
                       'Name'         ,'mpcTimer' , ...
                       'ExecutionMode','fixedrate'        );
    optTimer     = timer('TimerFcn'  ,@OnlineLoop, ...
                       'StartDelay'   , options.Topt/options.scaleTime, ...
                       'Period'       , options.Topt/options.scaleTime, ...
                       'Name'         ,'optTimer' , ...
                       'ExecutionMode','fixedrate'        );
    dispTimer     = timer('TimerFcn'  ,@DispatchLoop, ...
                       'StartDelay'   , options.Resolution*3600/options.scaleTime, ...
                       'Period'       , options.Resolution*3600/options.scaleTime, ...
                       'Name'         ,'dispTimer' , ...
                       'ExecutionMode','fixedrate'        );
    fanTimer     = timer('TimerFcn'  ,@writeThermalLoad, ...
                       'StartDelay'   , 0            , ...
                       'Period'       , options.Tmpc/options.scaleTime, ...
                       'Name'         ,'fanTimer' , ...
                       'ExecutionMode','fixedrate'        );
   %% start all timers
    start(mpcTimer);
    start(optTimer);
    start(dispTimer);
    if RealTime
        start(fanTimer);
    end
else %run quickly with dispatch calculated then MPC loop called
    while Virtual
        t1 = 0;
        while t1<dt1*3600%(3600*options.Resolution)
            t2 = 0;
            while t2<options.Topt
                tic
                MPCloop
                if t2==0 && Si==2 && timersOn
                    disp(['time for first MPC loop is ' num2str(toc)])
                end
                t1 = t1+options.Tmpc;
                t2 = t2+options.Tmpc;
            end
            if Virtual && t1<dt1*3600%(3600*options.Resolution)
                tic
                OnlineLoop
                if t1==options.Tmpc && Si==2 && timersOn
                    disp(['time for first Online loop is ' num2str(toc)])
                end
            end
        end
        if Virtual
            tic
            DispatchLoop
            if Si==2 && timersOn
                disp(['time for first complete dispatch loop is ' num2str(toc)])
            end
        end
    end
end
Input = Dispatch.Dispatch.GeneratorInput;
Outs = fieldnames(Last24hour);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Chiller') && ~ismember('E',Outs)%don't include cost if it shows up in generator demand
        Input(:,i) = 0;
    end
end
Time = (Dispatch.Dispatch.Timestamp(2:end)'-Dispatch.Dispatch.Timestamp(1));
Dispatch.NetCost = NetCostCalc(Dispatch.Dispatch.GeneratorState,Time,Input);

Dispatch.Baseline = RunBaseline; %finish simulation by running baseline
Plant.Dispatch = Dispatch;