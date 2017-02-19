function MPCloop(varargin)
global Plant RealTime Virtual scaleTime DateSim %supplied by GUI
global Operation  t_mpc  spGEN CurrentState OnOff GendX GenAvailTime RestartTime %determined in this loop
global SSmpc GenSSindex dischEff chargeEff selfDisch UB%loaded in load generators
global UBmpc LBmpc setGen ShutdownRamp MinThresh Threshold %determined by dispatch & online optimization

%Step 1: Identify the generator categories
%Step 2: load the current and previous states into the MPC state model, also load the derivatives (dX)
%Step 3: Realtime (collect actual state of generators), Virtual (use the MPC states to simulate response of generators 1 step forward)
%Step 4: Balance any instantaneous error between supply & demand using any grid or storage available 
%Step 5: Determine if any generators should turn on or off at this moment. If some turn off, set-up the ramp-down behavior until they are completely off
%step 6: Determine the restart time for any generators that have just turned off
%step 7: actually turn on/off generators determined in step 5
%step 8: enforce any grid constraints
%step 9: determine new set points for the dispatchable generators using the MPC
%step 10: Send these set points to the generators

%% variable descriptions
%RealTime: if connected to labview and running hardware
%Virtual: Running a simulation only, set to zero when the end of the test data set is reached
%scaletime: the ratio of emulated time to time in the test-data. For example 24 hours of test data can be run in 3 hours with a scaletime of 8. scaletime enlarges any energy storage and slows down the transient response of any generators
%DateSim: Current time in the simulation.

%Operation: data between subsequent dispatch optimiztions. note that energy storage is specified in terms of output kW not capacity kWh
%t_mpc: time step of MPC loop, overwites Operation each time the dispatch runs again
%spGEN: The current setting of each generator as determined by the MPC.
%CurrentState: Current state of all generators (kW) & storage (kWh) 
%OnOff: the desired generator state from the controller (may differ from actual on/off because it leaves generators on as they ramp down to min power, then shuts them off)
%GendX: The derivative of the generator state 
%GenAvailTime: The next time a generator is available to turn on, this is
%initialized in RunHMPC
%RestartTime: The amount of time necessary (in hours) to wait for restarting after a generator has turned off 

%SSmpc: the agregated state-space model seen by the MPC, the time may be scaled by both Tmpc and scaletime.
%GenSSindex: A matrix used to relate the primary and secondary states of each generator to their index in the states of the agregated state space model SSmpc. 
%dischEff: Discharge efficiency of energy storage
%chargeEff: charge efficiency of energy storage
%selfDisch: the amount of self dicharging per hour for a storage system

%UBmpc: upper bound of generators seen by MPC, can be different from the global upper bound due to ramping constraints.
%LBmpc: lower bound of generators seen by MPC, can be different from the global lower bound due to ramping constraints
%setGen the current setting of each generator as determined by the on-line dispatch.
%PrevDemand: The demand the last time MPC loop ran.
%ShutdownRamp: if a generator has just been shut down (with OnOff not actualOnOff) this describes the residual power that will be produced during shutdown so that it can be subtracted from the optimization of the remaining generators
%MinThresh: The minimum threshold of power purchased from the grid.
GetCurrentTime
RealData = GetCurrentData(DateSim);
Tmpc = Plant.optimoptions.Tmpc;
if isempty(spGEN)
    spGEN = setGen;%used the first time MPCloop runs
end
Operation.Timestamp(t_mpc) = DateSim;
Operation.Temperature(t_mpc) = RealData.Temperature;
S = fieldnames(RealData.Demand);
for i = 1:1:length(S)
    Operation.Demand.(S{i})(t_mpc) = RealData.Demand.(S{i});
end
%% set-up for MPC, weights
r = size(SSmpc.A);
s = size(SSmpc.C);
MPCweights = zeros(s(1)+r(1),s(1)+r(1));
MPCweights(r(1)+1:end,r(1)+1:end) =eye(s(1));
J = GenSSindex.Secondary(:,2);
MPCweights(r(1)+J,r(1)+J) = 0; %%currently do nothing with secondary states (CHP)% dont use MPC to control secondary setpoint, can change this later


nG = length(Plant.Generator);
Hratio = zeros(1,nG);
OpMats = fieldnames(Plant.OpMatB.QP);
Organize = Plant.OpMatB.Organize;
for s = 1:1:length(OpMats) %not currently set up well for 2 optimization matrices (chillers done sequentially)
    thisSeq = Organize.(OpMats{s}).thisSeq;
    stor = Organize.(OpMats{s}).stor;
    storC = Organize.(OpMats{s}).storC;
    storH = Organize.(OpMats{s}).storH;
    utility = Organize.(OpMats{s}).utility;
    utilC = Organize.(OpMats{s}).utilC;
    utilH = Organize.(OpMats{s}).utilH;
    chill = Organize.(OpMats{s}).chill;
    heater = Organize.(OpMats{s}).heater;
    renew = Organize.(OpMats{s}).renew;
    if strcmp(OpMats{s},'E')
        for j = 1:1:length(thisSeq)
            i = thisSeq(j);
            if isfield(Plant.Generator(i).OpMatB.output,'H')
                Hratio(i) = Plant.Generator(i).OpMatB.output.H;
            end
        end
    end
    allStor = sort([stor, storC, storH]);
    allGen = sort([thisSeq, chill, heater]);
    allUtility = [utility, utilC, utilH];
    CHPindex = nonzeros((1:nG).*(Hratio>0))';
    if ~isempty(renew)
        [RenPower,~] = RenewableOutput(DateSim,0,'Actual');
        Operation.GeneratorState(t_mpc,renew) = RenPower;
    end
    if ~isempty(chill) || ~isempty(utilC) || ~isempty(storC)
        DemandC = RealData.Demand.C;
    else DemandC = [];
    end
    if ~isempty(heater) || ~isempty(utilH) || ~isempty(storH) || ~isempty(CHPindex)
        DemandH = RealData.Demand.H;
    else DemandH = [];
    end
    %% set-up for MPC initial cond %update the previous state (and previous derivative) in the SS models
    X = zeros(r(1),2);
    X(GenSSindex.Primary(:,1),1) = CurrentState.Generators(allGen);
    X(GenSSindex.Primary(:,1)+1,1) = GendX(allGen); %derivatives of the generator state
    X(GenSSindex.Secondary(:,1),1) = CurrentState.Generators(CHPindex).*Hratio(CHPindex); %CHP is currently the only secondary output
    X(GenSSindex.Secondary(:,1)+1,1) = GendX(CHPindex).*Hratio(CHPindex); %CHP is currently the only secondary output

    if RealTime %get measured outputs (current state)
        %% %%% Get curent Data
        [CurrentState.Generators,Operation.Cogen(t_mpc,:),ActualOnOff,Operation.GeneratorInput(t_mpc,:)] = GetCurrentState;
        CurrentState.Generators(allGen) = max(0,CurrentState.Generators(allGen)); %eliminate negative readings on ICE & mGT
        CurrentState.Generators(allStor) = CurrentState.Generators(allStor)/(3600)*scaleTime; %convert kJ to kWh & scale storage  
        Operation.GeneratorState(t_mpc,allGen) = CurrentState.Generators(allGen);     
        %% --- Determine grid & energy storage balance
        Demand = RealData.Demand.(OpMats{s}) - sum(CurrentState.Generators(thisSeq));
        Operation.GeneratorState(t_mpc,[stor,utility]) = GRIDandStorage(Demand,setGen,stor,utility);
        if ~isempty(DemandC)
            DemandC = DemandC - sum(CurrentState.Generators(chill));
            Operation.GeneratorState(t_mpc,[storC,utilC]) = GRIDandStorage(DemandC,setGen,storC,utilC);
        end
        if ~isempty(DemandH)
            DemandH = DemandH - sum(CurrentState.Generators(heater)) - sum(CurrentState.Generators(CHPindex).*Hratio(CHPindex));
            Operation.GeneratorState(t_mpc,[storH,utilH]) = GRIDandStorage(DemandH,setGen,storH,utilH);
        end
        Operation.GeneratorInput(t_mpc,allUtility) = Operation.GeneratorState(t_mpc,allUtility);

        X(GenSSindex.Primary(:,1),2) = CurrentState.Generators(allGen)';
        X(GenSSindex.Primary(:,1)+1,2) = (CurrentState.Generators(allGen)'-X(GenSSindex.Primary(:,1),1))./(Tmpc/scaleTime); %derivative
        X(GenSSindex.Secondary(:,1),2) = Operation.Cogen(t_mpc,CHPindex)';
        X(GenSSindex.Secondary(:,1)+1,2) = (Operation.Cogen(t_mpc,CHPindex)'-X(GenSSindex.Secondary(:,1),1))./(Tmpc/scaleTime);%derivative
        Y(GenSSindex.Primary(:,2),1)=Operation.GeneratorState(t_mpc,allGen)';
        Y(GenSSindex.Secondary(:,2),1) = Operation.Cogen(t_mpc,CHPindex)';
    else
        %simulate generator output (current state)
        newX = SSmpc.A*X(:,1)+SSmpc.B*spGEN(allGen); %generator states
        X(:,2) = X(:,1)+(newX-X(:,1))/scaleTime;
        Y(:,1) = SSmpc.C*X(:,2); 
        z = Y(:,1)>1e-10;
        Y(:,1) = Y(:,1).*z; %eliminate numerical errors
        Operation.GeneratorState(t_mpc,allGen) = Y(GenSSindex.Primary(:,2),1)';%Primary outputs: Y may have additiona states, for example CHP has both electric and heat state
        Operation.Cogen(t_mpc,CHPindex) = Y(nonzeros(GenSSindex.Secondary(:,2)),1)';
        for i = 1:1:nG
            if ~isempty(Plant.Generator(i).Output)
                cap = Plant.Generator(i).Output.Capacity*UB(i);
            end
            eff = [];
            if strcmp(Plant.Generator(i).Type,'Electric Generator') || strcmp(Plant.Generator(i).Type,'CHP Generator')
                eff = Plant.Generator(i).Output.Electricity;
            elseif strcmp(Plant.Generator(i).Type,'Chiller')
                eff = Plant.Generator(i).Output.Cooling;
            elseif strcmp(Plant.Generator(i).Type,'Heater')
                eff = Plant.Generator(i).Output.Heat;    
            end
            if ~isempty(eff)
                Operation.GeneratorInput(t_mpc,i) = Operation.GeneratorState(t_mpc,i)/interp1(cap,eff,Operation.GeneratorState(t_mpc,i));
                if isnan(Operation.GeneratorInput(t_mpc,i))
                    Operation.GeneratorInput(t_mpc,i) = 0;
                end
            end
        end
        %% --- Determine grid & energy storage balance
        Demand = RealData.Demand.(OpMats{s}) - sum(CurrentState.Generators(thisSeq));
        Operation.GeneratorState(t_mpc,[stor,utility]) = GRIDandStorage(Demand,setGen,stor,utility);
        if ~isempty(DemandC)
            DemandC = DemandC - sum(CurrentState.Generators(chill));
            Operation.GeneratorState(t_mpc,[storC,utilC]) = GRIDandStorage(DemandC,setGen,storC,utilC);
        end
        if ~isempty(DemandH)
            DemandH = DemandH - sum(CurrentState.Generators(heater)) - sum(CurrentState.Generators(CHPindex).*Hratio(CHPindex));
            Operation.GeneratorState(t_mpc,[storH,utilH]) = GRIDandStorage(DemandH,setGen,storH,utilH);
        end
        Operation.GeneratorInput(t_mpc,allUtility) = Operation.GeneratorState(t_mpc,allUtility);
        CurrentState.Generators([allGen,allUtility]) = Operation.GeneratorState(t_mpc,[allGen,allUtility]);
        for j = 1:1:length(allStor) 
            i = allStor(j);
            loss = ((selfDisch(i)*UB(i))/(dischEff(i)*chargeEff(i)))*(Tmpc/3600)/scaleTime; %energy loss (power * time/ scaleTime)
            if Operation.GeneratorState(t_mpc,i)>0 %discharging (factor in discharging efficiency)
                Charge = -Operation.GeneratorState(t_mpc,i)*(1/dischEff(i))*(Tmpc/3600)/scaleTime; %energy leaving storage (power/eff * time/ scaleTime)
            else
                Charge = -Operation.GeneratorState(t_mpc,i)*(chargeEff(i))*(Tmpc/3600)/scaleTime; %energy leaving storage (power*eff * time/ scaleTime)
            end
            CurrentState.Generators(i) = CurrentState.Generators(i)-loss+Charge;
        end
    end
end
[TurnOn, TurnOff] = ThresholdEnforce(RealData.Demand);

if ~isempty(TurnOn)|| ~isempty(TurnOff)
    OnOff(TurnOn) = 1;
    OnOff(TurnOff) = 0;
    if ~isempty(TurnOff) %% (!!assume shut-down is 2x as fast as start-up)
         Threshold.E.t(TurnOff) = inf; %ensure generator does not turn back on imediately
         GenAvailTime(TurnOff) = DateSim+RestartTime(TurnOff)/24*scaleTime;%if you are turning it off, it is available later
         for k = 1:1:length(TurnOff)
            ramp_rate = Plant.Generator(TurnOff(k)).OpMatB.Ramp.b(1)*2; %ramp rate, (!!assume shut-down is 2x as fast as start-up)
            toff = (CurrentState.Generators(TurnOff(k))-0)./ramp_rate;%time to shutdown in hours
            toff = ceil(toff*3600); %rounded time to shutdown in seconds
            ShutdownRamp(TurnOff(k)).t = linspace(DateSim,DateSim+toff/3600/24,(toff/Tmpc)+1);
            ShutdownRamp(TurnOff(k)).Pow = linspace(CurrentState.Generators(TurnOff(k)),0,(toff/Tmpc)+1);
            r = round(Plant.optimoptions.Topt/Tmpc);
            ShutdownRamp(TurnOff(k)).t2 = [ShutdownRamp(TurnOff(k)).t(1:r:end) ShutdownRamp(TurnOff(k)).t(end)];
            ShutdownRamp(TurnOff(k)).Pow2 = [ShutdownRamp(TurnOff(k)).Pow(1:r:end) ShutdownRamp(TurnOff(k)).Pow(end)];
         end
    end
    Timestamp = zeros(1,length(Plant.Online));
    for t = 1:1:length(Plant.Online)
        Timestamp(t) = Plant.Online(t).Timestamp;
    end
    if nnz(DateSim<Timestamp)>0%if you've run the onlineloop for this dispatch loop already, re-run it
        OnlineLoop %re-run with new generator config
    else
        disp('Should never get here because dispatch loop should re-run: line 188 MPC loop. Maybe during realtime')
%         Plant.Online.Timestamp = Plant.Online.Timestamp+Plant.optimoptions.Resolution/24;%if you are running OnlineLoop for the first time during this timestep, re -run it after updating the timestamp
%         OnlineLoop
    end
end

for j = 1:1:length(allGen)
    i = allGen(j);
    if CurrentState.Generators(i)<LBmpc(i)
        if CurrentState.Generators(i)>(.05*LBmpc(i)) && OnOff(i)==0
            GenAvailTime(i) = DateSim+RestartTime(i)/24*scaleTime;% reset time for restart
        end
    end
end

if RealTime %actually turn on/off generators 
    TurnOffNow = [];
    TurnOnNow = [];
    for j = 1:1:length(allGen)
        i = allGen(j);
        if OnOff(i)>0 && ActualOnOff(i)==0 && DateSim>=GenAvailTime(i)
            TurnOnNow(end+1) = i;
        elseif OnOff(i)==0 && ActualOnOff(i)>0 && CurrentState.Generators(i)<=1.05*LBmpc(i)
            TurnOffNow(end+1) = i;
        end
    end %where generators actually turn on/off (must reach within 5% of LB during shutdown)
    if ~isempty(TurnOnNow)|| ~isempty(TurnOffNow)
        SetGeneratorOnOff(TurnOnNow,TurnOffNow);
    end
end
GendX(allGen) = X(GenSSindex.Primary(:,1)+1,2);

%% correct set-point if grid import threshold is exceeded
if ~isempty(MinThresh) && Operation.GeneratorState(t_mpc,(MinThresh>0))<MinThresh(MinThresh>0)
    genEactive =[];
    nEnow =[];
    for i = 1:1:nG
        if OnOff(i)>0 && isfield(Plant.Generator(i).OpMatA.output,'E')
            genEactive(end+1) = i;
            nEnow(end+1) = i;
        end
    end
    import = sum(MinThresh(MinThresh>0) - Operation.GeneratorState(t_mpc,(MinThresh>0)));
    while sum(setGen(nEnow))>sum(LBmpc(genEactive)) && import<0%change the generator set points to get back above this threshold
        p = setGen(nEnow)-LBmpc(genEactive)';
        np = nnz(p);
        lessPowerEach = (MinThresh-import)/np;
        a = min(p,lessPowerEach);
        setGen(nEnow) = setGen(nEnow)-a;
        import = import+sum(a);
    end
end

setGen = max(setGen,LBmpc');

%% use MPC to find new generator set-points
spMPC = zeros(length(GenSSindex.Primary(:,2))+length(GenSSindex.Secondary(:,2)),1);
spMPC(GenSSindex.Primary(:,2),1) = 0; 
spMPC(GenSSindex.Primary(:,2),1) = setGen(allGen);
spMPC(GenSSindex.Secondary(:,2),1) = Hratio(CHPindex)'.*setGen(CHPindex);%Y(J ,max(1,t-1));
Xf = [X(:,2)-X(:,1); Y(:,1)-spMPC(:,1);]; % Xf = [x(t)-x(t-1); y(t)- SP]
[spGEN(allGen)]=MPCsolve(spGEN(allGen),SSmpc.A,SSmpc.B,SSmpc.C,MPCweights,30,Xf);
% if nnz(spGEN<min(LBmpc',zeros(size(spGEN))))>0
%     disp('warning set points are too low in mpc loop, resetting to LB')
%     spGEN(spGEN<LBmpc') = LBmpc(spGEN<LBmpc');
% end

for j = 1:1:length(allGen)
    i = allGen(j);
    if OnOff(i)>0 && UBmpc(i)<=LBmpc(i)
        spGEN(i) = LBmpc(i); %force set-point to be lower bound as generator is brought on-line
    end
    if OnOff(i)==0 && CurrentState.Generators(i)>LBmpc(i)
        if ~isempty(ShutdownRamp)
            n = nnz(ShutdownRamp(i).t<=(DateSim+Tmpc/3600/24));
            spGEN(i) = ShutdownRamp(i).Pow(n);%change spgen if controlled shutdown is occuring
        end
    end
end
t_mpc = t_mpc+1;
DispLoopdT = buildTimeVector(Plant.optimoptions);%this is the first time step's dt in the dispatch loop
if DateSim>=Plant.Online(end).Timestamp || t_mpc>floor(DispLoopdT(1)*3600/Tmpc)%if you have completed one time interval, then start back at one for the next timestep
    t_mpc = 1;
end

if RealTime %send command setpoints
    SetGeneratorOutput(spGEN);
end
if ~RealTime && ~Virtual
    recordFromMPCloop %record final operation because dispatchloop will not run again
end