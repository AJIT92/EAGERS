function GenOutput = StepByStepDispatch(Demand,scaleCost,dt,IC,limit,StorageProfile)
% Time is the time from the current simulation time (DateSim), positive numbers are forward looking
global Plant UB LB
nG = length(Plant.Generator);
[nS,~] = size(scaleCost);
TestCombos = zeros(nS,1);

if ~isempty(IC)
    GenOutput = zeros(nS+1, nG);
    GenOutput(1,:) = IC;
else GenOutput = zeros(nS, nG);%nS should equal 1 in this case (finding IC)
end

if isempty(StorageProfile)
    StorageProfile = UB;
end

%% note: do not need to account for self discharging of storage. This was added to demand when doing initial GenDisp, and now we are using the results of that
Outs = Plant.optimoptions.Outputs;
if nnz(strcmp('H',Outs))>0 && nnz(strcmp('E',Outs))>0
    Outs = Outs(~(strcmp('H',Outs))); %combine heating optimization with electric due to CHP generators
end
if nnz(strcmp('C',Outs))>0 && nnz(strcmp('E',Outs))>0 && Plant.optimoptions.sequential == 0
    Outs = Outs(~(strcmp('C',Outs))); %combine cooling optimization with electric 
end
for t = 1:1:nS %for every timestep
    marginal = instantMarginalCost(StorageProfile(t,:),scaleCost(t,:));%update marginal cost
    Organize = Plant.OneStep.Organize;%indices (rows and columns) associated with each generator, allowing generators to be removed later
    EC = zeros(1,nG);
    for s = 1:1:length(Outs);%for each demand (electrical/heating or chilling)
        QP = Plant.OneStep.(Outs{s})(t); % quadratic programing matrices (H, f, Aeq, beq, A, b, ub, lb)
        QP = updateMatrices1Step(QP,scaleCost(t,:),marginal,IC,dt(t),Organize.(Outs{s}));
        thisSeq = Organize.(Outs{s}).thisSeq;
        stor = Organize.(Outs{s}).stor;
        storC = Organize.(Outs{s}).storC;
        storH = Organize.(Outs{s}).storH;
        utility = Organize.(Outs{s}).utility;
        utilC = Organize.(Outs{s}).utilC;
        utilH = Organize.(Outs{s}).utilH;
        chill = Organize.(Outs{s}).chill;
        heater = Organize.(Outs{s}).heater;
        allStor = [stor, storC, storH];
        allGen = [thisSeq, chill, heater];
        allUtility = [utility, utilC, utilH];
        
        mat = Organize.(Outs{s}).Demand{1};
        index = Organize.(Outs{s}).Demand{2};
        QP.(mat)(index) =Demand.(Outs{s})(t); %the first row of the equivalent matrix is the demand
        if strcmp(Outs{s},'E') && isfield(Demand,'H')
            Heating.Demand = Demand.H(t);
            Heating.stor = storH;
            Heating.gen = heater;
            Heating.util = utilH;
            Heating.chp = Organize.(Outs{s}).CHPindex;
            Heating.Hratio = Organize.(Outs{s}).Hratio;
            Heating.Organize = Organize.('H');
        else Heating = [];
        end
        if strcmp(Outs{s},'E') && isfield(Demand,'C') && Plant.optimoptions.sequential == 0
            Cooling.Demand = Demand.C(t);
            Cooling.stor = storC;
            Cooling.gen = chill;
            Cooling.util = utilC;
            Cooling.Organize = Organize.('C');
        else Cooling = [];
        end
        for j = 1:1:length(allUtility)
            i = allUtility(j);
            n_states = Organize.(Outs{s}).States{i};
            n_inf = isinf(QP.ub(n_states));
            if nnz(n_inf)>0 %remove infinity bounds
                QP.ub(n_states(n_inf>0)) = max(5*Demand.(Outs{s})(t),5*sum(UB(~isinf(UB)))); %upper bound of utility is 5 times demand or sum of other generation
            end
        end
        if strcmp(limit,'initially constrained')
            IC(allStor) = StorageProfile(t,allStor); %refer to initial dispatch for SOC of storage to calculate marginal value
        end
        MaxPower = zeros(1,nG);
        MinPower = zeros(1,nG);
        constCost = zeros(1,nG);
        for j = 1:1:length(allGen)
            i = allGen(j);
            if Plant.Generator(i).Enabled ==0
                thisSeq = thisSeq(thisSeq~=i);
                chill = chill(chill~=i);
                heater = heater(heater~=i);
            else
                if isempty(IC)
                    MaxPower(i) = UB(i); 
                    MinPower(i) = 0; 
                else
                    MaxPower(i) = min(UB(i),IC(i) + Plant.Generator(i).OpMatB.Ramp.b(1)*sum(dt(1:t)));
                    MinPower(i) = max(0,IC(i) - Plant.Generator(i).OpMatB.Ramp.b(2)*sum(dt(1:t)));
                end
                if isfield(Plant.Generator(i).OpMatB,'constCost') %all of these cost terms need to be scaled later on
                    constCost(i) = Plant.Generator(i).OpMatB.constCost.*scaleCost(t,i)*dt(t);
                end
            end
        end
        for j = 1:1:length(allUtility)
            i = allUtility(j);
            MaxPower(i) = UB(i); 
            MinPower(i) = LB(i); 
        end
        for j = 1:1:length(allStor)
            i = allStor(j);
            if isempty(IC)
                MaxPower(i) = Plant.Generator(i).OpMatB.Ramp.b(2);
                MinPower(i) = 0;%-Plant.Generator(i).OpMatB.Ramp.b(1); 
            else
                MaxPower(i) = min(IC(i)/dt(t),Plant.Generator(i).OpMatB.Ramp.b(2)); %minimum of peak discharge, and SOC/time (completely discharging storage in next step)
                MinPower(i) = max(-(UB(i)-IC(i))/dt(t),-Plant.Generator(i).OpMatB.Ramp.b(1)); % maximum of peak charging and (UB - SOC)/dt (completely charging in next time step
            end
        end

        [BestDispatch,TestCombos(t)] = eliminateCombinations(QP,Organize.(Outs{s}),Heating,Cooling,constCost,dt(t),IC,MaxPower,MinPower,thisSeq,stor,utility);%% combination elimination loop
        if isempty(BestDispatch)
            disp(['No feasible combination of generators at step' num2str(t)]);
            BestDispatch = IC;
            BestDispatch(allStor) = 0;
            TestCombos(t)=0;
        end
        EC([allGen, allUtility])= BestDispatch([allGen, allUtility]);
        if ~isempty(IC)
            EC(allStor) = IC(allStor) - BestDispatch(allStor)*dt(t); %change in storage for this power output
        end
    end
    if strcmp(limit, 'constrained')%if its constrained but not initially constrained then make the last output the initial condition
        IC = EC;
    end
    GenOutput(t+(~isempty(IC)),:) = EC;
end
% disp(['Average # of combinations tested is ' num2str(sum(TestCombos)/nS)]);