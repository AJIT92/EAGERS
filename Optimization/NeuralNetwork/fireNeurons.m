  function LockedNet = fireNeurons(FirstDisp,GenDemand,marginCost,scaleCost)
global Plant timersOn dischEff UB LB Si DateSim CurrentState
global lockedNet lockederror net %these are created in this function and used in each subsequent step
global Last24hour
%input: inputs for neuron dispatch and network training
%output: dispatch from network
Time = buildTimeVector(Plant.optimoptions);%% set up dt vector of time interval length
dt = Time - [0, Time(1:end-1)];
nG = length(Plant.Generator);
nS = length(Time);
allgens = 1:nG;
stor = allgens(dischEff>0);
gens = allgens(and(UB~=inf,dischEff==0));
gensandstor = [gens,stor];   
ngs = length(gensandstor);
ngens = length(gens);
Out = Plant.optimoptions.Outputs;
renews = []; %number of renewables
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Source, 'Renewable')
        renews(end+1) = i;
    end
end

if Si==2
    training = true;
else training = false;
end
%% create training data for 7 days at different start times
if training %train the first time through
    Tdays = min(100, floor(length(Plant.Data.Demand.E)/(24*4)));%number of training data sets (Training days). You can only train for the number of days you have data for
    DateSim1 = DateSim;
    histdataTime = datevec(Plant.Data.Timestamp(1));
    DateSimVec = datevec(DateSim);
    DateSim = datenum([histdataTime(1), DateSimVec(2:5), 0]);%make sure you are training over the same year you have data for, don't include seconds
    Last24hour1 = Last24hour;
    tic
    GenDisp = zeros(nS*Tdays,nG);
    scaleCostT = zeros(nS*Tdays,length(scaleCost(1,:)));
    FirstDispT = zeros(nS*Tdays,length(FirstDisp(1,:)));
    predictDisp = [];
    maxGenDem = 0;
    for i = 1:1:length(Out)
        DemandT.(Out{i}) = zeros(nS*Tdays, length(GenDemand.(Out{i})(1,:)));%this is demand Total: all the demand including what is provided by renewables over all the training days
    end
    for day = 0:1:Tdays-1%run for number of training days
        if day==0
            for i = 1:1:nG
                if isfield(Plant.Generator(i).OpMatA.output, 'C') %make sure chillers have positive initial conditions
                    if CurrentState.Generators(i)<0
                        CurrentState.Generators(i) = LBmpc(i);
                    end
                end
            end
            IC = CurrentState.Generators;
        else IC = predictDisp(1,:);
            DateSim = DateSim+1;
        end
        if ~isempty(dischEff) %must be a storage system
            IC(stor)=max(IC(stor).*(1-(1-dischEff(stor))*dt(1)),LB(stor));%scale storage value by discharge efficiency (this scaling is used in optimizations)
        end
        IC = min(UB,IC);
        scaleCost = updateGeneratorCost(Time); %% All feedstock costs were assumed to be 1 when building matrices 
        if day<1
            PredictDispatch = ones(length(Time)+1,1)*IC;
        else PredictDispatch = [IC;predictDisp(2:end,:)];
        end
        marginCost = updateMarginalCost(PredictDispatch,scaleCost,Time);%the dispatch is whatever has been dispatched so far, except for the initial condition.
        %change demand to be the historical data + some gaussian noise
        Organize = Plant.OpMatA.Organize;
        
        [QPall,~] = updateMatrices(Plant.OpMatA.QP,Organize,IC,Time,scaleCost,marginCost,[]);
        %the demand should be historical data + gaussian noise (-renewable
        %generation if it is added to the QP)
        if ~isempty(renews)%renewable generation
            [renewgen, ~]= RenewableOutput(DateSim,Time,'Predict');%Add forecast of renewable power @ resulution equal to the Time vector
            renewgen = sum(renewgen,2);%this assumes that all renewable generation goes to electric demand
        else
            renewgen = zeros(nS,1);
        end
        ibeq = 0;%index of demand within beq matrix
        for j = 1:1:length(Out)
            histData = Plant.Data.Demand.(Out{j})((DateSim-Plant.Data.Timestamp(1))*24*4+round((Time-Time(1))./0.25)+1);%assumes data is every 15 minutes
            sizeHistData = size(histData);
            if sizeHistData(1) == 1 %if format is in a single row, change it to a single collumn
                Last24hour.(Out{j}) = histData;
                histData = histData';
            else Last24hour.(Out{j}) = histData';
            end
            sigma = std(histData)/2;
            qpdemand = max(histData + randn(size(histData))*sqrt(sigma),0);%historical data plus gaussian noise
            DemandT.(Out{j})(1+nS*day:nS*(1+day),:) = qpdemand;
            if strcmp('H',Out{j}) && Plant.optimoptions.excessHeat
                QPall.(Out{1}).b(1:nS) = qpdemand;
            elseif ~Plant.optimoptions.sequential && strcmp('C',Out{j})
                QPall.(Out{1}).beq(ibeq+1:ibeq+nS) = qpdemand;
                ibeq = ibeq+nS;
            else
                QPall.(Out{1}).beq(ibeq+1:ibeq+nS) = qpdemand-renewgen;
                ibeq = ibeq+nS;
            end
        end
            
        %% Step 3 Determine initial dispatch
        [FirstDisp,~,Feasible] = DispatchQP(QPall,Organize,[]);

        %% Step 4 Refine Dispatch & Plot
        for j = 1:1:length(Out)
            StorPower.(Out{j}) = zeros(nS,1);
            GenDemand.(Out{j}) = zeros(nS,1);
            for i = 1:1:nG
                if isfield(Plant.Generator(i).OpMatA.output,Out{j})
                    if ~isfield(Plant.Generator(i).OpMatA,'Stor') %don't include the storage use/charging in this new profile;
                        GenDemand.(Out{j}) = GenDemand.(Out{j}) + FirstDisp(2:end,i)*Plant.Generator(i).OpMatB.output.(Out{j});
                    else  StorPower.(Out{j}) = StorPower.(Out{j}) + (FirstDisp(1:end-1,i)-FirstDisp(2:end,i))./dt'; 
                    end
                end
            end
        end
        OptimalState = StepByStepDispatch(GenDemand,scaleCost,dt,IC,'initially constrained',FirstDisp);

        %% check if any gen can be completley removed
        Locked = true(nS+1,nG);
        %% Rule 1 if off-line for entire optimal dispatch, then the generator should be off
        for i = 1:1:ngens
            if nnz(OptimalState(:,gens(i)))==0
                Locked(:,gens(i))=false;
            end
        end
        [FirstDisp, ~, Feasible] = DispatchQP(QPall,Organize,Locked);%this is the dispatch with fit A
        if Feasible ~=1
            Locked = true(nS+1,nG); %reset to all on default
        end

        Organize = Plant.OpMatB.Organize;
        marginCost = updateMarginalCost(FirstDisp,scaleCost,Time);
        [QPall,Forecast] = updateMatrices(Plant.OpMatB.QP,Organize,IC,Time,scaleCost,marginCost,[]); %update fit B matrices
        %make sure that demand is correct
        ibeq = 0;%index of demand within beq matrix
        for j = 1:1:length(Out)
            if strcmp('H',Out{j}) && Plant.optimoptions.excessHeat
                QPall.(Out{1}).b(1:nS) = DemandT.(Out{j})(1+nS*day:nS*(1+day),:);
            elseif ~Plant.optimoptions.sequential && strcmp('C',Out{j})
                QPall.(Out{1}).beq(ibeq+1:ibeq+nS) = DemandT.(Out{j})(1+nS*day:nS*(1+day),:);
                ibeq = ibeq+nS;
            else
                QPall.(Out{1}).beq(ibeq+1:ibeq+nS) = DemandT.(Out{j})(1+nS*day:nS*(1+day),:)-renewgen;
                ibeq = ibeq+nS;
            end
        end
        [DayDisp, dX] = FilterGenerators(QPall,Organize,IC,Forecast,FirstDisp,OptimalState,scaleCost,Locked);
        GenDisp(1+nS*day:nS*(day+1),:) = DayDisp(2:end,:);
        scaleCostT(1+nS*day:nS*(day+1),:) = scaleCost;
        FirstDispT(1+nS*day:nS*(1+day),:) = FirstDisp(2:end,:);
        predictDisp = DayDisp;%[DayDisp(ceil(.3*nS):end,:);DayDisp(1:ceil(.3*nS)-1,:)];
        maxGenDem = max(max(GenDemand.E),maxGenDem);
        Last24hour.Timestamp = DateSim+Time./24;
    end
    GenDisp = [GenDisp(1,:);GenDisp];
    scaleCost = scaleCostT;
    FirstDisp = [FirstDispT;FirstDispT(end,:)];
    for i = 1:1:length(Out)
        GenDemand.(Out{i}) = DemandT.(Out{i});
    end
    Time = repmat(Time,1,Tdays);
    dt = repmat(dt,1,Tdays);
    nS = length(Time);
    DateSim = DateSim1;
    Last24hour = Last24hour1;
end

%% create train and dispatch a locked network
lockInputs = zeros(nS, 5*ngens+2+length(Out));
for i = 1:1:length(Out)
    lockInputs(:,i) = GenDemand.(Out{i})./(max(GenDemand.(Out{i}))*1.5);%demands
end
%lockInputs(:,1) = GenDemand.E./(maxGenDem*1.5);%demand
lockInputs(:,2) = (scaleCost(:,1).*dt')./max(scaleCost(:,1).*dt');%cost of utility electrical
lockInputs(:,3) = FirstDisp(1:end-1,stor)./UB(stor);%SOC of battery at previous step
for i = 1:1:ngens
    j = gens(i);
    if ~isempty(renews)
        j1 = j-nnz(renews<j);%this index is used for f and H matrices
    else j1 = j;
    end
    if ~isfield(Plant.Generator(j).OpMatB.output,'C') || ~Plant.optimoptions.sequential %if its not a chiller or the chillers are run concurrently pull costs from OneStep.E
        lockInputs(:,3+i) = (scaleCost(:,2).*Plant.OneStep.E(1).f(j1).*dt')./max(scaleCost(:,2).*Plant.OneStep.E(1).f(j1).*dt');%linear cost
        lockInputs(:,3+ngens+i) = (scaleCost(:,2).*Plant.OneStep.E(1).H(j1,j1).*dt')./max(scaleCost(:,2).*Plant.OneStep.E(1).H(j1,j1).*dt');%quadratic cost
    else %if its a chiller run in sequence before generators, pull costs from OneStep.C
        lockInputs(:,3+i) = (scaleCost(:,2).*Plant.OneStep.C(1).f(j1).*dt')./max(scaleCost(:,2).*Plant.OneStep.C(1).f(j1).*dt');%linear cost
        lockInputs(:,3+ngens+i) = (scaleCost(:,2).*Plant.OneStep.C(1).H(j1,j1).*dt')./max(scaleCost(:,2).*Plant.OneStep.C(1).H(j1,j1).*dt');%quadratic cost
    end
    lockInputs(:,3+2*ngens+i) = (scaleCost(:,2).*Plant.Generator(j).OpMatB.constCost.*dt')./max(scaleCost(:,2).*Plant.Generator(j).OpMatB.constCost.*dt');%onCost
    lockInputs(:,3+3*ngens+i) = min(FirstDisp(1:end-1,j)+Plant.Generator(j).OpMatB.Ramp.b(1),UB(j))./UB(j);%upper bound
    lockInputs(:,3+4*ngens+i) = max(FirstDisp(1:end-1,j)+Plant.Generator(j).OpMatB.Ramp.b(1),LB(j))./UB(j);%lower bound
end

%prevent Nans
lockInputs(isnan(lockInputs)) = 0;

if training %train it the first time through, use it later
    lockedNet = Neural_Network(length(lockInputs(1,:)),ngens,'classify',1);
    %train locked network
    [lockedNet, lockederror] = trainNetwork(lockedNet, (GenDisp(2:end,gens)>0)',lockInputs);
    nS = nS/Tdays;
end
%use the network to find which generators to lock off or on
LockedNet = forward(lockedNet, lockInputs(1:nS,:));
LockedNet(LockedNet>=0.7) = 1;
LockedNet(LockedNet<0.7) = 0;

%% create train and dispatch dispatch network
% geninputs = zeros(nS, length(lockInputs(1,:))+3*length(stor)+length(LockedNet(:,1)));
% geninputs(:,1:length(lockInputs(1,:))) = lockInputs;
% index = length(lockInputs(1,:));
% %index = index+ngens*2;
% %add the discharge efficiency, and marginal costs
% for i = 1:1:length(stor)
%     j = stor(i);
%     geninputs(:,index+i) = dischEff(j);
%     %geninputs(:,index+length(stor)+i) = FirstDisp(1:end-1,j);%SOC at previous step
%     geninputs(:,index+length(stor)*2+i) = (marginCost.E.Min)./max(max(scaleCost));%margin cost at time t
%     geninputs(:,index+length(stor)*3+i) = marginCost.E.Max./max(max(scaleCost));
% end
% index = index+length(stor)*3;
% 
% allLocked = LockedNet;
% %if a generator is locked off, change all inputs to 0
% geninputs(:,index+1:end) = LockedNet';
% % if sum(sum(abs(lockederror)))>nS/10 %if incorrect for more than 5%, only use as an input
% %     allLocked = ones(ngens,nS);
% % end
%     
% if training
%     if sum(sum(abs(lockederror)))>nS/10%if incorrect for more than 5%, use the training Locked config
%         allLocked = (GenDisp(2:end,gens)>0)';
%     end
%     %rescale the generator dispatch for the upper/lower bound
%     net = Neural_Network(length(geninputs(1,:)),ngens);%ngs);
%     for i = 1:1:ngens%ngs
%         j = gens(i);
%         geninputs1 = geninputs.*(ones(length(geninputs(1,:)),1)*allLocked(i,:))';
%         componentnet = Neural_Network(length(geninputs(1,:)),1);
%         GenDispTrain = max((GenDisp(2:end,j)'-LB(j))./(UB(j)-LB(j)),0);
%         [componentnet, sqrerror] = trainNetwork(componentnet, GenDispTrain, geninputs1);
%         if nnz(sqrerror>.005)>nS*ngs/8 %if more than an 8th the results have high error try training again
%             shuffle = round(rand(nS,1).*nS);
%             shuffle(shuffle==0) = 1;
%             geninputsshuff = geninputs1(shuffle, :);
%             GenDispTrainshuff = GenDispTrain(:,shuffle);
%             [componentnet2,sqrerror2] = trainNetwork(componentnet, GenDispTrainshuff, geninputsshuff);
%             if mean(mean(sqrerror2))<mean(mean(sqrerror))
%                 componentnet = componentnet2;
%                 sqrerror = sqrerror2;
%             end
%         end
%         net.Wlayer1(:,i) = componentnet.Wlayer1;
%     end
% end
%  
% if training && timersOn
%     disp(['time for training neural network is ' num2str(toc)])
% end
% 
% % tic
% netDisp = zeros(nS,ngs);
% for i = 1:1:ngens %add in storage later
%     j = gens(i);
%     geninputs1 = geninputs.*(ones(length(geninputs(1,:)),1)*allLocked(i,:))';
%     netDisp(:,i) = (geninputs1*net.Wlayer1(:,i)).*(UB(j)-LB(j))+LB(j).*allLocked(i,:)';
%     netDisp(netDisp(:,i)<LB(j)/2,i) = 0;
% end
% % if training && timersOn
% %     disp(['time for dispatching neural network is ' num2str(toc)])
% % end
% if training && nnz(abs(netDisp-GenDisp(2:end,gensandstor))>0.1)>0 %if there is a .1kW difference between the dispatch from the net and dispatch from StepByStep
%     disp('Neural Network dispatch differs from OprimalState. training needs revision.')
% end

