function Baseline = RunBaseline(IC)
% the goal here is to run simple dynamic economic dispatch
% we can explore other alternatives for what the baseline is
global Plant MapWaitbarHandle TestData
Outs = fieldnames(TestData.Demand);
Date = TestData.Timestamp(1);
options = Plant.optimoptions;
options.tspacing = 'constant';
Time = buildTimeVector(options);
nG = length(Plant.Generator);
Renewable = zeros(length(Time),nG);
t = 1;
steps = nnz(TestData.Timestamp<TestData.Timestamp(1)+Plant.optimoptions.Horizon/24); %number of data points in horizon
Baseline.NetCost = 0;
MapWaitbarHandle=waitbar(0,'Running Baseline Dynamic Economic Dispatch');
while Date<TestData.Timestamp(1)+Plant.optimoptions.Interval
    Baseline.Timestamp = Date+[0,Time]./24;
    scaleCost = updateGeneratorCost(Baseline.Timestamp);
    margincost = updateMarginalCost(IC,scaleCost,Time);

    for S = 1:1:length(Outs)
        Demand.(Outs{S}) = TestData.Demand.(Outs{S})(t:t+nS-1);
    end
    for i = 1:1:nG
        if strcmp(Plant.Generator(i).Source,'Renewable') && Plant.Generator(i).Enabled
            Renewable(:,i) = RenewableOutput(Plant.Generator(i).VariableStruct,Date,Time,'Actual');
        end
    end

    QP = updateMatrices(Plant.OpMatB.QP,IC,Time,scaleCost,margincost,Demand,Renewable,[]); 
    
    Locked = ones(length(Time),nG); % will keep all generators on at all times. May not be feasible
    [GenDisp,~,Feasible] = DispatchQP(QP,Locked);%this is the dispatch with fit B, and all generators on
    
    if Feasible ==1
        Baseline.GeneratorState(t:t+nS,:) = GenDisp;
        Baseline.NetCost = Baseline.NetCost + NetCostCalc(GenDisp,Baseline.Timestamp,'Dispatch');
        t = t+ steps;
        Date = Date + Plant.optimoptions.Horizon/24;
        IC = GenDisp(end,:);
    else
        disp('baseline dispatch is infeasible with all generators on')
        Date = TestData.Timestamp(1)+Plant.optimoptions.Interval+1;
    end
end
%%--%%
close(MapWaitbarHandle)
MapWaitbarHandle =[];