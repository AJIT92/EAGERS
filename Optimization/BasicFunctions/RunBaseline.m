function Baseline = RunBaseline
%% Perform cost analysis & compare to baseline (control without forward looking optimization)
% the goal here is to run simple dydnamic economic dispatch
% we can explore other alternatives for what the baseline is
global Plant UB MapWaitbarHandle
Time = buildTimeVector(Plant.optimoptions);
Baseline.Timestamp = Plant.Data.Timestamp(1)+[0,Time]./24;
%% need to repeat the following if test is longer than horizon
scaleCost = updateGeneratorCost(Time);
margincost = updateMarginalCost(ones(length(Time)+1,1)*UB,scaleCost,Time);
IC = zeros(1,length(Plant.Generator));
stor =[];
for i = 1:1:length(Plant.Generator)
    if isfield(Plant.Generator(i).OpMatB,'Stor')
       stor(end+1) = i; 
    end
end
IC(stor) = .5*UB(stor);
scaleCost(:,stor)=0; %remove cst of storage
QPall = Plant.OpMatB.QP;
Outs = fieldnames(QPall);
Organize = Plant.OpMatB.Organize;
[QPall,~] = updateMatrices(QPall,Organize,IC,Time,scaleCost,margincost,[]); 
MapWaitbarHandle=waitbar(0,'Running Baseline Dynamic Economic Dispatch');
Locked = ones(length(Time),length(UB)); % will keep all generators on at all times. May not be feasible
[GenDisp,  ~, Feasible] = DispatchQP(QPall,Organize,Locked);%this is the dispatch with fit B, and all generators on
close(MapWaitbarHandle)
MapWaitbarHandle =[];
if Feasible ==1
    Input = 0*GenDisp;
    for i = 1:1:length(Plant.Generator)
        chiller = 0;
        if ~isempty(Plant.Generator(i).Output)
            cap = Plant.Generator(i).Output.Capacity*UB(i);
        end
        eff = [];
        if strcmp(Plant.Generator(i).Type,'Electric Generator') || strcmp(Plant.Generator(i).Type,'CHP Generator')
            eff = Plant.Generator(i).Output.Electricity;
        elseif strcmp(Plant.Generator(i).Type,'Chiller') && ~ismember('E',Outs)%don't include cost if it shows up in generator demand
            eff = Plant.Generator(i).Output.Cooling;
            if ~Plant.optimoptions.sequential
                chiller = 1;
            end
        elseif strcmp(Plant.Generator(i).Type,'Heater')
            eff = Plant.Generator(i).Output.Heat;    
        end
        if ~isempty(eff) && ~chiller %dont add the cost of a chiller if you ran E and C simultaneously, or you will double count the chiller demand
            Input(:,i) = GenDisp(:,i)./interp1(cap,eff,GenDisp(:,i));
        elseif strcmp(Plant.Generator(i).Type,'Utility')
            Input(:,i) = GenDisp(:,i);
        end
    end
    Baseline.GeneratorState = GenDisp;
    Baseline.NetCost = NetCostCalc(GenDisp,Time,Input);
else
    disp('baseline dispatch is infeasible with all generators on')
end