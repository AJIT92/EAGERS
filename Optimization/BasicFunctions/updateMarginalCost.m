function marginal = updateMarginalCost(Dispatch,scaleCost,Time)
global Plant
dt = Time-[0;Time(1:end-1)];
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
for i = 1:1:length(networkNames)
    storType.(networkNames{i}) = [];
    if strcmp(networkNames{i},'Electrical')
        Out.E = [];
    elseif strcmp(networkNames{i},'DistrictHeat')
        Out.H = [];
    elseif strcmp(networkNames{i},'DistrictCool')
        Out.C = [];
    elseif strcmp(networkNames{i},'Hydro')
        Out.W = [];
    end
end
marginal =[];
nG = length(Plant.Generator);
UB = zeros(nG,1);
for i = 1:1:nG
    states = Plant.Generator(i).OpMatB.states;
    for j = 1:1:length(states);
        UB(i) = UB(i) + Plant.Generator(i).OpMatB.(states{j}).ub;
    end
end
marginCost = zeros(4,nG);
stor = [];
CHP = [];
I = zeros(1,nG);
for i = 1:1:nG
    if Plant.Generator(i).Enabled
        s = Plant.Generator(i).OpMatB.states;
        if length(s)>1 && isfield(Plant.Generator(i).OpMatB, 'constCost') %all of these cost terms need to be scaled later on
            I(i) = Plant.Generator(i).OpMatB.(s{1}).ub;
            marginCost(1,i) = Plant.Generator(i).OpMatB.(s{1}).f;
            marginCost(2,i) = Plant.Generator(i).OpMatB.(s{2}).f;
            marginCost(3,i) = Plant.Generator(i).OpMatB.(s{2}).H;
            marginCost(4,i) = Plant.Generator(i).OpMatB.constCost/Plant.Generator(i).Size;%constant cost/upperbound
        elseif isfield(Plant.Generator(i).OpMatA,'Stor')
            stor(end+1) = i;
            if strcmp(Plant.Generator(i).Source,'Electricity')
                storType.Electrical(end+1) = i;
            elseif strcmp(Plant.Generator(i).Source,'Heat')
                storType.DistrictHeat(end+1) = i;
            elseif strcmp(Plant.Generator(i).Source,'Cooling')
                storType.DistrictCool(end+1) = i;
            elseif strcmp(Plant.Generator(i).Source,'Water')
                storType.Hydro(end+1) = i;
            end
        elseif~isempty(s) %utilities and single state generators (linear cost term)
            I(i) = UB(i);
            marginCost(1,i) = Plant.Generator(i).OpMatB.(s{1}).f;
            if isfield(Plant.Generator(i).OpMatB, 'constCost')
                marginCost(4,i) = Plant.Generator(i).OpMatB.constCost/Plant.Generator(i).Size;%constant cost/upperbound
            end
        end
        if I(i)>0
            S = fieldnames(Plant.Generator(i).OpMatA.output);
            if length(S)==2 && ismember('H',S) && ismember('E',S)
                Out.E(end+1) = i;
                CHP(end+1) = i;
            elseif isfield(Out, (S{1}))
                Out.(S{1})(end+1) = i;
            end
        end
    end
end
[m,n] = size(Dispatch);
if m>1
    m = m-1; %eliminate IC unless solving for IC
    GenOnDuringCharge = false(m,n);
end
minMarginCost = nan(m,n);
maxMarginCost = nan(m,n);

for i = 1:1:length(I)
    minMarginCost(:,i) = scaleCost(:,i).*(marginCost(1,i)+marginCost(4,i)); %marginal cost in linear segment+onCost
    if minMarginCost(:,i)<0 %negative slope at cost fit near LB (try 1/2 way to UB)
        half = UB(i)- (UB(i)-I(i))/2;
        minMarginCost(:,i) = scaleCost(:,i).*(marginCost(2,i) + (half-I(i))*marginCost(3,i)) + scaleCost(:,i).*marginCost(4,i); %marginal cost 1/2 way in quadratic segment
    end
    maxMarginCost(:,i) = minMarginCost(:,i);
    if I(i)>0 && I(i) ~= UB(i)
        maxMarginCost(:,i) = scaleCost(:,i).*(marginCost(2,i) + (UB(i)-I(i))*marginCost(3,i)) + scaleCost(:,i).*marginCost(4,i); %marginal cost in quadratic segment
    end
end
Type = fieldnames(storType);
for i = 1:1:length(Type)
    if ~isfield(marginal,(Type{i}))
        marginal.(Type{i}) = [];
    end
    stor = storType.(Type{i});
    if strcmp(Type{i},'Electrical')
        gen = Out.E;
    elseif strcmp(Type{i},'DistrictHeat')
        gen = Out.H;
    elseif strcmp(Type{i},'DistrictCool')
        gen = Out.C;
    elseif strcmp(Type{i},'Hydro')
        gen = Out.W;
    end
    if strcmp(Type{i},'DistrictCool') && Plant.optimoptions.sequential == 0 %chillers have no cost (show up as electric load)
        Egen = Out.E;
        for j = 1:1:length(gen)
            Cratio = Plant.Generator(gen(j)).Output.Cooling(end);
            minMarginCost(:,gen(j)) = min(minMarginCost(:,Egen),[],2)/Cratio;
            maxMarginCost(:,gen(j)) = max(maxMarginCost(:,Egen),[],2)/Cratio;
        end
    end
    
    if exist('GenOnDuringCharge','var') 
        ChargeIndex = false(m,1);
        for j = 1:1:length(stor)
            ChargeIndex = max(ChargeIndex,Dispatch(2:end,stor(j))-Dispatch(1:end-1,stor(j))>0);%[false;(Dispatch(2:end,stor(i))-Dispatch(1:end-1,stor(i)))>0]; %timesteps where charging
        end
        GenOnDuringCharge(ChargeIndex,gen) = Dispatch(ChargeIndex,gen)>0;
        if nnz(ChargeIndex)>0 && nnz(GenOnDuringCharge)>0
            maxMarginCost(GenOnDuringCharge==1) = NaN; %if the storage is charging, then its margin cost can be greater than the cost of the generators that are on
        end
    end
    
    MinThisType = minMarginCost(:,gen);
    MaxThisType = maxMarginCost(:,gen);
    if ~isempty(CHP)
        if strcmp(Type{i},'Electrical')
            for j = 1:1:length(CHP)
                [~,k] = ismember(CHP(j),gen);
                MinThisType(:,k) = 0.75*MinThisType(:,k); %assign 25% of the generator cost to the heat production
                MaxThisType(:,k) = 0.75*MaxThisType(:,k); %assign 25% of the generator cost to the heat production
            end
        end
        if strcmp(Type{i},'DistrictHeat')
            MinThisType(:,end+1:end+length(CHP)) = minMarginCost(:,CHP)*.25; %assign 25% of the generator cost to the heat production
            MaxThisType(:,end+1:end+length(CHP)) = maxMarginCost(:,CHP)*.25; %assign 25% of the generator cost to the heat production
        end
    end
    
    minOn_t = min(MinThisType,[],2); 
    maxOn_t = min(MaxThisType,[],2);
    
    timedivideMin = sum(dt(minOn_t~=0));
    timedivideMax = sum(dt(maxOn_t>0));
    minOn_t(isnan(minOn_t))=0;
    maxOn_t(isnan(maxOn_t))=0;
    marginal.(Type{i}).Min = sum(minOn_t(minOn_t~=0).*dt(minOn_t~=0))/timedivideMin;%make the cost proportional to amount of time at that cost
    marginal.(Type{i}).Max = sum(maxOn_t(maxOn_t>0).*dt(maxOn_t>0))/timedivideMax;
end