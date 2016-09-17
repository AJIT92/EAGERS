function marginal = updateMarginalCost(Dispatch,scaleCost,Time)
global Plant UB
% Time = buildTimeVector(Plant.optimoptions);
dt = Time-[0,Time(1:end-1)];
Outs = Plant.optimoptions.Outputs;
for i = 1:1:length(Outs)
    Out.(Outs{i}) = [];
    storType.(Outs{i}) = [];
end
marginal =[];
nG = length(Plant.Generator);
marginCost = zeros(4,nG);
stor = [];
Out.CHP = [];
I = zeros(1,nG);
for i = 1:1:nG
    s = Plant.Generator(i).OpMatB.states;
    if isfield(Plant.Generator(i).OpMatB,'constCost') %all of these cost terms need to be scaled later on
        I(i) = Plant.Generator(i).OpMatB.(s{1}).ub;
        marginCost(1,i) = Plant.Generator(i).OpMatB.(s{1}).f;
        marginCost(2,i) = Plant.Generator(i).OpMatB.(s{2}).f;
        marginCost(3,i) = Plant.Generator(i).OpMatB.(s{2}).H;
    elseif isfield(Plant.Generator(i).OpMatA,'Stor')
        stor(end+1) = i;
        S = fieldnames(Plant.Generator(i).OpMatA.output);
        storType.(S{1})(end+1) = i;
    elseif~isempty(s) %utilities and single state generators (linear cost term)
        marginCost(1,i) = Plant.Generator(i).OpMatB.(s{1}).f;
        I(i) = UB(i);
    end
    if I(i)>0
        S = fieldnames(Plant.Generator(i).OpMatA.output);
        if length(S)==2 && ismember('H',S) && ismember('E',S)
            Out.E(end+1) = i;
            Out.CHP(end+1) = i;
        else
            Out.(S{1})(end+1) = i;
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
    minMarginCost(:,i) = scaleCost(:,i).*marginCost(1,i); %marginal cost in linear segment
    if minMarginCost(:,i)<0 %negative slope at cost fit near LB (try 1/2 way to UB)
        half = UB(i)- (UB(i)-I(i))/2;
        minMarginCost(:,i) = scaleCost(:,i).*(marginCost(2,i) + (half-I(i))*marginCost(3,i)); %marginal cost 1/2 way in quadratic segment
    end
    maxMarginCost(:,i) = minMarginCost(:,i);
    if I(i)>0 && I(i) ~= UB(i)
        maxMarginCost(:,i) = scaleCost(:,i).*(marginCost(2,i) + (UB(i)-I(i))*marginCost(3,i)); %marginal cost in quadratic segment
    end
end
Type = fieldnames(storType);
for i = 1:1:length(Type)
    if ~isfield(marginal,(Type{i}))
        marginal.(Type{i}) = [];
    end
    gen = Out.(Type{i});
    stor = storType.(Type{i});
    if strcmp(Type{i},'C') && Plant.optimoptions.sequential == 0 %chillers have no cost (show up as electric load)
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
            minMarginCost(GenOnDuringCharge==1) = NaN; %if the storage is charging, then its margin cost can be greater than the cost of the generators that are on
        end
    end
    
    MinThisType = minMarginCost(:,gen);
    MaxThisType = maxMarginCost(:,gen);
    if ~isempty(Out.CHP)
        if strcmp(Type{i},'E')
            for j = 1:1:length(Out.CHP)
                [~,k] = ismember(Out.CHP(j),gen);
                MinThisType(:,k) = 0.75*MinThisType(:,k); %assign 25% of the generator cost to the heat production
                MaxThisType(:,k) = 0.75*MaxThisType(:,k); %assign 25% of the generator cost to the heat production
            end
        end
        if strcmp(Type{i},'H')
            MinThisType(:,end+1:end+length(Out.CHP)) = minMarginCost(:,Out.CHP)*.25; %assign 25% of the generator cost to the heat production
            MaxThisType(:,end+1:end+length(Out.CHP)) = maxMarginCost(:,Out.CHP)*.25; %assign 25% of the generator cost to the heat production
        end
    end
    
    minOn_t = min(MinThisType,[],2); 
    maxOn_t = min(MaxThisType,[],2);
    
    timedivideMin = sum(dt(minOn_t~=0));
    timedivideMax = sum(dt(maxOn_t>0));
    minOn_t(isnan(minOn_t))=0;
    maxOn_t(isnan(maxOn_t))=0;
    marginal.(Type{i}).Min = sum(minOn_t(minOn_t~=0).*dt(minOn_t~=0)')/timedivideMin;%make the cost proportional to amount of time at that cost
    marginal.(Type{i}).Max = sum(maxOn_t(maxOn_t>0).*dt(maxOn_t>0)')/timedivideMax;
end