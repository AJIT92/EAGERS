function marginal = instantMarginalCost(Dispatch,scaleCost)
global Plant
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
nG = length(Plant.Generator);  
marginCost = zeros(1,nG);
stor = [];
CHP = [];
I = zeros(1,nG);
for i = 1:1:nG
    if Plant.Generator(i).Enabled %only use enabled gens
        s = Plant.Generator(i).OpMatB.states;
        if isfield(Plant.Generator(i).OpMatB,'constCost') %all of these cost terms need to be scaled later on       
            I(i) = Plant.Generator(i).OpMatB.(s{1}).ub;
            if isempty(Dispatch) || Dispatch(i)<=I(i)
                marginCost(i) = Plant.Generator(i).OpMatB.(s{1}).f;
            elseif Dispatch(i)>I(i)
                marginCost(i) = Plant.Generator(i).OpMatB.(s{2}).f + (Dispatch(i)-I(i))*Plant.Generator(i).OpMatB.(s{2}).H;
            end
        elseif isfield(Plant.Generator(i).OpMatB,'Stor')
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
            marginCost(i) = Plant.Generator(i).OpMatB.(s{1}).f;
            for j = 1:1:length(s);
                I(i) = I(i) + Plant.Generator(i).OpMatB.(s{j}).ub;
            end
        end
        if I(i)>0
            S = fieldnames(Plant.Generator(i).OpMatB.output);
            if length(S)==2 && ismember('H',S) && ismember('E',S)
                Out.E(end+1) = i;
                CHP(end+1) = i;
            elseif isfield(Out,(S{1}))
                Out.(S{1})(end+1) = i;
            end
        end
    end
end
MarginCost = scaleCost.*marginCost; %marginal cost
Type = fieldnames(storType);
for i = 1:1:length(Type)
    if i==1 || ~isfield(marginal,Type{i})
        marginal.(Type{i}) = [];
    end
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
            MarginCost(gen(j)) = min(MarginCost(Egen))/Cratio;
        end
    end
    
    ThisType = MarginCost(gen);
    if ~isempty(CHP)
        if strcmp(Type{i},'Electrical')
            for j = 1:1:length(CHP)
                [~,k] = ismember(CHP(j),gen);
                ThisType(k) = 0.75*ThisType(k); %assign 25% of the generator cost to the heat production
            end
        end
        if strcmp(Type{i},'DistrictHeat')
            ThisType(end+1:end+length(CHP)) = MarginCost(:,CHP)*.25; %assign 25% of the generator cost to the heat production
        end
    end
    marginal.(Type{i}) = min(ThisType); 
end