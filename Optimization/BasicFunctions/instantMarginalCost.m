function marginal = instantMarginalCost(Dispatch,scaleCost)
global Plant UB
Outs = Plant.optimoptions.Outputs;
for i = 1:1:length(Outs)
    Out.(Outs{i}) = [];
    storType.(Outs{i}) = [];
end
nG = length(Plant.Generator);
marginCost = zeros(1,nG);
stor = [];
Out.CHP = [];
I = zeros(1,nG);
for i = 1:1:nG
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
        S = fieldnames(Plant.Generator(i).OpMatA.output);
        storType.(S{1})(end+1) = i;
    elseif~isempty(s) %utilities and single state generators (linear cost term)
        marginCost(i) = Plant.Generator(i).OpMatB.(s{1}).f;
        I(i) = UB(i);
    end
    if I(i)>0
        S = fieldnames(Plant.Generator(i).OpMatB.output);
        if length(S)==2 && ismember('H',S) && ismember('E',S)
            Out.E(end+1) = i;
            Out.CHP(end+1) = i;
        else
            Out.(S{1})(end+1) = i;
        end
    end
end
MarginCost = scaleCost.*marginCost; %marginal cost
Type = fieldnames(storType);
for i = 1:1:length(Type)
    if i==1 || ~isfield(marginal,Type{i})
        marginal.(Type{i}) = [];
    end
    gen = Out.(Type{i});
    
    if strcmp(Type{i},'C') && Plant.optimoptions.sequential == 0 %chillers have no cost (show up as electric load)
        Egen = Out.E;
        for j = 1:1:length(gen)
            Cratio = Plant.Generator(gen(j)).Output.Cooling(end);
            MarginCost(gen(j)) = min(MarginCost(Egen))/Cratio;
        end
    end
    
    ThisType = MarginCost(gen);
    if ~isempty(Out.CHP)
        if strcmp(Type{i},'E')
            for j = 1:1:length(Out.CHP)
                [~,k] = ismember(Out.CHP(j),gen);
                ThisType(k) = 0.75*ThisType(k); %assign 25% of the generator cost to the heat production
            end
        end
        if strcmp(Type{i},'H')
            ThisType(end+1:end+length(Out.CHP)) = MarginCost(:,Out.CHP)*.25; %assign 25% of the generator cost to the heat production
        end
    end
    marginal.(Type{i}) = min(ThisType); 
end