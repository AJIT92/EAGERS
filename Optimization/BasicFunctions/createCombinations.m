function [K, Nodes,Outs,include] = createCombinations(QP,netDemand)
global Plant
Outs = fieldnames(Plant.Data.Demand);
nG = length(Plant.Generator);
if any(strcmp('E',Outs)) && any(strcmp('H',Outs))
    %check if there is a CHP system
    isCHP = false;
    for i = 1:1:nG
        if strcmp(Plant.Generator(i).Type,'CHP Generator')
            isCHP = true;
            break
        end
    end
    if isCHP
        Outs = Outs(~strcmp('H',Outs)); %heaters are part of electricity case
    end
end
K = zeros(0,nG);
lines = 0;
for s = 1:1:length(Outs)
    inc = false(nG,1);
    for i = 1:1:nG
        if isfield(Plant.Generator(i).OpMatA.output,Outs{s})
            if isempty(strfind(Plant.Generator(i).Type,'Utility')) &&  isempty(strfind(Plant.Generator(i).Type,'Storage')) && Plant.Generator(i).Enabled
                inc(i) = true;
            end
        end
    end
    
    include.(Outs{s}) = find(inc);
    if (strcmp(Outs{s},'E') && isCHP) %combine heaters into CHP case
        for i = 1:1:nG
            if  isempty(strfind(Plant.Generator(i).Type,'Utility')) &&  isempty(strfind(Plant.Generator(i).Type,'Storage')) && Plant.Generator(i).Enabled && isfield(Plant.Generator(i).OpMatA.output,'H')
                inc(i) = true;
            end
        end
    end
    ninc = ~inc;
    inc = find(inc);
    ninc = find(ninc);
    n = length(inc);
    K = [K; zeros(2^n,nG)]; %K is a matrix of all the possible generator on/off combinations 
    if ~isempty(inc) && ~isempty(ninc)
        K(lines+1:end,ninc) = ones(2^n,1)*ninc'; % all systems that are always included
    end
    for j = 1:1:n %all combinations of generators are listed below
        z = 2^(n-j);
        r=0;
        while r+z<=2^n
            K(lines+r+1:lines+r+z,inc(j)) = inc(j);
            r=r+2*z;
        end
    end
    
    %% -- %%
    %% test if each combination is capable of meeting demand
    %if the generators lower bounds are higher than the demand, then the case is invalid and should be removed. 
    %if the generators upper bounds are lower than the demand, then this case does not produce enough and should be removed.
    %%could improve this section to check if feasible with network losses (not sure how)
    req = [];
    if strcmp(Outs{s},'E')
        req = QP.Organize.Balance.Electrical; %rows of Aeq associated with electric demand
    elseif strcmp(Outs{s},'H')
        req = QP.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
    elseif strcmp(Outs{s},'C')
        req = QP.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
    end
    Nodes.(Outs{s}) = length(req); %number of nodes in this output category
    limitU.(Outs{s}) = zeros(1,nG);
    limitL.(Outs{s}) = zeros(1,nG);
    for i = 1:1:nG
        states = QP.Organize.States{i};
        if~isempty(states)
            if any(isinf(QP.ub(states)))
                if any(QP.Aeq(req,states)~=0)
                    limitU.(Outs{s})(i) = inf;
                end
            else
                limitU.(Outs{s})(i) = limitU.(Outs{s})(i) + sum(QP.Aeq(req,states)*QP.ub(states));
            end
            if any(isinf(QP.lb(states)))
                if any(QP.Aeq(req,states)~=0)
                    limitL.(Outs{s})(i) = -inf;
                end
            else
                if strcmp(Outs{s},'H') && isfield(QP.Organize,'HeatVented') && any(QP.Organize.HeatVented(1,:))
                    lb = min(0,sum(QP.Aeq(req,states)*QP.lb(states)));
                else lb = sum(QP.Aeq(req,states)*QP.lb(states));
                end
                limitL.(Outs{s})(i) = limitL.(Outs{s})(i) + lb;
            end
        end
    end
    sumUB = sum((K(lines+1:end,inc)>0).*(ones(2^n,1)*limitU.(Outs{s})(inc)),2);
    sumLB = sum((K(lines+1:end,inc)>0).*(ones(2^n,1)*limitL.(Outs{s})(inc)),2);
    if ~isempty(ninc)
        if isempty(sumUB)
            sumUB = 0;
        end
        if isempty(sumLB)
            sumLB = 0;
        end
        sumUB = sumUB + sum((K(lines+1:end,ninc)>0).*(ones(2^n,1)*limitU.(Outs{s})(ninc)),2);
        sumLB = sumLB + sum((K(lines+1:end,ninc)>0).*(ones(2^n,1)*limitL.(Outs{s})(ninc)),2);
    end
    keep = [true(lines,1); ((sumUB>=netDemand.(Outs{s})) & (sumLB<=netDemand.(Outs{s})))];%keep the rows where the ub is capable of meeting demand and the lower bound is low enough to meet demand
    K = K(keep,:);
    lines = nnz(keep);
end
