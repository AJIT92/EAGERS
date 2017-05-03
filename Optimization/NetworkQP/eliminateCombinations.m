function [FeasibleDispatch,Cost,Binary,I,TestedCombos] = eliminateCombinations(QPsingle,netDemand,dt)
%This function identifies all feasible generator combinations that can meet
%the demand at this moment in time
%These feasible combinations are then tested, begining with the options
%that have the fewest number of generators
%Some combinations are avoided if they can be pre-emptively determined to
%be more costly
%The function returns the feasible generator dispatches, the cost of that 
%dispatch, and the on/off binary matrix that represents those combinations
options = optimset('Algorithm','interior-point-convex','MaxIter',100,'Display','none');
options2 = optimset('MaxIter',100,'Display','none');
nG = length(QPsingle.constCost);
[~,n] = size(QPsingle.organize);
nL = n-nG;
minRate = zeros(1,nG);
maxRate = zeros(1,nG);
I = zeros(1,nG);
for i = 1:1:nG
    k = QPsingle.organize{i};
    if ~isempty(k)
        minRate(i) = QPsingle.f(k(1));
        I(i) =  QPsingle.ub(k(1));
        if length(k)==2
            maxRate(i) = QPsingle.f(k(2));
        end
    end
end
%% create a matrix of all possible combinations
inc = find(and(QPsingle.Organize.Dispatchable==1, QPsingle.Organize.Enabled==1));
ninc = find(and(QPsingle.Organize.Dispatchable==0, QPsingle.Organize.Enabled==1));
n = length(inc);
lines = 2^n;
K = zeros(lines,nG); %K is a matrix of all the possible generator on/off combinations 
if ~isempty(ninc)
    K(:,ninc) = ones(lines,1)*(ninc); % all systems that are always included
end
for j = 1:1:n %all combinations of generators are listed below
    z = 2^(n-j);
    r=0;
    while r+z<=lines
        K(r+1:r+z,inc(j)) = inc(j);
        r=r+2*z;
    end
end

%% test if each combination is capable of meeting demand
%if the generators lower bounds are higher than the demand, then the case is invalid and should be removed. 
%if the generators upper bounds are lower than the demand, then this case does not produce enough and should be removed.
Outs = fieldnames(netDemand);
for s = 1:1:length(Outs)
    %could improve this section to check if feasible with network losses (not sure how)
    req = [];
    if strcmp(Outs{s},'E')
        req = QPsingle.Organize.Balance.Electrical; %rows of Aeq associated with electric demand
    elseif strcmp(Outs{s},'H')
        req = QPsingle.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
    elseif strcmp(Outs{s},'C')
        req = QPsingle.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
    end
    Nodes.(Outs{s}) = length(req); %number of nodes in this output category
    limitU.(Outs{s}) = zeros(1,nG);
    limitL.(Outs{s}) = zeros(1,nG);
    for i = 1:1:nG
        states = QPsingle.Organize.States{i};
        if~isempty(states)
            if any(isinf(QPsingle.ub(states)))
                if any(QPsingle.Aeq(req,states)~=0)
                    limitU.(Outs{s})(i) = inf;
                end
            else
                limitU.(Outs{s})(i) = limitU.(Outs{s})(i) + sum(QPsingle.Aeq(req,states)*QPsingle.ub(states));
            end
            if any(isinf(QPsingle.lb(states)))
                if any(QPsingle.Aeq(req,states)~=0)
                    limitL.(Outs{s})(i) = -inf;
                end
            else
                if strcmp(Outs{s},'H') && isfield(QPsingle.Organize,'HeatVented') && any(QPsingle.Organize.HeatVented(1,:))
                    lb = min(0,sum(QPsingle.Aeq(req,states)*QPsingle.lb(states)));
                else lb = sum(QPsingle.Aeq(req,states)*QPsingle.lb(states));
                end
                limitL.(Outs{s})(i) = limitL.(Outs{s})(i) + lb;
            end
        end
    end
    sumUB = sum((K(:,inc)>0).*(ones(lines,1)*limitU.(Outs{s})(inc)),2);
    sumLB = sum((K(:,inc)>0).*(ones(lines,1)*limitL.(Outs{s})(inc)),2);
    if ~isempty(ninc)
        if isempty(sumUB)
            sumUB = 0;
        end
        if isempty(sumLB)
            sumLB = 0;
        end
        sumUB = sumUB + sum((K(:,ninc)>0).*(ones(lines,1)*limitU.(Outs{s})(ninc)),2);
        sumLB = sumLB + sum((K(:,ninc)>0).*(ones(lines,1)*limitL.(Outs{s})(ninc)),2);
    end
    keep = (sumUB>=netDemand.(Outs{s})).*(sumLB<=netDemand.(Outs{s}));%keep the rows where the ub is capable of meeting demand and the lower bound is low enough to meet demand
    lines = nnz(keep);
    K = K(keep>0,:);
end

bestCost = inf;%initialize result
row = 1;
FeasibleDispatch = zeros(lines,nG+nL);
HeatVent = zeros(lines,1);
Cost = zeros(lines,1);
Binary = true(lines,nG);
feas = 0;
I = [];
if ~isempty(K) 
    %sort the combinations of generators by smallest number of generators
    nzK = sum(K>0,2);%this is the number of active generators per combination (nonzeros of K)
    [~, line] = sort(nzK); %sort the rows by number of generators that are on
    K = K(line,:);
    %% test the cases for cost
    while row<=length(K(:,1)) %run the quadprog/linprog for all cases with the least number of generators
        QP = disableGenerators(QPsingle,[],K(row,:));
        if nnz(QP.H) ~=0 %only do quadratic programming if there are quadratic terms
            [dispatch, result,flag1]  = quadprog(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options); 
        elseif ~isempty(QP.f) 
            [dispatch, result,flag1] = linprog(QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options2);
        elseif Demand==0
            flag1=1;
            result = 0;
            dispatch = [];%% need to fix this for when demand ==0
        end
        if flag1 == 1
            feas = feas+1; %count of the feasible options tested
            Binary(feas,:) = K(row,:)>0;
            Cost(feas) = result+sum(QP.constCost(Binary(feas,:))); %add the constant term to cost of generators that are on
            for j = 1:1:nG
                FeasibleDispatch(feas,j) = sum(dispatch(QP.organize{j})); %record this combination of outputs (SOC for storage)
            end
            for j = nG+1:1:nG+nL
                FeasibleDispatch(feas,j) = dispatch(QP.organize{j}); %record lines
            end
            if isfield(QP.Organize,'HeatVented')
                HeatVent(feas) = sum(dispatch(nonzeros(QP.Organize.HeatVented(1,:))));
            end
            if Cost(feas) < bestCost
                bestCost = Cost(feas);
                I = feas;
                if row<length(K(:,1))%no need to do for last row
                    %remove combinations that are the best combination and additional more expensive generators
                    for s = 1:1:length(Outs)
                        if Nodes.(Outs{s}) ==1 % The following definitely works with only 1 node
                            bestCostPerkWh = Cost(feas)/(netDemand.(Outs{s})*dt);
                            nRow = length(K(:,1))-row;%how many rows do you have left
                            include = inc(limitU.(Outs{s})(inc)>0);%only check generators for this type
                            if ~isempty(inc) && row<length(K(:,1))
                                posCheaperGen = include(logical((minRate(include)<bestCostPerkWh).*(1-(K(row,include)>0))));%possibly cheaper generators that are not on for this case, but could be on
                            else
                                posCheaperGen = [];
                            end
                            if ~isempty(posCheaperGen)
                                ckRows = ~any((ones(nRow,1)*K(row,include)-K(row+1:end,include))>0,2); %identify future rows that have the current set of active generators + more
                                noCheaper = ~any(K(row+1:end,posCheaperGen),2); %select columns of K with possibly cheaper generators, if row is all zeros it can be eliminated
                                rmRows = ckRows&noCheaper;
                                if nnz(rmRows)>0
                                    K = K([true(row,1);~rmRows],:);  
                                end
                            end
                            %remove combinations that swap the current
                            %combination with a more expensive generator
                            for m = 1:1:length(include)
                                if K(row,include(m))==1 && row<length(K(:,1))
                                    cheapGen = include(m);
                                    expensiveGens = minrate(include)>minrate(cheapGen);
                                    for q = 1:1:length(expensiveGens)
                                        addGen = zeros(1,nG);
                                        addGen(expensiveGens(q)) = 1;
                                        addGen(cheapGen) = -1;
                                        swapK = K(row,:)+addGen;
                                        K = K~=swapK;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        else
            %% not converging
        end
        row = row+1;
    end
    if feas == 0
        uiwait(msgbox('no feasible combination found in eliminateCombinations line 160','Error','modal'))
    else
        FeasibleDispatch(abs(FeasibleDispatch)<1e-3) = 0; %remove tiny outputs because they are most likely rounding errors
    end
end
Cost = Cost(1:feas);
Binary = Binary(1:feas,:);
FeasibleDispatch = FeasibleDispatch(1:feas,:);
TestedCombos = row-1;