function [BestDispatch,tests] = eliminateCombinations(QPsingle,netDemand)
%eliminates combinations of generators that are not feasible according to
%their constraints or because there is another combination that is cheaper
%it then returns the minimum cost dispatched and the generator outputs that
%correspond with that cost.
options = optimset('Algorithm','interior-point-convex','MaxIter',100,'Display','none');
options2 = optimset('MaxIter',100,'Display','none');
nG = length(QPsingle.constCost);
BestDispatch = [];
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
inc = find(QPsingle.Organize.Dispatchable==1);
ninc = find(QPsingle.Organize.Dispatchable==0);
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
    r =[];
    req = [];
    if strcmp(Outs{s},'E')
        req = QPsingle.Organize.Balance.Electrical; %rows of Aeq associated with electric demand
    elseif strcmp(Outs{s},'H')
        if QPsingle.excessHeat
            r = QPsingle.Organize.Imbalance.DistrictHeat; %rows of A associated with heat demand
        else
            req = QPsingle.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
        end
    elseif strcmp(Outs{s},'C')
        req = QPsingle.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
    end
    limitU.(Outs{s}) = zeros(1,nG);
    limitL.(Outs{s}) = zeros(1,nG);
    for i = 1:1:nG
        states = QPsingle.Organize.States{i};
        if~isempty(states)
            if any(isinf(QPsingle.ub(states)))
                if any(QPsingle.Aeq(req,states)~=0) || any(QPsingle.A(r,states)~=0)
                    limitU.(Outs{s})(i) = inf;
                end
            else
                if ~isempty(req)
                    limitU.(Outs{s})(i) = limitU.(Outs{s})(i) + sum(QPsingle.Aeq(req,states)*QPsingle.ub(states));
                end
                if ~isempty(r)
                    limitU.(Outs{s})(i) = limitU.(Outs{s})(i) + sum(QPsingle.A(r,states)*QPsingle.ub(states));
                end

            end
            if any(isinf(QPsingle.lb(states)))
                if any(QPsingle.Aeq(req,states)~=0) || any(QPsingle.A(r,states)~=0)
                    limitL.(Outs{s})(i) = -inf;
                end
            else
                if ~isempty(req)
                    limitL.(Outs{s})(i) = limitL.(Outs{s})(i) + sum(QPsingle.Aeq(req,states)*QPsingle.lb(states));
                end
                if ~isempty(r)
                    limitL.(Outs{s})(i) = limitL.(Outs{s})(i) + sum(QPsingle.A(r,states)*QPsingle.lb(states));
                end
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

result = inf;%initialize result
row = 1;
if ~isempty(K) 
    %sort the combinations of generators by smallest number of generators
    nzK = sum(K>0,2);%this is the number of active generators per combination (nonzeros of K)
    [~, line] = sort(nzK); %sort the rows by number of generators that are on
    K = K(line,:);
    %% test the cases for cost
    while row<=length(K(:,1)) %run the quadprog/linprog for all cases with the least number of generators
        QP = disableGenerators(QPsingle,[],K(row,:));
        if nnz(QP.H) ~=0 %only do quadratic programming if there are quadratic terms
            [dispatch, result2,flag1]  = quadprog(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options); 
        elseif ~isempty(QP.f) 
            [dispatch, result2,flag1] = linprog(QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options2);
        elseif Demand==0
            flag1=1;
            result2 = 0;
            dispatch = [];
        end
        if ~isempty(result2)
            result2 = result2+sum(QP.constCost(logical(K(row,:)>0))); %add the constant term to cost of generators that are on
            currentDispatch = zeros(1,nG);
            for j = 1:1:nG
                currentDispatch(j) = sum(dispatch(QP.organize{j})); %record this combination of outputs (SOC for storage)
            end
            if result2 < result
                if flag1 == 1
                    result = result2;
                    BestDispatch = currentDispatch;
                    %% need to improve the following to distinguish between generator categories
%                     if row<length(K(:,1))%no need to do for last row
%                         %remove combinations that are the best combination and additional more expensive generators
%                         for s = 1:1:length(Outs)
%                             bestCostPerkW = result2/(netDemand.(Outs{s})*dt);
%                             nRow = length(K(:,1))-row;
%                             if ~isempty(inc)
%                                 posCheaperGen = inc(logical((minRate(inc)<bestCostPerkW).*(1-(K(row,inc)>0))));%possibly cheaper generators that are not on for this case, but could be on
%                             else
%                                 posCheaperGen = [];
%                             end
%                             if ~isempty(posCheaperGen)
%                                 ckRows = ~any((ones(nRow,1)*K(row,inc)-K(row+1:end,inc))>0,2); %identify future rows that have the current set of active generators + more
%                                 noCheaper = ~any(K(row+1:end,posCheaperGen),2); %select columns of K with possibly cheaper generators, if row is all zeros it can be eliminated
%                                 rmRows = ckRows&noCheaper;
%                                 if nnz(rmRows)>0
%                                     K = K([true(row,1);~rmRows],:);  
%                                 end
%                             end
%                         end
%                     end
                elseif flag1 ==0
                    %% potential good solution, but not converging (should look into this)
                end
            elseif flag1==0
                %% not converging, but solution apears to be worse than the best, so ignore
            end
        else
            %should never get here, should investigate

        end
        row = row+1;
    end
    if isempty(BestDispatch)
        uiwait(msgbox('no feasible combination found in eliminateCombinations line 315','Error','modal'))
    else
        BestDispatch(abs(BestDispatch)<1e-3) = 0; %remove tiny outputs because they are most likely rounding errors
    end
end
tests = row-1;