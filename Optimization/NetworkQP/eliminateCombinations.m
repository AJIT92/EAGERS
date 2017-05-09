function [FeasibleDispatch,Cost,feasible,Iterate,HeatVent,K] = eliminateCombinations(QPsingle,netDemand,On,parallel,K,i,bestCost,dt)
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
FeasibleDispatch = zeros(1,nG+nL);
QP = disableGenerators(QPsingle,[],On);
%             tic
if nnz(QP.H) ~=0 %only do quadratic programming if there are quadratic terms
    [dispatch, result,flag1,Output]  = quadprog(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options); 
elseif ~isempty(QP.f) 
    [dispatch, result,flag1,Output] = linprog(QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options2);
elseif isempty(netDemand)%% need to fix this for when demand ==0
end
%             timeQP(i) = toc;
Iterate = Output.iterations;
if flag1==1
    for j = 1:1:nG
        FeasibleDispatch(j) = sum(dispatch(QP.organize{j})); %record this combination of outputs (SOC for storage)
    end
    for j = nG+1:1:nG+nL
        FeasibleDispatch(j) = dispatch(QP.organize{j}); %record lines
    end
    if isfield(QP.Organize,'HeatVented')
        HeatVent = sum(dispatch(nonzeros(QP.Organize.HeatVented(1,:))));
    end
    
    feasible = true;
    Cost = result+sum(QPsingle.constCost(On>0));
else
    feasible = false;
    Cost = inf;
    HeatVent = 0;
end
FeasibleDispatch(abs(FeasibleDispatch)<1e-3) = 0; %remove tiny outputs because they are most likely rounding errors

%%reduce size of K if you don't have parallel computing toolbox
if ~parallel
    if isempty(bestCost)
        bestCost = inf;
    end
    if Cost < bestCost
        %remove combinations that are the best combination and additional more expensive generators
        Outs = fieldnames(netDemand);
        for s = 1:1:length(Outs)
            req = [];
            if strcmp(Outs{s},'E')
                req = QP.Organize.Balance.Electrical; %rows of Aeq associated with electric demand
            elseif strcmp(Outs{s},'H')
                req = QP.Organize.Balance.DistrictHeat; %rows of Aeq associated with heat demand
            elseif strcmp(Outs{s},'C')
                req = QP.Organize.Balance.DistrictCool; %rows of Aeq associated with cool demand
            elseif strcmp(Outs{s},'W')
                req = QP.Organize.Balance.Hydro; %rows of Aeq associated with cool demand
            end

            if length(req) ==1 % The following definitely works with only 1 node
                bestCostPerkWh = Cost/(netDemand.(Outs{s})*dt);
                nRow = length(K(:,1))-i;%how many rows do you have left
                include = false(nG,1);
                for j = 1:1:nG
                    if QP.Organize.Dispatchable(j)
                        states = QP.Organize.States{j};
                        if any(QP.Aeq(req,states))
                            include(j) = true;
                        end
                    end
                end
                if ~isempty(include) && i<length(K(:,1))
                    posCheaperGen = include(logical((minRate(include)<bestCostPerkWh).*(1-(K(i,include)>0))));%possibly cheaper generators that are not on for this case, but could be on
                else
                    posCheaperGen = [];
                end
                if ~isempty(posCheaperGen)
                    ckRows = ~any((ones(nRow,1)*K(i,include)-K(i+1:end,include))>0,2); %identify future rows that have the current set of active generators + more
                    noCheaper = ~any(K(i+1:end,posCheaperGen),2); %select columns of K with possibly cheaper generators, if row is all zeros it can be eliminated
                    rmRows = ckRows&noCheaper;
                    if nnz(rmRows)>0
                        K = K([true(i,1);~rmRows],:);  
                    end
                end
                %remove combinations that swap the current
                %combination with a more expensive generator
                for m = 1:1:length(include)
                    if K(i,include(m))==1 && i<length(K(:,1))
                        cheapGen = include(m);
                        expensiveGens = minrate(include)>minrate(cheapGen);
                        for q = 1:1:length(expensiveGens)
                            addGen = zeros(1,nG);
                            addGen(expensiveGens(q)) = 1;
                            addGen(cheapGen) = -1;
                            swapK = K(i,:)+addGen;
                            K = K~=swapK;
                        end
                    end
                end
            end
        end
    end
end