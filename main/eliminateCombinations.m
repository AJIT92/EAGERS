function [BestDispatch,tests] = eliminateCombinations(QPsingle,Organize,Heating,Cooling,constCost,dt,IC,MaxPower,MinPower,gen,stor,utility)
%eliminates combinations of generators that are not feasible according to
%their constraints or because there is another combination that is cheaper
%it then returns the minimum cost dispatched and the generator outputs that
%correspond with that cost.
global UB LB
options = optimset('Algorithm','interior-point-convex','MaxIter',100,'Display','none');
options2 = optimset('MaxIter',100,'Display','none');
nG = length(MaxPower);
BestDispatch = [];
minRate = zeros(1,nG);
maxRate = zeros(1,nG);
maxRateQ = zeros(1,nG);
I = zeros(1,nG);
for j = 1:1:length(gen)
    i = gen(j);
    k = QPsingle.organize{i};
    if ~isempty(k)
        minRate(i) = QPsingle.f(k(1));
        I(i) =  QPsingle.ub(k(1));
        if length(k)==2
            maxRate(i) = QPsingle.f(k(2));
            maxRateQ(i) = QPsingle.H(k(2),k(2));
        end
    end
end
[~,index] = sort(minRate(gen));    %sort least expensive to most expensive
gen = gen(index);
ninc =[];
%% sort generator order
if ~isempty(Heating) %always enable heaters (negligable assume start-up costs and LB is usually ~0)
    mat = Heating.Organize.Demand{1};
    index = Heating.Organize.Demand{2};
    QPsingle.(mat)(index) = Heating.Demand;
    heater = Heating.gen;
    storH = Heating.stor;
    utilH = Heating.util;
    CHPindex = Heating.chp;
    Hratio = Heating.Hratio;
    ninc = [heater storH utilH]; %lock these on.
else heater = [];
end
if~isempty(Cooling) %find best cooling combo first
    chill = Cooling.gen;
    storC = Cooling.stor;
    utilC = Cooling.util;
    mat = Cooling.Organize.Demand{1};
    index = Cooling.Organize.Demand{2};
    QPsingle.(mat)(index) = Cooling.Demand;
    if strcmp(mat,'beq')
        mat2 = 'Aeq';
    else mat2 = 'A';
    end
    for j = 1:1:length(chill)
        i = chill(j);
        k = Organize.States{i};
        COP = QPsingle.(mat2)(index,k);
        minRate(i) = 1/COP*min(minRate(gen));
    end
    for j = 1:1:length(storC)
        i = storC(j);
        k = Organize.States{i};
        minRate(i) = QPsingle.f(k);
    end
    k = 1;
    while k<=length(chill)
        i = chill(k);
        if MinPower(i)>=LB(i)%if the initial condition is greater than the lower bound plus the ramp rate for that time step then lock it on
            chill = chill(chill~=i);%remove rom the list of generators that can be turned on/off
            ninc(end+1) = i; %force it to always be on
        else k = k+1;
        end
    end
    
    %% if it is here we are doing chilling with electric, and constant COP is assumed
    %find combination that meets demand with constraints enforced (easy since linear cost)
    UBcool = zeros(nG,1);
    LBcool = zeros(nG,1);
    for j = 1:1:length(chill)
        i = chill(j);
        UBcool(i) = MaxPower(i);
        LBcool(i) = MinPower(i);
    end
    for j = 1:1:length(Cooling.stor)
        i = Cooling.stor(j);
        UBcool(i) = MaxPower(i);
        LBcool(i) = MinPower(i);
    end
    UBcool(utilC) = UB(utilC);
    LBcool(utilC) = LB(utilC);
    %% maybe take this part out, but it should speed it up in some cases.
    if isempty(storC) %no storage, pre-emptively solve chillers
        [~,index] = sort(minRate(chill));    %sort least expensive to most expensive
        chill = chill(index);
        %% create a matrix of all possible combinations
        n = length(chill);
        lines = 2^n;
        K = zeros(lines,nG); %K is a matrix of all the possible generator on/off combinations 
        if ~isempty(Cooling.util)
            K(:,Cooling.util) = ones(lines,1)*(Cooling.util);
        end
        if ~isempty(Cooling.stor)
            K(:,Cooling.stor) = ones(lines,1)*(Cooling.stor);
        end
        for j = 1:1:n %all combinations of generators are listed below
            z = 2^(n-j);
            r=0;
            while r+z<=lines
                K(r+1:r+z,chill(j)) = chill(j);
                r=r+2*z;
            end
        end
        sumUB = (K>0)*UBcool;
        sumLB = (K>0)*LBcool;
        keep = (sumUB>=Cooling.Demand).*(sumLB<=Cooling.Demand);%keep the rows where the ub is capable of meeting demand and the lower bound is low enough to meet demand
        K = K(keep>0,:);
        %sort the combinations of generators by smallest number of generators
        nzK = sum(K>0,2);%this is the number of active generators per combination (nonzeros of K)
        [nzK, line] = sort(nzK); %sort the rows by number of generators that are on
        K = K(line,:);
        K = K(nzK==min(nzK),:);
        cost = K*minRate';
        [~,costIndex] = min(cost);
        chill = nonzeros(K(costIndex(1),chill))'; %first combo that fits UB and LB should be cheapest.
    else
        ninc = [ninc storC utilC];
    end
else
    chill = [];
    storC = [];
    utilC = [];
end

allGen = [gen, chill, heater];
Demand = QPsingle.beq(1);
for i = 1:1:nG %change upper and lower bounds
    k = Organize.States{i};
    if ismember(i,allGen)
        ublimit = MaxPower(i);
        lblimit = max(MinPower(i),QPsingle.lb(k(1)));
        for j = 1:1:length(k) %generators with multiple states are picewise functions starting from 0 and adding up to UB
            QPsingle.ub(k(j)) = min(QPsingle.ub(k(j)),ublimit);
            ublimit = ublimit -QPsingle.ub(k(j));
            QPsingle.lb(k(j)) = min(QPsingle.lb(k(j)),lblimit);
            lblimit = lblimit -QPsingle.lb(k(j));
        end
        if ismember(i,gen) && MinPower(i)>=LB(i)%if the initial condition is greater than the lower bound plus the ramp rate for that time step then lock it on
            gen = gen(gen~=i);%remove rom the list of generators that can be turned on/off
            ninc(end+1) = i; %force it to always be on
        end
    elseif MaxPower(i)>0 %storage
        if Organize.IC(i)>0 %only in threshold case
            QPsingle.beq(Organize.IC(i)) = IC(i);
        end
        QPsingle.ub(k(1)) = MaxPower(i);
        QPsingle.lb(k(1)) = MinPower(i);
    end
end

%% create a matrix of all possible combinations
GandC = [gen, chill];
ninc = [ninc stor utility];
n = length(GandC);
lines = 2^n;
K = zeros(lines,nG); %K is a matrix of all the possible generator on/off combinations 
if ~isempty(ninc)
    K(:,ninc) = ones(lines,1)*(ninc); % all systems that are always included
end
for j = 1:1:n %all combinations of generators are listed below
    z = 2^(n-j);
    r=0;
    while r+z<=lines
        K(r+1:r+z,GandC(j)) = GandC(j);
        r=r+2*z;
    end
end
if Demand~=0 && nnz(K(end,:))==0
    K = K(1:end-1,:);%eliminate zero generator zero case when no grid or storage and there is demand
end
%% test if each combination is capable of meeting demand
%if the generators lower bounds are higher than the demand, then the case is invalid and should be removed. 
%if the generators upper bounds are lower than the demand, then this case does not produce enough and should be removed.
sumUB = (K(:,gen)>0)*MaxPower(gen)';
sumLB = (K(:,(gen))>0)*LB(gen)';
if ~isempty(ninc)
    sumUB = sumUB + (K(:,ninc))*MaxPower(ninc)';
    sumLB = sumLB + (K(:,ninc))*MinPower(ninc)';
end
keep = (sumUB>=Demand).*(sumLB<=Demand);%keep the rows where the ub is capable of meeting demand and the lower bound is low enough to meet demand
K = K(keep>0,:);

%for the dual heating + electric case, remove combos that can't meet both demands
if ~isempty(Heating) %combined heating and electric generators
    UBheat = zeros(nG,1);
    LBheat = zeros(nG,1);
    UBheat([heater storH utilH]) = MaxPower([heater storH utilH]);
    LBheat([heater storH utilH]) = MinPower([heater storH utilH]);
    UBheat(CHPindex) = MaxPower(CHPindex).*Hratio;
    LBheat(CHPindex) = MinPower(CHPindex).*Hratio;
    sumUB = (K>0)*UBheat;
    sumLB = (K>0)*LBheat;
    keep = (sumUB>=Heating.Demand).*(sumLB<=Heating.Demand);%keep the rows where the ub is capable of meeting demand and the lower bound is low enough to meet demand
    K = K(keep>0,:);
end

% for the combine cooling and electric case, best chiller combo was already determined
if ~isempty(Cooling) %combined heating and electric generators
    if isempty(storC) %pre-determined chiller dispatch 
        lines = length(K(:,1));
        K(:,[chill storC utilC]) = ones(lines,1)*[chill storC utilC]; %lock on the chosen chillers
        chill = [];
    end
    sumUB = (K(:,[chill,storC,utilC])>0)*UBcool([chill,storC,utilC]);
    sumLB = (K(:,([chill,storC,utilC]))>0)*LBcool([chill,storC,utilC]);
    keep = (sumUB>=Cooling.Demand).*(sumLB<=Cooling.Demand);%keep the rows where the ub is capable of meeting demand and the lower bound is low enough to meet demand
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
        QP = disableGenerators(QPsingle,Organize,[],K(row,:));
        onNow = nonzeros(K(row,gen))';
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
            result2 = result2+sum(constCost(logical(K(row,:)>0))); %add the constant term to cost of generators that are on
            currentDispatch = zeros(1,nG);
            for j = 1:1:nG
                currentDispatch(j) = sum(dispatch(QP.organize{j})); %record this combination of outputs (SOC for storage)
            end
            if row<length(K(:,1)) % remove any combinations that swap one of the current generators with a more expensive one
                actRate = zeros(1,nG);
                costNow = 0;
                for j = 1:1:length(onNow)
                    i = onNow(j);
                    pow = currentDispatch(j);
                    if pow>I(i)
                        actRate(i) = maxRate(i) + (pow-I(i))*maxRateQ(i);
                        costNow = costNow +  dt*(maxRate(i)*I(i) + (pow-I(i))^2*maxRateQ(i)) + constCost(i);
                    else
                        actRate(i) = minRate(i);
                        costNow = costNow + minRate(i)*I(i)*dt + constCost(i);
                    end
                end
                % calculate chiller cost to rule out those?
                
                if ~isempty(utility)
                    for j=1:1:length(utility)
                        i = utility(j);
                        k = QP.organize{i};
                        costNow = costNow + dispatch(k(1))*QP.f(k(1));
                    end
                end
                bestCostPerkW = costNow/(sum(currentDispatch(onNow))*dt);
            end
            if result2 < result
                if flag1 == 1
                    result = result2;
%                     bestcombo = nonzeros(K(row,:))';
                    BestDispatch = currentDispatch;
                    if row<length(K(:,1))%no need to do for last row
                        %remove combinations that are the best combination and additional more expensive generators
                        nRow = length(K(:,1))-row;
                        posCheaperGen = gen(logical((minRate(gen)<bestCostPerkW).*(1-(K(row,gen)>0))));%possibly cheaper generators that are not on for this case, but could be on
                        if ~isempty(posCheaperGen)
                            ckRows = ~any((ones(nRow,1)*K(row,gen)-K(row+1:end,gen))>0,2); %identify future rows that have the current set of active generators + more
                            noCheaper = ~any(K(row+1:end,posCheaperGen),2); %select columns of K with possibly cheaper generators, if row is all zeros it can be eliminated
                            rmRows = ckRows&noCheaper;
                            if nnz(rmRows)>0
                                K = K([true(row,1);~rmRows],:);  
                            end
                        end
                    end
                elseif flag1 ==0
                    %% potential good solution, but not converging (should look into this)
                end
%             elseif flag1 ==1 && row<length(K(:,1)) %remove combinations that contain worse combinations %no need to do for last row
%                 badcombo = nonzeros(K(row,:))';%remove combinations that contain worse combinations
%                 rmRows = ismember(K(row+1:end,badcombo),badcombo,'rows');
%                 dontRmv = ismember(K(row+1:end,bestcombo),bestcombo,'rows');%rows with the current best set of gen + extra
%                 rmRows(rmRows) = rmRows(rmRows)-dontRmv(rmRows); %dont remove rows that are the current bad combination & the current best combination
%                 if nnz(rmRows)>0
%                     K = K([true(row,1);~rmRows],:);  
%                 end
            elseif flag1==0
                %% not converging, but solution apears to be worse than the best, so ignore
            end
            if row<length(K(:,1)) % remove any combinations that swap one of the current generators with a more expensive one
                for j = 1:1:length(onNow)
                    i = onNow(j);
                    posCheaperGen = nonzeros(gen.*(minRate(gen)<actRate(i)).*(1-(K(row,gen)>0)));%possibly cheaper generators that are not on for this case, but could be on
                    if ~isempty(posCheaperGen)
                        ckRows = all(K(row+1:end,onNow(onNow~=i)),2);%rows that have everything from onNow except j
                        noCheaper = ~any(K(row+1:end,posCheaperGen),2);%rows that dont add a cheaper gen
                        rmRows = ckRows&noCheaper;
                        if nnz(rmRows)>0
                            K = K([true(row,1);~rmRows],:);  
                        end
                    end
                end
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