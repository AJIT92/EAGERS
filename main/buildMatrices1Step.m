function QPall = buildMatrices1Step(dt)
%this function creates the a seperate set of QP matrices for each output
%(electric & heating together, cooling, steam etc)
% It makes a seperate QP matrix for each time interval in the optimization horizon if dt is not constant
% if stor is true then energy storage is considered
global Plant 
nG = length(Plant.Generator);
Outs = Plant.optimoptions.Outputs;
seq = 1;
while seq<=length(Outs)
    xL = 0;
    req = 1; % row index of the Aeq matrix and beq vector
    r = 0; % row index of the A matrix & b vector
    if strcmp('E',Outs{seq})
        Organize.E.Demand = {'beq',1};
        include = {'CHP Generator', 'Electric Generator'};
        if nnz(strcmp('H',Outs))>0
            include(end+1) = {'Heater'};
            Outs = Outs(~(strcmp('H',Outs))); %combine heating optimization with electric due to CHP generators
            if Plant.optimoptions.excessHeat == 1
                Organize.H.Demand = {'b',1};
                r = 1;
            else Organize.H.Demand = {'beq',2};
                req = 2;
            end 
        end
        if nnz(strcmp('C',Outs))>0 && Plant.optimoptions.sequential == 0 %%load the chillers separately if the user has selected the sequential option
            include(end+1) = {'Chiller'};
            Outs = Outs(~(strcmp('C',Outs))); 
            req = req+1;
            Organize.C.Demand = {'beq',req};
        end
    elseif strcmp('C',Outs{seq})>0 && Plant.optimoptions.sequential == 1
        include = {'Chiller'};
        Organize.C.Demand = {'beq',1};
    end

    thisSeq = zeros(nG,1);
    stor = zeros(1,nG);
    storC = zeros(1,nG);
    storH = zeros(1,nG);
    chill = zeros(1,nG);
    heater = zeros(1,nG);
    utility = zeros(1,nG);
    utilC = zeros(1,nG);
    utilH = zeros(1,nG);
    Hratio = zeros(1,nG);
    renew = zeros(1,nG);
    QP.organize = cell(1,nG);
    for i = 1:1:nG
        %% identify system type
        if ismember(Plant.Generator(i).Type,include)
            if Plant.optimoptions.sequential == 0 && strcmp(Plant.Generator(i).Type,'Chiller')
                chill(i) = i;
            elseif ~strcmp(Outs{seq},'H') && strcmp(Plant.Generator(i).Type,'Heater')
                heater(i) = i;
            else thisSeq(i) = i;
                if strcmp(Plant.Generator(i).Type,'CHP Generator')
                    Hratio(i) = Plant.Generator(i).OpMatB.output.H;
                end
            end
        elseif strcmp(Plant.Generator(i).Type,'Utility')
            if isfield(Plant.Generator(i).OpMatB.output,Outs{seq}) 
                utility(i) =i;
            elseif isfield(Plant.Generator(i).OpMatA.output,'C') && Plant.optimoptions.sequential == 0 
                utilC(i) = i;
            elseif isfield(Plant.Generator(i).OpMatA.output,'H')
                utilH(i) = i;
            end
        elseif isfield(Plant.Generator(i).OpMatB,'Stor')
            if Plant.optimoptions.sequential == 0  && strcmp('E',Outs{seq}) &&  isfield(Plant.Generator(i).OpMatB.output,'C')
                storC(i) = i;
            elseif strcmp('E',Outs{seq}) && (isfield(Plant.Generator(i).OpMatB.output,'H'))
                storH(i) = i;
            elseif isfield(Plant.Generator(i).OpMatB.output,Outs{seq})
                stor(i) = i; %storage in this seq
            end
        elseif strcmp('E',Outs{seq}) && strcmp(Plant.Generator(i).Source,'Renewable')
            renew(i) = i;
        end
        
        Organize.(Outs{seq}).IC(i) = 0;
        Organize.(Outs{seq}).States(i) ={''};
        Organize.(Outs{seq}).Equalities(i) = {[]};
        Organize.(Outs{seq}).Inequalities(i) = {[]};
        Gen = Plant.Generator(i).OpMatB;
        if thisSeq(i) || chill(i) || heater(i) || utility(i) || utilC(i) || utilH(i)
            nX = length(Gen.states);
            if isfield(Gen,'link') && isfield(Gen.link,'beq') %not really using this anymore
                neq = length(Gen.link.beq);
                Organize.(Outs{seq}).Equalities(i) = {req+1:req+neq};
                req = req + neq;
            end
        elseif stor(i) || storC(i) || storH(i)
            nX =1; %only 1 state for storage in the 1-step optimization
        else nX = [];
        end
        if ~isempty(nX)
            Organize.(Outs{seq}).States(i) = {xL+1:xL+nX};
            QP.organize{i} = linspace(xL+1,xL+nX,nX);%generator with one or multiple states
            xL = xL + nX;
        end
    end
    thisSeq = nonzeros(thisSeq)';
    chill = nonzeros(chill)';
    heater = nonzeros(heater)';
    utility = nonzeros(utility)';
    utilC = nonzeros(utilC)';
    utilH = nonzeros(utilH)';
    stor = nonzeros(stor)';
    storC = nonzeros(storC)';
    storH = nonzeros(storH)';
    renew = nonzeros(renew)';
    CHPindex = nonzeros((1:nG).*(Hratio>0))';
    Organize.(Outs{seq}).thisSeq = thisSeq;
    Organize.(Outs{seq}).stor = stor;
    Organize.(Outs{seq}).storC = storC;
    Organize.(Outs{seq}).storH = storH;
    Organize.(Outs{seq}).utility = utility;
    Organize.(Outs{seq}).utilC = utilC;
    Organize.(Outs{seq}).utilH = utilH;
    Organize.(Outs{seq}).chill = chill;
    Organize.(Outs{seq}).heater = heater;
    Organize.(Outs{seq}).renew = renew;
    Organize.(Outs{seq}).CHPindex = CHPindex;
    Organize.(Outs{seq}).Hratio = nonzeros(Hratio)';
    
    %% build matrices
    QP.f = zeros(xL,1);
    QP.Aeq = zeros(req,xL);
    QP.beq = zeros(req,1);
    QP.A = zeros(r,xL);
    QP.b = zeros(r,1);
    QP.lb = zeros(xL,1);
    QP.ub = inf*ones(xL,1);%if a state has now bound, then the upper bound is infinite
    for tS = 1:1:length(dt)
        H = zeros(xL,1);%QP.H needs to be in the loop so that is doesn't change from a diagonal to a vector every other timestep
        for i = 1:1:nG
            if ismember(i,[thisSeq, chill, heater, utility, utilH, utilC]) 
                Gen = Plant.Generator(i).OpMatB;
                states = Gen.states;
                k = Organize.(Outs{seq}).States{i};
                gOuts = fieldnames(Gen.output);
                for p = 1:1:length(gOuts) %load all outputs into the correct equality or inequality equation
                    mat = Organize.(gOuts{p}).Demand{1};
                    index = Organize.(gOuts{p}).Demand{2};
                    if strcmp(mat,'beq')
                        mat = 'Aeq';
                    else mat = 'A';
                    end
                    QP.(mat)(index,k) = Gen.output.(gOuts{p});
                    if strcmp(Plant.Generator(i).Type,'Utility') && length(k) ==2 %utility with sellback
                        QP.(mat)(index,k(2)) = -Gen.output.(gOuts{p});
                    end
                end
                %% Generally have done away with link states for non storage
                if isfield(Gen,'link')  %link is a field if there is more than one state
                    if isfield(Gen.link,'eq')
                        req = Organize.(Outs{seq}).Equalities{i};
                        QP.Aeq(req, k) = Gen.link.eq;
                        QP.beq(req) = Gen.link.beq;
                    end
                    if isfield(Gen.link,'ineq')
                        r = Organize.(Outs{seq}).Inequalities{i};
                        QP.A(r, k) = Gen.link.ineq;
                        QP.b(r) = Gen.link.bineq;
                    end
                end

                %% costs and bounds, each state has 1 value
                k = Organize.(Outs{seq}).States{i}(1);% first index for the states of the x vector in C = x'Hx+f'x
                for j = 1:1:length(states)
                    H(k+j-1) = Gen.(states{j}).H.*dt(tS); %costs scaled by length of time interval
                    QP.f(k+j-1) = Gen.(states{j}).f.*dt(tS); %costs scaled by length of time interval
                    QP.lb(k+j-1) = Gen.(states{j}).lb;
                    QP.ub(k+j-1) = Gen.(states{j}).ub;
                end
            elseif ismember(i,[stor, storC, storH]) %% add storage (as single state for power output, not SOC)
                Gen = Plant.Generator(i).OpMatB;
                k = Organize.(Outs{seq}).States{i}(1);% first index for the states of the x vector in C = x'Hx+f'x
                StorOut = fieldnames(Gen.output);
                for p = 1:1:length(StorOut)
                    mat = Organize.(StorOut{p}).Demand{1};
                    index = Organize.(StorOut{p}).Demand{2};
                    if strcmp(mat,'beq')
                        mat = 'Aeq';
                    else mat = 'A';
                    end
                    QP.(mat)(index,k) = 1;
                end
                %% bounds
                QP.ub(k) = Gen.Ramp.b(2); % peak power production (peak discharge)
                QP.lb(k) = -Gen.Ramp.b(1); % peak charging (negative power)
                %% Costs are figured out during the update stage
            end
        end
        QP.H = diag(H);
        QPall.(Outs{seq})(tS) = QP;
    end
    seq = seq+1;
end
QPall.Organize = Organize;      