function QPall = buildMatrices1Step(dt)

%this function creates the a seperate set of QP matrices for each output
%(electric & heating together, cooling, steam etc)
% It makes a seperate QP matrix for each time interval in the optimization horizon if dt is not constant
% if stor is true then energy storage is considered

global Model_dir Plant

nodes = Plant.nodalVar.nodes;
genNames = Plant.nodalVar.genNames;
genStates = Plant.nodalVar.genStates; 


nG = length(Plant.Generator);
Outs = Plant.optimoptions.Outputs;
seq = 1;

for i = 1:1:nG
    genNames(i,1) = {Plant.Generator(i).Name};
end
nodes = length(Plant.Network);


while seq<=length(Outs)
    xL = 0; %number of states
    req = 1; % row index of the Aeq matrix and beq vector
    r = 0; % row index of the A matrix & b vector
    if strcmp('E',Outs{seq})
        Organize.E.Demand = {'beq',1};
        include = {'CHP Generator', 'Electric Generator'};
    end

    thisSeq = zeros(nG,1);
    stor = zeros(1,nG);
    utility = zeros(1,nG);
    Hratio = zeros(1,nG);
    renew = zeros(1,nG);
    QP.organize = cell(1,nG);
    
    for i = 1:1:nodes
        gen = Plant.Network(i).gen;
        for j = 1:1:length(gen)
            s = strfind(gen{j},'.');
            genName = gen{j}(s+1:end);
            I = find(strcmp(genName,genNames),1,'first');
            states = Plant.Generator(I).OpMatB.states;
            statesIndex = genStates{i,j};
        %% identify system type
            if ismember(Plant.Generator(I).Type,include)
                if ismember(Plant.Generator(I).Type,include)
                    thisSeq(I) = I;
                    if strcmp(Plant.Generator(I).Type,'CHP Generator')
                        Hratio(I) = Plant.Generator(I).OpMatB.output.H;
                    end
                end
            elseif strcmp(Plant.Generator(I).Type,'Utility')
                if isfield(Plant.Generator(I).OpMatB.output,Outs{seq}) 
                    utility(I) = I;
                end
            elseif isfield(Plant.Generator(I).OpMatB,'Stor')
                if isfield(Plant.Generator(I).OpMatB.output,Outs{seq})
                    stor(I) = I; %storage in this seq
                end
            elseif strcmp('E',Outs{seq}) && strcmp(Plant.Generator(I).Source,'Renewable')
                renew(I) = I;
            end

            Organize.(Outs{seq}).IC(I) = 0;
            Organize.(Outs{seq}).States(I) ={''}; 
            Organize.(Outs{seq}).Equalities(I) = {[]};
            Organize.(Outs{seq}).Inequalities(I) = {[]};
            Gen = Plant.Generator(I).OpMatB;
            if thisSeq(I) || utility(I)
                nX = length(Gen.states);
                if isfield(Gen,'link') && isfield(Gen.link,'beq') %not really using this anymore                 
                    neq = length(Gen.link.beq);
                    Organize.(Outs{seq}).Equalities(I) = {[]}; %associated with multiple generators, so leaving blank
                    req = req + neq;
                end
            elseif stor(I)
                nX =1; %only 1 state for storage in the 1-step optimization
            else nX = []; 
            end
            if ~isempty(nX) 
                Organize.(Outs{seq}).States(I) = {xL+1:xL+nX}; 
                QP.organize{I} = linspace(xL+1,xL+nX,nX);%generator with one or multiple states; 
                xL = xL + nX;
            end
        end 
    end

    
    thisSeq = nonzeros(thisSeq)';
    utility = nonzeros(utility)';
    stor = nonzeros(stor)';
    renew = nonzeros(renew)';
    CHPindex = nonzeros((1:nG).*(Hratio>0))';
    Organize.(Outs{seq}).thisSeq = thisSeq;
    Organize.(Outs{seq}).stor = stor;
    Organize.(Outs{seq}).utility = utility;
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
        for i = 1:1:nodes
            gen = Plant.Network(i).gen;
            for j = 1:1:length(gen)
                s = strfind(gen{j},'.');
                genName = gen{j}(s+1:end);
                I = find(strcmp(genName,genNames),1,'first');
                statesIndex = genStates{i,j};
                if ismember(I,[thisSeq, utility]) 
                    Gen = Plant.Generator(I).OpMatB;
                    states = Gen.states;
                    k = Organize.(Outs{seq}).States{I};
                    %gOuts = fieldnames(Gen.output); %use this when set up for heating
                    gOuts = fieldnames(Plant.Generator(3).OpMatB.output); %Placeholder to prevent running 'H'
                    for p = 1:1:length(gOuts) %load all outputs into the correct equality or inequality equation
                        if p<=1
                            mat = Organize.(gOuts{p}).Demand{1};
                            index = Organize.(gOuts{p}).Demand{2};
                            if strcmp(mat,'beq')
                                mat = 'Aeq';
                            else mat = 'A';
                            end
                            QP.(mat)(index,k) = Gen.output.(gOuts{p});
                            if strcmp(Plant.Generator(I).Type,'Utility') && length(k) ==2 %utility with sellback
                                QP.(mat)(index,k(2)) = -Gen.output.(gOuts{p});
                            end
                        end 
                    end
                    %% Generally have done away with link states for non storage
                    if isfield(Gen,'link')  %link is a field if there is more than one state
                        if isfield(Gen.link,'eq')
                            req = Organize.(Outs{seq}).Equalities{I};
                            QP.Aeq(req, k) = Gen.link.eq;
                            QP.beq(req) = Gen.link.beq;
                        end
                        if isfield(Gen.link,'ineq')
                            r = Organize.(Outs{seq}).Inequalities{I};
                            QP.A(r, k) = Gen.link.ineq;
                            QP.b(r) = Gen.link.bineq;
                        end
                    end

                    %% costs and bounds, each state has 1 value
                    k = Organize.(Outs{seq}).States{I}(1);% first index for the states of the x vector in C = x'Hx+f'x
                    for j = 1:1:length(states)
                        H(k+j-1) = Gen.(states{j}).H.*dt(tS); %costs scaled by length of time interval
                        QP.f(k+j-1) = Gen.(states{j}).f.*dt(tS); %costs scaled by length of time interval
                        QP.lb(k+j-1) = Gen.(states{j}).lb;
                        QP.ub(k+j-1) = Gen.(states{j}).ub;
                    end
                elseif ismember(I,[stor]) %% add storage (as single state for power output, not SOC)
                    Gen = Plant.Generator(I).OpMatB;
                    k = Organize.(Outs{seq}).States{I}(1);% first index for the states of the x vector in C = x'Hx+f'x
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
        end 
        QP.H = diag(H);
        if isfield(Organize,'H') && Plant.optimoptions.excessHeat && ~strcmp(Outs{seq},'C')
            QP.A(Organize.H.Demand{2},:) = -QP.A(Organize.H.Demand{2},:);
        end
        QPall.(Outs{seq})(tS) = QP;
    end
    seq = seq+1;
end
QPall.Organize = Organize;      