function QP = disableGenerators(QP,Organize,Locked,Enabled)
% QP are the updated QP matrices
% Organize is the record of which indices are associated with each generator
% Locked is a matrix nSx nG with zeros when a generator is locked off
% Enabled is a variable to include each generator in the optimization or not.
global Plant
nodes = Plant.nodalVar.nodes;
genNames = Plant.nodalVar.genNames;


if ~isempty(Locked)
    [m,n] = size(Locked);
    for j = 1:1:n
        LB = sum(QP.lb(QP.organize{2,j}));
        
        %% modify lower bounds if ramp rates are to slow (or steps too small)
        r = Organize.Inequalities{j};
        if~isempty(r)
            l = length(r);
            r = r(l-2*(m-1)+1:l);
            index = (1:m-1)';
            starts = nonzeros(index.*((Locked(2:end,j)-Locked(1:m-1,j))>0));
            stops = nonzeros(index.*(Locked(1:m-1,j)-(Locked(2:end,j))>0));
            if ~isempty(starts)
                % example, the 4th index of locked is one and it takes 3 steps to turn on (0, .25LB, .7LB, >LB), then starts(1) = 3, n = 3, e = 3. We only need to change LB at 4 & 5 (33% and 67%), 3 was already set to zero
                for k = 1:1:length(starts)
                    e = starts(k);
                    n = 1;
                    RampUp = QP.b(r(2*e-1));
                    while RampUp<LB && (e+n)<(m-1)
                        RampUp = RampUp+QP.b(r(2*(e+n)-1));
                        n= n+1; %# of steps to turn on
                    end
                    if RampUp>=LB 
                        for t =1:1:n-1
                            nX = QP.organize{e+t,j};
                            RampUp = RampUp-QP.b(r(2*(e+t-1)-1));
                            QP.lb(nX) = (LB-RampUp)*(QP.lb(nX)/LB); %lb to ensure that it is off when locked =0
                        end
                    else
                        RampUp = 0; %don't have enough tome to fully come on, ramp up as much as possible
                        for t =1:1:n-1
                            nX = QP.organize{e+t,j};
                            RampUp = RampUp-QP.b(r(2*(e+t-1)-1));
                            QP.lb(nX) = RampUp*(QP.lb(nX)/LB); %lb to ensure that it is off when locked =0
                        end
                    end
                end
            end
            if ~isempty(stops)
                % example, the 6th index of locked is zero and it takes 3 steps to shut down (LB, .67LB, .33LB, 0), then stops(1) = 5, p = 3, e = 3. We only need to change LB at 4 & 5 (67% and 33%), 6 was already set to zero
                for k = 1:1:length(stops)
                    e = stops(k);
                    p=1;
                    RampDown = QP.b(r(2*e));
                    while RampDown<LB && (e-p)>0
                        RampDown = RampDown+QP.b(r(2*(e-p)));
                        p= p+1; %# of steps to turn off
                    end
                    e = e-p+1;
                    RampDown = QP.b(r(2*e));
                    for t =1:1:p-1
                        nX = QP.organize{e+t,j};
                        QP.lb(nX) = (LB-RampDown)*(QP.lb(nX)/LB); %lb to ensure that it is off when locked =0
                        RampDown = RampDown+QP.b(r(2*(e+t)));
                    end
                    if RampDown<LB
                        t=p;
                        while RampDown<LB && (e+t)<(m-1)
                            Locked(e+t,j) = true; %need to force gen on, cant ramp down soon enough
                            nX = QP.organize{e+t,j};
                            QP.lb(nX) = (LB-RampDown)*(QP.lb(nX)/LB); %lb to ensure that it is off when locked =0
                            RampDown = RampDown+QP.b(r(2*(e+t)));
                            t = t+1;
                        end
                    end
                end
            end
        end
        for i = 1:1:m
            if ~Locked(i,j)
                nX = QP.organize{i,j};
                QP.lb(nX) = 0;
                QP.ub(nX) = 0;
            end
        end
    end
end

%remove disabled generators from optimization matrices
if nnz(Enabled)<length(Enabled) %may be some disabled generators
    [r,~] = size(QP.A);
    [req,xL] = size(QP.Aeq);
    rkeep = linspace(1,r,r)';
    reqkeep = linspace(1,req,req)';
    xkeep = linspace(1,xL,xL)';
    rmStates = 0;
    [m,n] = size(QP.organize);
    for t = 1:1:m
        for i = 1:1:nodes
            gen = Plant.Network(i).gen;
            for j = 1:1:length(gen)
                s = strfind(gen{j},'.');
                genName = gen{j}(s+1:end);
                I = find(strcmp(genName,genNames),1,'first');
                if Enabled(I) ==0
                    states = Organize.States{I}';
                    xkeep(states)=0;
                    if ~isempty(Organize.Equalities{I})
                        reqkeep(Organize.Equalities{I}')=0;
                    end
                    if ~isempty(Organize.Inequalities{I})
                        rkeep(Organize.Inequalities{I}')=0;
                    end
                    if ~isempty(states)
                        cell = length(QP.organize{t,I});
                    else cell = length(states);
                    end 
                       rmStates = rmStates + cell;
                       QP.organize{t,I} = [];
                else
                    if ~isempty(QP.organize{t,I})
                        QP.organize{t,I} = QP.organize{t,I} - rmStates; %reduce the coordinate of other generators if generators ahead of it were removed
                    end
                end
            end 
        end 
    end 
    rkeep = nonzeros(rkeep);
    reqkeep = nonzeros(reqkeep);
    xkeep = nonzeros(xkeep);
    QP.H = QP.H(xkeep,xkeep);
    QP.f = QP.f(xkeep);
    QP.Aeq = QP.Aeq(reqkeep,xkeep);
    QP.beq = QP.beq(reqkeep);
    if r>0 && nnz(rkeep)>0
        QP.A = QP.A(rkeep,xkeep);
        QP.b = QP.b(rkeep);
    else QP.A = [];
        QP.b = [];
    end
    QP.lb = QP.lb(xkeep);
    QP.ub = QP.ub(xkeep);
end