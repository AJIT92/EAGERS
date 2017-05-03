function QP = disableGenerators(QP,Locked,Enabled)
% QP are the updated QP matrices
% Organize is the record of which indices are associated with each generator
% Locked is a matrix nSx nG with zeros when a generator is locked off
% Enabled is a variable to include each generator in the optimization or not.
[m,n] = size(QP.organize);
nS = m-1;
nG = length(QP.Organize.Dispatchable);
nL = n - nG;
if ~isempty(Locked)
    Enabled = ones(nG,1);
    for j = 1:1:nG
        if all(Locked(:,j)==0)
                Enabled(j) = 0;
        elseif any(Locked(:,j)==0) %is generator off at any point
            ineqRows = [];
            if isfield(QP.Organize,'SpinReserve')
                ineqRows = QP.Organize.SpinReserve{j};
                SRstates = QP.Organize.SpinReserveStates(:,j);
            end
            for t = 1:1:nS
                if ~Locked(t+1,j)
                    QP.lb(QP.organize{t+1,j}) = 0;
                    QP.ub(QP.organize{t+1,j}) = 0;
                    if ~isempty(ineqRows)
                        QP.A(ineqRows(t),SRstates(t)) = 0;%%remove from calculation of spinning reserve
                    end
                end
            end
            %% modify lower bounds if ramp rates are to slow (or steps too small)
            Pmin = QP.lb(QP.organize{2,j});
            rows = QP.Organize.Ramping{j};
            states = QP.Organize.States{j};
            nt = length(states)/nS; %number of states per timestep
            starts = nonzeros((1:nS)'.*((Locked(2:end,j)-Locked(1:nS,j))>0));
            if ~isempty(starts)
                % example, the 4th index of locked is one and it takes 3 steps to turn on (0, .25LB, .7LB, >LB), then starts(1) = 3, n = 3, e = 3. We only need to change LB at 4 & 5 (33% and 67%), 3 was already set to zero
                for k = 1:1:length(starts)
                    e = starts(k);
                    if e<nS %cant change lb after nS, so nothing to do if starting at last step
                        n = 0;
                        Pmax = QP.b(rows(2*e-1));
                        for f = 1:1:nt %starting with lb of 1st state in generator, then moving on
                            while Pmax<sum(Pmin(1:f)) && (e+n)<nS
                                QP.lb(states(nt*(e+n-1)+f)) = Pmax - sum(QP.lb(states(nt*(e+n-1)+(1:(f-1))))); %lb to ensure that it is off when locked =0
                                QP.lb(states(nt*(e+n-1)+(f+1):nt)) = 0; %other states of this generator during start-up
                                n= n+1; %# of steps to turn on
                                Pmax = Pmax+QP.b(rows(2*(e+n-1)-1));
                            end
                        end
                    end
                end
            end
            stops = nonzeros((1:nS)'.*(Locked(1:nS,j)-(Locked(2:end,j))>0));
            if ~isempty(stops)
                % example, the 6th index of locked is zero and it takes 3 steps to shut down (LB, .67LB, .33LB, 0), then stops(1) = 5, p = 3, e = 3. We only need to change LB at 4 & 5 (67% and 33%), 6 was already set to zero
                for k = 1:1:length(stops)
                    e = stops(k);
                    if e>1 %cant change lb of IC, so nothing to do if shutting down on this step
                        n = 1;
                        Pmax = QP.b(rows(2*(e-1)));
                        for f = 1:1:nt %starting with lb of 1st state in generator, then moving on
                            while Pmax<sum(Pmin(1:f))
                                QP.lb(states(nt*(e-n)+f)) = Pmax - sum(QP.lb(states(nt*(e-n)+1:(f-1)))); %lb to ensure that it is off when locked =0
                                QP.lb(states(nt*(e-n)+(f+1):nt)) = 0; %other states of this generator 
                                Pmax = Pmax+QP.b(rows(2*(e-n-1)));
                                n= n+1; %# of steps to turn off
                            end
                        end
                    end
                end
            end
        end
    end
end

%remove disabled generators from optimization matrices
rmv = false;
for i = 1:1:length(Enabled)
    if Enabled(i) ==0  && ~isempty(QP.Organize.States{i})
        rmv = true;  %there is a disabled generator with states to remove
    end
end
if rmv
    [r,~] = size(QP.A);
    [req,xL] = size(QP.Aeq);
    rkeep = linspace(1,r,r)';
    reqkeep = linspace(1,req,req)';
    xkeep = linspace(1,xL,xL)';
    for i = 1:1:length(Enabled)
        if Enabled(i) ==0  && ~isempty(QP.Organize.States{i})%remove disabled generator
            if isfield(QP.Organize,'IC') && QP.Organize.IC(i)>0
                xkeep(QP.Organize.IC(i)) = 0;
            end
            states = QP.Organize.States{i};
            xkeep(states) = 0; %removes states associated with this generator.
            eq_rows = QP.Organize.Equalities{i};
            reqkeep(eq_rows) = 0; %equality constraints
            if isfield(QP.Organize,'SpinRow')
                rkeep(QP.Organize.SpinRow{i}) = 0; % spinning reserve inequality constraints
            end
            %% fix QP.organize to be alligned with the reduced # of states
            if isfield(QP.Organize,'Ramping') %multi-time step
                s = length(QP.organize{1,i});
                QP.organize{1,i} = [];
                ineq_rows = QP.Organize.Ramping{i};
                rkeep(ineq_rows) = 0; %ramping inequality constraint
                nt = nnz((QP.Organize.States{i}(2:end)-QP.Organize.States{i}(1:end-1))==1)/nS+1;%number of states per timestep
                for j = i+1:1:length(Enabled)
                    if ~isempty(QP.organize{1,j})
                        QP.organize{1,j} = QP.organize{1,j} -s; %lower other IC by 1
                    end
                end
                for t = 1:1:nS
                    for j = 1:1:(nG+nL)
                        if j==i
                            QP.organize{t+1,j} = [];
                            s = s+nt; %cumulative # of states removed
                        elseif ~isempty(QP.organize{t+1,j})
                            QP.organize{t+1,j} = QP.organize{t+1,j} - s; %lower other gen and line state numbers by s
                        end
                    end
                    if ~isempty(QP.Organize.HeatVented)
                        QP.Organize.HeatVented(t,:) = max(0,QP.Organize.HeatVented(t,:) - s);
                    end
                end
            else %single time step case
                QP.organize{1,i} = [];
                s = length(QP.Organize.States{i});
                for j = i+1:1:(nG+nL)
                    if ~isempty(QP.organize{1,j})
                        QP.organize{1,j} = QP.organize{1,j} -s; %lower other states by s
                    end
                end
                if ~isempty(QP.Organize.HeatVented)
                    QP.Organize.HeatVented = max(0,QP.Organize.HeatVented - s);
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