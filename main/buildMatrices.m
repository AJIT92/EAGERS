function QPall = buildMatrices(fit,Time)
%builds constant matrices for multi-time-step optimization
%Fit A includes energy storage and uses the fit with zero y-intercept
%Fit B does not include energy storage and uses non-zero intercept
%Demands, initial conditions, and utility costs must updated prior to optimization
global Plant 
Op = strcat('OpMat',fit);
nG = length(Plant.Generator);

dt = Time - [0, Time(1:end-1)];
nS = length(Time);

Outs = Plant.optimoptions.Outputs;
Organize.States={};
Organize.Equalities = {};
Organize.Inequalities = {};
Organize.IC = zeros(nG,1);
seq = 1;
while seq<=length(Outs);
    xL = 0;
    req = 0; % row index of the Aeq matrix and beq vector
    r = 0; % row index of the A matrix & b vector
    if strcmp('E',Outs{seq})>0
        Organize.E.Stor = zeros(1,nG);
        Organize.E.Demand = {'beq',(1:nS)};
        req = req+nS;
        include = {'CHP Generator', 'Electric Generator'};
        if nnz(strcmp('H',Outs))>0
            include(end+1) = {'Heater'};
            Outs = Outs(~(strcmp('H',Outs))); %combine heating optimization with electric due to CHP generators
            if Plant.optimoptions.excessHeat == 1
                Organize.H.Demand = {'b',(r+1:r+nS)};
                r = r+nS;
            else Organize.H.Demand = {'beq',(req+1:req+nS)};
                req = req+nS;
            end 
        end
        if nnz(strcmp('C',Outs))>0 && Plant.optimoptions.sequential == 0 %%load the chillers separately if the user has selected the sequential option
            include(end+1) = {'Chiller'};
            Outs = Outs(~(strcmp('C',Outs))); 
            Organize.C.Demand = {'beq',(req+1:req+nS)};
            req = req+nS;
        end
    elseif strcmp('C',Outs{seq})>0
        Organize.C.Stor = zeros(1,nG);
        Organize.C.Demand = {'beq',(1:nS)};
        include = {'Chiller'};
        req = req+nS;
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
    
    QP.organize = cell(nS+1,nG);
    for i = 1:1:nG
        Gen = Plant.Generator(i).(Op);
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
        if thisSeq(i) || chill(i) || heater(i) || stor(i) || utility(i) || storC(i) || storH(i) || utilC(i) || utilH(i)
            if isfield(Gen,'Ramp') 
                Organize.IC(i) = req+1;
                QP.organize{1,i} = xL+1; %output state organized into matrix of time vs. generator (IC)
                ic = 1;
            else ic =0; %initial condition
            end
            s = length(Gen.states);%generator with multiple states
            xLend = xL + s*nS+ic;
            Organize.States(i) = {xL+1:xLend};
            if stor(i) || storC(i) || storH(i)
                s = 1; % storage only outputs the 1st state
            end
            for t = 1:1:nS
                QP.organize{t+1,i} = xL+ic+t:nS:xL+(s-1)*nS+ic+t; %output state organized into matrix of time vs. generator 
            end
            neq = eval(Gen.req);
            Organize.Equalities(i) = {req+1:req+neq};
            req = req + neq;
            ineq = eval(Gen.r);
            Organize.Inequalities(i) = {r+1:r+ineq};
            r = r + ineq;
            xL = xLend;
        else %if this is not a member of this sequence, but in case it is the last generator
            Organize.States(i) = {[]};
            Organize.Equalities(i) = {[]};
            Organize.Inequalities(i) = {[]};
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
    allStor = [stor, storC, storH];
    allGen = [thisSeq, chill, heater];
    allUtility = [utility, utilC, utilH];
        
    QP.H = zeros(xL,1);
    QP.f = zeros(xL,1);
    QP.Aeq = zeros(req,xL);
    QP.beq = zeros(req,1);
    QP.A = zeros(r,xL);
    QP.b = zeros(r,1);
    QP.lb = zeros(xL,1);
    QP.ub = inf*ones(xL,1);%if a state has no upper bound, then the upper bound is infinite

    %% go through list of generators & build matrices
    for i = 1:1:nG
        k = Organize.States{i};% first index for the states of the x vector in C = x'Hx+f'x
        req = Organize.Equalities{i}; %first row of equality matrix (not associated with output)
        r = Organize.Inequalities{i}; %first row of inequality matrix (not associated with output)
        Gen = Plant.Generator(i).(Op);
        states = Gen.states;
        if~isempty(k)
            k = k(1);
        end
        if~isempty(r)
            r = r(1);
        end
        if ismember(i,[allGen, allStor, allUtility])
            % Initial condition (only 1 per generator output, not per state: must correspond to 1st state, utility does not have an IC
            if isfield(Gen,'Ramp') %has ramping constraint, thus it has an initial condition in order to enforce this
                QP.Aeq(req(1),k) = 1; %initial value in beq will be updated prior to running optimization
                k = k+1; % move to state at t=1
                if length(req)>1
                    req = req(2:end); %move to next equality constraint
                end
            end

            % Output {'E', 'C', 'H', 'S', 'H2'...}
            %this creates the first rows of the A and Aeq matrixes which
            %specify the outputs equality for each demand type at each interval
            gOuts = fieldnames(Gen.output);
            for p = 1:1:length(gOuts) %load all outputs into the correct equality or inequality equation
                mat = Organize.(gOuts{p}).Demand{1};
                index = Organize.(gOuts{p}).Demand{2}';
                if strcmp(mat,'beq')
                    mat = 'Aeq';
                else mat = 'A';
                end
                %% add heat smoothing here
                %%
                if ismember(i,allStor) %energy storage -- 2nd term (difference is discharging)
                    for t = 1:1:length(index)
                        QP.(mat)(index(t),k+t-1) = -1/dt(t); %final state at end of time step
                        QP.(mat)(index(t),k+t-2) = 1/dt(t); %Initial state at start of time step (always =1) final state can be affected by self-discharging (value in Gen.output = -1/(1-selfDisch))
                        if ismember('Y',states) 
                            QP.(mat)(index(t),k+nS+t-1) = -1/dt(t); %charging state (negative load)
                        end
                    end
                elseif ismember(i,allUtility)
                    QP.(mat)(index,k:k+nS-1) = Gen.output.(gOuts{p})*eye(nS);
                    if length(states)==2 %sellback state is 2nd
                        QP.(mat)(index,k+nS:k+2*nS-1) = Gen.output.(gOuts{p})*(-eye(nS));
                    end
                else
                    for s = 1:1:length(states)
                        QP.(mat)(index,k+(s-1)*nS:k+s*nS-1) = Gen.output.(gOuts{p})*eye(nS);
                    end
                end
            end

            % link states
            if isfield(Gen,'link')  %link is a field if there is more than one state and the states are linked by an inequality or an equality
                if isfield(Gen.link,'eq')
                    for s = 1:1:length(states)
                        QP.Aeq(req, k+nS*(s-1):k+nS*s-1) = Gen.link.eq(s)*eye(nS);
                    end
                    QP.beq(req) = Gen.link.beq;
                end
                if isfield(Gen.link,'ineq') 
                    [m,p] = size(Gen.link.ineq);
                    for j = 1:1:m
                        if j==1 && m ==3 %first inequality constraint sets up charging state (power): requires -SOC(t-1)/dt + SOC(t)/dt - Charging < 0
                            for t=1:1:nS% can't overwrite whole sub-matrix like before, so must do one at a time
                                QP.A(r+t-1, k+t-2) = -1/dt(t);%SOC(t-1)/dt, current state is given scalar 1
                                QP.A(r+t-1, k+t-1) = Gen.link.ineq(1,1)/dt(t);%SOC(t)/dt
                                QP.A(r+t-1, k+t+nS-1) = Gen.link.ineq(1,2)/dt(t);%charging state (power, don't divide by dt, because already in kW)
                            end
                        else
                            for s = 1:1:p
                                QP.A(r:r+nS-1, k+nS*(s-1):k+nS*s-1) = Gen.link.ineq(j,s)*eye(nS);
                            end
                        end
                        QP.b(r:r+nS-1) = Gen.link.bineq(j);
                        r = r+nS;
                    end
                end
            end

            % ramp rates
            if isfield(Gen, 'Ramp')
                for t = 1:1:nS
                    QP.A(r:r+1, k+t-2:k+t-1) = Gen.Ramp.A;% 2x2 matrix  {-SOC(t-1), SOC(t); SOC(t-1), -SOC(t);} 
                    MaxRamp = Gen.(states{1}).ub;
                    if ismember(i,allGen) && length(states)>1 %generator with multiple states
                        for s = 2:1:length(states)
                            MaxRamp = MaxRamp + Gen.(states{s}).ub;
                            if t==1 %IC (only do 2nd column, secondary states do not have IC)
                                QP.A(r:r+1, k+nS*(s-1)+t-1) = Gen.Ramp.A(:,2);% 1x2 matrix {X(t);-X(t);} 
                            else
                                QP.A(r:r+1, k+nS*(s-1)+t-2:k+nS*(s-1)+t-1) = Gen.Ramp.A;% 2x2 matrix {-X(t-1), X(t); X(t-1), -X(t);} 
                            end
                        end
                    end
                    QP.b(r:r+1) = min(MaxRamp,Gen.Ramp.b*dt(t)); %scale the ramping constraint for variable time steps
                    r = r+2;
                end
            end

            %% costs and bounds, each state has 1 value at each time interval
            %bounds don't matter at IC (constrained by equality)
            for j = 1:1:length(states)
                QP.H(k:k+nS-1) = Gen.(states{j}).H.*dt; %costs scaled by length of time interval
                QP.f(k:k+nS-1) = Gen.(states{j}).f.*dt; %costs scaled by length of time interval
                QP.lb(k:k+nS-1) = Gen.(states{j}).lb;
                QP.ub(k:k+nS-1) = Gen.(states{j}).ub;
                k = k+nS; %add nS states 
            end
        end
    end
    QP.H = diag(QP.H);
    QPall.QP.(Outs{seq})= QP;% quadratic programing matrices (H, f, Aeq, beq, A, b, ub, lb) -- parts may need to be updated later 
    seq = seq+1;
end
QPall.Organize = Organize; %indices (rows and columns) associated with each generator, allowing generators to be removed later
QPall.Timestamp = Time(1);