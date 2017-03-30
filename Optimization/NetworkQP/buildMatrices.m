function QP = buildMatrices(fit,dt)
%builds constant matrices for multi-time-step optimization
%Fit A includes energy storage and uses the fit with zero y-intercept
%Fit B does not include energy storage and uses non-zero intercept
%Demands, initial conditions, and utility costs must updated prior to optimization
global Plant 
Plant.optimoptions.SpinReserve = 0;
Op = strcat('OpMat',fit);
nG = length(Plant.Generator);
nL = length(Plant.subNet.lineNames);


nS = length(dt);

Organize.States=cell(1,nG+nL);
Organize.Equalities = cell(1,nG);
Organize.IC = zeros(nG,1);
Organize.Ramping = cell(1,nG);
Organize.StorageInequalities = cell(1,nG);
Organize.Transmission = cell(1,nL);
Organize.Dispatchable = zeros(1,nG);

networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
%% First organize the order of states (set cost, and bounds: H, f, lb, ub)
% IC for each generator & storage
% states for each generator/storage at t = 1
% spinning reserve for generator/storage at t = 1 if option is selected
% states for each transmission line/pipe at t = 1
% repeat order of generators and lines for t = 2:nS
ic = 0; % row index of the Aeq matrix and beq vector
H = []; f = []; lb =[]; ub = [];
QP.organize = cell(nS+1,nG+nL);
QP.constCost = zeros(1,nG);
Organize.SpinReserveStates = zeros(nS,nG+1);
for i = 1:1:nG
    if isfield(Plant.Generator(i).(Op),'Ramp') 
        ic = ic+1;%initial condition state
        QP.organize{1,i} = ic; %output state organized into matrix of time vs. generator (IC) 
        H(end+1) = 0;
        f(end+1) = 0;
        lb(end+1) = 0;
        ub(end+1) = Plant.Generator(i).Size*10; %oversized uper bound in case it actually overproduces at some moment.
    else
        QP.organize{1,i} = []; %output state organized into matrix of time vs. generator (IC)   
    end
end
xL = ic;
for i = 1:1:nG
    Gen = Plant.Generator(i).(Op);
    s = length(Gen.states);%generator with multiple states
    if s>0
        if ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
            QP.organize{2,i} = xL+1; %output state for storage is only SOC
        else
            QP.organize{2,i} = xL+1:xL+s; %output is sum of multiple states at each time step
            if isempty(strfind(Plant.Generator(i).Type,'Utility'))
                Organize.Dispatchable(i) = 1;
                if isfield(Plant.Generator(i).(Op),'constCost') 
                     QP.constCost(i) = Plant.Generator(i).OpMatB.constCost;
                end
            end
        end
        for k = 1:1:s
            H(end+1) = Gen.(Gen.states{k}).H;
            f(end+1) = Gen.(Gen.states{k}).f;
            lb(end+1) = Gen.(Gen.states{k}).lb;
            ub(end+1) = Gen.(Gen.states{k}).ub;
        end
        if Plant.optimoptions.SpinReserve && isfield(Gen.OpMatA.output,'E')
            if ~isempty(strfind(Plant.Generator(i).Type,'Storage')) || Organize.Dispatchable(i) %storage and dispatchable generators have spinning reserve capacity
                s = s+1;
                Organize.SpinReserveStates(1,i) = xL +s; %state of spinning reserve at time 1
                H(end+1) = 0;
                f(end+1) = 0;
                lb(end+1) = 0;
                ub(end+1) = Gen.Size;
            end
        end
        Organize.States(i)= {xL+1:xL+s};
        xL = xL + s;
    end
end
if Plant.optimoptions.SpinReserve
    Organize.SpinReserveStates(1,end) = xL+1;
    xL = xL + 1; % Add single state for reserve shortfall state which is equal to reserve target + reserve sold as ancillary service - actual spinning reserve
    H(end+1) = 0;
    f(end+1) = 0;
    lb(end+1) = 0;
    ub(end+1) = sum(ub(nonzeros(Organize.SpinReserveStates(1,1:nG))));
end

for i = 1:1:nL % 3 states for each line (state of the line and penalty term in each direction)
    QP.organize{1,nG+i} = []; %no initial condition for lines
    QP.organize{2,nG+i} = xL+1; %line state organized into matrix of time vs. line state
    Organize.States(nG+i)= {[xL+1, xL+2, xL+3]};
    H(end+1:end+3) = [0 0 0];
    f(end+1:end+3) = [0 0 0];
    lb(end+1:end+3) = [-Plant.subNet.lineLimit(i),0,0];
    ub(end+1:end+3) = [Plant.subNet.lineLimit(i), Plant.subNet.lineLimit(i)*(1-Plant.subNet.lineEff(i)), Plant.subNet.lineLimit(i)*(1-Plant.subNet.lineEff(i))];
    xL = xL + 3;
end


%% Next organize equality equations (Aeq, and Demand to locate beq later)
% IC for each generator and storage device
% Electric energy balance @ each Electric subNet node  at t = 1
% Heat balance @ each DistricHeat subNet node at t = 1
% Cooling balance @ each DistrictCool subNet node at t = 1
% Any generator link equalities (linking states within a generator)
% Repeat power balance equalities and link equalities at t = 2:nS
req = 0; % row index of the Aeq matrix and beq vector
beq = [];
for i = 1:1:nG
    if isfield(Plant.Generator(i).(Op),'Ramp') 
        req = req+1;%initial condition state will be identity matrix at start of Aeq
        Organize.IC(i) = req;
    else Organize.IC(i) = 0;
    end
end

% The following puts together the energy balance equations
% 1 equation for each subNet node
% Nodes were agregated if their line efficiencies were 1
storIndex = [];
storRow = [];
storStates = [];
for net = 1:1:length(networkNames)
    n = length(Plant.subNet.(networkNames{net}));
    Organize.Balance.(networkNames{net}) = [];
    for i = 1:1:n
        if strcmp(networkNames{net},'DistrictHeat') && Plant.optimoptions.excessHeat == 1
            %skip and do in inequality section
        else
            req = req+1;%there is an energy balance at this node
            Organize.Balance.(networkNames{net})(end+1) = req;
            %%identify generators at this node
            genI = Plant.subNet.(networkNames{net})(i).Equipment;
            for j = 1:1:length(genI)
                states = Organize.States{genI(j)};%states associated with gerator i
                if strcmp(networkNames{net},'Electrical') && isfield(Plant.Generator(genI(j)).(Op).output,'E')
                    if strcmp(Plant.Generator(genI(j)).Type,'Electric Storage')
                        %record where you put values, so you can add t-1 state and divide by dt(t)
                        storIndex(end+1) = genI(j);
                        storRow(end+1) = req;
                        storStates(end+1) = states(1);
                    else
                        Aeq(req,states) = Plant.Generator(genI(j)).(Op).output.E;
                    end
                elseif strcmp(networkNames{net},'DistrictHeat') && isfield(Plant.Generator(genI(j)).(Op).output,'H')
                    if strcmp(Plant.Generator(genI(j)).Type,'Thermal Storage')
                        %record where you put values, so you can add t-1 state and divide by dt(t)
                        storIndex(end+1) = genI(j);
                        storRow(end+1) = req;
                        storStates(end+1) = states(1);
                    else
                        Aeq(req,states) = Plant.Generator(genI(j)).(Op).output.H;
                    end
                elseif strcmp(networkNames{net},'DistrictCool') && isfield(Plant.Generator(genI(j)).(Op).output,'C')
                    if strcmp(Plant.Generator(genI(j)).Type,'Thermal Storage')
                        %record where you put values, so you can add t-1 state and divide by dt(t)
                        storIndex(end+1) = genI(j);
                        storRow(end+1) = req;
                        storStates(end+1) = states(1);
                    else
                        Aeq(req,states) = Plant.Generator(genI(j)).(Op).output.C;
                    end
                end
            end
            %%identify lines coming in and out
            connect = Plant.subNet.(networkNames{net})(i).connections;
            for j = 1:1:length(connect)
                I = [];
                while isempty(I)
                    nName = Plant.subNet.(networkNames{net})(i).nodes(1); %name of current subnet node
                    I = find(strcmp(strcat(nName,'_',networkNames{net},'_',connect{j}),Plant.subNet.lineNames),1,'first');
                    if isempty(I)
                        I = find(strcmp(strcat(connect{j},'_',networkNames{net},'_',nName),Plant.subNet.lineNames),1,'first');
                        if~isempty(I)
                            dir = -1; %reverse direction
                        else dis('error: line does not exist')
                        end
                    else
                        dir = 1; %forward direction
                    end
                end
                linestates = Organize.States{nG+I};
                if dir ==1
                    Aeq(req,linestates) = [-1,0,-1]; %forward direction A-->B, is positive, thus positive transmission is power leaving the node, the penalty from b->a is power not seen at a
                else
                    Aeq(req,linestates) = [1,-1,0];%reverse direction B-->A, is positive, thus positive power is power entering the node, the penalty from a->b is power not seen at b
                end
            end
            %%note any demands at this node
            load = Plant.subNet.(networkNames{net})(i).Load;
            for j = 1:1:length(load)
                Demand.(networkNames{net})(load(j)) = req;
            end
        end
    end
end

            
%link is a field if there is more than one state and the states are linked by an inequality or an equality
for i = 1:1:nG
    if isfield(Plant.Generator(i).(Op),'link') && isfield(Plant.Generator(i).(Op).link,'eq')
        [m,n] = size(Plant.Generator(i).(Op).link.eq);
        states = Organize.States{i};%states associated with gerator i
        for k = 1:1:m
            Aeq(req+k,states) = Plant.Generator(i).(Op).link.eq(k,:);
            beq(req+k) = Plant.Generator(i).(Op).link.beq(k);
        end
        if m==1
            Organize.Equalities(i) = {req+1}; 
        else
            Organize.Equalities(i) = {req+1:req+m}; 
        end
        req = req+m;
    end
end

%% Organize inequalities
% Distric heating energy inequalities (inequality because heat can be rejected)
% 2 ramping constraints for each generator at each time step
% constraints for each Energy storage 
% 2 constraints for each transmission line
% repeat all of the above for each time 2:nS

r = 0; % row index of the A matrix & b vector
A = [];
%Distric heating energy inequalities (inequality because heat can be rejected)
storIndexT = [];
storRowT = [];
storStatesT = [];
for net = 1:1:length(networkNames)
    n = length(Plant.subNet.(networkNames{net}));
    Organize.Imbalance.(networkNames{net}) = [];
    for i = 1:1:n
        if strcmp(networkNames{net},'DistrictHeat') && Plant.optimoptions.excessHeat == 1
            r = r+1;%there is an energy balance at this node
            Organize.Imbalance.(networkNames{net})(end+1) = r;
            %%identify generators at this node
            genI = Plant.subNet.(networkNames{net})(i).Equipment;
            for j = 1:1:length(genI)
                states = Organize.States{genI(j)}-ic;%states associated with gerator i
                if strcmp(networkNames{net},'DistrictHeat') && isfield(Plant.Generator(genI(j)).(Op).output,'H')
                    if strcmp(Plant.Generator(genI(j)).Type,'Thermal Storage')
                        %record where you put values, so you can add t-1 state and divide by dt(t)
                        storIndexT(end+1) = genI(j);
                        storRowT(end+1) = r;
                        storStatesT(end+1) = states(1);
                    else
                        A(r,states) = Plant.Generator(genI(j)).(Op).output.H;
                    end
                end
            end
            %%identify lines coming in and out
            connect = Plant.subNet.(networkNames{net})(i).connections;
            for j = 1:1:length(connect)
                I = [];
                while isempty(I)
                    nName = Plant.subNet.(networkNames{net})(i).nodes(1); %name of current subnet node
                    I = find(strcmp(strcat(nName,'_',networkNames{net},'_',connect{j}),Plant.subNet.lineNames),1,'first');
                    if isempty(I)
                        I = find(strcmp(strcat(connect{j},'_',networkNames{net},'_',nName),Plant.subNet.lineNames),1,'first');
                        if~isempty(I)
                            dir = -1; %reverse direction
                        else dis('error: line does not exist')
                        end
                    else
                        dir = 1; %forward direction
                    end
                end
                linestates = Organize.States{nG+I};
                if dir ==1
                    A(r,linestates) = [1,-1,0]; %forward direction A-->B, is positive
                else
                    A(r,linestates) = [-1,0,-1];%reverse direction B-->A, is positive
                end
            end
            %%note any demands at this node
            load = Plant.subNet.(networkNames{net})(i).Load;
            for j = 1:1:length(load)
                Demand.(networkNames{net})(load(j)) = r;
            end
        end
    end
end

%Ramping Inequalities
rampRow = [];
rampStates = {};
for i = 1:1:nG
    if isfield(Plant.Generator(i).(Op),'Ramp') 
        states = Organize.States{i}-ic;%states associated with generator i at t = 1
        %%if storage, ramping only affects 1st state)
        if ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
            states = states(1);
        end
        rampRow(end+1) = r+1;
        rampStates(end+1) = {states};
        b(r+1:r+2,1) = Plant.Generator(i).(Op).Ramp.b;
        Organize.Ramping(i) = {[r+1,r+2]}; 
        r = r+2;
    end
end
%Energy Storage inequalities
storIndex2 = [];
storRow2 = [];
storStates2 = [];
for i = 1:1:nG 
    if isfield(Plant.Generator(i).(Op),'link') && isfield(Plant.Generator(i).(Op).link,'ineq')%link is a field if there is more than one state and the states are linked by an inequality or an equality
        states = Organize.States{i}-ic;%states associated with gerator i
        [m,n] = size(Plant.Generator(i).(Op).link.ineq);
        b(r+1:r+m) = Plant.Generator(i).(Op).link.bineq;
        storIndex2(end+1) = i;
        storRow2(end+1) = r+1; %row of first constraint
        storStates2(end+1) = states(1); %first storage state
        Organize.StorageInequalities(i) = {r+1:r+m}; 
        r = r+m;
    end
end
%Transmission line inequalites (penalty terms)
for i = 1:1:nL
    linestates = Organize.States{nG+i}-ic;
    A(r+1,linestates) = [(1-Plant.subNet.lineEff(i)), -1, 0];% Pab*(1-efficiency) < penalty a to b
    A(r+2,linestates) = [-(1-Plant.subNet.lineEff(i)), 0, -1];% -Pab*(1-efficiency) < penalty b to a
    Organize.Transmission(i) = {r+1:r+2}; 
    r = r+2;
end

%spinning reserve inequalities (sum all spinning reserves & individual spinning reserves) 
%currently only implemented for electric power
SpinRow = cell(nG,1);
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Electrical') && Plant.optimoptions.SpinReserve
        Organize.SpinReserve = r;
        A(r+1,nonzeros(Organize.SpinReserveStates(1,:))) = -1; %Inequality for spinning reserve shortfall
        r = r+1;
        for i = 1:1:nG
            if isfield(Plant.Generator(i).OpMatA.output,'E') && (~isempty(strfind(Plant.Generator(i).Type,'Storage')) || Organize.Dispatchable(i))%electric storage & dispatchable electric generators have spinning reserve capacity
                SpinRow(i) = {r+1:r+2};
                r = r+2;
            end
        end
    end
end
        
%% Finally, expand by the number of time steps
% # of states per time step is xL - ic
% Order of states will repeat for nS time steps
t1States = xL - ic;
totalStates = nS*t1States+ic;
for t= 2:nS
    for n = 1:1:length(QP.organize(2,:))
        QP.organize{t+1,n} = QP.organize{t,n} + t1States;
    end
end
QP.H = zeros(totalStates,1);
QP.f = zeros(totalStates,1);
QP.lb = zeros(totalStates,1);
QP.ub = zeros(totalStates,1);

QP.H(1:ic+t1States)  = H';
QP.f(1:ic+t1States)  = f';
QP.lb(1:ic+t1States)  = lb';
QP.ub(1:ic+t1States)  = ub';
H = H(ic+1:end);
f = f(ic+1:end);
lb = lb(ic+1:end);
ub = ub(ic+1:end);
for t= 2:nS
    QP.H(ic+(t-1)*t1States+1:ic+t*t1States)  = H';
    QP.f(ic+(t-1)*t1States+1:ic+t*t1States)  = f';
    QP.lb(ic+(t-1)*t1States+1:ic+t*t1States)  = lb';
    QP.ub(ic+(t-1)*t1States+1:ic+t*t1States)  = ub';
end
QP.H = diag(QP.H);

% number of generator equality constraints & energy balances at each time step will be req - ic
t1Balances = req - ic;
totalEqualities = nS*t1Balances + ic;
QP.Aeq = zeros(totalEqualities,totalStates);
QP.beq = zeros(totalEqualities,1);
[m,n] = size(Aeq);
if req>m || xL>n
    Aeq(req,xL) = 0;%make sure Aeq is the right size
end
if req>length(beq)
    beq(req,1) = 0; %make sure beq is the right length
end
Aeq = Aeq(ic+1:end,ic+1:end);%remove initial condition equalities and states
beq = beq(ic+1:end,1);%remove initial condition equalities
QP.Aeq(1:ic,1:ic) = eye(ic); %initial condition identity matrix
for t= 1:nS
    QP.Aeq(ic+(t-1)*t1Balances+1:ic+t*t1Balances,ic+(t-1)*t1States+1:ic+t*t1States) = Aeq;
    QP.beq(ic+(t-1)*t1Balances+1:ic+t*t1Balances)  = beq;
end
%storage
if ~isempty(storRow)
    for t= 1:nS
        for i = 1:1:length(storRow)
            QP.Aeq((t-1)*t1Balances+storRow(i),(t-1)*t1States+storStates(i)) = -1/dt(t); %SOC at t
            if ismember('Y',Plant.Generator(storIndex(i)).(Op).states)
                QP.Aeq((t-1)*t1Balances+storRow(i),(t-1)*t1States+storStates(i)+1) = -1/dt(t); %charging penalty at t
            end
            if t==1
                storIC = nnz(Organize.IC(1:storIndex(i))); %order in IC
                QP.Aeq(storRow(i),storIC) = 1/dt(t);%SOC at IC
            else
                QP.Aeq((t-1)*t1Balances+storRow(i),(t-2)*t1States+storStates(i)) = 1/dt(t); %SOC at t-1
            end
        end
    end
end

% number of generator inequality constraints & energy imbalances at each time step will be r 
% there are 2 ramping constraints on each generator/storage
t1ineq = r;
totalInequalities = nS*t1ineq;
QP.A = zeros(totalInequalities,totalStates);
QP.b = zeros(totalInequalities,1);

if ~isempty(A)
    [m,n] = size(A);
    if m<t1ineq || n<t1States
        A(t1ineq,t1States) = 0;
    end
else A = zeros(t1ineq,t1States);
end
if length(b)<t1ineq
    b(t1ineq) = 0;
end
for t= 1:nS
    A2 = A; %A2 will be A(t) to capture changing step sizes, but will keep parts of A that don't change with step size
    A_tPrev = zeros(t1ineq,t1States); %matrix of the same size, will collect terms associated with states at previous time step, but same inequality rows:  QP.A = [A_tPrev, A2];
    A_IC = zeros(t1ineq,ic);
    %Storage in energy balance (district heating with heat rejection)
    if ~isempty(storRowT)
        for i = 1:1:length(storRowT)
            A2(storRowT(i),storStatesT(i)) = -1/dt(t); %SOC at t
            if ismember('Y',Plant.Generator(storIndexT(i)).(Op).states)
                A2(storRowT(i),storStatesT(i)+1) = -1/dt(t); %charging penalty at t
            end
            if t==1
                storIC = nnz(Organize.IC(1:storIndexT(i))); %order in IC
                A_IC(storRowT(i),storIC) = 1/dt(t);%SOC at IC
            else
                A_tPrev(storRowT(i),storStatesT(i)) = 1/dt(t); %SOC at t-1
            end
        end
    end
    %Ramping
    if ~isempty(rampRow)
        for i = 1:1:length(rampRow)
            Rstates = rampStates{i};
            A2(rampRow(i),Rstates) = 1/dt(t);%ramp up 
            A2(rampRow(i)+1,Rstates) = -1/dt(t);%ramp down
            if t ==1 %ramping from initial condition
                A_IC(rampRow(i),i) = -1/dt(t); %ramp up 
                A_IC(rampRow(i)+1,i) = 1/dt(t); %ramp down 
            else %condition at previous time step
                A_tPrev(rampRow(i),Rstates) = -1/dt(t); %ramp up 
                A_tPrev(rampRow(i)+1,Rstates) = 1/dt(t); %ramp down 
            end
        end
    end
    %Storage inequalities constraints
    if ~isempty(storRow2)
        for i = 1:1:length(storRow2)
            link = Plant.Generator(storIndex2(i)).(Op).link.ineq;
            [m2,n2] = size(link);
            for k = 1:1:m2
                if k ==1 && ismember('Y',Plant.Generator(storIndex2(i)).(Op).states) % (SOC(t) - SOC(t-1))*(1-efficiency) < charging penalty
                    A2(storRow2(i),storStates2(i)) = link(1,1)/dt(t);  % SOC at t  
                    A2(storRow2(i),storStates2(i)+1) = link(1,2)/dt(t);  % charging state at t: value is -1/(1-efficiency)
                    if t ==1 %SOC change from IC
                        storIC = nnz(Organize.IC(1:storIndex2(i))); %order in IC
                        A_IC(storRow2(i),storIC) = -1/dt(t); % SOC at t-1
                    else
                        A_tPrev(storRow2(i),storStates2(i)) = -1/dt(t);  % SOC at t-1
                    end
                else
                    states = storStates2(i)+(1:n2) - 1;
                    A2(storRow2(i)+k-1,states) = link(k,:);
                end
            end
        end
    end
    
    %Transmission
    %%no connection to previous or later time steps, and no dependence on step size. Everything already done in A
    
    %spinning reserve inequalities (sum all spinning reserves & individual spinning reserves) 
    %currently only implemented for electric power
    if Plant.optimoptions.SpinReserve
        for i = 1:1:nG
            row = SpinRow{i};
            if Organize.Dispatchable(i) && isfield(Plant.Generator(i).OpMatA.output,'E')% dispatchable generators have spinning reserve capacity
                A2(row(1),Organize.SpinReserveStates(1,i)) = 1/dt(t); %SR + power(t) - power(t-1)<= ramprate*dt
                A2(row(1),Organize.States{i}) = 1/dt(t);
                if t ==1 %ramping from IC
                    genIC = nnz(Organize.IC(1:i)); %order in IC
                    A_IC(row(1),genIC) = -1/dt(t); % Power at t-1
                else
                    A_tPrev(row(1),Organize.States{i}) = -1/dt(t); % Power at t-1
                end
                b(row(1)) = Plant.Generator(i).(Op).Ramp.b(1);%ramp up constraint
                
                A2(row(2),Organize.SpinReserveStates(1,i)) = 1; %SR + power <= Size
                A2(row(2),Organize.States{i}) = 1;
                b(row(2)) = Plant.Generator(i).Size; %max capacity constraint
            elseif ~isempty(strfind(Plant.Generator(i).Type,'Storage')) && isfield(Plant.Generator(i).OpMatA.output,'E') %electric storage 
                A2(row(1),Organize.SpinReserveStates(1,i)) = 1; %SR + (SOC(t-1) - SOC(t))/dt <= peak discharge
                A2(row(1),QP.organize{2,i}) = -1/dt(t);
                A2(row(2),Organize.SpinReserveStates(1,i)) = 1; %SR - SOC(t-1)/dt <= 0
                if t ==1 %SOC change from IC
                    storIC = nnz(Organize.IC(1:i)); %order in IC
                    A_IC(row(1),storIC) = 1/dt(t); % SOC at t-1
                    A_IC(row(2),storIC) = -1/dt(t); % SOC at t-1
                else
                    A_tPrev(row(1),QP.organize{2,i}) = 1/dt(t); % SOC at t-1
                    A_tPrev(row(2),QP.organize{2,i}) = -1/dt(t); % SOC at t-1
                end
                b(row(1)) = Plant.Generator(i).(Op).Ramp.b(2);%peak discharge constraint
            end
        end
    end
            
    %put into QP.A
    QP.A((t-1)*t1ineq+1:t*t1ineq,ic+(t-1)*t1States+1:ic+t*t1States) = A2;
    if t==1
        QP.A(1:t1ineq,1:ic) = A_IC;
    else
        QP.A((t-1)*t1ineq+1:t*t1ineq,ic+(t-2)*t1States+1:ic+(t-1)*t1States) = A_tPrev;
    end
    QP.b((t-1)*t1ineq+1:t*t1ineq,1) = b;
end

%list of states associated with each generator at all time steps
for i = 1:1:length(QP.organize(1,:))
    if ~isempty(Organize.States{i})
        s = Organize.States{i};
        n = length(s);
        states = [];
        for t = 1:1:nS
            states(end+1:end+n) = s + (t-1)*t1States;
        end
        Organize.States(i) = {states};
    end
end

%setup indexing for beq (and b) affiliated with each demand
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'DistrictHeat') && Plant.optimoptions.excessHeat == 1
        for i = 1:1:length(Demand.(networkNames{net}))
            ineq1 = Demand.(networkNames{net})(i);%this is the inequality constraint associated with this demand at t = 1
            Organize.Demand.(networkNames{net})(i) = {ineq1:t1ineq:totalInequalities};%evaluating this string produces a vector of all the equality rows associated with this demand at t = 1 to nS
        end
    else
        for i = 1:1:length(Demand.(networkNames{net}))
            eq1 = Demand.(networkNames{net})(i); %this is the equality constraint associated with this demand at t = 1
            Organize.Demand.(networkNames{net})(i) = {eq1:t1Balances:totalEqualities}; %evaluating this string produces a vector of all the equality rows associated with this demand at t = 1 to nS
        end
    end
    bal = Organize.Balance.(networkNames{net});
    n = length(bal);
    Organize.Balance.(networkNames{net}) = zeros(n,nS);
    for t = 1:1:nS
        Organize.Balance.(networkNames{net})(1:n,t) = bal + (t-1)*t1Balances;
    end
    imbal = Organize.Imbalance.(networkNames{net});
    n = length(imbal);
    Organize.Imbalance.(networkNames{net}) = zeros(n,nS);
    for t = 1:1:nS
        Organize.Imbalance.(networkNames{net})(1:n,t) = imbal + (t-1)*t1ineq;
    end
end

% Equalities for each generator
for i = 1:1:nG
    if isfield(Plant.Generator(i).(Op),'link') && isfield(Plant.Generator(i).(Op).link,'eq')
        s = Organize.Equalities{i};
        n = length(s);
        rows = [];
        for t = 1:1:nS
            rows(end+1:end+n) = s + (t-1)*t1Balances;
        end
    else
        rows = [];
    end
    Organize.Equalities(i) = {rows};
end


% Ramping each generator/storage
for i = 1:1:nG
    if isfield(Plant.Generator(i).(Op),'Ramp') 
        s = Organize.Ramping{i};
        ramp = [];
        for t = 1:1:nS
            ramp(end+1:end+2) = s + (t-1)*t1ineq;
        end
        Organize.Ramping(i) = {ramp};
    end
end

% Storage inequalities for each generator
for i = 1:1:nG
    if (isfield(Plant.Generator(i).(Op),'link') && isfield(Plant.Generator(i).(Op).link,'ineq'))
        s = Organize.StorageInequalities{i};
        n = length(s);
        stor = [];
        for t = 1:1:nS
            stor(end+1:end+n) = s + (t-1)*t1ineq;
        end
        Organize.StorageInequalities(i) = {stor};
    end
end

%transmission inequalities
for i = 1:1:nL
    s = Organize.Transmission{i};
    line = [];
    for t = 1:1:nS
        line(end+1:end+2) = s + (t-1)*t1ineq;
    end
    Organize.Transmission(i) = {line};
end
QP.Organize = Organize; %indices (rows and columns) associated with each generator, allowing generators to be removed later