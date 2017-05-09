function QP = buildMatrices(Op,dt)
%builds constant matrices for multi-time-step optimization
%Demands, initial conditions, and utility costs must updated prior to optimization
global Plant 
nG = length(Plant.Generator);
nS = length(dt);

networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
for net = 1:1:length(networkNames)
    nLinet(net) = length(Plant.subNet.lineNames.(networkNames{net}));
end
nL = sum(nLinet);

Organize.States=cell(1,nG+nL);
Organize.Equalities = cell(1,nG);
Organize.IC = zeros(nG+nL,1);
Organize.Ramping = cell(1,nG);
Organize.Inequalities = cell(1,nG);
Organize.Transmission = cell(1,nL);
Organize.Dispatchable = zeros(1,nG);
%% First organize the order of states (set cost, and bounds: H, f, lb, ub)
% IC for each generator & storage
% states for each generator/storage at t = 1
% spinning reserve for generator/storage at t = 1 if option is selected
% states for each transmission line/pipe at t = 1
% state for P_loss in heating energy balance
% state for spinning reserve shortfall (target - cumulative of all active generators) and SR provided as ancillary service
% repeat order of generators and lines for t = 2:nS
ic = 0; % row index of the Aeq matrix and beq vector
H = []; f = []; lb =[]; ub = [];
QP.organize = cell(nS+1,nG+nL);
QP.constCost = zeros(1,nG);
Organize.SpinReserveStates = zeros(nS,nG+2);
for i = 1:1:nG
    if isfield(Plant.Generator(i).(Op),'Ramp') 
        ic = ic+1;%initial condition state
        QP.organize{1,i} = ic; %output state organized into matrix of time vs. generator (IC)   
        Organize.IC(i) = ic;
    end
end
nLcum = 0; %cumulative line #
for net = 1:1:length(networkNames)
    if ~strcmp(networkNames{net},'Hydro')
        nLcum = nLcum+nLinet(net); %all hydro lines have initial state
    else
        for i = 1:1:nLinet(net) 
            ic = ic+1;%initial condition state
            QP.organize{1,nG + nLcum + i} = ic; %output state organized into matrix of time vs. generator (IC)   
            Organize.IC(nG + nLcum + i) = ic;
        end
    end
end
    
xL = 0;
GenNames = cell(nG,1);
for i = 1:1:nG
    Gen = Plant.Generator(i).(Op);
    GenNames{i} = Plant.Generator(i).Name;
    s = length(Gen.states);%generator with multiple states
    if s>0
        if ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
            QP.organize{2,i} = xL+1+ic; %output state for storage is only SOC
        else
            QP.organize{2,i} = xL+1+ic:xL+s+ic; %output is sum of multiple states at each time step
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
        Organize.States(i)= {xL+1:xL+s};
        if Plant.optimoptions.SpinReserve && isfield(Gen.output,'E') && ~isfield(Gen.output,'C') %electric systems that are not electric chillers
            if ~isempty(strfind(Plant.Generator(i).Type,'Storage')) || Organize.Dispatchable(i) %storage and dispatchable generators have spinning reserve capacity
                Organize.SpinReserveStates(1,i) = xL +s+1; %state of spinning reserve at time 1
                H(end+1) = 0;
                f(end+1) = 0;
                if isempty(strfind(Plant.Generator(i).Type,'Storage'))
                    lb(end+1) = 0;
                    ub(end+1) = Plant.Generator(i).Size;
                else %ub for storage discharge
                    if strcmp(Plant.Generator(i).Type,'Hydro Storage')
                        lb(end+1) =  0;
                        ub(end+1) = Plant.Generator(i).VariableStruct.MaxGenCapacity;
                    else
                        lb(end+1) = -Plant.Generator(i).OpMatA.Stor.PeakCharge;
                        ub(end+1) = Plant.Generator(i).OpMatA.Stor.PeakDisch;
                    end
                end
            end
            xL = xL + s + 1;
        else
            xL = xL + s;
        end
    end
end

%states for transmission lines
nLcum = 0; %cumulative line #
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Hydro')
        eff = [];
        limit = Plant.subNet.lineLimit.(networkNames{net});
        minimum = Plant.subNet.lineMinimum.(networkNames{net}); 
    else
        eff = Plant.subNet.lineEff.(networkNames{net});
        limit = Plant.subNet.lineLimit.(networkNames{net}); 
        minimum = zeros(nLinet(net),1);
    end
    for i = 1:1:nLinet(net) 
        nLcum = nLcum+1;
        QP.organize{1,nG+nLcum} = []; %no initial condition for lines
        QP.organize{2,nG+nLcum} = xL+1+ic; %line state organized into matrix of time vs. line state
        if isempty(eff) || length(eff(i,:))==1 || eff(i,2)==0 %uni-directional transfer, 1 state for each line
            Organize.States(nG+nLcum)= {xL+1};
            H(end+1) = 0;
            f(end+1) = 0;
            lb(end+1) = minimum(i);
            ub(end+1) = limit(i,1);
            xL = xL + 1;
        else% bi-directional transfer, 3 states for each line (state of the line and penalty term in each direction)
            Organize.States(nG+nLcum)= {[xL+1, xL+2, xL+3]};
            H(end+1:end+3) = [0 0 0];
            f(end+1:end+3) = [0 0 0];
            lb(end+1:end+3) = [-limit(i,2),0,0];
            ub(end+1:end+3) = [limit(i,1), limit(i,1)*(1-eff(i,1)), limit(i,2)*(1-eff(i,2))];
            xL = xL + 3;
        end
    end
end

if any(strcmp('DistrictHeat',networkNames)) && Plant.optimoptions.excessHeat == 1
    %%assume heat can be lost any any node in the network that has a device producing heat
    n = length(Plant.subNet.('DistrictHeat'));
    Organize.HeatVented =zeros(nS,n); %matrix for the state associated with venting heat at each district heating node, at each time step
    for i = 1:1:n
        genI = Plant.subNet.('DistrictHeat')(i).Equipment;%%identify generators at this node
        hasHeater = 0;
        for j = 1:1:length(genI)
            if isfield(Plant.Generator(genI(j)).(Op).output,'H')
                hasHeater = 1;
                break
            end
        end
        if hasHeater
            Organize.HeatVented(1,i) = (xL+1);
            xL = xL + 1; % Add single state for heat that is ventd to make energy equality true
            H(end+1) = 0;
            f(end+1) = 0;
            lb(end+1) = 0;
            ub(end+1) = inf;
        end
    end
end

if Plant.optimoptions.SpinReserve && any(Organize.SpinReserveStates(1,1:nG)>0)
    Organize.SpinReserveStates(1,nG+1) = xL+1; %spinning reserve shortfall at t=1
    Organize.SpinReserveStates(1,nG+2) = xL+2; %SR provided as ancillary service at t=1
    xL = xL + 2; % Add states for reserve shortfall and ancillary service: shortfall state is equal to reserve target + reserve sold as ancillary service - actual spinning reserve
    H(end+1:end+2) = 0;
    f(end+1:end+2) = 0;
    lb(end+1:end+2) = 0;
    ub(end+1:end+2) = sum(ub(nonzeros(Organize.SpinReserveStates(1,1:nG))));
end
%% Next organize equality equations (Aeq, and Demand to locate beq later)
% IC for each generator and storage device
% Electric energy balance @ each Electric subNet node  at t = 1, including transmission lines
% Heat balance @ each DistricHeat subNet node at t = 1...
% Cooling balance @ each DistrictCool subNet node at t = 1
% Any generator link equalities (linking states within a generator)
% Repeat power balance equalities and link equalities at t = 2:nS
req = 0; % row index of the Aeq matrix and beq vector
beq = [];
Aeq = [];

% The following puts together the energy balance equations
% 1 energy/mass balance equation for each subNet node
% Nodes were agregated if their line efficiencies were 1
storRow = zeros(nG,1);
hydroRowMass = zeros(nG,1);
hydroRowEnergy = zeros(nG,1);
nLcum = 0; %cumulative line #
for net = 1:1:length(networkNames)
    n = length(Plant.subNet.(networkNames{net}));
    Organize.Balance.(networkNames{net}) = [];
    Organize.Demand.(networkNames{net}) = [];
    for i = 1:1:n
        req = req+1;%there is an energy/mass balance at this node
        Organize.Balance.(networkNames{net})(end+1) = req;
        %%identify generators at this node
        genI = Plant.subNet.(networkNames{net})(i).Equipment;
        for j = 1:1:length(genI)
            states = Organize.States{genI(j)};%states associated with gerator i
            if any(strcmp(Plant.Generator(genI(j)).Type,{'Electric Storage';'Thermal Storage';}))
                storRow(genI(j)) = req;%record where you put values, so you can add t-1 state and divide by dt(t)
            elseif strcmp(Plant.Generator(genI(j)).Type,'Hydro Storage')
                if strcmp(networkNames{net},'Electrical')
                    hydroRowEnergy(genI(j)) = req; %stor the mass balance equality row
                elseif strcmp(networkNames{net},'Hydro')
                    hydroRowMass(genI(j)) = req; %stor the mass balance equality row
                end
            elseif strcmp(networkNames{net},'Electrical') && isfield(Plant.Generator(genI(j)).(Op).output,'E')
                Aeq(req,states) = Plant.Generator(genI(j)).(Op).output.E;
            elseif strcmp(networkNames{net},'DistrictHeat') && isfield(Plant.Generator(genI(j)).(Op).output,'H')
                Aeq(req,states) = Plant.Generator(genI(j)).(Op).output.H;
            elseif strcmp(networkNames{net},'DistrictCool') && isfield(Plant.Generator(genI(j)).(Op).output,'C')
                Aeq(req,states) = Plant.Generator(genI(j)).(Op).output.C;
            end
        end
        %%identify lines coming in and out of the node
        if ~strcmp(networkNames{net},'Hydro') %hydro is done later because of the time of transfer of the lines
            connect = Plant.subNet.(networkNames{net})(i).connections;
            for j = 1:1:length(connect)
                nName = Plant.subNet.(networkNames{net})(i).nodes(1); %name of current subnet node
                I = find(strcmp(strcat(nName,'_',networkNames{net},'_',connect{j}),Plant.subNet.lineNames.(networkNames{net})),1,'first');
                dir = 1; %forward direction
                if isempty(I)
                    I = find(strcmp(strcat(connect{j},'_',networkNames{net},'_',nName),Plant.subNet.lineNames.(networkNames{net})),1,'first');
                    dir = -1; %reverse direction
                end
                linestates = Organize.States{nG+nLcum+I};
                if length(linestates)==1 %uni-directional transfer
                    if dir == 1 %power/water is leaving this node
                        Aeq(req,linestates) = -1;
                    else %power/water is entering this node
                        Aeq(req,linestates) = Plant.subNet.lineEff.(networkNames{net})(I);
                    end
                else %bi-directional power transfer
                    if dir ==1
                        Aeq(req,linestates) = [-1,0,-1]; %forward direction (this node, i)-->(connected node), is positive, thus positive transmission is power leaving the node, the penalty from b->a is power not seen at a
                    else
                        Aeq(req,linestates) = [1,-1,0];%reverse direction (connected node)-->(this node, i), is positive, thus positive power is power entering the node, the penalty from a->b is power not seen at b
                    end
                end
            end
        end
        %%any heat loss term to balance equality
        if strcmp('DistrictHeat',networkNames{net}) && Plant.optimoptions.excessHeat == 1 && Organize.HeatVented(1,i)~=0
            Aeq(req,Organize.HeatVented(1,i)) = -1;
        end
        %%note any demands at this node
        if ~isempty(Plant.subNet.(networkNames{net})(i).Load)
            Organize.Demand.(networkNames{net})(i) = Plant.subNet.(networkNames{net})(i).Load;
        end
    end
    nLcum = nLcum + nLinet(net);
end

for i = 1:1:nG
    %link is a field if there is more than one state and the states are linked by an inequality or an equality
    if isfield(Plant.Generator(i).(Op),'link') && isfield(Plant.Generator(i).(Op).link,'eq')
        [m,~] = size(Plant.Generator(i).(Op).link.eq);
        states = Organize.States{i};%states associated with gerator i
        for k = 1:1:m
            Aeq(req+k,states) = Plant.Generator(i).(Op).link.eq(k,:);
            beq(req+k) = Plant.Generator(i).(Op).link.beq(k);
        end
        Organize.Equalities(i) = {req+1:req+m}; 
        req = req+m;
    end
end

%% Organize inequalities
% 2 ramping constraints for each generator at each time step
% constraints for each Energy storage 
% 2 constraints for each bi-directional transmission line
% repeat all of the above for each time 2:nS

r = 0; % row index of the A matrix & b vector
%Ramping Inequalities
Ramping = zeros(nG,1);
HydroRow = zeros(nG,1);
for i = 1:1:nG
    if strcmp(Plant.Generator(i).Type,'Hydro Storage')
        HydroRow(i) = r+1;
        r = r+3; %3 constraints for a hydro plant
    elseif isfield(Plant.Generator(i).(Op),'Ramp') 
        Ramping(i) = r+1; 
        r = r+2;
    end
end
%Generator inequalities
for i = 1:1:nG 
    if isfield(Plant.Generator(i).(Op),'link') && isfield(Plant.Generator(i).(Op).link,'ineq')%link is a field if there is more than one state and the states are linked by an inequality or an equality
        [m,~] = size(Plant.Generator(i).(Op).link.ineq);
        Organize.Inequalities(i) = {r+1:r+m}; 
        r = r+m;
    end
end

%Transmission line inequalities (penalty terms)
nLcum = 0; %cumulative line #
for net = 1:1:length(networkNames)
    if ~strcmp(networkNames{net},'Hydro')
        eff = Plant.subNet.lineEff.(networkNames{net});
        for i = 1:1:nLinet(net)
            if length(eff(i,:))==1 || eff(i,2)==0 %uni-directional transfer, 1 state for each line
                %do nothing, no inequalities linking penalty states
            else%bi-directional power transfer
                Organize.Transmission(nLcum+i) = {r+1:r+2}; 
                r = r+2;
            end
        end
    end
    nLcum = nLcum + nLinet(net);
end

%spinning reserve inequalities (sum all spinning reserves & individual spinning reserves) 
%currently only implemented for electric power
SpinRow = zeros(nG,1);
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Electrical') && Plant.optimoptions.SpinReserve
        Organize.SpinReserve = r+1;
        r = r+1;
        for i = 1:1:nG
            if isfield(Plant.Generator(i).OpMatA.output,'E') && ~isfield(Gen.output,'C') && (~isempty(strfind(Plant.Generator(i).Type,'Storage')) || Organize.Dispatchable(i))%electric storage & dispatchable electric generators have spinning reserve capacity
                SpinRow(i) = r+1;
                r = r+2;
            end
        end
    end
end
        
%% Finally, expand by the number of time steps
% # of states per time step is xL 
% Order of states will repeat for nS time steps
t1States = xL;
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

QP.lb(1:ic)  = -inf;
QP.ub(1:ic)  = inf;
for t= 1:nS
    QP.H(ic+(t-1)*t1States+1:ic+t*t1States)  = H';
    QP.f(ic+(t-1)*t1States+1:ic+t*t1States)  = f';
    QP.lb(ic+(t-1)*t1States+1:ic+t*t1States)  = lb';
    QP.ub(ic+(t-1)*t1States+1:ic+t*t1States)  = ub';
end
QP.H = diag(QP.H);

% number of generator equality constraints & energy balances at each time step will be req
t1Balances = req;
totalEqualities = nS*t1Balances + ic;
QP.Aeq = zeros(totalEqualities,totalStates);
QP.beq = zeros(totalEqualities,1);
[m,n] = size(Aeq);
if req>m || t1States>n
    Aeq(req,t1States) = 0;%make sure Aeq is the right size
end
if req>length(beq)
    beq(req,1) = 0; %make sure beq is the right length
end
QP.Aeq(1:ic,1:ic) = eye(ic); %initial condition identity matrix
for t= 1:nS
    QP.Aeq(ic+(t-1)*t1Balances+1:ic+t*t1Balances,ic+(t-1)*t1States+1:ic+t*t1States) = Aeq;
    QP.beq(ic+(t-1)*t1Balances+1:ic+t*t1Balances)  = beq;
end
%storage
for i = 1:1:nG
    if any(strcmp(Plant.Generator(i).Type,{'Electric Storage';'Thermal Storage';}))
        states = Organize.States{i};
        eff = Plant.Generator(i).(Op).Stor.DischEff;
        for t= 1:nS
            QP.Aeq((t-1)*t1Balances+storRow(i)+ic,(t-1)*t1States+states(1)+ic) = -eff/dt(t); %SOC at t
            if ismember('Y',Plant.Generator(i).(Op).states)
                QP.Aeq((t-1)*t1Balances+storRow(i)+ic,(t-1)*t1States+states(2)+ic) = -1/dt(t); %charging penalty at t
            end
            if t==1
                storIC = nnz(Organize.IC(1:i)); %order in IC
                QP.Aeq(storRow(i)+ic,storIC) = eff/dt(t);%SOC at IC
            else
                QP.Aeq((t-1)*t1Balances+storRow(i)+ic,(t-2)*t1States+states(1)+ic) = eff/dt(t); %SOC at t-1
            end
        end
    end
end

%Hydro Equalities
if any(hydroRowMass)
    nLcum = 0; %cumulative line #
    for net = 1:1:length(networkNames)
        if strcmp(networkNames{net},'Hydro')
            upriver = {};
            downriver = {};
            spill = {};
            upLines = [];
            downLines = [];
            spillLines = [];
            for i = 1:1:length(Plant.subNet.lineNames.Hydro)
                name = Plant.subNet.lineNames.Hydro{i};
                k = strfind(name,'_');
                if strcmp(name(k(1)+1:k(2)-1),'Spill')
                    spill(end+1) = {name(1:k(1)-1)};
                    spillLines(end+1) = i;
                else
                    upriver(end+1) = {name(1:k(1)-1)}; %node names of the upriver node (origin of line segment)
                    upLines(end+1) = i;
                    downriver(end+1) = {name(k(2)+1:end)};% node names of the downriver node
                    downLines(end+1) = i;
                end
            end
            for m = 1:1:length(Plant.subNet.Hydro)
                equip = Plant.subNet.Hydro(m).Equipment;
                J = upLines(strcmp(Plant.subNet.Hydro(m).nodes,upriver));%lines leaving this node, i.e. this node is the upriver node (should be at most 1)
                K = downLines(strcmp(Plant.subNet.Hydro(m).nodes,downriver));%lines entering this node, i.e. this node is the downriver node
                S = spillLines(strcmp(Plant.subNet.Hydro(m).nodes,spill));% spill flow of this node, subtracts from the power generation  (should be at most 1)
                for t = 1:1:nS
                    %river segments flowing into this node
                    for j = 1:1:length(K)
                        T = Plant.subNet.lineTime.Hydro(K(j));
                        tt = sum(dt(1:t));
                        if tt<T                            
                            %Do nothing; the upriver flow rate will be updated in updateMatrices
                        elseif tt>T && tt<=T+dt(1)%between initial condition & first step
                            frac = (tt-T)/dt(1);%portion of flow from step 1, remainder from line initial condition
                            QP.Aeq((t-1)*t1Balances+hydroRowMass(m)+ic,Organize.States{nG+nLcum+K(j)}+ic)= frac/(12.1*dt(t)); % Qupriver at step 1, conversion factor is from 1000 ft^3/s to 1000 acre ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
                            QP.Aeq((t-1)*t1Balances+hydroRowMass(m)+ic,nnz(Organize.IC(1:nG+nLcum+K(j))))= (1-frac)/(12.1*dt(t)); % Qupriver initial condition, conversion factor is from 1000 ft^3/s to 1000 acre ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
                        else
                            step = 2;
                            while tt>(T+sum(dt(1:step)))
                                step = step+1;
                            end
                            frac = (tt-(T+sum(dt(1:step))))/dt(step);%portion of flow from step, remainder from previous step
                            QP.Aeq((t-1)*t1Balances+hydroRowMass(m)+ic,(step-1)*t1States+Organize.States{nG+nLcum+K(j)}+ic)= frac/(12.1*dt(t)); % Qupriver, conversion factor is from 1000 ft^3/s to 1000 acre ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
                            QP.Aeq((t-1)*t1Balances+hydroRowMass(m)+ic,(step-2)*t1States+Organize.States{nG+nLcum+K(j)}+ic)= (1-frac)/(12.1*dt(t)); % Qupriver, conversion factor is from 1000 ft^3/s to 1000 acre ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
                        end 
                    end 
                    %water flow out of the node
                    QP.Aeq((t-1)*t1Balances+hydroRowMass(m)+ic,(t-1)*t1States+Organize.States{nG+nLcum+J}+ic)= -1/(12.1*dt(t)); %Qdownriver, conversion factor is from 1000 ft^3/s to 1000 acre ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
                    
                end 
                for k = 1:1:length(equip)
                    I = equip(k); %Index in generator list
                    if strcmp(Plant.Generator(I).Type,'Hydro Storage')
                        for t = 1:1:nS
                            states = Organize.States{I};%states associated with generator i
                            %SOC of the reservior in 1000 acre ft.
                            QP.Aeq((t-1)*t1Balances+hydroRowMass(m)+ic,(t-1)*t1States+states(1)+ic) = 1/dt(t); %SOC at t
                            if t==1
                                hydroIC = nnz(Organize.IC(1:i)); %order in IC
                                QP.Aeq(hydroRowMass(m),hydroIC) = -1/dt(t);%SOC at IC
                            else
                                QP.Aeq((t-1)*t1Balances+hydroRowMass(m)+ic,(t-2)*t1States+states(1)+ic) = -1/dt(t); %SOC at t-1
                            end
                            %Converting water flow to power
                            QP.Aeq((t-1)*t1Balances+hydroRowEnergy(m)+ic,(t-1)*t1States+Organize.States{nG+nLcum+J}+ic)= Plant.Generator(I).OpMatA.output.E/dt(t); %Qdownriver * energy conversion factor (1000 ft^3/s to kW)
                            if ~isempty(S)
                                QP.Aeq((t-1)*t1Balances+hydroRowEnergy(m)+ic,(t-1)*t1States+Organize.States{nG+nLcum+S}+ic) = -Plant.Generator(I).OpMatA.output.E/dt(t); %Spillway flow (subtracted from downriver flow)
                            end
                        end
                    else
                        %% add water district here
                    end
                end
            end
        end
        nLcum = nLcum + nLinet(net);
    end  
end  

% number of generator inequality constraints & energy imbalances at each time step will be r 
% there are 2 ramping constraints on each generator/storage
t1ineq = r;
totalInequalities = nS*t1ineq;
QP.A = zeros(totalInequalities,totalStates);
QP.b = zeros(totalInequalities,1);
r = 0;
s = ic;
for t= 1:nS
    for i = 1:1:nG
        states = Organize.States{i};
        %Ramping (Hydro ramping handled later)
        if Ramping(i)>0
            %%if storage, ramping only affects 1st state)
            if ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
                rampStates = states(1);
            else rampStates = states;
            end
            QP.A(r+Ramping(i),s+rampStates) = 1/dt(t);%ramp up 
            QP.A(r+Ramping(i)+1,s+rampStates) = -1/dt(t);%ramp down
            if t ==1 %ramping from initial condition
                QP.A(Ramping(i),nnz(Organize.IC(1:i))) = -1/dt(t); %ramp up 
                QP.A(Ramping(i)+1,nnz(Organize.IC(1:i))) = 1/dt(t); %ramp down 
            else %condition at previous time step
                QP.A(r+Ramping(i),s-t1States+rampStates) = -1/dt(t); %ramp up 
                QP.A(r+Ramping(i)+1,s-t1States+rampStates) = 1/dt(t); %ramp down 
            end
            QP.b(r+[Ramping(i),Ramping(i)+1],1) = Plant.Generator(i).(Op).Ramp.b;
        end
        %Inequalities constraints
        ineqRow = Organize.Inequalities{i};
        for k = 1:1:length(ineqRow)
            if k == 1 && ~isempty(strfind(Plant.Generator(i).Type,'Storage')) && ismember('Y',Plant.Generator(i).(Op).states)
                QP.A(r+ineqRow(k),s+states(1)) = Plant.Generator(i).(Op).link.ineq(1,1)/dt(t);  % SOC at t  
                QP.A(r+ineqRow(k),s+states(2)) = Plant.Generator(i).(Op).link.ineq(1,2)/dt(t);  % charging state at t: value is -1/(1-efficiency)
                if t ==1 %SOC change from IC
                    QP.A(ineqRow(k),nnz(Organize.IC(1:i))) = -Plant.Generator(i).(Op).link.ineq(1,1)/dt(t); % SOC at t-1
                else
                    QP.A(r+ineqRow(k),s-t1States+states(1)) = -Plant.Generator(i).(Op).link.ineq(1,1)/dt(t);  % SOC at t-1
                end
            else
                QP.A(r+ineqRow(k),s+states) = Plant.Generator(i).(Op).link.ineq(k,:);
            end
            QP.b(r+ineqRow(k)) = Plant.Generator(i).(Op).link.bineq(k);
        end
    end
    %Transmission
    %%no connection to previous or later time steps, and no dependence on step size. 
    for i = 1:1:nL
        lineRow = Organize.Transmission{i};
        if~isempty(lineRow)
            QP.A(r+lineRow(1),s+Organize.States{nG+i}) = [(1-eff(i,1)), -1, 0];% Pab*(1-efficiency) < penalty a to b
            QP.A(r+lineRow(2),s+Organize.States{nG+i}) = [-(1-eff(i,2)), 0, -1];% -Pab*(1-efficiency) < penalty b to a
        end
    end
    
    %spinning reserve inequalities (sum all spinning reserves & individual spinning reserves) 
    %currently only implemented for electric power
    if Plant.optimoptions.SpinReserve
        SRstates = nonzeros(Organize.SpinReserveStates(1,1:nG+1));
        QP.A(r+Organize.SpinReserve,s+SRstates) = -1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
        SRancillary = Organize.SpinReserveStates(1,nG+2);
        QP.A(r+Organize.SpinReserve,s+SRancillary) = 1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
        for i = 1:1:nG
            SRstate = Organize.SpinReserveStates(1,i);
            states = Organize.States{i};
            if SpinRow(i)~=0 && isempty(strfind(Plant.Generator(i).Type,'Storage'))% dispatchable generators have spinning reserve capacity
                QP.A(r+SpinRow(i),s+SRstate) = 1/dt(t); %SR + power(t) - power(t-1)<= ramprate*dt
                QP.A(r+SpinRow(i),s+states) = 1/dt(t);
                if t ==1 %ramping from IC
                    QP.A(SpinRow(i),nnz(Organize.IC(1:i))) = -1/dt(t); % Power at t-1
                else
                    QP.A(r+SpinRow(i),s-t1States+states) = -1/dt(t); % Power at t-1
                end
                QP.b(r+SpinRow(i)) = Plant.Generator(i).(Op).Ramp.b(1);%ramp up constraint
                
                QP.A(r+SpinRow(i)+1,s+SRstates) = 1; %SR + power <= Size
                QP.A(r+SpinRow(i)+1,s+states) = 1;
                QP.b(r+SpinRow(i)+1) = Plant.Generator(i).Size; %max capacity constraint
            elseif SpinRow(i)~=0 
                if ~strcmp(Plant.Generator(i).Type,'Hydro Storage')%electric storage 
                    eff = Plant.Generator(i).(Op).Stor.DischEff;
                    QP.A(r+SpinRow(i),s+SRstate) = 1; %SR + eff*(SOC(t-1) - SOC(t))/dt <= peak discharge
                    QP.A(r+SpinRow(i),s+states(1)) = -eff/dt(t);
                    QP.A(r+SpinRow(i)+1,s+SRstate) = 1; %SR - SOC(t-1)/dt <= 0
                    if t ==1 %SOC change from IC
                        QP.A(SpinRow(i),nnz(Organize.IC(1:i))) = eff/dt(t); % SOC at t-1
                        QP.A(SpinRow(i)+1,nnz(Organize.IC(1:i))) = -eff/dt(t); % SOC at t-1
                    else
                        QP.A(r+SpinRow(i),s-t1States+states) = eff/dt(t); % SOC at t-1
                        QP.A(r+SpinRow(i)+1,s-t1States+states) = -eff/dt(t); % SOC at t-1
                    end
                    QP.b(r+SpinRow(i)) = Plant.Generator(i).(Op).Ramp.b(2);%peak discharge constraint
                else
                    %%Hydro spinning reserve
                    
                end
            end
        end
    end
    r = r+t1ineq;
    s = s+t1States;
end

%% Hydro Inequalities (2 ramping constraints and max generator flow (all applied to the flow through the generators)
nLcum = 0; %cumulative line #
for net = 1:1:length(networkNames)
    if ~strcmp(networkNames{net},'Hydro')
        nLcum = nLcum + nLinet(net);
    else
        for m = 1:1:length(Plant.subNet.Hydro)
            equip = Plant.subNet.Hydro(m).Equipment;
            J = upLines(strcmp(Plant.subNet.Hydro(m).nodes,upriver));%lines leaving this node, i.e. this node is the upriver node (should be at most 1)
            S = spillLines(strcmp(Plant.subNet.Hydro(m).nodes,spill));% spill flow of this node, subtracts from the power generation  (should be at most 1)
            for k = 1:1:length(equip)
                I = equip(k); %Index in generator list
                if strcmp(Plant.Generator(I).Type,'Hydro Storage')
                    r = 0;
                    s = ic;
                    for t = 1:1:nS
                        QP.A(r+HydroRow(I),s+Organize.States{nG+nLcum+J}) = 1; %max flow constraint
                        QP.A(r+HydroRow(I)+1,s+Organize.States{nG+nLcum+J}) = Plant.Generator(I).(Op).output.E;%ramp up constraint total flow at t
                        QP.A(r+HydroRow(I)+2,s+Organize.States{nG+nLcum+J}) = -Plant.Generator(I).(Op).output.E;%ramp down constraint total flow at t
                        
                        QP.b(r+HydroRow(I)) = Plant.Generator(I).VariableStruct.MaxGenFlow;%max flow constrain
                        QP.b(r+HydroRow(I)+1) = Plant.Generator(I).VariableStruct.RampUp;%max increase in power
                        QP.b(r+HydroRow(I)+2) = Plant.Generator(I).VariableStruct.RampDown;%max decrease in power
                        if t == 1
                            QP.A(r+HydroRow(I)+1,nnz(Organize.IC(1:nG+nLcum+J))) = -Plant.Generator(I).(Op).output.E;%ramp up constraint total flow at t-1
                            QP.A(r+HydroRow(I)+2,nnz(Organize.IC(1:nG+nLcum+J))) = Plant.Generator(I).(Op).output.E;%ramp down constraint total flow at t-1
                        else
                            QP.A(r+HydroRow(I)+1,s-t1States+Organize.States{nG+nLcum+J}) = -Plant.Generator(I).(Op).output.E;%ramp up constraint total flow at t-1
                            QP.A(r+HydroRow(I)+2,s-t1States+Organize.States{nG+nLcum+J}) = Plant.Generator(I).(Op).output.E;%ramp down constraint total flow at t-1
                        end
                        if ~isempty(S) %spill flow (power flow = total flow - spill flow)
                            QP.A(r+HydroRow(I),s+Organize.States{nG+nLcum+S}) = -1; %max flow constraint
                            QP.A(r+HydroRow(I)+1,s+Organize.States{nG+nLcum+S}) = -1; %ramp up constraint spill flow at t
                            QP.A(r+HydroRow(I)+2,s+Organize.States{nG+nLcum+S}) = 1; %ramp down constraint spill flow at t
                            if t == 1
                                QP.A(r+HydroRow(I)+1,nnz(Organize.IC(1:nG+nLcum+S))) = Plant.Generator(I).(Op).output.E;%ramp up constraint spill flow at t-1
                                QP.A(r+HydroRow(I)+2,nnz(Organize.IC(1:nG+nLcum+S))) = -Plant.Generator(I).(Op).output.E;%ramp down constraint spill flow at t-1
                            else
                                QP.A(r+HydroRow(I)+1,s-t1States+Organize.States{nG+nLcum+S}) = Plant.Generator(I).(Op).output.E;%ramp up constraint spill flow at t-1
                                QP.A(r+HydroRow(I)+2,s-t1States+Organize.States{nG+nLcum+S}) = -Plant.Generator(I).(Op).output.E;%ramp down constraint spill flow at t-1
                            end
                        end
                        r = r+t1ineq;
                        s = s+t1States;
                    end
                end
            end
        end
    end
end


%% Put together indexing of matrix equation order and state order for use later
%list of states associated with each generator at all time steps
for i = 1:1:nG+nL
    if ~isempty(Organize.States{i})
        s = Organize.States{i};
        n = length(s);
        states = [];
        for t = 1:1:nS
            states(end+1:end+n) = s + (t-1)*t1States + ic;
        end
        Organize.States(i) = {states};
    end
end

%setup indexing for beq (and b) affiliated with each demand
for net = 1:1:length(networkNames)
    bal = Organize.Balance.(networkNames{net});
    n = length(bal);
    Organize.Balance.(networkNames{net}) = zeros(n,nS);
    for t = 1:1:nS
        Organize.Balance.(networkNames{net})(1:n,t) = bal + (t-1)*t1Balances+ic;
    end
end

% Equalities for each generator
for i = 1:1:nG
    if isfield(Plant.Generator(i).(Op),'link') && isfield(Plant.Generator(i).(Op).link,'eq')
        s = Organize.Equalities{i};
        n = length(s);
        rows = [];
        for t = 1:1:nS
            rows(end+1:end+n) = s + (t-1)*t1Balances+ic;
        end
    else
        rows = [];
    end
    Organize.Equalities(i) = {rows};
end

% Ramping each generator/storage
for i = 1:1:nG
    if isfield(Plant.Generator(i).(Op),'Ramp') 
        ramp = [];
        for t = 1:1:nS
            ramp(end+1:end+2) =(t-1)*t1ineq+[Ramping(i), Ramping(i)+1];
        end
        Organize.Ramping(i) = {ramp};
    end
end

% Inequalities for each generator
for i = 1:1:nG
    if (isfield(Plant.Generator(i).(Op),'link') && isfield(Plant.Generator(i).(Op).link,'ineq'))
        s = Organize.Inequalities{i};
        n = length(s);
        ineqRow = [];
        for t = 1:1:nS
            ineqRow(end+1:end+n) = s + (t-1)*t1ineq;
        end
        Organize.Inequalities(i) = {ineqRow};
    end
end

%transmission inequalities
for i = 1:1:nL
    s = Organize.Transmission{i};
    line = [];
    if ~isempty(s) %only bi-directional lines
        for t = 1:1:nS
            line(end+1:end+2) = s + (t-1)*t1ineq;
        end
    end
    Organize.Transmission(i) = {line};
end

%heat loss terms
if any(strcmp('DistrictHeat',networkNames)) && Plant.optimoptions.excessHeat == 1
    for i = 1:1:length(Plant.subNet.('DistrictHeat'))
        if Organize.HeatVented(1,i)~=0
            Organize.HeatVented(:,i) = (Organize.HeatVented(1,i)+ic:t1States:totalStates)'; %matrix of the states associated with venting heat, organized by the nodes in the district heating network
        end
    end
end

%spin reserve
if Plant.optimoptions.SpinReserve
    %spin reserve states
    for i = 1:1:nG+2
        if Organize.SpinReserveStates(1,i)~=0
            Organize.SpinReserveStates(1:end,i) = (Organize.SpinReserveStates(1,i):t1States:totalStates)'; %shortened version works because there is 1 state per time step
        end
    end
    %spin reserve inequalities (2 per generator)
    Organize.SpinRow = cell(nG,1);
    for i = 1:1:nG
        if SpinRow(i)~=0
            row = [];
            for t = 1:1:nS
                row(end+1:end+2) =(t-1)*t1ineq + [SpinRow(i),SpinRow(i)+1];
            end
            Organize.SpinRow(i) = {row};
        end
    end
    %cumulative spin reserve rows in inequality
    Organize.SpinReserve = (Organize.SpinReserve:t1ineq:nS*t1ineq)';
end

QP.Organize = Organize; %indices (rows and columns) associated with each generator, allowing generators to be removed later