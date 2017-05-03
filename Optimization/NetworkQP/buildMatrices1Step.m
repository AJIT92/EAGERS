function QP = buildMatrices1Step
%builds constant matrices for step-by-step optimization
%similar to multi-step except:
%energy storage looks like a generator with 1 state (can be pos or neg)
%no ramping constraints
%Fit A includes energy storage and uses the fit with zero y-intercept
%Fit B does not include energy storage and uses non-zero intercept
%Demands, upper/lower bounds, and utility costs must updated prior to optimization
global Plant 
Op = 'OpMatB';
nG = length(Plant.Generator);

networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
for net = 1:1:length(networkNames)
    nLinet(net) = length(Plant.subNet.lineNames.(networkNames{net}));
end
nL = sum(nLinet);

Organize.States = cell(1,nG);
Organize.Equalities = cell(1,nG);
Organize.Inequalities = cell(1,nG);
Organize.StorageEquality = zeros(1,nG); %row of Aeq coresponding to energy balance that includes this storage device (expected storage output shows up as negative in beq)
Organize.Dispatchable = zeros(1,nG);
Organize.SpinReserveStates = zeros(1,nG+2);
Organize.Transmission = cell(1,nL);
QP.organize = cell(1,nG+nL);
QP.constCost = zeros(1,nG);
%% First organize the order of states (set cost, and bounds: H, f, lb, ub)
% states for each generator/storage
% spinning reserve for generator/storage if option is selected
% states for each transmission line/pipe 
% state for P_loss in heating energy balance
% state for spinning reserve shortfall (target - cumulative of all active generators) and SR provided as ancillary service
H = []; f = []; lb =[]; ub = [];
xL = 0;
for i = 1:1:nG
    Gen = Plant.Generator(i).(Op);
    s = length(Gen.states);%generator with multiple states
    if s>0
        if ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
            s = 1;% Storage treated as generator with 1 state
            H(end+1) = 0;
            f(end+1) = 0;
            lb(end+1) = -Gen.Ramp.b(1);
            ub(end+1) = Gen.Ramp.b(2);
        else %dispatchable generator or utility
            for k = 1:1:s
                H(end+1) = Gen.(Gen.states{k}).H;
                f(end+1) = Gen.(Gen.states{k}).f;
                lb(end+1) = Gen.(Gen.states{k}).lb;
                ub(end+1) = Gen.(Gen.states{k}).ub;
            end
            if isempty(strfind(Plant.Generator(i).Type,'Utility'))
                Organize.Dispatchable(i) = 1;
                if isfield(Plant.Generator(i).(Op),'constCost') 
                     QP.constCost(i) = Plant.Generator(i).OpMatB.constCost;
                end
            end
        end
        QP.organize{1,i} = xL+1:xL+s; %output is sum of multiple states at each time step
        Organize.States(i)= {xL+1:xL+s};
        if Plant.optimoptions.SpinReserve && isfield(Gen.output,'E')
            if ~isempty(strfind(Plant.Generator(i).Type,'Storage')) || Organize.Dispatchable(i) %storage and dispatchable generators have spinning reserve capacity
                Organize.SpinReserveStates(1,i) = xL + s +1; %state of spinning reserve at time 1
                H(end+1) = 0;
                f(end+1) = 0;
                lb(end+1) = 0;
                ub(end+1) = Plant.Generator(i).Size;
            end
            xL = xL + s + 1;
        else xL = xL + s;
        end
        
        
    end
end
%states for transmission lines
nLcum = 0; %cumulative line #
for net = 1:1:length(networkNames)
    for i = 1:1:nLinet(net) 
        nLcum = nLcum+1;
        if strcmp(networkNames{net},'Hydro')
            
        else 
            eff = Plant.subNet.lineEff.(networkNames{net});
            limit = Plant.subNet.lineLimit.(networkNames{net}); 
            QP.organize{1,nG+nLcum} = xL+1; %line state 
            if length(eff(i,:))==1 || eff(i,2)==0 %uni-directional transfer, 1 state for each line
                Organize.States(nG+nLcum)= {xL+1};
                H(end+1:end+3) = 0;
                f(end+1:end+3) = 0;
                lb(end+1:end+3) = 0;
                ub(end+1:end+3) = limit(i,1);
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
end

if any(strcmp('DistrictHeat',networkNames)) && Plant.optimoptions.excessHeat == 1
    %%assume heat can be lost any any node in the network that has a device producing heat
    n = length(Plant.subNet.('DistrictHeat'));
    Organize.HeatVented = zeros(1,n);
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
    xL = xL + 2; % Add single state for reserve shortfall state which is equal to reserve target + reserve sold as ancillary service - actual spinning reserve
    H(end+1:end+2) = 0;
    f(end+1:end+2) = 0;
    lb(end+1:end+2) = 0;
    ub(end+1:end+2) = sum(ub(nonzeros(Organize.SpinReserveStates(1,1:nG))));
end

%% Next organize equality equations (Aeq, and Demand to locate beq later)
% Electric energy balance @ each Electric subNet node  at t = 1
% Heat balance @ each DistricHeat subNet node at t = 1
% Cooling balance @ each DistrictCool subNet node at t = 1
% Any generator link equalities (linking states within a generator)
req = 0; % row index of the Aeq matrix and beq vector
beq = [];

% The following puts together the energy balance equations
% 1 equation for each subNet node
% Nodes were agregated if their line efficiencies were 1
QP.excessHeat = Plant.optimoptions.excessHeat;
for net = 1:1:length(networkNames)
    n = length(Plant.subNet.(networkNames{net}));
    Organize.Balance.(networkNames{net}) = [];
    Organize.Demand.(networkNames{net}) = [];
    if strcmp(networkNames{net},'Hydro')
        
    else
        for i = 1:1:n
            req = req+1;%there is an energy balance at this node
            Organize.Balance.(networkNames{net})(i) = req;
            %%identify generators at this node
            genI = Plant.subNet.(networkNames{net})(i).Equipment;
            for j = 1:1:length(genI)
                states = Organize.States{genI(j)};%states associated with gerator i
                if any(strcmp(Plant.Generator(genI(j)).Type,{'Electric Storage';'Thermal Storage';}))
                    Aeq(req,states(1)) = 1; %storage converted to "generator"
                    Organize.StorageEquality(genI(j)) = req;
                elseif strcmp(networkNames{net},'Electrical') && isfield(Plant.Generator(genI(j)).(Op).output,'E')
                    Aeq(req,states) = Plant.Generator(genI(j)).(Op).output.E;
                elseif strcmp(networkNames{net},'DistrictHeat') && isfield(Plant.Generator(genI(j)).(Op).output,'H')
                    Aeq(req,states) = Plant.Generator(genI(j)).(Op).output.H;
                elseif strcmp(networkNames{net},'DistrictCool') && isfield(Plant.Generator(genI(j)).(Op).output,'C')
                    Aeq(req,states) = Plant.Generator(genI(j)).(Op).output.C;
                end
            end
            %%identify lines coming in and out
            connect = Plant.subNet.(networkNames{net})(i).connections;
            for j = 1:1:length(connect)
                nName = Plant.subNet.(networkNames{net})(i).nodes(1); %name of current subnet node
                I = find(strcmp(strcat(nName,'_',networkNames{net},'_',connect{j}),Plant.subNet.lineNames.(networkNames{net})),1,'first');
                dir = 1; %forward direction
                if isempty(I)
                    I = find(strcmp(strcat(connect{j},'_',networkNames{net},'_',nName),Plant.subNet.lineNames.(networkNames{net})),1,'first');
                    dir = -1; %reverse direction
                end
                linestates = Organize.States{nG+I};
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
            %%any heat loss term to balance equality
            if strcmp('DistrictHeat',networkNames{net}) && Plant.optimoptions.excessHeat == 1 && Organize.HeatVented(1,i)~=0
                Aeq(req,Organize.HeatVented(1,i)) = -1;
            end
            %%note any demands at this node
            Organize.Demand.(networkNames{net})(i) = Plant.subNet.(networkNames{net})(i).Load;
        end
    end
end

hydroIndex = [];
hydroRow = [];
hydroSOC = [];
for i = 1:1:nG
    %link is a field if there is more than one state and the states are linked by an inequality or an equality
    if isfield(Plant.Generator(i).(Op),'link') && isfield(Plant.Generator(i).(Op).link,'eq')
        [m,n] = size(Plant.Generator(i).(Op).link.eq);
        states = Organize.States{i};%states associated with gerator i
        for k = 1:1:m
            Aeq(req+k,states) = Plant.Generator(i).(Op).link.eq(k,:);
            beq(req+k) = Plant.Generator(i).(Op).link.beq(k);
        end
        Organize.Equalities(i) = {req+1:req+m}; 
        req = req+m;
    elseif strcmp(Plant.Generator(i).Type,'Hydro')
        states = Organize.States{i};%states associated with generator i
        Aeq(req+1,states(1:2)) = 1; %Flow for power and spill;
        hydroIndex(end+1) = i;
        hydroRow(end+1) = req+1;
        hydroSOC(end+1) = states(3);
        req = req + 1;
    end
end

%spinning reserve: no ramping constraint, so SR can be an equalities (sum all spinning reserves & individual spinning reserves) 
%currently only implemented for electric power
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Electrical') && Plant.optimoptions.SpinReserve
        Organize.SpinRow = cell(1,nG);
        for i = 1:1:nG
            if isfield(Plant.Generator(i).OpMatA.output,'E') && ~isfield(Plant.Generator(i).OpMatA.output,'C') && (~isempty(strfind(Plant.Generator(i).Type,'Storage')) || Organize.Dispatchable(i))%electric storage & dispatchable electric generators have spinning reserve capacity
                states = Organize.States{i};%states associated with gerator i
                Organize.SpinRow(i) = {req+1};
                Aeq(req+1,states) = 1;
                beq(req+1) = sum(ub(states(1:end-1)));%max power
                req = req+1;
            end
        end
    end
end
%% Organize inequalities
% No ramping constraints or energy storage inequalities (energy storage is a generator)
% 2 constraints for each bi-directional transmission line

r = 0; % row index of the A matrix & b vector
A = [];

%Transmission line inequalites (penalty terms)
nLcum = 0; %cumulative line #
for net = 1:1:length(networkNames)
    for i = 1:1:nLinet(net)
        nLcum = nLcum + 1;
        if strcmp(networkNames{net},'Hydro')

        else
            eff = Plant.subNet.lineEff.(networkNames{net});
            if length(eff(i,:))==1 || eff(i,2)==0 %uni-directional transfer, 1 state for each line
                %do nothing, no inequalities linking penalty states
            else%bi-directional power transfer
                linestates = Organize.States{nG+i};
                A(r+1,linestates) = [(1-eff(i,1)), -1, 0];% Pab*(1-efficiency) < penalty a to b
                A(r+2,linestates) = [-(1-eff(i,2)), 0, -1];% -Pab*(1-efficiency) < penalty b to a
                Organize.Transmission(nLcum) = {r+1:r+2}; 
                r = r+2;
            end
        end
    end
end
%cumulative spinning reserve shortfall: currently only implemented for electric power
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Electrical') && Plant.optimoptions.SpinReserve
        Organize.SpinReserve = r+1;
        A(r+1,nonzeros(Organize.SpinReserveStates(1,1:nG+1))) = -1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
        A(r+1,nonzeros(Organize.SpinReserveStates(1,nG+2))) = 1; %Inequality for spinning reserve shortfall:  -(shortfall) - sum(SR(i)) + SR ancillary <= -SR target
        r = r+1;
    end
end

%% Finally, build out QP matrices
QP.H = diag(H);
QP.f =f';
QP.lb = lb';
QP.ub = ub';

[m,n] = size(Aeq);
if req>m || xL>n
    Aeq(req,xL) = 0;%make sure Aeq is the right size
end
if req>length(beq)
    beq(req,1) = 0; %make sure beq is the right length
end
QP.Aeq = Aeq;
QP.beq = beq;

if ~isempty(A)
    [m,n] = size(A);
    if m<r|| n<xL
        A(r,xL) = 0;
    end
else A = zeros(r,xL);
end
QP.A = A;
QP.b = zeros(r,1);
QP.Organize = Organize; %indices (rows and columns) associated with each generator, allowing generators to be removed later