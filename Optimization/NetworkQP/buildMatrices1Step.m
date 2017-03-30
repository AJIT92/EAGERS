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
nL = length(Plant.subNet.lineNames);

Organize.States = cell(1,nG);
Organize.Equalities = cell(1,nG);
Organize.Inequalities = cell(1,nG);
Organize.StorageEquality = zeros(1,nG); %row of Aeq coresponding to energy balance that includes this storage device (expected storage output shows up as negative in beq)
Organize.Dispatchable = zeros(1,nG);
QP.organize = cell(1,nG+nL);
QP.constCost = zeros(1,nG);

networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
%% First organize the order of states (set cost, and bounds: H, f, lb, ub)
% states for each generator/storage
% states for each transmission line/pipe 
H = []; f = []; lb =[]; ub = [];
xL = 0;
for i = 1:1:nG
    Gen = Plant.Generator(i).(Op);
    s = length(Gen.states);%generator with multiple states
    if s>0
        if ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
            QP.organize{1,i} = xL+1; % Storage treated as generator with 1 state
            Organize.States(i)= {xL+1};
            H(end+1) = 0;
            f(end+1) = 0;
            lb(end+1) = -Gen.Ramp.b(1);
            ub(end+1) = Gen.Ramp.b(2);
            xL = xL + 1;
        else
            QP.organize{1,i} = xL+1:xL+s; %output is sum of multiple states at ech time step
            Organize.States(i)= {xL+1:xL+s};
            for k = 1:1:s
                H(end+1) = Gen.(Gen.states{k}).H;
                f(end+1) = Gen.(Gen.states{k}).f;
                lb(end+1) = Gen.(Gen.states{k}).lb;
                ub(end+1) = Gen.(Gen.states{k}).ub;
            end
            xL = xL + s;
            if isempty(strfind(Plant.Generator(i).Type,'Utility'))
                Organize.Dispatchable(i) = 1;
                if isfield(Plant.Generator(i).(Op),'constCost') 
                     QP.constCost(i) = Plant.Generator(i).OpMatB.constCost;
                end
            end
        end
    end
end
for i = 1:1:nL % 3 states for each line (state of the line and penalty term in each direction)
    QP.organize{1,nG+i} = xL+1; %line state 
    Organize.States(nG+i)= {[xL+1, xL+2, xL+3]};
    H(end+1:end+3) = [0 0 0];
    f(end+1:end+3) = [0 0 0];
    lb(end+1:end+3) = [-Plant.subNet.lineLimit(i),0,0];
    ub(end+1:end+3) = [Plant.subNet.lineLimit(i), Plant.subNet.lineLimit(i)*(1-Plant.subNet.lineEff(i)), Plant.subNet.lineLimit(i)*(1-Plant.subNet.lineEff(i))];
    xL = xL + 3;
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
    for i = 1:1:n
        if strcmp(networkNames{net},'DistrictHeat') && Plant.optimoptions.excessHeat == 1
            QP.excessHeat =1;
            %skip and do in inequality section
        else
            req = req+1;%there is an energy balance at this node
            Organize.Balance.(networkNames{net})(n) = req;
            %%identify generators at this node
            genI = Plant.subNet.(networkNames{net})(i).Equipment;
            for j = 1:1:length(genI)
                states = Organize.States{genI(j)};%states associated with gerator i
                if strcmp(networkNames{net},'Electrical') && isfield(Plant.Generator(genI(j)).(Op).output,'E')
                    if strcmp(Plant.Generator(genI(j)).Type,'Electric Storage')
                        Aeq(req,states) = 1; %storage converted to "generator"
                        Organize.StorageEquality(genI(j)) = req;
                    else
                        Aeq(req,states) = Plant.Generator(genI(j)).(Op).output.E;
                    end
                elseif strcmp(networkNames{net},'DistrictHeat') && isfield(Plant.Generator(genI(j)).(Op).output,'H')
                    if strcmp(Plant.Generator(genI(j)).Type,'Thermal Storage')
                        Aeq(req,states) = 1; %storage converted to "generator"
                        Organize.StorageEquality(genI(j)) = req;
                    else
                        Aeq(req,states) = Plant.Generator(genI(j)).(Op).output.H;
                    end
                elseif strcmp(networkNames{net},'DistrictCool') && isfield(Plant.Generator(genI(j)).(Op).output,'C')
                    if strcmp(Plant.Generator(genI(j)).Type,'Thermal Storage')
                        Aeq(req,states) = 1; %storage converted to "generator"
                        Organize.StorageEquality(genI(j)) = req;
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
                Organize.Demand.(networkNames{net})(load(j)) = req;
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
% 2 constraints for each transmission line

r = 0; % row index of the A matrix & b vector
A = [];
%Distric heating energy inequalities (inequality because heat can be rejected)
for net = 1:1:length(networkNames)
    n = length(Plant.subNet.(networkNames{net}));
    Organize.Imbalance.(networkNames{net}) = [];
    for i = 1:1:n
        if strcmp(networkNames{net},'DistrictHeat') && Plant.optimoptions.excessHeat == 1
            r = r+1;%there is an energy balance at this node
            Organize.Imbalance.(networkNames{net})(n) = r;
            %%identify generators at this node
            genI = Plant.subNet.(networkNames{net})(i).Equipment;
            for j = 1:1:length(genI)
                states = Organize.States{genI(j)};%states associated with gerator i
                if strcmp(networkNames{net},'DistrictHeat') && isfield(Plant.Generator(genI(j)).(Op).output,'H')
                    if strcmp(Plant.Generator(genI(j)).Type,'Thermal Storage')
                        A(r,states) = 1;
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
                Organize.Demand.(networkNames{net})(load(j)) = r;
            end
        end
    end
end

%Transmission line inequalites (penalty terms)
for i = 1:1:length(Plant.subNet.lineNames)
    linestates = Organize.States{nG+i};
    A(r+1,linestates) = [(1-Plant.subNet.lineEff(i)), -1, 0];% Pab*(1-efficiency) < penalty a to b
    A(r+2,linestates) = [-(1-Plant.subNet.lineEff(i)), 0, -1];% -Pab*(1-efficiency) < penalty b to a
    Organize.Transmission(i) = {strcat('[',num2str(r+1),':',num2str(r+2),']')}; 
    r = r+2;
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