function QP = updateMatrices1Step(QP,Demand,Renewable,scaleCost,dt,EC,StorPower,MinPower,MaxPower)
% update the equalities with the correct demand, and scale fuel and electric costs
% EC is the expected end condition at this time stamp (can be empty)
% StorPower is the expected output/input of any energy storage at this timestep (can be empty)
% MinPower and MaxPower define the range of this generator at this timestep
global Plant
nG = length(Plant.Generator);
marginal = instantMarginalCost(EC,scaleCost);%update marginal cost
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));

%% update demands and storage self-discharge
for net = 1:1:length(networkNames)
    if strcmp(networkNames{net},'Electrical')
        out = 'E';
        if Plant.optimoptions.SpinReserve
            QP.b(QP.Organize.SpinReserve) = -Plant.optimoptions.SpinReservePerc/100*sum(Demand.E,2);% -shortfall + SRancillary - SR generators - SR storage <= -SR target
        end
    elseif strcmp(networkNames{net},'DistrictHeat')
        out = 'H';
    elseif strcmp(networkNames{net},'DistrictCool')
        out = 'C';
    end
    for i = 1:1:length(Plant.subNet.(networkNames{net}))
        equip = Plant.subNet.(networkNames{net})(i).Equipment;
        eq = QP.Organize.Balance.(networkNames{net})(i);
        loads = QP.Organize.Demand.(networkNames{net})(i);
        QP.beq(eq) = sum(Demand.(out)(:,loads),2); %multile demands can be at the same node
        for j = 1:1:length(equip)
            k = equip(j);
            if strcmp(networkNames{net},'Electrical') && strcmp(Plant.Generator(k).Source,'Renewable')% subtract renewable generation 
                QP.beq(eq) = QP.beq(eq) - Renewable(k); %put renewable generation into energy balance at correct node
            end
            if ~isempty(strfind(Plant.Generator(k).Type,'Storage'))
                if (strcmp(Plant.Generator(k).Source,'Electricity') && strcmp(networkNames{net},'Electrical')) || (strcmp(Plant.Generator(k).Source,'Heat') && strcmp(networkNames{net},'DistrictHeat')) || (strcmp(Plant.Generator(k).Source,'Cooling') && strcmp(networkNames{net},'DistrictCool'))
                    loss = dt*(Plant.Generator(k).OpMatA.Stor.SelfDischarge*Plant.Generator(k).OpMatA.Stor.UsableSize);
                    QP.beq(eq) = QP.beq(eq) - loss; %account for self-discharge losses
                end
            end
        end
    end
end  

%% Update upper and lower bounds based on ramping constraint (if applicable)
if ~isempty(MinPower)
    for i = 1:1:nG
        states = QP.Organize.States{i};
        if isempty(strfind(Plant.Generator(i).Type,'Storage'))%all generators and utilities
            if MinPower(i)>sum(QP.lb(states))
                QP.Organize.Dispatchable(i) = 0; %can't shut off
                %raise the lower bounds
                Padd = MinPower(i) - sum(QP.lb(states));
                for f = 1:1:length(states) %starting with lb of 1st state in generator, then moving on
                    gap = QP.ub(states(f)) - QP.lb(states(f));
                    if gap>0
                        gap = min(gap,Padd);
                        QP.lb(states(f)) = QP.lb(states(f)) + gap;
                        Padd = Padd-gap;
                    end
                end
            end
            if MaxPower(i)<Plant.Generator(i).Size %lower the upper bounds
                Psub = sum(QP.ub(states)) - MaxPower(i);
                for f = length(states):-1:1 %starting with lb of 1st state in generator, then moving on
                    if QP.ub(states(f))>0
                        gap = min(QP.ub(states(f)),Psub);
                        QP.ub(states(f)) = QP.ub(states(f)) - gap;
                        Psub = Psub-gap;
                    end
                end
            end
        else
            %update storage output range to account for what is already scheduled
            QP.lb(states) = QP.lb(states) - StorPower(i);
            QP.ub(states) = QP.ub(states) - StorPower(i);
            %update spinning reserve (max additional power from storage
            if isfield(QP.Organize,'SpinRow') && QP.Organize.SpinRow{i}~=0
                QP.beq(QP.Organize.SpinRow{i}) = sum(QP.ub(states(1:end-1))) - StorPower(i);
            end
        end
    end
end

%% update costs
H = diag(QP.H);%convert to vector
for i = 1:1:nG
    if ~isempty(QP.Organize.States{i})
        states = QP.Organize.States{i}; %states of this generator
        if isempty(strfind(Plant.Generator(i).Type,'Storage'))%all generators and utilities
            H(states) = H(states)*scaleCost(i)*dt;
            QP.f(states) = QP.f(states)*scaleCost(i)*dt;
        else
            if strcmp(Plant.Generator(i).Source,'Electricity')
                type = 'Electrical';
            elseif strcmp(Plant.Generator(i).Source,'Heat')
                type = 'DistrictHeat';
            elseif strcmp(Plant.Generator(i).Source,'Cooling')
                type = 'DistrictCool';
            elseif strcmp(Plant.Generator(i).Source,'Water')
                type = 'Hydro';
            end
            if ~isempty(StorPower)%remove expected storage output from beq
                QP.beq(QP.Organize.StorageEquality(i)) = QP.beq(QP.Organize.StorageEquality(i)) - StorPower(i);
            end
            %penalize deviations from the expected storage output
            %don't penalize spinning reserve state
            if isempty(EC)
                QP.f(states(1)) = marginal.(type);
                if any(strcmp(networkNames,'Electrical')) && strcmp(type,'H') % first initialization give arbitrarily high cost to storage (but not thermal storage if in conjunction with electric dispatch)
                    H(states(1)) = 0;
                else
                    H(states(1)) = 1e8*marginal.(type);
                end
            else
                a = 4; %sets severity of quadratic penalty
                MaxCharge =min((Plant.Generator(i).OpMatA.Stor.UsableSize-EC(i))/dt,Plant.Generator(i).OpMatB.Ramp.b(1)); %minimum of peak Charge, and SOC/time (completely charging storage in next step)
                QP.f(states(1)) = marginal.(type);
                if MaxCharge == 0
                    H(states(1)) = 0;
                else
                    H(states(1)) = a*2*QP.f(states(1))/MaxCharge;  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                end
            end
        end
    end
end 
QP.constCost = QP.constCost.*scaleCost*dt;
if Plant.optimoptions.SpinReserve
    SRshort = QP.Organize.SpinReserveStates(1,nG+1);%cumulative spinning reserve shortfall at t = 1 --> nS
    SRancillary = QP.Organize.SpinReserveStates(1,nG+2);%Ancillary spinning reserve value (negative cost) at t = 1 --> nS
    if Plant.optimoptions.SpinReservePerc>5 %more than 5% spinning reserve
        SpinCost = 2*dt./(Plant.optimoptions.SpinReservePerc/100*sum(Demand.(Outs{s}),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    else
        SpinCost = 2*dt./(0.05*sum(Demand.(Outs{s}),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    end
    H(SRshort) = SpinCost;%effectively $2 per kWh at point where shortfall = spin reserve percent*demand or $2 at 5%
    QP.f(SRshort) = 0.05*dt; % $0.05 per kWh
end
QP.H = diag(H);%convert back to matrix