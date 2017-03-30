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
%% update demands: both beq & b (heating w/ rejection);
Outs = fieldnames(Demand);
for n = 1:1:length(Outs)
    if strcmp(Outs{n},'H') && Plant.optimoptions.excessHeat == 1
        for i = 1:1:length(Demand.H)
            ineq = QP.Organize.Demand.DistrictHeat(i);%this is the inequality constraint associated with this demand at t = 1:nS
            QP.b(ineq) = QP.b(ineq) + Demand.H(i); %multile demands can be at the same node
        end
    else
        for i = 1:1:length(Demand.(Outs{n}))
            if strcmp(Outs{n},'E')
                eq = QP.Organize.Demand.Electrical(i);%this is the inequality constraint associated with this demand at t = 1:nS
            elseif strcmp(Outs{n},'H')
                eq = QP.Organize.Demand.DistrictHeat(i);%this is the inequality constraint associated with this demand at t = 1:nS
            elseif strcmp(Outs{n},'C')
                eq = QP.Organize.Demand.DistrictCool(i);%this is the inequality constraint associated with this demand at t = 1:nS
            end
            QP.beq(eq) = QP.beq(eq) + Demand.(Outs{n})(i); %multile demands can be at the same node
        end
    end
end

% subtract renewable generation and storage self-discharge
for net = 1:1:length(networkNames)
    for n = 1:1:length(Plant.subNet.(networkNames{net}))
        equip = Plant.subNet.(networkNames{net})(n).Equipment;
        for j = 1:1:length(equip)
            i = equip(j);
            gen = Plant.Generator(i);
            if strcmp(networkNames{net},'Electrical') && strcmp(gen.Source,'Renewable')
                eq = QP.Organize.Balance.Electrical(n);
                QP.beq(eq) = QP.beq(eq) - Renewable(i); %put renewable generation into energy balance at correct node
            end
            if ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
                stor = Plant.Generator(i).OpMatA.Stor;
                loss = dt*(stor.SelfDischarge*stor.UsableSize)/(stor.DischEff*stor.ChargeEff);
                if strcmp(Plant.Generator(i).Source,'Electricity')
                    eq = QP.Organize.Balance.Electrical(n);
                    QP.beq(eq) = QP.beq(eq) - loss; %account for self-discharge losses
                elseif strcmp(Plant.Generator(i).Source,'Heat')
                    if Plant.optimoptions.excessHeat == 1
                        ineq = QP.Organize.Imbalance.DistrictHeat(n);
                        QP.b(ineq) = QP.b(ineq) - loss; %account for self-discharge losses
                    else
                        eq = QP.Organize.Balance.DistrictHeat(n);
                        QP.beq(eq) = QP.beq(eq) - loss; %account for self-discharge losses
                    end
                elseif strcmp(Plant.Generator(i).Source,'Cooling')
                    eq = QP.Organize.Balance.DistrictCool(n);
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
            if ~isempty(StorPower)%remove expected storage output from beq
                QP.beq(QP.Organize.StorageEquality(i)) = QP.beq(QP.Organize.StorageEquality(i)) - StorPower(i);
            end
            %penalize deviations from the expected storage output
            Out = fieldnames(Plant.Generator(i).OpMatB.output);
            if isempty(EC) && ~strcmp(Out{1},'H') % first initialization give arbitrarily high cost to storage (but not thermal storage if in conjunction with electric dispatch)
                QP.f(states) = marginal.(Out{1});
                H(states) = 1e8*marginal.(Out{1});
            elseif isempty(EC)%thermal storage 1st time
                QP.f(states) = marginal.(Out{1});
                H(states) = 0;
            else
                a = 4; %sets severity of quadratic penalty
                MaxCharge =min((Plant.Generator(i).OpMatA.Stor.UsableSize-EC(i))/dt,Plant.Generator(i).OpMatB.Ramp.b(1)); %minimum of peak Charge, and SOC/time (completely charging storage in next step)
                QP.f(states) = marginal.(Out{1});
                if MaxCharge == 0
                    H(states) = 0;
                else
                    H(states) = a*2*QP.f(states)/MaxCharge;  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                end
            end
        end
    end
    QP.constCost = QP.constCost.*scaleCost*dt;
end 
QP.H = diag(H);%convert back to matrix