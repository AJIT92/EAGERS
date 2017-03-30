function [QP,Demand,Renewable] = updateMatrices(QP,IC,Time,scaleCost,marginCost,Demand,Renewable,EC)
% QP is the set of optimization matrices 
% IC is the intial condition
% Time is the vector of time steps
% scale cost is the matrix multiplier of cost for each generator at each time step
% MarginCost is the marginal cost of generation (used in the final value of the energy storage)
% EC is the end condition for the threshold optimization
global Plant

nG = length(Plant.Generator);
nS = length(Time);
dt = Time' - [0; Time(1:end-1)'];

networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));

%% update demands: both beq & b (heating w/ rejection);
Outs = fieldnames(Demand);
Outs = Outs(~strcmp(Outs,'T'));
for n = 1:1:length(Outs)
    if strcmp(Outs{n},'H') && Plant.optimoptions.excessHeat == 1
        for i = 1:1:length(Demand.H(1,:))
            ineq = QP.Organize.Demand.DistrictHeat{i};%this is the inequality constraint associated with this demand at t = 1:nS
            QP.b(ineq) = QP.b(ineq) + Demand.H(:,i); %multile demands can be at the same node
        end
    else
        for i = 1:1:length(Demand.(Outs{n})(1,:))
            if strcmp(Outs{n},'E')
                eq = QP.Organize.Demand.Electrical{i};%this is the inequality constraint associated with this demand at t = 1:nS
            elseif strcmp(Outs{n},'H')
                eq = QP.Organize.Demand.DistrictHeat{i};%this is the inequality constraint associated with this demand at t = 1:nS
            elseif strcmp(Outs{n},'C')
                eq = QP.Organize.Demand.DistrictCool{i};%this is the inequality constraint associated with this demand at t = 1:nS
            end
            QP.beq(eq) = QP.beq(eq) + Demand.(Outs{n})(:,i); %multile demands can be at the same node
        end
    end
end

% subtract renewable generation and storage self-discharge
for net = 1:1:length(networkNames)
    for n = 1:1:length(Plant.subNet.(networkNames{net}))
        equip = Plant.subNet.(networkNames{net})(n).Equipment;
        for j = 1:1:length(equip)
            i = equip(j);
            if strcmp(networkNames{net},'Electrical') && strcmp(Plant.Generator(i).Source,'Renewable')
                eq = QP.Organize.Balance.Electrical(n,:);
                QP.beq(eq) = QP.beq(eq) - Renewable(:,i); %put renewable generation into energy balance at correct node
            end
            if ~isempty(strfind(Plant.Generator(i).Type,'Storage'))
                stor = Plant.Generator(i).OpMatA.Stor;
                loss = dt*(stor.SelfDischarge*stor.UsableSize)/(stor.DischEff*stor.ChargeEff);
                if strcmp(Plant.Generator(i).Source,'Electricity') && strcmp(networkNames{net},'Electrical')
                    eq = QP.Organize.Balance.Electrical(n,:);
                    QP.beq(eq) = QP.beq(eq) - loss; %account for self-discharge losses
                elseif strcmp(Plant.Generator(i).Source,'Heat') && strcmp(networkNames{net},'DistrictHeat')
                    if Plant.optimoptions.excessHeat == 1
                        ineq = QP.Organize.Imbalance.DistrictHeat(n,:);
                        QP.b(ineq) = QP.b(ineq) - loss; %account for self-discharge losses
                    else
                        eq = QP.Organize.Balance.DistrictHeat(n,:);
                        QP.beq(eq) = QP.beq(eq) - loss; %account for self-discharge losses
                    end
                elseif strcmp(Plant.Generator(i).Source,'Cooling') && strcmp(networkNames{net},'DistrictCool')
                    eq = QP.Organize.Balance.DistrictCool(n,:);
                    QP.beq(eq) = QP.beq(eq) - loss; %account for self-discharge losses
                end
            end
        end
    end
end         

%update initial conditions
for i = 1:1:nG
    if QP.Organize.IC(i)>0 %note #1: IC of storage is already sclaed by discharge losses
        nIC = nnz(QP.Organize.IC(1:i)); %order in IC
        QP.beq(nIC) = IC(i);
    end
 end

%update costs
H = diag(QP.H);
for i = 1:1:nG
    if ~isempty(QP.Organize.States{i})
        states = QP.Organize.States{i};
        nt = length(states)/nS; %number of states per timestep
        if isempty(strfind(Plant.Generator(i).Type,'Storage'))%all generators and utilities
            for t = 1:1:nS
                Xn = states((t-1)*nt+1:t*nt);
                H(Xn) = H(Xn)*scaleCost(t,i); 
                QP.f(Xn) = QP.f(Xn)*scaleCost(t,i); 
            end          
        else % storage systems
             %% update storage costs
            Out = fieldnames(Plant.Generator(i).OpMatA.output);
            s_end = states(end)-nt+1; %final state of charge
            StorSize = Plant.Generator(i).OpMatA.X.ub;
            BuffSize = Plant.Generator(i).OpMatA.W.ub;
            if isempty(EC)
                if strcmp(Out{1},'H')
                    Max = 0.8*marginCost.H.Max;
                    Min = 0;
                elseif strcmp(Out{1},'E')
                    Max = 1.25*marginCost.E.Max;
                    Min = 0.95*marginCost.E.Min;
                elseif strcmp(Out{1},'C')
                    Max = 1.25*marginCost.C.Max;
                    Min = 0.65*marginCost.C.Min;
                end
                a1 = -Max; % fitting C = a1*SOC + a2*SOC^2 so that dC/dSOC @ 0 = -max & dC/dSOC @ UB = -min
                a2 = (Max - Min)/(2*StorSize);
                H(s_end) = 2*a2;%quadratic final value term loaded into SOC(t=nS)  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                QP.f(s_end) = a1;%linear final value term loaded into SOC(t=nS)
                if BuffSize>0
                    for t = 1:1:nS
                        Xn = states(t*nt-1:t*nt);
                        QP.f(Xn) = Min;%this is the linear buffer term loaded into the lower & upper buffer
                        H(Xn) = 2*(2*Max-Min)/(2*BuffSize);%this is the quadratic buffer term loaded into the lower & upper buffer
                    end
                end
            else %used in Online loop when there is a target EC determined by dispatch loop
                %final SOC deviation cost (quadratic cost for error = SOC(end) - EC
                rows = QP.Organize.StorageInequalities{i};
                nr = length(states)/nS;
                PeakChargePower = Plant.Generator(i).OpMatA.Ramp.b(1);
                dSOC_10perc = .1*PeakChargePower*Time(end); %energy (kWh) if charging at 10%
                H(s_end) = -2*marginCost.(Out{1}).Min/dSOC_10perc;%quadratic final value term loaded into SOC(t=nS)  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                QP.f(s_end) = -marginCost.(Out{1}).Min;%linear final value term loaded into SOC(t=nS)
                nIC = nnz(QP.Organize.IC(1:i)); %order in IC
                QP.beq(nIC) = IC(I)-EC(I);
                for t = 1:1:nS
                    Xn = states((t-1)*nt+1);
                    Rn = rows(t*nr);
                    QP.lb(Xn) = -EC(i);%change lb so that SOC = 0 coresponds to EC
                    QP.ub(Xn) = StorSize - EC(i);%change ub so that SOC = 0 coresponds to EC
                    QP.b(Rn-1) = -BuffSize + EC(i);%change lb so that SOC = 0 coresponds to EC (adding EC because there is a -1 in front of SOC in this inequality)
                    QP.b(Rn) = StorSize-BuffSize - EC(i);%change lb so that SOC = 0 coresponds to EC
                end       
            end
        end
    end
end     
QP.H = diag(H);