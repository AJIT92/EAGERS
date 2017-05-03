function [QP,Demand,Renewable] = updateMatrices(QP,IC,Time,scaleCost,marginCost,Demand,Renewable,EC)
% QP is the set of optimization matrices 
% IC is the intial condition
% Time is the vector of time steps
% scale cost is the matrix multiplier of cost for each generator at each time step
% MarginCost is the marginal cost of generation (used in the final value of the energy storage)
% Demand is the forecasted loads
% Renewable is the forecasted uncontrollable generation
% EC is the end condition for the threshold optimization
global Plant
nG = length(Plant.Generator);
nS = length(Time);
dt = Time - [0; Time(1:end-1)];

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
        eq = QP.Organize.Balance.(networkNames{net})(i,:);
        loads = QP.Organize.Demand.(networkNames{net})(i);
        QP.beq(eq) = sum(Demand.(out)(:,loads),2); %multile demands can be at the same node
        for j = 1:1:length(equip)
            k = equip(j);
            if strcmp(networkNames{net},'Electrical') && strcmp(Plant.Generator(k).Source,'Renewable')% subtract renewable generation 
                QP.beq(eq) = QP.beq(eq) - Renewable(:,k); %put renewable generation into energy balance at correct node
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
            %Out = fieldnames(Plant.Generator(i).OpMatA.output);
            Out = Plant.optimoptions.Outputs;
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
                a1 = -Max; % fitting C = a1*SOC + 0.5*a2*SOC^2 so that dC/dSOC @ 0 = -max & dC/dSOC @ UB = -min
                a2 = (Max - Min)/(StorSize);
                H(s_end) = a2;%quadratic final value term loaded into SOC(t=nS)  %quadprog its solving C = 0.5*x'*H*x + f'*x
                QP.f(s_end) = a1;%linear final value term loaded into SOC(t=nS)
                if BuffSize>0
                    for t = 1:1:nS
                        Xn = states(t*nt-1:t*nt);
                        QP.f(Xn) = Min;%this is the linear buffer term loaded into the lower & upper buffer
                        H(Xn) = (2*Max-Min)/BuffSize;%this is the quadratic buffer term loaded into the lower & upper buffer
                    end
                end
            else %used in Online loop when there is a target EC determined by dispatch loop
                %final SOC deviation cost (quadratic cost for error = SOC(end) - EC
                rows = QP.Organize.Inequalities{i};
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
if Plant.optimoptions.SpinReserve
    SRshort = QP.Organize.SpinReserveStates(:,nG+1);%cumulative spinning reserve shortfall at t = 1 --> nS
    SRancillary = QP.Organize.SpinReserveStates(:,nG+2);%Ancillary spinning reserve value (negative cost) at t = 1 --> nS
    if Plant.optimoptions.SpinReservePerc>5 %more than 5% spinning reserve
        SpinCost = 2*dt./(Plant.optimoptions.SpinReservePerc/100*sum(Demand.(Outs{s}),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    else
        SpinCost = 2*dt./(0.05*sum(Demand.(Outs{s}),2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    end
    H(SRshort) = SpinCost;%effectively $2 per kWh at point where shortfall = spin reserve percent*demand or $2 at 5%
    QP.f(SRshort) = 0.05*dt; % $0.05 per kWh
end
QP.H = diag(H);
QP.Renewable = Renewable;