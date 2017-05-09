function QP = updateMatrices(QP,IC,Date,Time,scaleCost,marginCost,Demand,EC)
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
QP.Renewable = Demand.Renewable;
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
for net = 1:1:length(networkNames)
    nLinet(net) = length(Plant.subNet.lineNames.(networkNames{net}));
end
nL = sum(nLinet);
%% update demands and storage self-discharge
nLcum = 0; %cumulative line #
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
    elseif strcmp(networkNames{net},'Hydro')
        out = 'W';
    end
    for i = 1:1:length(Plant.subNet.(networkNames{net})) %run through all the nodes in this network
        equip = Plant.subNet.(networkNames{net})(i).Equipment; %equipment at this node
        eq = QP.Organize.Balance.(networkNames{net})(i,:);%balance at this node
        loads = QP.Organize.Demand.(networkNames{net})(i); %load at this node
        QP.beq(eq) = sum(Demand.(out)(:,loads),2); %multiple demands can be at the same node, or none
        for j = 1:1:length(equip)
            k = equip(j);
            if strcmp(networkNames{net},'Electrical') && any(QP.Renewable(:,k))% subtract renewable generation 
                QP.beq(eq) = QP.beq(eq) - QP.Renewable(:,k); %put renewable generation into energy balance at correct node
            end
            if ~isempty(strfind(Plant.Generator(k).Type,'Storage'))
                if (strcmp(Plant.Generator(k).Source,'Electricity') && strcmp(networkNames{net},'Electrical')) || (strcmp(Plant.Generator(k).Source,'Heat') && strcmp(networkNames{net},'DistrictHeat')) || (strcmp(Plant.Generator(k).Source,'Cooling') && strcmp(networkNames{net},'DistrictCool')) || (strcmp(Plant.Generator(k).Source,'Water') && strcmp(networkNames{net},'Hydro'))
                    loss = dt*(Plant.Generator(k).OpMatA.Stor.SelfDischarge*Plant.Generator(k).OpMatA.Stor.UsableSize);
                    QP.beq(eq) = QP.beq(eq) - loss; %account for self-discharge losses
                end
            end
        end
    end
    if strcmp(networkNames{net},'Hydro')
        downriver = {};
        downLines = [];
        for i = 1:1:length(Plant.subNet.lineNames.Hydro)
            name = Plant.subNet.lineNames.Hydro{i};
            k = strfind(name,'_');
            if strcmp(name(k(1)+1:k(2)-1),'Hydro')
                downriver(end+1) = {name(k(2)+1:end)};% node names of the downriver node
                downLines(end+1) = i;
            end
        end
        for i = 1:1:length(Plant.subNet.Hydro) %run through all the nodes in this network
            %% Node inflows (source/sink terms and upstream flow at time t-T ago if t<T)
            %need to be able to forecast sink/source
            eq = QP.Organize.Balance.Hydro(i,:);%mass balance at this node
            SourceSink = getHydroSourceSink(Date+Time,node);
            QP.beq(eq)= QP.beq(eq) - SourceSink./(12.1*dt); %source/sink term converted from 1000 ft^3/s to acre-ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
            K = downLines(strcmp(Plant.subNet.Hydro(i).nodes,downriver));%lines entering this node, i.e. this node is the downriver node
            for j = 1:1:length(K)
                Date2 = Date+Time-Plant.subNet.lineTime.Hydro(K(j));
                n = nnz(Date2<Date);
                InFlow = getHydroFlows(Date2(1:n),K(j));%river segments flowing into this node
                QP.beq(eq(1:n))= QP.beq(eq(1:n)) - InFlow./(12.1*dt(1:n)); % Qupriver, conversion factor is from 1000 ft^3/s to 1000 acre ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
            end
        end
    end
    nLcum = nLcum + nLinet(net);
end         

%update initial conditions
for i = 1:1:nG+nL
    if QP.Organize.IC(i)>0 
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
            if strcmp(Plant.Generator(i).Source,'Electricity')
                type = 'Electrical';
            elseif strcmp(Plant.Generator(i).Source,'Heat')
                type = 'DistrictHeat';
            elseif strcmp(Plant.Generator(i).Source,'Cooling')
                type = 'DistrictCool';
            elseif strcmp(Plant.Generator(i).Source,'Water')
                type = 'Hydro';
            end
             %% update storage costs
            s_end = states(end)-nt+1; %final state of charge
            StorSize = Plant.Generator(i).OpMatA.X.ub;
            BuffSize = Plant.Generator(i).OpMatA.W.ub;
            if isempty(EC)
                if strcmp(Plant.Generator(i).Source,'Heat')
                    Max = 0.8*marginCost.DistrictHeat.Max;
                    Min = 0;
                elseif strcmp(Plant.Generator(i).Source,'Electricity')
                    Max = 1.25*marginCost.Electrical.Max;
                    Min = 0.95*marginCost.Electrical.Min;
                elseif strcmp(Plant.Generator(i).Source,'Cooling')
                    Max = 1.25*marginCost.DistrictCool.Max;
                    Min = 0.65*marginCost.DistrictCool.Min;
                elseif strcmp(Plant.Generator(i).Source,'Water')
                    Max = 1.25*marginCost.Hydro.Max;
                    Min = 0.65*marginCost.Hydro.Min;
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
                H(s_end) = -2*marginCost.(type).Min/dSOC_10perc;%quadratic final value term loaded into SOC(t=nS)  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                QP.f(s_end) = -marginCost.(type).Min;%linear final value term loaded into SOC(t=nS)
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
        SpinCost = 2*dt./(Plant.optimoptions.SpinReservePerc/100*sum(Demand.E,2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    else
        SpinCost = 2*dt./(0.05*sum(Demand.E,2));% -shortfall + SRancillary - SR generators - SR storage <= -SR target
    end
    H(SRshort) = SpinCost;%effectively $2 per kWh at point where shortfall = spin reserve percent*demand or $2 at 5%
    QP.f(SRshort) = 0.05*dt; % $0.05 per kWh
end
QP.H = diag(H);
