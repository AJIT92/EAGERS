function [QPall,Demand] = updateMatrices(QPall,Organize,IC,Time,scaleCost,marginCost,EC)
% IC is the intial condition
% fit refers to quadratic cost fit A or B
% Locked is a matrix nSx nG with zeros when a generator is locked off
% Stor is a variable to include storage in the optimization or not.
global Plant DateSim UB
nG = length(Plant.Generator);
nS = length(Time);

Demand = updateForecast(DateSim,Time);%% function that creates demand vector with time intervals coresponding to those selected
Demand = AccountForSelfDischarge(Demand);

Outs = fieldnames(QPall);    
for seq = 1:1:length(Outs)
    QP = QPall.(Outs{seq});% quadratic programing matrices (H, f, Aeq, beq, A, b, ub, lb)
    thisSeq = Organize.(Outs{seq}).thisSeq;
    stor = Organize.(Outs{seq}).stor;
    storC = Organize.(Outs{seq}).storC;
    storH = Organize.(Outs{seq}).storH;
    utility = Organize.(Outs{seq}).utility;
    utilC = Organize.(Outs{seq}).utilC;
    utilH = Organize.(Outs{seq}).utilH;
    chill = Organize.(Outs{seq}).chill;
    heater = Organize.(Outs{seq}).heater;
    allStor = [stor, storC, storH];
    allGen = [thisSeq, chill, heater];
    allUtility = [utility, utilC, utilH];

    %update demands: both beq & b (heating);
    for p = 1:1:length(Plant.optimoptions.Outputs) %reserve nS steps for each output type
        mat = Organize.(Plant.optimoptions.Outputs{p}).Demand{1};
        index = Organize.(Plant.optimoptions.Outputs{p}).Demand{2};
        if strcmp(mat,'b')
            QP.b(index) = -Demand.(Plant.optimoptions.Outputs{p});%Ax<beq, so b and A must be negative so you are producing more than enough
        else
            QP.beq(index) = Demand.(Plant.optimoptions.Outputs{p});
        end
    end

    %update initial condition
    for i = 1:1:nG
        if ismember(i,[allGen, allStor]) && Organize.IC(i)~=0
            QP.beq(Organize.IC(i)) = IC(i);
            if ~isempty(Organize.States{i})
                QP.ub(Organize.States{i}(1)) = IC(i);
            end
        end
        if ismember(i, allUtility)
            k = Organize.States{i}';
            if nnz(isinf(QP.ub(k)))>0
                p = fieldnames(Plant.Generator(i).OpMatA.output); %utilities with infinite upper bound
                QP.ub(k) = 100*max(sum(UB(~isinf(UB))),max(Demand.(p{1})));
            end
        end
    end

    %update costs
    H = diag(QP.H);
    for i = 1:1:nG
        if ismember(i,[allGen, allUtility])
            scaleC = [];
            k = Organize.States{i}';
            if ~isempty(k)
                if isfield(Plant.Generator(i).OpMatA,'Ramp') %has ramping constraint, thus it has an initial condition in order to enforce this
                    k = k(2:end);%%first state corresponds to IC
                end
                for j = 1:1:floor(length(k)/nS)
                    scaleC(end+1:end+nS,1) = scaleCost(:,i);
                end
                H(k) = H(k).*scaleC;
                QP.f(k) = QP.f(k).*scaleC;
            end
        elseif ismember(i,allStor)
            %% update storage costs
            %Storage States %SOC(t+1), charging power, upper buffer, lower buffer, self discharge 
            nX = Organize.States{i}; %first state corresponds to IC
            nSi = nS+1; %final state of charge
            StorSize = Plant.Generator(i).OpMatA.X.ub;
            BuffSize = Plant.Generator(i).OpMatA.W.ub;
            if length(nX) == (4*nS+1)
                XlowBuff = 2*nS+2:3*nS+1; %lower bound (W state) state #
                XhighBuff = 3*nS+2:4*nS+1; %lower bound (W state) state #
                IneqLowBuff = nS+1:2*nS; %lower bound (W state) inequality row
                IneqHighBuff = 2*nS+1:3*nS; %lower bound (W state) inequality row
            elseif length(nX) == (3*nS+1)
                XlowBuff = nS+2:2*nS+1;
                XhighBuff = 2*nS+2:3*nS+1;
                IneqLowBuff = 1:nS;
                IneqHighBuff = nS+1:2*nS;
            end
            if ~isempty(EC)
                %final SOC deviation cost (quadratic cost for error = SOC(end) - EC
                Out = fieldnames(marginCost);
                for j = 1:1:length(Out)
                    if isfield(Plant.Generator(i).OpMatA,'Stor') && isfield(Plant.Generator(i).OpMatA.output,Out{j})
                        %% may need to increase the scale of the quadratic penalty to keep closer to original soln
                        %penalize more than 10% deviation in power output
                        PeakChargePower = Plant.Generator(i).OpMatA.Ramp.b(1);
                        dSOC_10perc = .1*PeakChargePower*Time(end); %energy (kWh) if charging at 10%
                        H(nX(nSi)) = -2*marginCost.(Out{j}).Min/dSOC_10perc;%quadratic final value term loaded into SOC(t=nS)  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                        QP.f(nX(nSi)) = -marginCost.(Out{j}).Min;%linear final value term loaded into SOC(t=nS)
%                         QP.lb(nX(1:nSi)) = -EC(i);%change lb so that SOC = 0 coresponds to EC
%                         QP.ub(nX(1:nSi)) = StorSize - EC(i);%change ub so that SOC = 0 coresponds to EC
                        %need to find the inequality rows and adjust b term to account for this offset
                        r = Organize.Inequalities{i};
%                         QP.b(r(IneqLowBuff)) = -BuffSize + EC(i);%change lb so that SOC = 0 coresponds to EC (adding EC because there is a -1 in front of SOC in this inequality)
%                         QP.b(r(IneqHighBuff)) = StorSize-BuffSize - EC(i);%change lb so that SOC = 0 coresponds to EC
                        %need to change IC to IC-EC
%                         QP.beq(Organize.IC(i)) = IC(i)-EC(i);
                    end
                end
            else % dispatch optimization with final state & buffer states
                Out = fieldnames(Plant.Generator(i).OpMatA.output);
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
                H(nX(nSi)) = 2*a2;%quadratic final value term loaded into SOC(t=nS)  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                QP.f(nX(nSi)) = a1;%linear final value term loaded into SOC(t=nS)
                if BuffSize>0
                    QP.f(nX(XlowBuff)) = Min;%this is the linear buffer term loaded into the upper buffer
                    QP.f(nX(XhighBuff)) = Min;%this is the linear buffer term loaded into the lower buffer
                    H(nX(XlowBuff)) = 2*(2*Max-Min)/(2*BuffSize);%this is the quadratic buffer term loaded into the upper buffer
                    H(nX(XhighBuff)) = 2*(2*Max-Min)/(2*BuffSize);%this is the quadratic buffer term loaded into the lower buffer    
                end
            end
        end
    end
    QP.H = diag(H);
    QPall.(Outs{seq})= QP;% quadratic programing matrices (H, f, Aeq, beq, A, b, ub, lb) 
end