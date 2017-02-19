function QP = updateMatrices1Step(QP,scaleCost,marginal,IC,dt,Organize)
global Plant UB 
nG = length(UB);
allStor = [Organize.stor, Organize.storC, Organize.storH];
allGen = [Organize.thisSeq, Organize.chill, Organize.heater, Organize.utility, Organize.utilC, Organize.utilH];

%% update costs
H = diag(QP.H);%convert to vector
for i = 1:1:nG
    k = Organize.States{i}; %states of this generator
    if ismember(i,allGen)
        H(k) = H(k)*scaleCost(i);
        QP.f(k) = QP.f(k)*scaleCost(i);
    elseif ismember(i,allStor)
        Output = fieldnames(Plant.Generator(i).OpMatB.output);
        if isempty(IC) && ismember(i,Organize.stor) % first initialization give arbitrarily high cost to storage (but not thermal storage if in conjunction with electric dispatch)
            QP.f(k) = marginal.(Output{1});
            H(k) = 1e8*marginal.(Output{1});
        elseif isempty(IC)%thermal storage 1st time
            QP.f(k) = marginal.(Output{1});
            H(k) = 0;
        else
            MaxCharge =min((UB(i)-IC(i))/dt,Plant.Generator(i).OpMatB.Ramp.b(1)); %minimum of peak Charge, and SOC/time (completely charging storage in next step)
            QP.f(k) = marginal.(Output{1});
            if MaxCharge == 0
                H(k) = 0;
            else
                H(k) = 2*QP.f(k)/MaxCharge;  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
            end
        end
    end
end
QP.H = diag(H);%convert back to matrix