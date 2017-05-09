function [GenDisp,cost,Feasible] = DispatchQP(QP,Locked)
[m,n] = size(QP.organize);
nG = length(QP.Organize.IC);
nS = m-1;
QP = disableGenerators(QP,Locked,[]);%Disable generators here
if nnz(QP.H)==0
    options = optimset('Algorithm','interior-point','MaxIter',100,'Display','none');
    [GenSetting,cost,Feasible] = linprog(QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options); 
else
    options2 = optimset('Algorithm','interior-point-convex','MaxIter',100,'Display','none');
    [GenSetting,cost,Feasible] = quadprog(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options2);
end
GenDisp = zeros(nS+1,n);
if Feasible ~=1
    Feasible = false;
%     disp('Infeasible in DispatchQP');
else
    Feasible = true;
    for i = 1:1:nG
        if any(QP.Renewable(:,i)~=0)
            GenDisp(2:end,i) = QP.Renewable(:,i);
        else
            for t = 1:1:nS+1
                if ~isempty(QP.organize{t,i})
                    GenDisp(t,i) = sum(GenSetting(QP.organize{t,i}));%%put solution back into recognizable format
                end
            end
        end
    end
    for i = nG+1:1:n %Transmission lines
        for t = 2:1:nS+1
            GenDisp(t,i) = GenSetting(QP.organize{t,i});%%put solution back into recognizable format
        end
    end
    GenDisp(abs(GenDisp)<1e-1) = 0;
end
% % Check charging state
% global Plant
% for i = 1:1:nG
%     n = 0;
%     if ~isempty(QP.Organize.Inequalities(i))
%         n = n+1;
%         states = eval(QP.Organize.States(i));
%         chargeState(:,n) = zeros(nS,1);
%         for t = 1:1:nS
%             chargeState(t,1) = GenSetting(QP.organize{t+1,i}+1);%%charge state (if it exists) is adjacent to SOC
%         end
%         eff = Plant.Generator(i).OpMatA.Stor.DischEff*Plant.Generator(i).OpMatA.Stor.ChargeEff;
%         StorPower(:,n) = GenDisp(1:nS,i)-GenDisp(2:nS+1,i);
%         Charge = -StorPower(:,n).*(StorPower(:,n)<0);
%         correctChargePower = Charge*(1-eff);
%         error(:,n) = chargeState(:,n)- correctChargePower;
%     end
% end