function [GenDisp,cost,Feasible] = DispatchQP(QPall,Organize,Locked)
global Plant %chargeEff dischEff
Outs = fieldnames(QPall);
options = optimset('Algorithm','interior-point','MaxIter',100,'Display','none');
options2 = optimset('Algorithm','interior-point-convex','MaxIter',100,'Display','none');
[m,nG] = size(QPall.(Outs{1}).organize);
nS = m-1;
if ismember('C',Outs) && ismember('E',Outs)
    QP = QPall.C; %do chillers first
    Disp.C = zeros(nS+1,nG);
    Enabled = zeros(nG,1);
    
    thisSeq = Organize.C.thisSeq;
    stor = Organize.C.stor;
    storC = Organize.C.storC;
    storH = Organize.C.storH;
    utility = Organize.C.utility;
    utilC = Organize.C.utilC;
    utilH = Organize.C.utilH;
    chill = Organize.C.chill;
    heater = Organize.C.heater;
    allStor = [stor, storC, storH];
    allGen = [thisSeq, chill, heater];
    allUtility = [utility, utilC, utilH];
    for i = 1:1:length(allUtility)
        if Plant.Generator(allUtility(i)).Enabled
            Enabled(allUtility(i)) =1;
        end
    end
    for i = 1:1:length(allStor)
        if Plant.Generator(allStor(i)).Enabled
            Enabled(allStor(i)) =1;
        end
    end
    for i = 1:1:length(allGen)
        if Plant.Generator(allGen(i)).Enabled && (~isempty(Locked)&& any(Locked(:,allGen(i))))
            Enabled(allGen(i)) =1;
        end
    end
    QP = disableGenerators(QP,Organize,Locked,Enabled);%Disable generators here
    if nnz(QP.H)==0
        [GenSetting,cost,Feasible] = linprog(QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options); 
    else [GenSetting,cost,Feasible] = quadprog(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options2); 
    end
    for i = 1:1:nG
        for t = 1:1:nS+1
            Disp.C(t,i) = sum(GenSetting(QP.organize{t,i}));%%put solution back into recognizable format
        end
    end
    Disp.C(abs(Disp.C)<1e-2) = 0;
    Eload =0;%% add chilling demand to electric
    for i = 1:1:nG
        if nnz(Disp.C(:,i))>0
            if strcmp(Plant.Generator(i).Type,'Chiller') && strcmp(Plant.Generator(i).Source,'Electricity')% electric chillers
                eff = Plant.Generator(i).Output.Cooling;
                cap = Plant.Generator(i).Output.Capacity*Plant.Generator(i).Size;
                Eload = Eload + Disp.C(:,i)./interp1(cap,eff,Disp.C(:,i));
            end
        end
    end
    QP.E.beq(1:length(Eload)) = QP.E.beq+Eload; %add load to electric demand
    Outs = Outs(~(strcmp('C',Outs)));
end
for j = 1:1:length(Outs)
    QP = QPall.(Outs{j});
    Enabled = zeros(nG,1);
    Disp.(Outs{j}) = zeros(nS+1,nG);
    
    thisSeq = Organize.(Outs{j}).thisSeq;
    stor = Organize.(Outs{j}).stor;
%     storC = Organize.(Outs{j}).storC;
%     storH = Organize.(Outs{j}).storH;
    utility = Organize.(Outs{j}).utility;
%     utilC = Organize.(Outs{j}).utilC;
%     utilH = Organize.(Outs{j}).utilH;
%     chill = Organize.(Outs{j}).chill;
%     heater = Organize.(Outs{j}).heater;
    allStor = [stor];%, storC, storH];
    allGen = [thisSeq];%, chill, heater];
    allUtility = [utility];%, utilC, utilH];
%     stor = ~ismember([1:1:nG],[allStor,allGen,allUtility]);
%     Enabled(stor) = 1;
    for i = 1:1:length(allUtility)
        if Plant.Generator(allUtility(i)).Enabled
            Enabled(allUtility(i)) =1;
        end
    end
    for i = 1:1:length(allStor)
        if Plant.Generator(allStor(i)).Enabled
            Enabled(allStor(i)) =1;
        end
    end
    for i = 1:1:length(allGen)
        if Plant.Generator(allGen(i)).Enabled && (isempty(Locked) || any(Locked(:,allGen(i))))
            Enabled(allGen(i)) =1;
        end
    end
    QP = disableGenerators(QP,Organize,Locked,Enabled);%Disable generators here
    if nnz(QP.H)==0
        [GenSetting,cost,Feasible] = linprog(QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options); 
    else [GenSetting,cost,Feasible] = quadprog(QP.H,QP.f,QP.A,QP.b,QP.Aeq,QP.beq,QP.lb,QP.ub,[],options2); 
    end
    if Feasible==1
        for i = 1:1:nG
            for t = 1:1:nS+1
                Disp.(Outs{j})(t,i) = sum(GenSetting(QP.organize{t,i}));%%put solution back into recognizable format
            end
        end
    end
    Disp.(Outs{j})(abs(Disp.(Outs{j}))<1e-1) = 0;
        %% Check charging state
%         if isfield(Organize,'E') && isfield(Organize.E,'Stor') && Organize.E.Stor(i)==1
%             chargeState = zeros(nS,1);
%             for t = 1:1:nS
%                 chargeState(t,1) = GenSetting(QP.organize{t+1,i}+nS);%%put solution back into recognizable format
%             end
%             StorPower = GenDisp(1:nS,i)-GenDisp(2:nS+1,i);
%             Charge = -StorPower.*(StorPower<0);
%             correctChargePower = Charge*(1-(chargeEff(i)*dischEff(i)));
%             error = chargeState- correctChargePower;
%         end
end
Outs = fieldnames(QPall);
GenDisp = zeros(nS+1,nG);
for j = 1:1:length(Outs)
    for i = 1:1:nG
        if nnz(Disp.(Outs{j})(:,i))>0
            GenDisp(:,i) = Disp.(Outs{j})(:,i);
        end
    end
end