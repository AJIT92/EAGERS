function [GenDisp, dX] = FilterGenerators(QP,Organize,IC,Demand,FirstDisp,SecondDisp,scaleCost,preLock)
global Plant OnOff UB LB
Time = buildTimeVector(Plant.optimoptions);
dt = Time' - [0 Time(1:end-1)]';
[nS,nG] = size(SecondDisp);

%% load parameters
dX = zeros(nS-1,nG);
dXtotal = zeros(nS-1,nG);
OnCost = zeros(nS,nG);
Input = zeros(nS,nG);
CostPerkW = zeros(nG,1);
for i = 1:1:nG
    if isfield(Plant.Generator(i).OpMatB,'Ramp') 
        dX(:,i) = Plant.Generator(i).OpMatB.Ramp.b(1)*dt;
        dXtotal(:,i) = Plant.Generator(i).OpMatB.Ramp.b(1)*Time;
    end
    if isfield(Plant.Generator(i).OpMatB,'constCost') %all of these cost terms need to be scaled later on
        OnCost(2:end,i) = Plant.Generator(i).OpMatB.constCost*scaleCost(:,i).*dt;
    end
    if ~isempty(Plant.Generator(i).Output)
        cap = Plant.Generator(i).Output.Capacity*UB(i);
    end
    eff = [];
    if strcmp(Plant.Generator(i).Type,'Electric Generator') || strcmp(Plant.Generator(i).Type,'CHP Generator')
        eff = Plant.Generator(i).Output.Electricity;
    elseif strcmp(Plant.Generator(i).Type,'Chiller')
        eff = Plant.Generator(i).Output.Cooling;
    elseif strcmp(Plant.Generator(i).Type,'Heater')
        eff = Plant.Generator(i).Output.Heat;    
    end
    if ~isempty(eff)
        Input(:,i) = FirstDisp(:,i)./interp1(cap,eff,FirstDisp(:,i));
        Input(isnan(Input(:,i)),i)=0;
        CostPerkW(i) = sum(Input(2:end,i).*scaleCost(:,i))/sum(FirstDisp(2:end,i));
    end
end


cheapest=false;
while ~cheapest
    Locked = preLock;
    OriginSoln = true(nS,nG);
    MaxCharge = zeros(nS-1,nG);
    PrevCharge = zeros(nS-1,nG);
    Outs = fieldnames(QP);
    %% Identify relevant generators & sort by cost
    for seq=1:1:length(Outs)
        thisSeq = Organize.(Outs{seq}).thisSeq;
        stor = Organize.(Outs{seq}).stor;
        storC = Organize.(Outs{seq}).storC;
        storH = Organize.(Outs{seq}).storH;
        utility = Organize.(Outs{seq}).utility;
        utilC = Organize.(Outs{seq}).utilC;
        utilH = Organize.(Outs{seq}).utilH;
        chill = Organize.(Outs{seq}).chill;
        heater = Organize.(Outs{seq}).heater;
        Hratio = Organize.(Outs{seq}).Hratio;
        CHPindex = Organize.(Outs{seq}).CHPindex;
        allutil = [utility,utilC,utilH];

        [~,index] = sort(CostPerkW(thisSeq),1,'descend');
        thisSeq = thisSeq(index); %sorted from largest constant cost to smallest (general goal is to turn off more expensive ones first
        if~isempty(chill)
            [~,index] = sort(CostPerkW(chill),1,'descend');
            chill = chill(index); %sorted from largest constant cost to smallest (general goal is to turn off more expensive ones first
        end
        if~isempty(heater)
            [~,index] = sort(CostPerkW(heater),1,'descend');
            heater = heater(index); %sorted from largest constant cost to smallest (general goal is to turn off more expensive ones first
        end
        allGen = [thisSeq, chill, heater];
        allStor = [stor, storC, storH];
        OriginSoln(:,allGen) = (FirstDisp(:,allGen)~=0);

        for j = 1:1:length(allStor)
            i = allStor(j);
            MaxCharge(:,i) = max((UB(i)-FirstDisp(2:end,i))./dt,Plant.Generator(i).OpMatB.Ramp.b(1));
            PrevCharge(:,i) = (FirstDisp(2:end,i)-FirstDisp(1:end-1,i))./dt;
        end
    
        %% Rule 2 while 2nd dispatch corresponds to IC, extend until first time it needs to turn on (ramp rate to optimal setting the first time its on)   
        for i = 1:1:nG
            if ~OnOff(i) && ~ismember(i,allutil) && IC(i)==0%if it is set off and is already off, lock it off for the initial condition, don't lock off if there is an initial condition
                Locked(1,i) = false;
%             elseif ~OnOff(i) %if it is set off but is still ramping down, set it off at the soonest possible step
%                 t_zeroOutput = sum((IC(1,i)./dXtotal(:,i))>=1) + 1;%the number of steps until the generator reaches 0 
%                 Locked(t_zeroOutput+1,i) = false;
            end
        end
        %% check for feasibiltiy
        [GenDisp,Cost,Feasible] = DispatchQP(QP,Organize,Locked);
        if Feasible~=1
            %% Turn off enough to make the problem feasible with lower bound constraints (start with least used generators)
            %if there is a utility with inf UB and -inf LB, don't need to do one of these
            Locked = forceOffForFeasibility(IC,Locked,thisSeq,dX,Demand.(Outs{seq}),SecondDisp,sum(MaxCharge(:,stor),2),sum(PrevCharge(:,stor),2),utility,stor,dt,[],[]);
            if~isempty(chill)
                Locked = forceOffForFeasibility(IC,Locked,chill,dX,Demand.C,SecondDisp,sum(MaxCharge(:,storC),2),sum(PrevCharge(:,storC),2),utilC,storC,dt,[],[]);
            end
            if (~isempty(heater)||~isempty(CHPindex)) && Plant.optimoptions.excessHeat==0 %if you can reject heat, don't worry about this
                Locked = forceOffForFeasibility(IC,Locked,heater,dX,Demand.H,SecondDisp,sum(MaxCharge(:,storH),2),sum(PrevCharge(:,storH),2),utilH,storH,dt,CHPindex,Hratio);
            end
        end
    end
    % re-check for feasibiltiy
    [GenDisp,Cost,Feasible] = DispatchQP(QP,Organize,Locked);
    if Feasible~=1
        L2 = Locked;
        for j = 1:1:length(allGen)
            i = allGen(j);
            L2(2:end,i) = (FirstDisp(2:end,i)>=LB(i))&L2(2:end,i);%don't change initial condition
        end
        % re-check for feasibiltiy
        [GenDisp,Cost,Feasible] = DispatchQP(QP,Organize,L2);
        if Feasible==1
            Locked = L2;
        else
            disp('Cant find a feasible solution with lower bounds enforced (FilterGenerators line 114)')
    %                 uiwait(msgbox('Cant find a feasible solution with lower bounds enforced (FilterGenerators line 115)','Error','modal'))
            GenDisp = SecondDisp;
            cheapest = true;
        end
    end
    if Feasible==1
    %% Proceed with rules
        Cost = Cost + sum(sum(Locked.*OnCost));
        index = (1:nS-1)';
        for j = 1:1:length(allGen)
            %% Find when starts and stops occur in optimal dispatch
            i = allGen(j);
            starts = nonzeros(index.*(((SecondDisp(2:end,i)>0)-(SecondDisp(1:nS-1,i)>0))>0)); % if on at t = 3, start = 3, where IC is t=0
            if ~OnOff(i)&& SecondDisp(1,i)>0
                starts = [1;starts];
            end
            stops = nonzeros(index.*(((SecondDisp(1:nS-1,i)>0)-(SecondDisp(2:end,i)>0))>0)); % if off at t = 12, stop = 12, where IC is t=0
            if ~Locked(1,i) || ~OnOff(i)
                %% Rule 3: turn off for longer time at start if possible
                if~isempty(starts)
                    p = 0;
                    RampUp = dX(starts(1),i);
                    while RampUp<SecondDisp(starts(1)+1,i) && (starts(1)-p>0)
                        RampUp = RampUp+dX(starts(1)-p,i);
                        p = p+1;
                    end
                    if starts(1)-p>0
                        L2 = Locked;
                        L2(1:(starts(1)-p+1),i) = false;
                        [~,newCost,Feasible] = DispatchQP(QP,Organize,L2);
                        newCost = newCost + sum(sum(L2.*OnCost));
                        if Feasible==1 && newCost<Cost
                            Locked = L2;
                            Cost = newCost;
                        end
                    end
                    if length(starts)>1
                        starts = starts(2:end);
                    else starts =[];
                    end
                end
            end
            %% Rule 4: If off for a long enough segment in optimal dispatch, try turning off for as much of that as possible given ramp rates
            for k = 1:1:length(starts)
                if ~isempty(stops) && length(stops)>=k
                    if sum(dX(stops(k):starts(k),i))>(SecondDisp(stops(k),i) + SecondDisp(starts(k)+1,i) + 2*LB(i)) %can ramp all the way down, and back up, and a little more
                        L2 = Locked;
                        %find step when it can hit zero given setting at Disp(stops(k)
                        n=1;
                        RampDown = dX(stops(k),i);
                        while RampDown<SecondDisp(stops(k),i)
                            RampDown = RampDown+dX(stops(k)+n,i);
                            n = n+1;
                        end
                        p = 1;
                        RampUp = dX(starts(k),i);
                        while RampUp<SecondDisp(starts(k)+1,i)
                            RampUp = RampUp+dX(starts(k)-p,i);
                            p = p+1;
                        end
                        L2((stops(k)+n):(starts(k)-p+1),i) = false;
                        [~,newCost,Feasible] = DispatchQP(QP,Organize,L2);
                        newCost = newCost + sum(sum(L2.*OnCost));
                        if Feasible==1 && newCost<Cost
                            Locked = L2;
                            Cost = newCost;
                        end
                    end
                end
            end
            %% Rule 5: try turning off @ end
            if length(stops)>length(starts)
                n=stops(end);
                RampDown = dX(n,i);
                while RampDown<SecondDisp(stops(end),i) && n<(nS-1) && n>0
                    RampDown = RampDown+dX(n,i);
                    n = n-1;
                end
                if n<(nS-1)
                    L2 = Locked;
                    L2((n+1):nS,i) = false;
                    [~,newCost,Feasible] = DispatchQP(QP,Organize,L2);
                    newCost = newCost + sum(sum(L2.*OnCost));
                    if Feasible==1 && newCost<Cost
                        Locked = L2;
                        Cost = newCost;
                    end
                end
            end
        end

        for j = 1:1:length(allGen)
            %% Find when starts and stops occur in current Locked matrix
            i = allGen(j);
            starts = nonzeros(index.*((Locked(2:end,i)-Locked(1:nS-1,i))>0)); % if on at t = 3, start = 3, where IC is t=0
            stops = nonzeros(index.*((Locked(1:nS-1,i)-Locked(2:end,i))>0)); % if off at t = 12, stop = 12, where IC is t=0
            
            %% Rule 6: if on for a short time and sufficient capacity in remaining active generators, turn gen off completely
            if Locked(1,i)
                if length(stops)>1
                    stops = stops(2:end);
                else stops =[];
                end
            end
            for k = 1:1:length(stops)
                if sum(dX(starts(k):stops(k),i))>(UB(i) + LB(i)) && sum(Locked(:,i))<floor(nS/4)%can only ramp to 1/2 power and less than 1/4 of horizon
                    L2 = Locked;
                    L2(starts(k):(stops(k)+1),i)= false;
                    [~,newCost,Feasible] = DispatchQP(QP,Organize,L2);
                    newCost = newCost + sum(sum(L2.*OnCost));
                    if Feasible==1 && newCost<Cost
                        Locked = L2;
                        Cost = newCost;
                    end
                end
            end
            
        end
        %% prevent flicker
        for j = 1:1:length(allGen)
            i = allGen(j);
            starts = nonzeros(index.*((Locked(2:end,i)-Locked(1:nS-1,i))>0));
            stops = nonzeros(index.*((Locked(1:nS-1,i)-Locked(2:end,i))>0));
            if ~isempty(starts) && ~isempty(stops)
                stops = stops(stops<starts(end));
                Cost = Cost + sum(OnCost(starts,i)./dt(starts)).*5;
            end
            if ~isempty(stops) && ~isempty(starts)
                starts = starts(starts>stops(1));
                for k = 1:1:length(starts)
                    if starts(k)==stops(k)+1 %if you are starting immediately after stopping
                        L2 = Locked;
                        L2(stops(k)+1,i) = true;
                        [~,newCost,Feasible] = DispatchQP(QP,Organize,L2);
                        newCost = newCost + sum(sum(L2.*OnCost)) + sum(OnCost(starts,i)./dt(starts)).*5;%OnCost(i).*4 is the cost of shutting down and back on
                        if Feasible==1 && newCost<Cost
                            Locked = L2;
                            Cost = newCost;
                        end
                    end
                end
            end
        end
%         %% Rule 7: If very close to LB for a few steps, try shutting off
%         for j = 1:1:length(allGen)
%             i = allGen(j);
%             if Locked(1,i)==0 && max(Locked(:,i))>0%only mess with generators that are not already on, but turning on
%                 r = find(Locked(:,i)>0,1,'first'); %find first time generator is on 
%                 z=1;
%                 while r+z<nS && MarginalDispatch(r+z,i)<2*LB(i) && MarginalDispatch(r+z,i)<(UB(i)-LB(i)) %if it comes on-line and stays at low output
%                     z=z+1;
%                 end
%                 if z>2
%                     t=r;
%                     Feasible =1;
%                     L2 = Locked;
%                     while Feasible==1 && t<(r+z)
%                         L2(t,i)=0;
%                         [~,newCost,Feasible] = DispatchQP(QP,Organize,L2);
%                         newCost = newCost + sum(sum(L2.*OnCost));
%                         if Feasible==1 && newCost<Cost
%                             Locked = L2;
%                             Cost = newCost;
%                         else Feasible =0;
%                         end
%                         t=t+1;
%                     end
%                 end
%                 %% Rule 8: If can shut off for a few more steps, try shutting off
%                 r = find(Locked(:,i)>0,1,'first'); %find first time generator is on 
%                 Feasible =1;
%                 while Feasible==1 && ~isempty(r)
%                     L2 = Locked;
%                     L2(r,i)=0;
%                     [~,newCost,Feasible] = DispatchQP(QP,Organize,L2);
%                     newCost = newCost + sum(sum(L2.*OnCost));
%                     if Feasible==1 && newCost<Cost
%                         Locked = L2;
%                         Cost = newCost;
%                         r = find(Locked(:,i)>0,1,'first'); %find first time generator is on 
%                     else Feasible =0;
%                     end
%                 end
%                 r = find(Locked(:,i)>0,1,'last'); %find last time generator is on 
%                 Feasible =1;
%                 while Feasible==1 && ~isempty(r)
%                     L2 = Locked;
%                     L2(r,i)=0;
%                     [~,newCost,Feasible] = DispatchQP(QP,Organize,L2);
%                     newCost = newCost + sum(sum(L2.*OnCost));
%                     if Feasible==1 && newCost<Cost
%                         Locked = L2;
%                         Cost = newCost;
%                         r = find(Locked(:,i)>0,1,'last'); %find last time generator is on 
%                     else Feasible =0;
%                     end
%                 end
%             end
%         end
%         %% Rule 9: if it is switching generators, try switching 1 step earlier or later
%         %need to write this rule


        [GenDisp,Cost,~] = DispatchQP(QP,Organize,Locked);
        Cost = Cost + sum(sum(Locked.*OnCost));
        cheapest = true;

        %% check against original dispatch (currently thinks it saves money, but doesn't)
%             [GenDisp,originCost,Feasible] = DispatchQP(QP,Organize,OriginSoln>0);
%             if Feasible==1
%                 originCost = originCost + sum(sum(OriginSoln.*OnCost));
%             end
%             if ~(Feasible==1)||originCost>Cost
%                 GenDisp = DispatchQP(QP,Organize,Locked);
%                 cheapest = true; %if the new solution is cheaper, then stick with the new one
%             elseif nnz(Disp == MarginalDispatch)>0
%                 Disp = MarginalDispatch;
%             else cheapest = true;
%             end
    end
end