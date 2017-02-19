function Locked = forceOffForFeasibility(IC,Locked,thisSeq,dX,load,SecondDisp,MaxCharge,PrevCharge,util,stor,dt,CHPindex,Hratio)
global UB LB OnOff
if ~isempty(CHPindex)
    LBglobal = LB;
    LB(CHPindex) = LB(CHPindex).*Hratio;
    thisSeq = [thisSeq CHPindex];
end
[nS,~] = size(Locked);
if ~isempty(util)
    GridMin = sum(LB(util));
    GridThresh = max(0,GridMin);
else GridMin = 0;
    GridThresh = 0;
end
if length(MaxCharge) ==1 %no storage, need to make into vector
    MaxCharge = zeros(nS,1);
    PrevCharge = zeros(nS,1);
end
if ~isempty(stor)
    Smax = sum(UB(stor));
    SminNow = sum(IC(stor));
end
%% %sort generators by their cost / frequency of appearance in Optimal Dispatch
nOn = sum((SecondDisp(:,thisSeq)>0),1);
for j = 1:1:length(thisSeq)
    i = thisSeq(j);
    if IC(i)>0
        stop1 = nonzeros((1:(nS-1))'.*((SecondDisp(1:(nS-1),i)>0)-(SecondDisp(2:nS,i)>0)));
        if ~isempty(stop1)
            nOn(j) = nOn(j) + 3*stop1(1);%give more weight if generator is already on & has been on
        else nOn(j) = nOn(j) + 3*nS;
        end
    end
end
[~,I] = sort(nOn);
offSeq = thisSeq(I);
%initialize
LBnow = 0*LB;
LBnow(thisSeq) = IC(thisSeq);
LBnet = sum(LBnow);
Starting = thisSeq((LBnow(thisSeq)<LB(thisSeq))&OnOff(thisSeq));
Ending = thisSeq((LBnow(thisSeq)<LB(thisSeq))&(LBnow(thisSeq)>0)&~OnOff(thisSeq));
PrevOff =thisSeq(LBnow(thisSeq)<LB(thisSeq));
% go through steps, and shut things down when necessary
for t = 2:1:nS
    forceOff=0;
    RampUp = dX(t-1,:);
    RampDown = dX(t-1,:);
    if ~isempty(util)
        gridDown = load(t-1)-LBnet-GridMin;
    else gridDown = 0;
    end
    LBold = LBnet;
    for j = 1:1:length(Ending) %if already ending when this fcn starts, make sure locked =0 while ending, othrwise locked = 1 while ending
        i = Ending(j);
        if ~Locked(t-1,i)
            Locked(t,i) = false;
        end
    end
    PrevOn = thisSeq(LBnow(thisSeq)>=LB(thisSeq));
    i = 1;
    while i<=length(PrevOn)
        j = PrevOn(i);
        if LBnow(j)>LB(j)
            a = min(LBnow(j) - LB(j),RampDown(j));
            RampDown(j) = RampDown(j) - a;
            LBnow(j) = LBnow(j) - a;
        end
        if ismember(j,PrevOff)
            PrevOn= PrevOn(PrevOn~=j); %make a generator be on for a full step before counting on (needs to ramp up to be able to ramp down)
            i = i-1;
        end
        i = i+1;
    end
    PrevOff =thisSeq(LBnow(thisSeq)==0);
    newOn = thisSeq(Locked(t,thisSeq)>Locked(t-1,thisSeq));
    newOff = thisSeq(Locked(t,thisSeq)<Locked(t-1,thisSeq));
    i = 1;
    while i<=length(Starting)
        if LBnow(Starting(i))<LB(Starting(i))
            PrevOff(end+1) = Starting(i);
            a = min(LB(Starting(i))-LBnow(Starting(i)),RampUp(Starting(i)));
            RampUp(Starting(i)) = RampUp(Starting(i)) -a;
            LBnow(Starting(i)) = min(LBnow(Starting(i))+a,LB(Starting(i)));
            i = i+1;
        else
            Starting = Starting(Starting~=Starting(i));
        end
    end
    i = 1;
    while i<=length(Ending)
        if LBnow(Ending(i))>0
            PrevOn(end+1) = Ending(i);
            RampDown(Ending(i)) = min(LBnow(Ending(i)),RampDown(Ending(i)));
            LBnow(Ending(i)) = LBnow(Ending(i))-RampDown(Ending(i));
            i = i+1;
        else
            Ending = Ending(Ending~=Ending(i));
        end
    end
    for i =1:1:length(newOn)
        LBnow(newOn(i)) = min(LB(newOn(i)),RampUp(newOn(i)));
    end
    for i =1:1:length(newOff)
        LBnow(newOff(i)) = max(0,LB(newOff(i))-RampDown(newOff(i)));
    end
    if ~isempty(stor)
        if (SminNow + (sum(LBnow(thisSeq)) - load(t-1))*dt(t-1))>Smax 
            %turn something off previously?
            forceOff=1;
        elseif (SminNow + (sum(LBnow(thisSeq)) - load(t-1))*dt(t-1))<0
            %turn something on
            %re-calculate SminNow
            % should this ever happen, if starting with all things on?
        end
    end
    LBnet = sum(LBnow(thisSeq).*Locked(t,thisSeq));
    if t>2
        dLoad = load(t-1)-load(t-2);
    else dLoad = load(t-1)-sum(SecondDisp(1,thisSeq))-GridThresh;
    end
    if (LBnet+GridThresh-MaxCharge(t-1))>(load(t-1)+PrevCharge(t-1))%the lowest possible produced is greater than the load
        forceOff=1;
    elseif LBnet>LBold && (LBnet-LBold-dLoad)>(sum(RampDown(PrevOn))+MaxCharge(t-1)-PrevCharge(t-1)+gridDown)  %New generators are trying to come online, and the current ones must accomodate, the change in LB is greater than ramp rates + dLoad + what utility can accomodate
        forceOff=1;
    end
%             %% Create an opposite rule to keep generators on? (can get extra power from grid usually)
%             if -(LBnet-LBold-dLoad)>(sum(RampUp(PrevOn))-MaxCharge(t-1)+PrevCharge(t-1))  %the change in LB is greater than ramp rates + dLoad
%                 forceOn=1;
%             end
    
    %% %while forceOff is still true and there are generators that can turned off, turn something else off
    while forceOff && any(Locked(t,thisSeq))
        n =[];
       k=1;
       if ~isempty(newOn) %first try turning off generators trying to come on.
           while isempty(n) && k<=length(thisSeq)
               if ismember(offSeq(k),newOn)
                   n=offSeq(k);
                   Locked(t,n)= false; %Don't turn on least commonly used generator
                   LBnow(n) = max(0,LBnow(n)-RampUp(n)); %undo ramping up
                   newOn = newOn(newOn~=n);
               end
               k = k+1;
           end
       end
       k=0;
       while isempty(n) && k<length(thisSeq) %now turn off generators early
           k=k+1;
           if Locked(t,offSeq(k))
                n=offSeq(k);
                q = t-1;
                dXsumflip = zeros(t-1,1);
                for m=1:1:t-1
                    dXsumflip(m) = sum(dX(q:t-1,n));
                    q = q-1;
                end
                p = sum(LBnow(n)>dXsumflip)+1;%ceil(LBnow(n)/RampDown(n));
%                 p = p+((LBnow(n)-dX(t-p,n))>0);
                if all(Locked(t-p+1:t,n))
                    Locked(t-p+1:t,n) = false; %pre-emptively ramp down
                    LBnow(n) = 0;
                else
                    Locked(t,n)= false; %turns off least frequently called upon generator
                    LBnow(n) = max(0,LBnow(n)-RampDown(n)); %start ramping down
                    Ending(end+1) = n;
                end
           end
       end
       LBnet = sum(LBnow);
       if (LBnet+GridThresh-MaxCharge(t-1))<=load(t-1) &&  (LBnet-LBold-dLoad)<=(sum(RampDown(PrevOn))+MaxCharge(t-1)-PrevCharge(t-1)+gridDown) && (isempty(stor) || (SminNow + (sum(LBnow(thisSeq)) - load(t-1))*dt(t-1))<Smax)
           %if you can now meet both the previous demand, and can ramp
           %down fast enough, then exit this loop
           forceOff=0;
       end
    end
    if ~isempty(stor)
        SminNow = SminNow + (sum(LBnow(thisSeq)) - load(t-1))*dt(t-1);
    end
    Starting = [Starting, newOn(Locked(t,newOn))];
end
if ~isempty(CHPindex)
    LB = LBglobal; %return electric gen to original LB setting
end