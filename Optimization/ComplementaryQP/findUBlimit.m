function UBnow = findUBlimit(Locked,IC,dX,MaxDischarge,thisSeq,util,stor)
global UB
[nS,nG] = size(Locked);
UBnow = zeros(nS,nG);
UBnow(1,thisSeq) = IC(thisSeq);
UBnow(:,util) = UB(util);
UBnow(:,stor) = MaxDischarge(stor);
for j = 1:1:length(thisSeq)
    i = thisSeq(j);
    if ~any(Locked(2:nS,i))
        UBnow(2:nS) = 0;
    elseif all(Locked(2:nS,i))
        for t = 2:1:nS
            RampUp = dX(t-1,i);
            UBnow(t,i) = min(UB(i),UBnow(i)+RampUp);
        end
    else %at least 1 start or stop
        index = (1:nS-1)';
        starts = nonzeros(index.*((Locked(2:end,i)-Locked(1:nS-1,i))>0)); % if on at t = 3, start = 3, where IC is t=0
        stops = nonzeros(index.*((Locked(1:nS-1,i)-Locked(2:end,i))>0)); % if off at t = 12, stop = 12, where IC is t=0
        if isempty(starts) || starts(end)<stops(end) %end condiion is off
            starts(end+1) = nS+1;
            stops(end+1) = nS+2;
        elseif isempty(stops) || stops(end)<starts(end)%end condition is on
            stops(end+1) = nS+100; 
            starts(end+1) = nS+200;
        end
        for t = 2:1:nS
            RampUp = dX(t-1,i);
            if t== stops(1)
                UBnow(t,i) = 0;
                stops = stops(2:end);
            elseif t==starts(1)
                UBnow(t,i) = min(UB(i),RampUp);
                starts = starts(2:end);
            elseif starts(1)<stops(1)
                UBnow(t,i) = 0;
            else
                UBnow(t,i) = min(UB(i),UBnow(i)+RampUp,sum(dX(t+1:stop(1),i)));
            end
        end
    end
end