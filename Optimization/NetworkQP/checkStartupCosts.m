function [GenOutput,Binary] = checkStartupCosts(Alt,Binary,StartCost,dt,I,nL)
%This function finds the n shortest segments of time the generator is on or
%off, and compares the start-up cost to the alternative configuration that
%avoids the start or re-start
n = 4; % # of segments to check
nG = length(StartCost);
nS = length(Alt.Disp);
SkipOn = [];
SkipOff = [];
%find I(t), index that results in lowest cost once start-up is considered.
for k = 1:1:n
    %%Try and remove the shortest generator on segment
    [onSeg, ~] = segmentLength(Binary,StartCost,dt,SkipOn,SkipOff);
    if isempty(onSeg)
        break %no  segments to check
    else
        [~,Ion] = min(onSeg(:,4)+StartCost(onSeg(:,1))'/(2*max(StartCost))); %shortest segment, and if equal lengths, highest start cost
        %only allow generators that are on at begining or end to be involved (or that have smaller start-up cost)
        Ialt = I;
        i = onSeg(Ion,1); %generator index
        t1 = onSeg(Ion,2); %index when generator turns on
        t2 = onSeg(Ion,3);%index when generator shuts off
        notAllowed = (~Binary(t1,:) & ~Binary(t2,:)); %generators that are off at the beginning and end
        notAllowed(i) = true;
        dCost = zeros(t2-t1,1);
        for t = t1:1:t2-1
            Opt = Alt.Binary{t}; %all feasible options tested at this time
            cost2 = Alt.Cost{t}; %cost for these options
            cost2(any(Opt(:,notAllowed),2)) = inf;%make the cost infinite for all feasible combinations with notAllowed generators
            if any(~isinf(cost2))
                [dCost(t-t1+1),Ialt(t)] = min(cost2);
            else
                dCost = inf;%no feasible combinations without turning on a differnt generator (don't make changes to this segment)
                break
            end
        end
        if sum(dCost)<StartCost(i) %sum of the marginal increase in cost is less than the start-up cost
            I = Ialt; %use the alternative index
            for t = t1:1:t2-1
                Binary(t+1,:) = Alt.Binary{t}(I(t),:); %best alternative option tested at this time
            end
        else
            SkipOn(end+1,1:4) = onSeg(Ion,:); %add to list of segments to avoid
        end
    end
    
    %%Try and remove the shortest generator off segment
    [~, offSeg] = segmentLength(Binary,StartCost,dt,SkipOn,SkipOff);
    if isempty(offSeg)
        break %no  segments to check
    else
        [~,Ioff] = min(offSeg(:,4)+StartCost(offSeg(:,1))'/(2*max(StartCost))); %shortest segment, and if equal lengths, highest start cost
        %only allow generators that are on at begining or end to be involved (or that have smaller start-up cost)
        Ialt = I;
        i = offSeg(Ioff,1); %generator index
        t1 = offSeg(Ioff,2); %index when generator turns on
        t2 = offSeg(Ioff,3);%index when generator shuts off
        dCost = zeros(t2-t1,1);
        for t = t1:1:t2-1
            Opt = Alt.Binary{t}; %all feasible options tested at this time
            cost2 = Alt.Cost{t}; %cost for these options
            cost2(~Opt(:,i)) = inf;%make the cost infinite for all combinations without generator i
            if any(~isinf(cost2))
                [dCost(t-t1+1),Ialt(t)] = min(cost2);
            else
                dCost = inf;%no feasible combinations keeping this generator active (don't make changes to this segment)
                break
            end
        end
        if sum(dCost)<StartCost(i) %sum of the marginal increase in cost is less than the start-up cost
            I = Ialt; %use the alternative index
            for t = t1:1:t2-1
                Binary(t+1,:) = Alt.Binary{t}(I(t),:); %best alternative option tested at this time
            end
        else
            SkipOff(end+1,1:4) = offSeg(Ioff,:); %add to list of segments to avoid
        end
    end
    
    
    %b) don't create any smaller on/off segments (how?)
    
end

%pull the corresponding best dispatches with the start-cost considered
GenOutput = zeros(nS,nG+nL);
for t = 1:1:nS
    GenOutput(t,:) = Alt.Disp{t}(I(t),:);
end

function [onSeg, offSeg] = segmentLength(Binary,StartCost,dt,SkipOn,SkipOff)
onSeg = [];
offSeg = [];
nG = length(StartCost);
nS = length(dt);
Horizon = sum(dt);
% find length of segments that a generator is on or off
for i = 1:1:nG
    if StartCost(i)>0 && any(~Binary(:,i))
        starts = nonzeros((1:nS)'.*(Binary(2:end,i) & ~Binary(1:nS,i)));
        stops = nonzeros((1:nS)'.*(~Binary(2:end,i) & Binary(1:nS,i)));
        nOn = length(starts);
        nOff = length(stops);
        if nOn>0 && nOff>0 %only look at generators that turn both on and off during the window
            if stops(1)<starts(1) %generator is off for a segment
                Seg = [i, stops(1), starts(1), sum(dt(stops(1):starts(1)-1))];
                if (isempty(SkipOff) || ~any(ismember(SkipOff,Seg,'rows'))) && Seg(4)<Horizon/4 %segment hasn't previously been ruled out, and is less than 1/4 of the horizon
                    offSeg(end+1,1:4) = Seg;
                end
                stops = stops(2:end);
                nOff  = nOff - 1;
            end
            j = 0;
            while j < nOff
                j = j + 1;
                Seg = [i, starts(j), stops(j), sum(dt(starts(j):stops(j)-1))]; %index of generator, start index, stop index, duration of segment
                if (isempty(SkipOn) || ~any(ismember(SkipOn,Seg,'rows'))) && Seg(4)<Horizon/4 %segment hasn't previously been ruled out, and is less than 1/4 of the horizon
                    onSeg(end+1,1:4) = Seg;
                end
                if j<nOn
                    Seg = [i, stops(j), starts(j+1), sum(dt(stops(j):starts(j+1)-1))];
                    if (isempty(SkipOff) || ~any(ismember(SkipOff,Seg,'rows'))) && Seg(4)<Horizon/4 %segment hasn't previously been ruled out, and is less than 1/4 of the horizon
                        offSeg(end+1,1:4) = Seg;
                    end
                end
            end
        end
    end
end