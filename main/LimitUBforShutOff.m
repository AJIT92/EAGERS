function LimitUBforShutOff(GenDisp,dX)
%this function may be necessary now that on-line optimization has ramp rates included
global Plant OnOff UB LB UBopt LBopt
%OnOff: the desired generator state from the controller (may differ from actual on/off because it leaves generators on as they ramp down to min power, then shuts them off)
%UB: the upper limit (capacity) of each generator
%LB: the lower limit (capacity) of each generator when on
%UBmpc: upper bound of generators seen by MPC, can be different from the global upper bound due to ramping constraints.
%LBmpc: lower bound of generators seen by MPC, can be different from the global lower bound due to ramping constraints
include = {'CHP Generator', 'Electric Generator', 'Chiller'};
UBopt = UB;
LBopt = LB;
for i = 1:1:length(Plant.Generator)
    if ismember(Plant.Generator(i).Type,include)
        if GenDisp(1,i)<LB(i) && GenDisp(2,i)>GenDisp(1,i)%generator in process of starting
            UBopt(i) = min(UB(i),GenDisp(1,i)+dX(1,i));
            LBopt(i) = GenDisp(1,i); %in process of startin up, so cant go below current state
        elseif OnOff(i)>0 
            LBopt(i) = max(LB(i),GenDisp(1,i)-dX(1,i));
            r = find(GenDisp(2:end,i)>=LB(i),1,'first');
            if ~isempty(r)
                nextOff = find(GenDisp(r+2:end,i)<LB(i),1,'first')+r;
                if ~isempty(nextOff)
                    set = GenDisp(nextOff+1,i);
                    UBopt(i) = min(UB(i),set+sum(dX(1:nextOff,i)));
                end
            end
        elseif OnOff(i)==0 
            LBopt(i) = 0;
        end
    end
end