function [dX_dt, SS_1] = RampRateCalc(SS,LB,UB,Hratio)
global scaleTime
SS = SS(1,:);
if isfield(SS,'Dt')
    Dt = SS.Dt;
else Dt =1;
    SS.Dt = 1;
end
if Dt~=1 %convert to 1 second sampling time
    convSS = ss(SS.A,SS.B,SS.C,SS.D,Dt);
    newSS = d2d(convSS,1);
    r = length(newSS);
    SS_1.A = newSS(1).A;
    SS_1.B = newSS(1).B;
    SS_1.C = newSS(1).C;
    SS_1.D = newSS(1).D;
    for k = 2:1:r
        SS_1.C(end+1,:) = newSS(k).C;
        SS_1.D(end+1,:) = newSS(k).D;
    end
else SS_1 = SS;
end
x0 = LB;
if ~isempty(Hratio)
    x0(2) = LB*Hratio; %CHP heat produed per unit electricity
end
nS = round(4*3600/Dt)+1; % assume ramping is less than 4 hours (i forget why I made this limit)
t = linspace(0, Dt*(nS-1),nS);
u = UB*linspace(1,1,nS)';
[z,z2] = size(SS.C);
X0 = zeros(z2,1);
for k = 1:1:z
    X0(find(SS.C(k,:),1,'first'))=x0(k);%finds the first non-zero element in SS.C and makes X0 at that index = x0
end
SS = ss(SS.A,SS.B,SS.C,SS.D,Dt);
[y,t] = lsim(SS,u,t,X0);
tRise = t(find(y(:,1)>(.95*u(1)),1,'first'))/3600; %rise time in hours
if isempty(tRise)
    tRise = 4;
end
dX_dt = (UB.*(0.95)./tRise)./scaleTime;
