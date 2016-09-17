function GetCurrentTime
global RealTime Virtual DateSim Plant
if Virtual
    Tmpc = Plant.optimoptions.Tmpc;
    D = datevec(DateSim);
    D(6) = D(6)+Tmpc;
elseif RealTime
    D = clock;
end
DateSim = datenum(D);