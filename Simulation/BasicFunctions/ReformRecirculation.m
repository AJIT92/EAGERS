function [Rref_net,Rnet,FC,Reformer,AnInlet,AnOutlet] = ReformRecirculation(FuelMix,Reformer,Hin,FC)
% this function works if we want to find the recirculation which achieves a
% desired temperature in an adibatic reformer (requires more than standard iterations)
% It is commented out in line 873
errorR = 100;
errorOld=0;
count2 = 0;
Rnet.CH4 = [];
X0guess = .99;
while abs(errorR)>1
    Reformer.Inlet = FuelMix;
    [Reformer.Outlet,Rref_net,~,RefPerc] = equilibReform(Reformer.Inlet,FC.FuelPinit,0,Reformer.Outlet.T,1,'Q');
    for i = 1:1:length(FC.AnSpec)
        AnInlet.(FC.AnSpec{i}) = Reformer.Outlet.(FC.AnSpec{i})/FC.Cells;
    end
    if~isempty(Rnet.CH4)
        CH4max = min(AnInlet.CH4,AnInlet.H2O);
        CH4min = -min(AnInlet.CO,((AnInlet.H2 + AnInlet.CO - sum(FC.Current)/(2*F*1000))/3)*3/4);
        X0guess = (Rnet.CH4 - CH4min)/(CH4max -CH4min); 
    end
    [AnOutlet,Rnet] = equilib2D(AnInlet,AnOutlet.T,FC.FuelPinit,sum(FC.Current)/(2*F*1000),FC.FCtype,FC.AnPercEquilib,X0guess);
    if FC.ReformT > Reformer.Outlet.T
        errorR = log(1+FC.ReformT - Reformer.Outlet.T);
    else errorR = -real(log(-1+FC.ReformT - Reformer.Outlet.T));
    end
    %note: total energy of the fuel mix needs to be total
    %energy of the reformer outlet (adiabatic)
    Rold = FC.Recirc.Anode;
    if count2>=first
        if abs(errorR)>2
            FC.Recirc.Anode = max(.3,min(.99,FC.Recirc.Anode + errorR/100 - errorOld/300));
        else FC.Recirc.Anode = max(.3,min(.99,FC.Recirc.Anode + errorR/200 - errorOld/600));
        end
    end
    dr = FC.Recirc.Anode - Rold;
    nScale = (1-Rold)/((1-Rold)-dr);
    errorOld = errorR;
    %% find new fuel mix 
    for i = 1:1:length(FC.AnSpec)
        FuelMix.(FC.AnSpec{i}) = (StackAnIn.(FC.AnSpec{i}) + FC.Recirc.Anode*AnOutlet.(FC.AnSpec{i})*FC.Cells*nScale);
    end
    FC.Steam2Carbon = FuelMix.H2O/(FuelMix.CH4+.5*FuelMix.CO);
    if FC.Steam2Carbon<2
        FC.Recirc.Anode = FC.Recirc.Anode*2/FC.Steam2Carbon;
        for i = 1:1:length(FC.AnSpec)
            FuelMix.(FC.AnSpec{i}) = StackAnIn.(FC.AnSpec{i}) + FC.Recirc.Anode*AnOutlet.(FC.AnSpec{i})*FC.Cells*nScale;
        end
        FC.Steam2Carbon = FuelMix.H2O/(FuelMix.CH4+.5*FuelMix.CO);
    end
    %%find resulting temperature of mixture
    errorT = 100;
    Hnet = (Hin + FC.Recirc.Anode*enthalpy(AnOutlet)*FC.Cells*nScale);
    Cp = SpecHeat(AnOutlet);
    NetFlowMix = NetFlow(FuelMix);
    while abs(errorT)>1
        errorT = (Hnet-enthalpy(FuelMix))/(Cp*NetFlowMix);
        FuelMix.T = FuelMix.T + errorT;
    end    
    count2 = count2+1;
end