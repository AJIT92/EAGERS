function [block,FuelMix, Flow1,Flow3] = KineticCoef(block,Inlet,first)
%% find the kinetic coefficient which results in the net reforming determined by equilibrium
global Ru F
if ~isstruct(block.Current)
    a = block.Current;
    block.Current = [];
    block.Current.H2 = a;
    block.Current.CO = 0*a;
end
H2consume = sum(block.Current.H2)/(2*F*1000);
COconsume = sum(block.Current.CO)/(2*F*1000);
FuelMix = Inlet.Mixed;
Rnet.CH4 = sum(block.R_CH4);
RCH4old =0;
switch block.Reformer
    case {'internal','adiabatic'}
        CH4max = min(FuelMix.CH4,FuelMix.H2O)/block.Cells*block.RefSpacing;
        CH4min = -min(FuelMix.CO,((FuelMix.H2 + FuelMix.CO)/3)*3/4)/block.Cells*block.RefSpacing;
        X0guessRef = (sum(block.R_CH4ref)/block.RefPerc - CH4min)/(CH4max -CH4min);
        X0guessRef = max(min(X0guessRef,(1-1e-5)),1e-5);
        X0guess = max(min(block.AnPercEquilib,(1-1e-5)),1e-5);
    case 'direct'
        CH4max = min(FuelMix.CH4,FuelMix.H2O);
        CH4min = -min(FuelMix.CO-COconsume,((FuelMix.H2 + FuelMix.CO - (H2consume+COconsume))/3)*3/4);
        X0guess = (sum(block.R_CH4)/block.AnPercEquilib - CH4min)/(CH4max -CH4min);
        X0guess = max(min(X0guess,(1-1e-5)),1e-5);
    case 'external'
        CH4max = min((1-block.RefPerc)*FuelMix.CH4,FuelMix.H2O-1.8*block.RefPerc*FuelMix.CH4);
        CH4min = -min(FuelMix.CO+.2*block.RefPerc*FuelMix.CH4,((FuelMix.H2 + FuelMix.CO + 4*block.RefPerc*FuelMix.CH4 - (H2consume+COconsume))/3)*3/4);
        X0guess = (sum(block.R_CH4)/block.AnPercEquilib - CH4min)/(CH4max -CH4min);
        X0guess = max(min(X0guess,(1-1e-5)),1e-5);
end

%% first find equilibrium at outlet
count = 0;
Tol = 1e-3;
while abs((Rnet.CH4-RCH4old)/Rnet.CH4)>Tol% && (FC.Recirc.Anode>0 || count==0)
    RCH4old = Rnet.CH4;
    AnOutlet.T = mean(block.T.Flow1(block.Flow1Dir(:,end)));
    switch block.Reformer
        case 'external'
            AnInlet = FuelMix;
            AnInlet.CH4 = (1-block.RefPerc)*FuelMix.CH4;
            AnInlet.CO = (FuelMix.CO + 0.2*block.RefPerc*FuelMix.CH4);
            AnInlet.CO2 = (FuelMix.CO2 + 0.8*block.RefPerc*FuelMix.CH4);
            AnInlet.H2 = (FuelMix.H2 + 3.8*block.RefPerc*FuelMix.CH4);
            AnInlet.H2O = (FuelMix.H2O - 1.8*block.RefPerc*FuelMix.CH4);
            Flow3 =[];
        case 'internal'
            block.RefPlates = block.Cells/block.RefSpacing;
            RefInlet = FuelMix;
            RefOutlet.T = mean(block.T.Flow3(block.Flow3Dir(:,end)));
            [RefOutlet,Rref_net] = equilib2D(RefInlet,RefOutlet.T,block.Flow1_Pinit,0,0,block.FCtype,block.RefPerc,X0guessRef);
            AnInlet = RefOutlet;
        case 'direct'
            AnInlet = FuelMix;
            Flow3 =[];
        case 'adiabatic' %% find recirculation that achieves desired reformer temp
            Tol = 1e-2;
            Flow3.Outlet.T = block.ReformT;
            %% fixed recirculation
            Flow3.Inlet = FuelMix;
            [Flow3.Outlet,Rref_net,~,RefPerc] = equilibReform(Flow3.Inlet,block.Flow1_Pinit,0,Flow3.Outlet.T,1,'Q');
            
            K_CH4eq = (Flow3.Outlet.H2.^3.*Flow3.Outlet.CO)./(Flow3.Outlet.CH4.*Flow3.Outlet.H2O).*(block.Flow1_Pinit./NetFlow(Flow3.Outlet)).^2;
            K_WGSeq = Flow3.Outlet.CO2.*Flow3.Outlet.H2./(Flow3.Outlet.CO.*Flow3.Outlet.H2O);
            a = 4352.2./Flow3.Outlet.T - 3.99;
            K_WGS = exp(a);% Water gas shift equilibrium constant
            K_CH4 = 2459000*exp(-6.187*a);
            block.scaleK_CH4 = K_CH4eq./K_CH4;
            block.scaleK_WGS = K_WGSeq./K_WGS;
            
            block.ReformT = Flow3.Outlet.T; %update to have a better guess next time
            AnInlet = Flow3.Outlet; 
    end
    [AnOutlet,Rnet] = equilib2D(AnInlet,AnOutlet.T,block.Flow1_Pinit,H2consume*block.Cells,COconsume*block.Cells,block.FCtype,block.AnPercEquilib,X0guess);
    if block.Recirc.Anode >0 %only first time through
        errorR = 1;
        while abs(errorR)>1e-5 %% loop to find anode recirculation that meets steam2carbon design
            for i = 1:1:length(block.Spec1)
                FuelMix.(block.Spec1{i}) = Inlet.Flow1.(block.Spec1{i}) + block.Recirc.Anode*AnOutlet.(block.Spec1{i});
            end
            S2Ccheck = FuelMix.H2O/(FuelMix.CH4+.5*FuelMix.CO);
            errorR = (block.Steam2Carbon-S2Ccheck)/block.Steam2Carbon;
            block.Recirc.Anode = block.Recirc.Anode*(1+.9*errorR);
        end
        %%find resulting temperature of mixture
        errorT = 1;
        [~,Hin] = enthalpy(Inlet.Flow1);
        [~,Hout] = enthalpy(AnOutlet);
        Hnet = Hin + block.Recirc.Anode*Hout;
        Cp = SpecHeat(AnOutlet);
        NetFlowMix = NetFlow(FuelMix);
        while abs(errorT)>1e-3
            [~,Hmix] = enthalpy(FuelMix);
            errorT = (Hnet-Hmix)/(Cp*NetFlowMix);
            FuelMix.T = FuelMix.T + errorT;
        end  
    end
    count = count+1;
end 
% %% From running equilibrium function I can calculate exponential fit to WGS equilibrium: K_WGS = exp(4189.8./T -3.8242) : slightly different than in Paradis paper

%% Indirect Reformer
switch block.Reformer
    case 'internal'
        if first
            X_CH4in = RefInlet.CH4/NetFlow(RefInlet);
            X_CH4out = (RefInlet.CH4-Rref_net.CH4)/(NetFlow(RefInlet)+2*Rref_net.CH4);
            lambda = log(X_CH4out/X_CH4in)/(-block.columns); %exponential decay in CH4
            R_cumulative = zeros(length(block.Flow3Dir(:,1)),1);
            XCH4 = zeros(block.nodes,1);
            for i= 1:1:block.columns
                k = block.Flow3Dir(:,i);
                XCH4(k) = X_CH4in*exp(-i*lambda);
                if i == 1 % R = (in flow - outflow) = (Xin*flowin + Xout*(flowin +2*R)) solved for R
                    R.CH4(k,1) = (RefInlet.CH4/length(k) - XCH4(k).*(NetFlow(RefInlet)/length(k) +2*R_cumulative))./(1+2*XCH4(k));
                    R_cumulative = R_cumulative+R.CH4(k);
                else
                    R.CH4(k,1) = (XCH4(kold) - XCH4(k)).*(NetFlow(RefInlet)/length(k)+2*R_cumulative)./(1+2*XCH4(k));
                    R_cumulative = R_cumulative+R.CH4(k);
                end
                kold = k;
            end
            R.WGS = Rref_net.WGS/Rref_net.CH4*R.CH4; %assume same initial distribution for WGS reaction
        else
            for r = 1:1:block.rows 
                k_r = block.Flow3Dir(r,:);
                R.CH4(k_r,1) = block.R_CH4ref(k_r)*(Rref_net.CH4/block.rows)/sum(block.R_CH4ref(k_r));%make sure the total reforming is correct.
                R.WGS(k_r,1) = block.R_WGSref(k_r)*(Rref_net.WGS/block.rows)/sum(block.R_WGSref(k_r)); %assume same initial distribution for WGS reaction
            end
        end
        R.CH4 = R.CH4/block.RefPlates;
        R.WGS = R.WGS/block.RefPlates;
        RefCurrent.H2 = zeros(length(R.CH4),1);
        RefCurrent.CO = zeros(length(R.CH4),1);
        [R, Flow3, K] = FindKineticCoef(RefInlet,block.T.Flow3,R,block.Flow3Dir,Rref_net.CH4/block.RefPlates,RefCurrent,block.Flow1_Pinit,block.FCtype,block.RefPlates,block.method,1e-5);
        block.KineticCoeff3 = K;
        block.R_CH4ref = R.CH4;
        block.R_WGSref = R.WGS;  
    case 'adiabatic'
        %% solve for recircualtion that gives desired reformer Toutlet;
%         [Hin,~] = enthalpy(Inlet.AnodeIn); %total energy (not just sensible)
%         [Rref_net,Rnet,FC,Reformer,AnInlet,AnOutlet] = ReformRecirculation(FuelMix,Reformer,Hin,FC);
        %%---%
        R.CH4ref = Rref_net.CH4;
        R.WGSref = Rref_net.WGS;
        nout = NetFlow(Flow3.Outlet);
        X_CH4 = Flow3.Outlet.CH4./nout*block.Flow1_Pinit*1000; %partial pressures in Pa
        X_H2O = Flow3.Outlet.H2O./nout*block.Flow1_Pinit*1000; %partial pressures in Pa
        if strcmp(block.method,'Achenbach')
            block.KineticCoeff3 = R.CH4ref/(X_CH4.*exp(-8.2e4./(Ru*Flow3.Outlet.T)));
        elseif strcmp(block.method,'Leinfelder')
            block.KineticCoeff3 = R.CH4ref/(30.8e10*X_CH4.*X_H2O.*exp(-2.05e5./(Ru*Flow3.Outlet.T)));
        elseif strcmp(block.method,'Drescher')
            block.KineticCoeff3 = R.CH4ref/(288.52*X_CH4.*X_H2O.*exp(-1.1e4./(Ru*Flow3.Outlet.T))./(1+16*X_CH4+0.143*X_H2O.*exp(3.9e4./(Ru*Flow3.Outlet.T))));
        end
        block.R_CH4ref = R.CH4ref;
        block.R_WGSref = R.WGSref;
end
%% Anode Reforming
if first
    X_CH4in = AnInlet.CH4/NetFlow(AnInlet);
    X_CH4out = (AnInlet.CH4-Rnet.CH4)/(NetFlow(AnInlet)+2*Rnet.CH4);
    lambda = log(X_CH4out/X_CH4in)/(-block.columns); %exponential decay in CH4
    R_cumulative = zeros(length(block.Flow1Dir(:,1)),1);
    XCH4 = zeros(block.nodes,1);
    for i= 1:1:block.columns
        k = block.Flow1Dir(:,i);
        XCH4(k) = X_CH4in*exp(-i*lambda);
        if i == 1 % R = (in flow - outflow) = (Xin*flowin + Xout*(flowin +2*R)) solved for R
            R.CH4(k,1) = (AnInlet.CH4/length(k) - XCH4(k).*(NetFlow(AnInlet)/length(k) +2*R_cumulative))./(1+2*XCH4(k));
            R_cumulative = R_cumulative+R.CH4(k);
        else
            R.CH4(k,1) = (XCH4(kold) - XCH4(k)).*(NetFlow(AnInlet)/length(k) + 2*R_cumulative)./(1+2*XCH4(k));
            R_cumulative = R_cumulative+R.CH4(k);
        end
        kold = k;
    end
    R.WGS = Rnet.WGS/Rnet.CH4*R.CH4; %assume same initial distribution for WGS reaction
else
    R.CH4 = block.R_CH4*Rnet.CH4/sum(block.R_CH4); %keep same disribution as last time, but make the sum equal to the global calculation
    R.WGS = block.R_WGS*Rnet.WGS/sum(block.R_WGS); %assume same initial distribution for WGS reaction
end
R.CH4 = R.CH4/block.Cells;
R.WGS = R.WGS/block.Cells;
[R, Flow1, K] = FindKineticCoef(AnInlet,block.T.Flow1,R,block.Flow1Dir,Rnet.CH4/block.Cells,block.Current,block.Flow1_Pinit,block.FCtype,block.Cells,block.method,1e-5);
block.KineticCoeff1 = K;
block.R_CH4 = R.CH4;
block.R_WGS = R.WGS;

function [R, Flow, KinCoef] = FindKineticCoef(Inlet,T_Out,R, Dir, referenceR_CH4,Current,Pressure,Type,Cells,method,Tol)
global Ru F
specInterest = {'CH4','CO','CO2','H2','H2O'};
Flow = FCin2Out(T_Out,Inlet,Dir,Type,Cells,Current,R,'anode');
nout = NetFlow(Flow.Outlet);
X_CH4 = Flow.Outlet.CH4./nout*Pressure*1000; %partial pressures in Pa
X_H2O = Flow.Outlet.H2O./nout*Pressure*1000; %partial pressures in Pa
r = length(Dir(:,1));%rows
if strcmp(method,'Achenbach')
    KinCoef = R.CH4./(X_CH4.*exp(-8.2e4./(Ru*T_Out))); %best guess of KinCoef
elseif strcmp(method,'Leinfelder')
    KinCoef = R.CH4./(X_CH4.*30.8e10*X_H2O.*exp(-2.05e5./(Ru*T_Out)));
elseif strcmp(method,'Drescher')
    KinCoef = R.CH4./(X_CH4.*288.52*X_H2O.*exp(-1.1e4./(Ru*T_Out(Dir(:,1))))/(1+16*X_CH4+0.143*X_H2O.*exp(3.9e4./(Ru*T_Out))));
end
KinCoef = sum(KinCoef.*R.CH4/referenceR_CH4);

spec = fieldnames(Inlet);
count =0;
error = 1;   
while abs(error)>Tol %iterate to converge on a kinetic coefficients (if less CH4 in exhaust than equilibrium, smaller coefficient)
    R_CH4a = loopConverge(Flow,R,T_Out,Pressure,KinCoef,Dir,method);
    eK = 1e-6*KinCoef;
    R_CH4b = loopConverge(Flow,R,T_Out,Pressure,KinCoef+eK,Dir,method);
    dR_dK = (sum(R_CH4b) - sum(R_CH4a))/eK;
    error = (referenceR_CH4/sum(R_CH4a)-1);
    KinCoef = KinCoef + (referenceR_CH4 - sum(R_CH4a))/dR_dK; %adjust kinetic coefficient
    R.CH4 = loopConverge(Flow,R,T_Out,Pressure,KinCoef,Dir,method);

    R_CH4max = (1-1e-10)*min(Flow.Inlet.CH4,Flow.Inlet.H2O);
    R.CH4 = min(R.CH4,R_CH4max);
    for j= 1:1:length(Dir(1,:))
        k = Dir(:,j);
        if j == 1
            X.CH4 = Inlet.CH4/Cells/r - R.CH4(k);
            X.CO = Inlet.CO/Cells/r + R.CH4(k) - Current.CO(k)/(2*F*1000);
            X.CO2 = Inlet.CO2/Cells/r + Current.CO(k)/(2*F*1000);
            X.H2 = Inlet.H2/Cells/r + 3*R.CH4(k) - Current.H2(k)/(2*F*1000);%hydrogen consumed
            X.H2O = Inlet.H2O/Cells/r - R.CH4(k) + Current.H2(k)/(2*F*1000);% water produced
            if strcmp(Type,'MCFC')
                X.CO2 = X.CO2 + (Current.H2(k)+Current.CO(k))/(2*F*1000); % CO2 brought over
            end
            for i = 1:1:length(spec)
                if ~ismember(spec{i},specInterest)
                    X.(spec{i}) = Inlet.(spec{i})/Cells/r;
                end
            end
        else
            X.CH4 = X.CH4 - R.CH4(k);
            X.CO = X.CO + R.CH4(k) - Current.CO(k)/(2*F*1000);
            X.CO2 = X.CO2 + Current.CO(k)/(2*F*1000);
            X.H2 = X.H2 + 3*R.CH4(k) - Current.H2(k)/(2*F*1000);%hydrogen consumed
            X.H2O = X.H2O - R.CH4(k) + Current.H2(k)/(2*F*1000);% water produced
            if strcmp(Type,'MCFC')
                X.CO2 = X.CO2 + (Current.H2(k)+Current.CO(k))/(2*F*1000); % CO2 brought over
            end
        end

        R_COmin = -min(X.CO2,X.H2);
        R_COmax = min(X.H2O,X.CO);
        y0 = max(0+1e-5,min(1-1e-5,(R.WGS(k)-R_COmin)./(R_COmax-R_COmin)));
        
        y0 = Newton1D(y0,X,R_COmin,R_COmax,T_Out(k),Pressure,specInterest,1e-6,'GibbVal');

        R.WGS(k) = R_COmin+y0.*(R_COmax-R_COmin);
        X.CO = X.CO-R.WGS(k);
        X.CO2 = X.CO2+R.WGS(k);
        X.H2 = X.H2+R.WGS(k);
        X.H2O = X.H2O-R.WGS(k);
    end
    count = count+1;
end
% disp(strcat('FindKineticCoef count is:',num2str(count)))

function R_CH4 = loopConverge(Flow,R,T_Out,Pressure,KinCoef,Dir,method)
global Ru
%% find new reforming reaction rates
k = Dir(:,1);
n_in = NetFlow(Flow.Inlet);
n_in = n_in(k);
H2O_in = Flow.Inlet.H2O(k);
X_CH4in = Flow.Inlet.CH4(k)./n_in;
nout = NetFlow(Flow.Outlet);
X_CH4 = Flow.Outlet.CH4./nout;%initial guess of X_CH4
for j= 1:1:length(Dir(1,:))
    k = Dir(:,j);
    if strcmp(method,'Achenbach')
        C = exp(-8.2e4./(Ru*T_Out(k)));
    elseif strcmp(method,'Leinfelder')
        X_H2O = (H2O_in - R.CH4(k) - R.WGS(k))./(n_in+2*R.CH4(k))*Pressure*1000; %partial pressures in Pa
        C = 30.8e10*X_H2O.*exp(-2.05e5./(Ru*T_Out(k)));
    elseif strcmp(method,'Drescher')
        X_H2O = (H2O_in - R.CH4(k) - R.WGS(k))./(n_in+2*R.CH4(k))*Pressure*1000; %partial pressures in Pa
        C = 288.52*X_H2O.*exp(-1.1e4./(Ru*T_Out))./(1+16*X_CH4(k)*Pressure*1000+0.143*X_H2O.*exp(3.9e4./(Ru*T_Out(k))));
    end

    % use newton method to find X_CH4_out that makes R.CH4 so that R = X_CH4 in *flow in - X_CH4out*flow out
    dX1 = (X_CH4in.*n_in - KinCoef.*X_CH4(k)*Pressure*1000.*C)./(n_in + 2*KinCoef*X_CH4(k)*Pressure*1000.*C) - X_CH4(k);
    eX = 1e-5*X_CH4in;
    dX2 = (X_CH4in.*n_in - KinCoef.*(X_CH4(k)+eX)*Pressure*1000.*C)./(n_in + 2*KinCoef*(X_CH4(k)+eX)*Pressure*1000.*C) - (X_CH4(k)+eX);
    m = (dX2-dX1)./eX;
    b = dX1 - m.*X_CH4(k);
    X_CH4(k) = - b./m; 
    
    R.CH4(k) = KinCoef*X_CH4(k)*Pressure*1000.*C;
    %inlet to the next column
    X_CH4in = (X_CH4in.*n_in - R.CH4(k))./(n_in + 2*R.CH4(k));
    H2O_in = H2O_in - R.CH4(k) - R.WGS(k);
    n_in = n_in + 2*R.CH4(k);
%     errorJ = (X_CH4(k) - X_CH4in)./X_CH4in; %remaining error after 1 newton step
end
R_CH4 = R.CH4;