function Out = Reformer(t,Y, Inlet,block,string1)
%Nodal reformer model adiabatic or with heat exchange from second fluid
% Two or Four (2 or 4) inlets: {'Primary Flow', 'Primary Pout', 'Secondary Flow','Secondary Pout'}
global Ru Tags
newSpec = fieldnames(Inlet.Flow1);
for i = 1:1:length(block.spec1)
    if ~ismember(block.spec1{i},newSpec)
        Inlet.Flow1.(block.spec1{i}) = 0;
    end
end
S2C = Inlet.Flow1.H2O/(Inlet.Flow1.CH4 + .5*Inlet.Flow1.CO);
%% seperate out temperatures
nodes = block.nodes;
n = 0;
Primary.Outlet.T = Y(n+1:n+nodes);n = n+nodes;
Tsolid = Y(n+1:n+nodes);n = n+nodes;
s = 0;
N = 2;
if isfield(Inlet,'Flow2')
    newSpec = fieldnames(Inlet.Flow2);
    for i = 1:1:length(block.spec2)
        if ~ismember(block.spec2{i},newSpec)
            Inlet.Flow2.(block.spec2{i}) = 0;
        end
    end
    
    SecondaryOut.T = Y(n+1:n+nodes);n = n+nodes;
    PsecondaryIn = Y(end);
    NsecondaryOut = block.PfactorSecondary*(PsecondaryIn-Inlet.Flow2Pout);%total secondary flow out
    s = 1;
    N = 3;
end
PprimaryIn = Y(end-s);
NprimaryOut = block.PfactorPrimary*(PprimaryIn-Inlet.Flow1Pout);%total cold flow out

%% Cold flow
for i = 1:1:length(block.spec1)
    Primary.Outlet.(block.spec1{i}) = max(0,Y(n+1:n+nodes));n = n+nodes;
end
for j = 1:1:length(block.Flow1Dir(1,:));%1:columns
    k = block.Flow1Dir(:,j);
    r = length(k);
    if j==1
        Primary.Inlet.T(k,1) = Inlet.Flow1.T;
        for i = 1:1:length(block.spec1)
            Primary.Inlet.(block.spec1{i})(k,1) = Inlet.Flow1.(block.spec1{i})/r;
        end
    else
        Primary.Inlet.T(k,1) = Primary.Outlet.T(kprev,1);
        for i = 1:1:length(block.spec1)
            Primary.Inlet.(block.spec1{i})(k,1) = Primary.Outlet.(block.spec1{i})(kprev,1);
        end
    end
    kprev = k;
end
if isfield(Inlet,'Flow2')
    %% Hot flow
    for i = 1:1:length(block.spec2)
        SecondaryOut.(block.spec2{i}) = max(0,Y(n+1:n+nodes));n = n+nodes;
    end
    for j = 1:1:length(block.Flow2Dir(1,:));%1:columns
        k = block.Flow2Dir(:,j);
        r = length(k);
        if j==1
            SecondaryIn.T(k,1) = Inlet.Flow2.T;
            for i = 1:1:length(block.spec2)
                SecondaryIn.(block.spec2{i})(k,1) = Inlet.Flow2.(block.spec2{i})/r;
            end
        else
            SecondaryIn.T(k,1) = SecondaryOut.T(kprev);
            for i = 1:1:length(block.spec2)
                SecondaryIn.(block.spec2{i})(k,1) = SecondaryOut.(block.spec2{i})(kprev,1);
            end
        end
        kprev = k;
    end
    CooledOut.T  = mean(SecondaryOut.T(block.Flow2Dir(:,end))); %temperature 
    for i = 1:1:length(block.spec2)
        CooledOut.(block.spec2{i}) = sum(SecondaryOut.(block.spec2{i})(block.Flow2Dir(:,end)));
    end
    Flow2 = NetFlow(CooledOut);
end

ReformedOut.T  = mean(Primary.Outlet.T(block.Flow1Dir(:,end))); %temperature 
for i = 1:1:length(block.spec1)
    ReformedOut.(block.spec1{i}) = sum(Primary.Outlet.(block.spec1{i})(block.Flow1Dir(:,end)));
end
Flow1 = NetFlow(ReformedOut);

if strcmp(string1,'Outlet')
    %% Outlet Ports
    Out.Reformed = ReformedOut;
    Out.ReformedPin = PprimaryIn;
    Out.MeasureS2C = S2C;
    Out.MeasureReformT = ReformedOut.T;
    Tags.(block.name).H2flow = ReformedOut.H2;
    Tags.(block.name).MeasureS2C = S2C;
    Tags.(block.name).ReformedT = ReformedOut.T;
    Tags.(block.name).CH4 = ReformedOut.CH4;
    if isfield(Inlet,'Flow2')
        Out.Cooled = CooledOut;
        Out.CooledPin = PsecondaryIn;
        Out.MeasureCooledT = CooledOut.T;
        Tags.(block.name).CooledT = CooledOut.T;
    end
elseif strcmp(string1,'dY')
%     K_WGS = exp(4189.8./ReformerOut.T -3.8242);% Water gas shift equilibrium constant
    a = 4352.2./Primary.Outlet.T - 3.99;
    K_WGS = block.scaleK_WGS.*exp(a);% Water gas shift equilibrium constant
    K_CH4 = block.scaleK_CH4.*2459000.*exp(-6.187*a);
    CH4_eq = (Primary.Outlet.H2.^3.*Primary.Outlet.CO)./(K_CH4.*Primary.Outlet.H2O).*(PprimaryIn./NetFlow(Primary.Outlet)).^2;
    R.CH4 = (Primary.Inlet.CH4)-CH4_eq; %inlet CO + CO from reforming - outlet CO
    CO_eq = Primary.Outlet.CO2.*Primary.Outlet.H2./(K_WGS.*Primary.Outlet.H2O);
    R.WGS = (Primary.Inlet.CO+R.CH4)-CO_eq; %inlet CO + CO from reforming - outlet CO
    
    %confirm we don't violate anything casuing negative species
    R.CH4 = min([R.CH4,Primary.Inlet.CH4,Primary.Inlet.H2O],[],2);
    R.CH4 = max([R.CH4,-Primary.Inlet.CO,-Primary.Inlet.H2/3],[],2);
    R.WGS = min([R.WGS,Primary.Inlet.CO+R.CH4,Primary.Inlet.H2O-R.CH4],[],2);
    R.WGS = max([R.WGS,-Primary.Inlet.CO2,-(Primary.Inlet.H2+3*R.CH4)],[],2);
    %identify the target outflow that dY will converge to
    RefOut.T = Primary.Outlet.T;
    for i = 1:1:length(block.spec1)
        if strcmp(block.spec1{i},'CO2') 
            RefOut.CO2 = (Primary.Inlet.CO2 + R.WGS); 
        elseif strcmp(block.spec1{i},'H2') 
            RefOut.H2 = (Primary.Inlet.H2 + 3*R.CH4 + R.WGS);
        elseif strcmp(block.spec1{i},'H2O') 
            RefOut.H2O = (Primary.Inlet.H2O - R.CH4 - R.WGS); 
        elseif strcmp(block.spec1{i},'CH4') 
            RefOut.CH4 = (Primary.Inlet.CH4 - R.CH4); 
        elseif strcmp(block.spec1{i},'CO') 
            RefOut.CO = (Primary.Inlet.CO + R.CH4 - R.WGS);  
        else
            RefOut.(block.spec1{i}) = Primary.Inlet.(block.spec1{i});  
        end
    end 

    HoutCold = enthalpy(Primary.Outlet);
    scale = NetFlow(Primary.Outlet)./NetFlow(RefOut);
    HinCold = enthalpy(Primary.Inlet).*scale;
    Cp_cold = SpecHeat(Primary.Outlet);
    
    if isfield(Inlet,'Flow2')
        HoutHot = enthalpy(SecondaryOut);
        scale = NetFlow(SecondaryOut)./NetFlow(SecondaryIn);
        HinHot = enthalpy(SecondaryIn).*scale;
        Cp_hot = SpecHeat(SecondaryOut);
        QT = block.HTconv*Y(1:N*nodes) + block.HTcond*Y(1:N*nodes);
    else
        QT = block.HTconv*Y(1:N*nodes) + block.HTcond*Y(1:N*nodes);
    end

    dY = 0*Y;
    % time constants for states
    tC1 = (block.Vol_1*Cp_cold*PprimaryIn./(Ru*Primary.Outlet.T));
    tC2 = (block.Vol_Solid*block.Solid_Density*block.Solid_SpecHeat);
    for i=1:1:length(block.Flow1Dir(1,:))
        k = block.Flow1Dir(:,i);
        dY(k)= (QT(k) + HinCold(k) - HoutCold(k))./tC1(k); %Cold flow
        if i>1
            dY(k) = dY(k)+dY(kprev);
        end
        
        kprev = k;
    end
    dY(nodes+1:2*nodes)= QT(nodes+1:2*nodes)./tC2;  % Plate
    for j = 1:1:length(block.spec1)
        dY((1+j+s)*nodes+1:(2+j+s)*nodes)= (RefOut.(block.spec1{j})- Primary.Outlet.(block.spec1{j})) ./block.tC((1+j+s)*nodes+1:(2+j+s)*nodes); %all species concentration
    end 

    if isfield(Inlet,'Flow2')
        tC3 = (block.Vol_2*Cp_hot*PsecondaryIn./(Ru*SecondaryOut.T));
        for i=1:1:length(block.Flow2Dir(1,:))
            k = block.Flow2Dir(:,i);
            dY(2*nodes+k)= (QT(2*nodes+k) + HinHot(k) - HoutHot(k))./tC3(k); %Hot flow
            if i>1
                dY(2*nodes+k) = dY(2*nodes+k)+dY(2*nodes+kprev);
            end
            
            kprev = k;
        end
        n = (2+j);
        for j = 1:1:length(block.spec2)
            dY((n+j)*nodes+1:(n+1+j)*nodes) = (SecondaryIn.(block.spec2{j}) - SecondaryOut.(block.spec2{j}))./block.tC((n+j)*nodes+1:(n+1+j)*nodes);  %all  species concentration
        end
        dY(end) = (Flow2-NsecondaryOut)*Ru*Inlet.Flow2.T/block.Vol_2;
    end

    dY(end-s) = (Flow1-NprimaryOut)*Ru*Inlet.Flow1.T/block.Vol_1;
    Out = dY;
end