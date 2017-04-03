function FuelCellNernst(Flow1,Flow2,Current,T,P,block)
global F Ru Tags
if isstruct(Current)
    netCurrent = (Current.H2+Current.CO);
else netCurrent = Current;
end
%calclae local concentrations
n_flow1_in = NetFlow(Flow1.Inlet);
n_flow1_out = NetFlow(Flow1.Outlet);
n_flow2_in = NetFlow(Flow2.Inlet);
n_flow2_out = NetFlow(Flow2.Outlet);
switch block.FCtype
    case {'SOFC';'SOEC';'oxySOFC';}
        if isfield(block,'ClosedCathode') && block.ClosedCathode
            X.O2 = ones(block.nodes,1);
        else
            X.O2 = (Flow2.Outlet.O2+Flow2.Inlet.O2)./(n_flow2_in+n_flow2_out);
        end
    case {'MCFC';'MCEC';'oxyMCFC';}
        if isfield(block,'ClosedCathode') && block.ClosedCathode
            X.O2 = ones(block.nodes,1)/3;
            X.CO2c = 2*ones(block.nodes,1)/3;
        else
            X.O2 = (Flow2.Outlet.O2+Flow2.Inlet.O2)./(n_flow2_in+n_flow2_out);
            X.CO2c = (Flow2.Outlet.CO2+Flow2.Inlet.CO2)./(n_flow2_in+n_flow2_out);
        end
end
X.H2 = (Flow1.Outlet.H2+Flow1.Inlet.H2)./(n_flow1_in+n_flow1_out);
X.H2O = (Flow1.Outlet.H2O+Flow1.Inlet.H2O)./(n_flow1_in+n_flow1_out);
if isfield(Flow1.Outlet,'CO')
    X.CO = (Flow1.Outlet.CO+Flow1.Inlet.CO)./(n_flow1_in+n_flow1_out);
end
if isfield(Flow1.Outlet,'CO2')
    X.CO2a = (Flow1.Outlet.CO2+Flow1.Inlet.CO2)./(n_flow1_in+n_flow1_out);
end

k = block.Flow1Dir(:,1);
if min(Flow1.Inlet.H2(k))==0 && isfield(Flow1.Outlet,'CH4')%gives some reformed methane as anode inlet
    R.CH4(k) = Flow1.Inlet.CH4(k) - Flow1.Outlet.CH4(k); %Rate of methane reforming R.CH4
    K_WGS(k) = exp(4189.8./Flow1.Outlet.T(k) -3.8242);% Water gas shift equilibrium constant
    CO_eq(k) = Flow1.Outlet.CO2(k).*Flow1.Outlet.H2(k)./(K_WGS(k).*Flow1.Outlet.H2O(k));
    R.WGS(k) = (Flow1.Inlet.CO(k)+R.CH4(k))-CO_eq(k); %inlet CO + CO from reforming - outlet CO
    X.H2(k) = (Flow1.Outlet.H2(k) + 0.5*(3*R.CH4(k)+R.WGS(k)))./(n_flow1_in(k)+n_flow1_out(k));
    X.H2O(k) = (Flow1.Outlet.H2O(k) - 0.5*(R.CH4(k) + R.WGS(k))+ Flow1.Inlet.H2O(k))./(n_flow1_in(k)+n_flow1_out(k));
end

%% Calculate local voltages
currentDen = abs(netCurrent)/block.A_Node;% A/m^2
h = enthalpy(T,{'H2';'H2O';'O2'; 'CO'; 'CO2'});
s = entropy(T,{'H2';'H2O';'O2';'CO'; 'CO2'});
if any(strcmp(block.FCtype,{'MCFC';'oxyMCFC';}))
    CO2_frac = X.CO2c./X.CO2a;
else
    CO2_frac = 1;
end
E0_H2 = -((h.H2O-s.H2O.*T)-(h.H2-s.H2.*T)-.5*(h.O2-s.O2.*T))/(2*F);
frac_H2 = sqrt(X.O2).*X.H2./X.H2O.*CO2_frac;% reaction 1/2*O2 + H2 --> H2O    or   1/2*O2 + CO2 --> CO2 + H2O
V_Nernst_H2 = E0_H2 + Ru*T/(2*F).*log(abs(sqrt(P/101.325)*frac_H2));

if isfield(X,'CO') && isfield(X,'CO2a')
    E0_CO = -((h.CO2-s.CO2.*T)-(h.CO-s.CO.*T)-.5*(h.O2-s.O2.*T))/(2*F); %reference voltage for CO
    frac_CO = sqrt(X.O2).*X.CO./X.CO2a.*CO2_frac;% reaction 1/2*O2 + CO --> CO2    or   1/2*O2 + CO2 + CO --> 2*CO2
    V_Nernst_CO = E0_CO + Ru*T/(2*F).*log(abs(sqrt(P/101.325)*frac_CO)); %nernst voltage for CO reaction
end

if mean(netCurrent)<0 %electrolyzer or Fuel cell mode
    pow = -1;
else pow = 1;
end
switch block.FCtype
    case {'SOFC';'SOEC';'oxySOFC';} % note sign convention of current is negative for electrolyzers !!
        %removed the activation energy, because reversible SOFC shows linear through OCV
        V_concentration =0*netCurrent;
        V_activation =0*netCurrent;
    %     frac2 = X.O2-Ru*T/(4*F).*currentDen*block.t_Cath/block.DeffO2/(P*1000);
    %     V_activation = Ru*T/(4*F*block.alpha).*log(currentDen/block.Io/(P/101.325)./frac2);
    %     if min(currentDen/block.Io/(P/101.325)./frac2)<exp(1) %in very low current areas, make activation linear function of current density, lets activation loss linearly aproach 0 as current approaches 0
    %         k = find((currentDen/block.Io/(P/101.325)./frac2)<exp(1));
    %         V_activation(k) = Ru*T(k)/(4*F*block.alpha).*currentDen(k)/block.Io/exp(1);
    %     end        
        V_ohmic = currentDen*block.t_Membrane/block.ElecConst.*(exp(block.deltaG./T));
    
%         if isfield(X,'CO') && isfield(X,'CO2a')
%             %%calculate current attributable to direct CO reaction
%             m =.25; %value given in literature
%             E_A = 110000; % kJ/kmol // activation energy for anode (value from literature)
%             E_C = 160000; %kJ/kmol // activation energy for cathode (value from literature)
%             k_anode_h2 = 2.2e6; %pre-exponential factor for H2
%             k_anode_CO = 5.5e5; %pre-exponential factor for CO
%             k_cathode_O2= 1.5e8; %pre-exponential factor for O2
%             RA_H2 = 1./((2*F./(Ru*T)).*k_anode_h2.*((X.H2*P/101.325).^m).*exp(-E_A./(Ru*T))); % this is expression for resistance from hydrogen conversion
%             RA_CO = 1./((2*F./(Ru*T)).*k_anode_CO.*((X.CO*P/101.325).^m).*exp(-E_A./(Ru*T))); % this is expression for resistance from CO conversion
%             RC_O2 = 1./((4*F./(Ru*T)).*k_cathode_O2.*((X.O2*P/101.325).^m).*exp(-E_C./(Ru*T)));% this is expression for resistance from cathode
%             I_CO =max(0,(V_Nernst_CO - V_Nernst_H2 + RA_H2.*netCurrent)./(RA_H2+RA_CO)); %solving for current produced by CO
%             I_H2 = netCurrent-I_CO;   %solving for current from H2
%             Vloss = V_activation + V_concentration + V_ohmic + RA_H2.*I_H2 + RC_O2.*netCurrent;
%         else
            %% calculate current as only H2 reactions
            I_CO = 0*netCurrent;
            I_H2 = netCurrent;
            Vloss = V_activation + V_concentration + V_ohmic;
%         end
    case {'MCFC';'MCEC';'oxyMCFC';}
        V_activation = pow*Ru*T/(4*F*block.alpha).*log(abs(currentDen)/block.Io);
        if min(abs(currentDen)/block.Io)<exp(1) %in very low current areas, make activation linear function of current density, lets activation loss linearly aproach 0 as current approaches 0
            k = find((abs(currentDen)/block.Io)<exp(1));
            V_activation(k) = pow*Ru*T(k)/(4*F*block.alpha).*abs(currentDen(k))/block.Io/exp(1);
        end
        V_concentration = pow*Ru*T/(4*F)*(1+1/block.alpha).*log(block.J_L./(block.J_L-abs(currentDen)));
        I_CO = 0*netCurrent;
        I_H2 = netCurrent-I_CO;   %solving for current from H2
        
        V_ohmic = currentDen.*(block.Cr0+block.Cr1*(T-273));
        Vloss = V_activation+V_concentration+V_ohmic;
end
nVoltage = V_Nernst_H2 - pow*Vloss;%this is the voltage from the fuel cell
ASR = mean(Vloss./(currentDen/(100^2))); %should be in neighborhood of 0.25 for SOFC
Tags.(block.name).ASR = ASR;
Tags.(block.name).LocalCurrentDensity = netCurrent/block.A_Node/1e4;%current in A/cm2
Tags.(block.name).nCurrent = netCurrent;
Tags.(block.name).nVoltage = nVoltage;
Tags.(block.name).nPower = abs(netCurrent).*Tags.(block.name).nVoltage;% power per node in W
Tags.(block.name).Power = sum(Tags.(block.name).nPower)*block.Cells/1000;%power for stack in kW
Tags.(block.name).LocalNernst = V_Nernst_H2;
Tags.(block.name).LocalActivation = V_activation;
Tags.(block.name).LocalConcentration = V_concentration;
Tags.(block.name).LocalOhmic = V_ohmic;
Tags.(block.name).I_CO = I_CO;
Tags.(block.name).I_H2 = I_H2;