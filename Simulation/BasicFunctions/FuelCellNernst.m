function FuelCellNernst(nCurrent,T,P,X,block)
global F Ru Tags
[h,~] = enthalpy(T,{'H2';'H2O';'O2';});
s = entropy(T,{'H2';'H2O';'O2';});
E0 = -((h.H2O-s.H2O.*T)-(h.H2-s.H2.*T)-.5*(h.O2-s.O2.*T))/(2*F);
switch block.FCtype
    case {'SOFC';'SOEC';'oxySOFC';} % note sign convention of current is negative for electrolyzers !!
        frac = sqrt(X.O2).*X.H2./X.H2O;
        V_Nernst = E0 + Ru*T/(2*F).*log(abs(sqrt(P/101.325)*frac));
        currentDen = nCurrent/block.A_Node;

        %removed the activation energy, because reversible SOFC shos linear through OCV
    %     frac2 = X.O2-Ru*T/(4*F).*currentDen*block.t_Cath/block.DeffO2/(P*1000);
    %     V_activation = Ru*T/(4*F*block.alpha).*log(currentDen/block.Io/(P/101.325)./frac2);
    %     if min(currentDen/block.Io/(P/101.325)./frac2)<exp(1) %in very low current areas, make activation linear function of current density, lets activation loss linearly aproach 0 as current approaches 0
    %         k = find((currentDen/block.Io/(P/101.325)./frac2)<exp(1));
    %         V_activation(k) = Ru*T(k)/(4*F*block.alpha).*currentDen(k)/block.Io/exp(1);
    %     end
    %     V_ohmic = currentDen*block.t_Membrane.*T/block.ElecConst./(exp(-block.deltaG./(Ru*T)));
        V_activation =0*nCurrent;
        V_concentration =0*nCurrent;
        V_ohmic = currentDen*block.t_Membrane/block.ElecConst.*(exp(block.deltaG./T));
        Vloss = V_activation+V_concentration+V_ohmic;
        nVoltage = V_Nernst-Vloss;
        ASR = mean(Vloss./(currentDen/(100^2))); %should be in neighborhood of 0.25
    case {'MCFC';'MCEC';'oxyMCFC';}
        frac = sqrt(X.O2).*X.CO2c.*X.H2./X.H2O./X.CO2a;
        V_Nernst = E0 + Ru*T/(2*F).*log(abs(sqrt(P/101.325).*frac));
        currentDen = nCurrent/block.A_Node;
        if mean(currentDen)<0
            pow = -1;
        else pow = 1;
        end
        V_activation = pow*Ru*T/(4*F*block.alpha).*log(abs(currentDen)/block.Io);
        if min(abs(currentDen)/block.Io)<exp(1) %in very low current areas, make activation linear function of current density, lets activation loss linearly aproach 0 as current approaches 0
            k = find((abs(currentDen)/block.Io)<exp(1));
            V_activation(k) = pow*Ru*T(k)/(4*F*block.alpha).*abs(currentDen(k))/block.Io/exp(1);
        end
        V_concentration = pow*Ru*T/(4*F)*(1+1/block.alpha).*log(block.J_L./(block.J_L-abs(currentDen)));
        V_ohmic = currentDen.*(block.Cr0+block.Cr1*(T-273));
        Vloss = V_activation+V_concentration+V_ohmic;
        nVoltage = V_Nernst-Vloss;
        ASR = mean(Vloss./(currentDen/(100^2)));
end
Tags.(block.name).ASR = ASR;
Tags.(block.name).LocalCurrentDensity = nCurrent/block.A_Node/1e4;%current in A/cm2
Tags.(block.name).nCurrent = nCurrent;
Tags.(block.name).nVoltage = nVoltage;
Tags.(block.name).nPower = abs(nCurrent).*Tags.(block.name).nVoltage;% power per node in W
Tags.(block.name).Power = sum(Tags.(block.name).nPower)*block.Cells/1000;%power for stack in kW
Tags.(block.name).LocalNernst = V_Nernst;
Tags.(block.name).LocalActivation = V_activation;
Tags.(block.name).LocalConcentration = V_concentration;
Tags.(block.name).LocalOhmic = V_ohmic;
