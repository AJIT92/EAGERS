function MF = MassFlow(Flow)
MolarMass.CH4 = 16;
MolarMass.CO = 28;
MolarMass.CO2 = 44;
MolarMass.H2 = 2;
MolarMass.H2O = 18;
MolarMass.N2 = 28;
MolarMass.O2 = 32;

speciesName = fieldnames(Flow);
MF = 0;
for i = 1:1:length(speciesName)
    if ~strcmp(speciesName{i},'T')
        MF = MF + Flow.(speciesName{i})*MolarMass.(speciesName{i});
    end
end