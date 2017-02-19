stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID'; ... 
                'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';...
                'MO';'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';...
                'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';'USavg';};
data = struct;

for i=1:1:length(stateAbrev)
    state = char(stateAbrev(i));
    data.(state).meanElecAnnual = mean(ElecRate.(state).CommercialAnnual);
    data.(state).meanGas = mean(GasRate.(state).Commercial);
    data.(state).meanElec = mean(ElecRate.(state).Commercial);
    data.(state).ratioElecOverGas=data.(state).meanElec/data.(state).meanGas;
    Ratios(1,i) = data.(state).ratioElecOverGas;
    Elec(i,1) = data.(state).meanElec;
    ElecAnnual(i,1) = data.(state).meanElecAnnual;
    CommercialGas(i,1) = data.(state).meanGas
    ForGenGas(i,1) = mean(GasRate.(state).ForElectricGeneration);
end