function [BaselineEmission DispatchEmission] = EmissionsCalculated(RESULT,State,GridMix,Boiler)

stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};         
if length(char(State))>2
    Num = find(strcmp(State,stateName));
    State = char(stateAbrev(Num));
end
[CO2, NOx, SO2] = EmissionProfile(State,0,GridMix);

BoilerCO2 = [0.3988 0.5565 0.7066 0.4742];%lb CO2 per kWh
BoilerCO2 = BoilerCO2(Boiler);
BoilerNOx = 7.7e-4; %lbNOx/kWh
BoilerSO2 = 2.5e-3; %lb SO2/kWh

BaselineEmission = [sum(RESULT.Baseline.Elec.*CO2)/2e3 , 0 ,sum(BoilerCO2*RESULT.Baseline.Heat)/2e3;
                    sum(RESULT.Baseline.Elec.*NOx) , 0 ,sum(BoilerNOx*RESULT.Baseline.Heat);
                    sum(RESULT.Baseline.Elec.*SO2) , 0 ,sum(BoilerSO2*RESULT.Baseline.Heat);];
NewGridCO2 = RESULT.Dispatch.Elec.*CO2;
NewGridNOx = RESULT.Dispatch.Elec.*NOx;
NewGridSO2 = RESULT.Dispatch.Elec.*SO2;
NewGridCO2(NewGridCO2<0)=0;
NewGridNOx(NewGridNOx<0)=0;
NewGridSO2(NewGridSO2<0)=0;
DispatchEmission = [sum(NewGridCO2)/2e3 , sum(RESULT.eOut.CO2)/2e3 , sum(BoilerCO2*RESULT.eOut.BoilerHeatHour)/2e3;
                    sum(NewGridNOx) , sum(RESULT.eOut.NOx) , sum(BoilerNOx*RESULT.eOut.BoilerHeatHour);
                    sum(NewGridSO2) , sum(RESULT.eOut.SO2) , sum(BoilerSO2*RESULT.eOut.BoilerHeatHour);];