global Model_dir

sheetName = 'CA';
projdir=fullfile(Model_dir, 'project');
files=dir(fullfile(projdir,'*.mat'));
list=strrep({files.name},'.mat','');

projName=fullfile(Model_dir,'project',list{1});
load(projName);
costs = Project.Result.costOut;    
baselineCosts = [costs.NPVbaselineUseCharges ...
        costs.NPVbaselineDemCharges costs.NPVbaselineFuelCost ...
        costs.NPVbaselineOandM costs.NPVbaselineFinance;];
label1 = [{'NPC'},{'Energy'},{'Demand'},{'Fuel'},{'O & M'},{'Finance'}];

toWrite = fullfile(Model_dir,'results','NY Emission Results');
xlswrite(toWrite,label1,sheetName,'A1:F1');
xlswrite(toWrite,{'Baseline'},sheetName,'A2');
xlswrite(toWrite,baselineCosts,sheetName,'B2:F2');



%Emissions Results
label2 = [{'CO2'},{'Grid'},{'CHP'},{'Boiler'}];
label3 = [{'NOx'},{'Grid'},{'CHP'},{'Boiler'}];
label4 = [{'SO2'},{'Grid'},{'CHP'},{'Boiler'}];
xlswrite(toWrite,label2,sheetName,'I1:L1');
xlswrite(toWrite,label3,sheetName,'M1:P1');
xlswrite(toWrite,label4,sheetName,'Q1:T1');

[CO2, NOx, SO2] = EmissionProfile('CA',0,1);
DemandE = Project.Building.DemandE;
DemandH = Project.Building.DemandH;
Ts = 8760/length(DemandE);
BaselinePurchaseHour = zeros(1,8760);
BaselineHeatHour = zeros(1,8760);
for j = 1:1:8760
    BaselinePurchaseHour(j) = sum(DemandE(1+1/Ts*(j-1):1/Ts*j)*Ts);
    BaselineHeatHour(j) = sum(DemandH(1+1/Ts*(j-1):1/Ts*j)*Ts);
end

Boiler = 1;
BoilerCO2 = [0.3988 0.5565 0.7066 0.4742];%lb CO2 per kWh
BoilerCO2 = BoilerCO2(Boiler);
BoilerNOx = 7.7e-4; %lbNOx/kWh
BoilerSO2 = 2.5e-3; %lb SO2/kWh

BaselineCO2 = [sum(BaselinePurchaseHour.*CO2')/1e3 , 0 ,sum(BoilerCO2*BaselineHeatHour)/1e3;];
BaselineNOx = [sum(BaselinePurchaseHour.*NOx') , 0 ,sum(BoilerNOx*BaselineHeatHour);];
BaselineSO2 = [sum(BaselinePurchaseHour.*SO2') , 0 ,sum(BoilerSO2*BaselineHeatHour);];
xlswrite(toWrite,BaselineCO2,sheetName,'J2:L2');
xlswrite(toWrite,BaselineNOx,sheetName,'N2:P2');
xlswrite(toWrite,BaselineSO2,sheetName,'R2:T2');

for i=1:1:length(list)
    projName=fullfile(Model_dir,'project',list{i});
    load(projName);

    costs = Project.Result.costOut;
    costCHP = [costs.NPVnewUseCharges ...
              costs.NPVnewDemCharges costs.NPVnewFuelCost costs.NPVnewOandM ...
              costs.NPVnewFinance;];
    
    row = num2str(i+2);
    xlswrite(toWrite,{Project.Name},sheetName,strcat('A',row));
    xlswrite(toWrite,costCHP,sheetName,strcat('B',row,':','F',row));
    
    %Emissions Results
    DemandE = Project.Building.DemandE;
    DemandH = Project.Building.DemandH;
    
    RESULTS = Project.Result;
    DispatchCO2 = [sum(Project.Result.eOut.Grid_purchases_hourly_kWh.*CO2)/1e3 , sum(RESULTS.eOut.CO2)/1e3 , sum(BoilerCO2*RESULTS.eOut.BoilerHeatHour)/1e3;];
    DispatchNOx = [sum(Project.Result.eOut.Grid_purchases_hourly_kWh.*NOx) , sum(RESULTS.eOut.NOx) , sum(BoilerNOx*RESULTS.eOut.BoilerHeatHour);];
    DispatchSO2 = [sum(Project.Result.eOut.Grid_purchases_hourly_kWh.*SO2) , sum(RESULTS.eOut.SO2) , sum(BoilerSO2*RESULTS.eOut.BoilerHeatHour);];
    
    xlswrite(toWrite,DispatchCO2,sheetName,strcat('J',row,':','L',row));
    xlswrite(toWrite,DispatchNOx,sheetName,strcat('N',row,':','P',row));
    xlswrite(toWrite,DispatchSO2,sheetName,strcat('R',row,':','T',row));
    
end
  