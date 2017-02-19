%% Runs Projects
global Model_dir
projdir=fullfile(Model_dir, 'DesignProjects');
files=dir(fullfile(projdir,'*.mat'));
list=strrep({files.name},'.mat','');

for i=1:1:length(list)
    projName=fullfile(Model_dir,'project',list{i});
    load(projName);

    costs = Project.Result.costOut;    
    costs = [costs.NPVbaselineUseCharges ...
            costs.NPVbaselineDemCharges costs.NPVbaselineFuelCost ...
            costs.NPVbaselineOandM costs.NPVbaselineFinance; ...
            costs.NPVnewUseCharges ...
            costs.NPVnewDemCharges costs.NPVnewFuelCost costs.NPVnewOandM ...
            costs.NPVnewFinance;];
    label1 = [{'NPC'},{'Energy'},{'Demand'},{'Fuel'},{'O & M'},{'Finance'}];
    label2 = [{'Baseline'}; {'with CHP'}];
   
    
    toWrite = fullfile(Model_dir,'results','Test');      
    xlswrite(toWrite,costs,list{i},'B2:F3');
    xlswrite(toWrite,label1,list{i},'A1:F1');
    xlswrite(toWrite,label2,list{i},'A2:A3');
    
    
    %Emissions Results
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
    RESULTS = Project.Result;
    BaselineEmission = [sum(BaselinePurchaseHour.*CO2')/1e3 , 0 ,sum(BoilerCO2*BaselineHeatHour)/1e3;
                        sum(BaselinePurchaseHour.*NOx') , 0 ,sum(BoilerNOx*BaselineHeatHour);
                        sum(BaselinePurchaseHour.*SO2') , 0 ,sum(BoilerSO2*BaselineHeatHour);];
    DispatchEmission = [sum(Project.Result.eOut.Grid_purchases_hourly_kWh.*CO2)/1e3 , sum(RESULTS.eOut.CO2)/1e3 , sum(BoilerCO2*RESULTS.eOut.BoilerHeatHour)/1e3;
                        sum(Project.Result.eOut.Grid_purchases_hourly_kWh.*NOx) , sum(RESULTS.eOut.NOx) , sum(BoilerNOx*RESULTS.eOut.BoilerHeatHour);
                        sum(Project.Result.eOut.Grid_purchases_hourly_kWh.*SO2) , sum(RESULTS.eOut.SO2) , sum(BoilerSO2*RESULTS.eOut.BoilerHeatHour);];
    label3 = [{'Baseline'};{'CO2'};{'NOx'};{'SO2'}];
    label4 = [{'Grid'},{'CHP'},{'Boiler'}];    
    label5 = [{'with CHP'};{'CO2'};{'NOx'};{'SO2'}];
    xlswrite(toWrite,BaselineEmission,list{i},'B5:D7');
    xlswrite(toWrite,DispatchEmission,list{i},'B9:D11');
    xlswrite(toWrite,label3,list{i},'A4:A7')
    xlswrite(toWrite,label4,list{i},'B4:D4')
    xlswrite(toWrite,label4,list{i},'B8:D8')
    xlswrite(toWrite,label5,list{i},'A8:A11')

end