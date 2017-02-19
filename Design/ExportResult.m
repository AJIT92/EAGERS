function ExportResult(CallFcn, varargin)
global Project Model_dir
    stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
                 'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
                 'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
                 'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
                 'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
    stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
                  'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};    
    BuildType = {'Restaurant: full-service (sit down)';'Restaurant: quick-service (fast food)';'School: primary school';'School: secondary school';'Office: large office';'Office: medium office';'Office: small office';'Mid-rise apartment building';'Hospitality: large hotel';'Hospitality: small hotel/motel';'Health care: large hospital';'Health care: outpatient facility';'Retail: big-box, standalone retail store';'Retail: retail store located in a strip mall';'Retail: supermarket';'Unrefrigerated warehouse';'Multiple';'Real Building'};
MSG = msgbox('Exporting...');
if strcmp(CallFcn,'ViewResults')
     Boiler = varargin{1};
    name = char(Project.Name);
    [f,p]=uiputfile('*.xls','Save Project As...',fullfile(Model_dir,'results', name));
    if f==0; return; end
    filename = fullfile(p,f);

    %% System Description
    ColumnLabelSheet1 = {'Property','Value','Units'};
    xlswrite(filename,ColumnLabelSheet1,'System','A1:C1');
    RowLabelSheet1 = {'Project name';'Control Strategy';'Economic Case';'Building Name';'Electric utility';'Gas Utility';'User Type';'System Size';...
                      'Electric Chiller Size';'Absorption Chiller Size';'Thermal Storage Size';'Battery Size';'Wind Size';'Solar Size';'CHP Turndown';};
    xlswrite(filename,RowLabelSheet1,'System','A2:A16');

    data1 = {Project.Name; Project.Control.Name; Project.Economic.Name; Project.Building.Name; Project.Utilities.Grid.Name; Project.Utilities.NatGas.Name;};
    xlswrite(filename,data1,'System','B2:B8');

    CHPSize = 0;
    minSize = 0;
    for i = 1:1:length(Project.System.CHP)
        CHPSize = CHPSize+Project.System.CHP(i).SysSize(1);
        minSize = minSize+Project.System.CHP(i).SysSize(2);
    end
    Turndown = CHPSize/minSize;
    data2 = [Project.Result.Dispatch.SysSize; Project.Result.Dispatch.ChillerSize(1)/3.517; Project.Result.Dispatch.ChillerSize(2)/3.517; Project.Result.Dispatch.TESsize/3.517; Project.Result.Dispatch.BatterySize; Project.Result.Dispatch.RenewableSize(2); Project.Result.Dispatch.RenewableSize(1); Turndown;];
    xlswrite(filename,data2,'System','B9:B16');

    data3 = {'kW','Tons','Tons','Ton-hr','kW-hr','kW','kW','N/A'};
    xlswrite(filename,data3,'System','C9:C16');
    %% Building/load
    ColumnLabelSheet2 = {'Hour','Electricity (kWh)','Cooling (kWh)','Heating (kWh)','Generated Power (kW)', 'TES Shift', 'Battery Shift'};
    xlswrite(filename,ColumnLabelSheet2,'Building','A1:G1');
    data4 = [linspace(.25,8760,4*8760)' Project.Building.DemandE Project.Building.DemandC Project.Building.DemandH];
    data4(:,5) = Project.Result.eOut.GenPower;
    data4(:,6) = Project.Result.eOut.DemandEshift;
    data4(:,7) = Project.Result.eOut.DemandEshift2;
    xlswrite(filename,data4,'Building','A2:G35041');

    %% Dispatch
    ColumnLabelSheet3 = {'Hour','Baseline Demand','CHP Elec Generation (kWh)','Grid Power (kWh)','Heating (not met by CHP) (kWh)'};
    xlswrite(filename,ColumnLabelSheet3,'CHP','A1:E1');
    data5 = [linspace(1,8760,8760)' Project.Result.Baseline.Elec (Project.Result.Baseline.Elec-Project.Result.Dispatch.Elec) Project.Result.Dispatch.Elec Project.Result.eOut.BoilerHeatHour];
    xlswrite(filename,data5,'CHP','A2:E8761');
    ColumnLabelSheet4 = ['Property','Value','Units'];
    
    %% Finances
    xlswrite(filename,ColumnLabelSheet4,'Result','A1:C1');
    cost = Project.Result.costOut;
    NPV = [cost.NPVbaselineUseCharges; cost.NPVbaselineDemCharges; cost.NPVbaselineFuelCost; cost.NPVbaselineOandM; cost.NPVbaselineFinance; cost.NPVnewUseCharges; cost.NPVnewDemCharges; cost.NPVnewFuelCost; cost.NPVnewOandM; cost.NPVnewFinance;];

    data6 = {'CHP size';'CHP capacity factor'; 'Self Generation proportion'; 'Proportion of heating met by CHP';'Annual  baseline grid cost'; 'Annual cost with CHP ';...
             'Grid Cost with CHP';'Fuel Cost with CHP'; 'Avg Grid Electricity Price (Baseline)';'Avg Grid Electricty price (with CHP)'; 'Avg value of CHP Electricity';...
             'NPV Use Charges (baseline)';  'NPV Demand Charges (baseline)';  'NPV Fuel Cost (baseline)';  'NPV OandM (baseline)';  'NPV Finance (baseline)';  ...
             'NPV Use Charges (w/ CHP)';  'NPV Demand Charges (w/ CHP)';  'NPV Fuel Cost (w/ CHP)';  'NPV OandM (w/ CHP)';  'NPV Finance (w/ CHP)'; };
    xlswrite(filename,data6,'Result','A2:A22');
    data7 = [Project.Result.eOut.SysSize; Project.Result.eOut.Total_Electricity_produced_kWh/(Project.Result.eOut.SysSize*8760); Project.Result.eOut.SelfGen;...
             Project.Result.eOut.Total_DG_Heat_out_kWh/(Project.Result.eOut.Total_DG_Heat_out_kWh+ Project.Result.eOut.Total_peak_burner_heat_out_kWh); ...
             mean(cost.BaselineGridBill+cost.BaselineFuel); mean(cost.CostPerYear); mean(cost.NewGridBill); mean(cost.FuelCosts);...
             Project.Result.costOut.Baseline.TotalGridBill/sum(Project.Result.Baseline.Elec);...
             Project.Result.costOut.Dispatch.TotalGridBill/sum(Project.Result.Dispatch.Elec);...
             (Project.Result.costOut.Baseline.TotalGridBill-Project.Result.costOut.Dispatch.TotalGridBill)/Project.Result.Dispatch.ElecTotProd; NPV];
    xlswrite(filename,data7,'Result','B2:B22');

    %% Emissions
    ColumnLabelSheet5 = {'Year','Month','Electric Rate ($/kWh)','Gas Rate ($/mmBTU','Hour','Baseline CO2 (tons)','Baseline SO2 (lbs)','Baseline NOx (lbs)','w/ CHP CO2 (tons)','w/ CHP SO2 (lbs)','w/ CHP NOx (lbs)'};
    xlswrite(filename,ColumnLabelSheet5,'Forecast Rates & Emissions','A1:K1');
    for i = 1:1:20
        Years(1+12*(i-1):12*i) = 2013+i;
        Months(1+12*(i-1):12*i) = linspace(1,12,12);
    end
    data9 = [Years' Months' Project.Result.costOut.FuelPriceMonthly' Project.Result.costOut.GridPriceMonthly'];
    xlswrite(filename,data9,'Forecast Rates & Emissions','A2:D241');

    BoilerCO2 = [0.3988 0.5565 0.7066 0.4742];%lb CO2 per kWh
    BoilerCO2 = BoilerCO2(Boiler);
    BoilerNOx = 7.7e-4; %lbNOx/kWh
    BoilerSO2 = 2.5e-3; %lb SO2/kWh
   
    Ts = 8760/length(Project.Building.DemandE);
    DemandE = Project.Building.DemandE;
    DemandH = Project.Building.DemandH;
    BaselinePurchaseHour = zeros(1,8760);
    BaselineHeatHour = zeros(1,8760);
    for i = 1:1:8760
        BaselinePurchaseHour(i) = sum(DemandE(1+1/Ts*(i-1):1/Ts*i)*Ts);
        BaselineHeatHour(i) = sum(DemandH(1+1/Ts*(i-1):1/Ts*i)*Ts);
    end
    
    StateNum = find(strcmp(Project.State,stateName),1);

    State = stateAbrev(StateNum);
    [CO2, NOx, SO2] = EmissionProfile(State); %retrive emission profile
    data10 = linspace(1,8760,8760)';
    data10(:,2) = (BaselinePurchaseHour.*CO2' + BoilerCO2*BaselineHeatHour)/1e3;
    data10(:,3) = (BaselinePurchaseHour.*NOx' + BoilerNOx*BaselineHeatHour);
    data10(:,4) = (BaselinePurchaseHour.*SO2' + BoilerSO2*BaselineHeatHour);
    data10(:,5) = (Project.Result.eOut.Grid_purchases_hourly_kWh.*CO2 + Project.Result.eOut.CO2 + BoilerCO2*Project.Result.eOut.BoilerHeatHour)/1e3;
    data10(:,6) = (Project.Result.eOut.Grid_purchases_hourly_kWh.*NOx + Project.Result.eOut.NOx + BoilerNOx*Project.Result.eOut.BoilerHeatHour);
    data10(:,7) = (Project.Result.eOut.Grid_purchases_hourly_kWh.*SO2 + Project.Result.eOut.SO2 + BoilerSO2*Project.Result.eOut.BoilerHeatHour);
    xlswrite(filename,data10,'Forecast Rates & Emissions','E2:K8761');
elseif strcmp(CallFcn,'NationalSurvey')
    result = varargin{1};
    name = char(result.StudyName);
    [f,p]=uiputfile('*.xls','Save Project As...',fullfile(Model_dir,'results', name));
    if f==0; return; end
    filename = fullfile(p,f);
    FieldNames = fields(result);
    
    %% NPC results
    stateList5x = ['Demand Charge';stateName;'Energy Charges'; stateName; 'Fuel Charges'; stateName; 'O&M Charge';stateName;'Finance Charges';stateName;];
    xlswrite(filename,BuildType','NPC Baseline','B2:Q2');
    xlswrite(filename,stateList5x,'NPC Baseline','A2:A256');
    if ~strcmp(result.BuildName','All')
        NumBuild = 1;
    else
        NumBuild = 16;
    end
    Data = zeros(254,NumBuild);
    Data(1:50,1:NumBuild) = squeeze(result.NPCbaseline(1,:,:));
    Data(52:101,1:NumBuild) = squeeze(result.NPCbaseline(2,:,:));
    Data(103:152,1:NumBuild) = squeeze(result.NPCbaseline(3,:,:));
    Data(154:203,1:NumBuild) = squeeze(result.NPCbaseline(4,:,:));
    Data(205:254,1:NumBuild) = squeeze(result.NPCbaseline(5,:,:));
    xlswrite(filename,Data,'NPC Baseline','B3:Q256');
    
    xlswrite(filename,BuildType','NPC Dispatch','B2:Q2');
    xlswrite(filename,stateList5x,'NPC Dispatch','A2:A256');
    Data = zeros(254,NumBuild);
    Data(1:50,1:NumBuild) = squeeze(result.NPCdispatch(1,:,:));
    Data(52:101,1:NumBuild) = squeeze(result.NPCdispatch(2,:,:));
    Data(103:152,1:NumBuild) = squeeze(result.NPCdispatch(3,:,:));
    Data(154:203,1:NumBuild) = squeeze(result.NPCdispatch(4,:,:));
    Data(205:254,1:NumBuild) = squeeze(result.NPCdispatch(5,:,:));
    xlswrite(filename,Data,'NPC Dispatch','B3:Q256');
    
    if strcmp(result.BuildName,'all')
        for i = 1:1:length(FieldNames)
            B =size(result.(char(FieldNames{i})));
            if B(2) ==16 && ~isstruct(result.(char(FieldNames{i})))
                xlswrite(filename,BuildType',char(FieldNames{i}),'B1:Q1');
                xlswrite(filename,stateName,char(FieldNames{i}),'A2:A51');
                xlswrite(filename,result.(char(FieldNames{i})),char(FieldNames{i}),'B2:Q51');
            end
        end
    else %% write all to 1 sheet
    end
    stateList3x = ['CO2';stateName;'NOx'; stateName; 'SO2'; stateName;];
    xlswrite(filename,BuildType','Baseline Emissions','B2:Q2');
    xlswrite(filename,stateList3x,'Baseline Emissions','A2:A154');
    Data = zeros(152,NumBuild);
    Data2 = zeros(152,NumBuild);
    for j = 1:50
        for k = 1:NumBuild
            Data(j,k) = result.Baseline(j,k).CO2;
            Data(j+51,k) = result.Baseline(j,k).NOx;
            Data(j+102,k) = result.Baseline(j,k).SO2;
            Data2(j,k) = result.Dispatch(j,k).CO2;
            Data2(j+51,k) = result.Dispatch(j,k).NOx;
            Data2(j+152,k) = result.Dispatch(j,k).SO2;
        end
    end
    xlswrite(filename,Data,'Baseline Emissions','B3:Q154');
    xlswrite(filename,BuildType','Dispatch Emissions','B2:Q2');
    xlswrite(filename,stateList3x,'Dispatch Emissions','A2:A154');
    xlswrite(filename,Data2,'Dispatch Emissions','B3:Q154');
            
end
close(MSG)
