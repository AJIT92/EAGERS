global Plant
Pant.Name = 'ColumbiaRiverDams';
%%Columbia river: 'GC';'CJ';'W';'RR';'RI';'WNPM';'PR';'MN';'JD';'TD';'BNVL';
%%Snake River: 'BRWNL';'OX';'HC';'LWRGRNT';'LTLGSE';'LWRMONUM';'ICEHRBR';


%%rows 1->366 (hourly): October 2nd, 2015 -> October 1st, 2016; 
Plant.Data.Timestamp = linspace(datenum([2015 10 2 1 0 0]),datenum([2016 10 2 0 0 0]),8784);
%% Need temperatures for WA on these dates
% Plant.Data.Temperature = [];
%%columns 1->18 (kcfs or KW): Grand Coulee, Chief Joseph, Wells, Rocky Ridge, Rock Island, Wanapum, PriestRiver, McNary, John Day, The Dalles, Bonneville, Brownlee, Oxbow, Hells Canyon, Lower Granite, Little Goose, Lower Monumental, Ice Harbor;
%SourcesandSinks = Mass balance; Sinks and Sources in river segments; Negative values = sinks; Positive values = sources;
load('SourcesandSinks'); load('Powerflow');load('Spillflow');load('Outflow');load('Inflow');load('PowerGen');
Plant.Data.Hydro.SourcesandSinks = SourcesandSinks; 
%Powerflow = Flow used to produce power
Plant.Data.Hydro.Powerflow = Powerflow;
%Spillflow = Flow that is spilled not used for power
Plant.Data.Hydro.Spillflow = Spillflow;
%Outflow = Full amount of discharge flow
Plant.Data.Hydro.Outflow = Outflow;
%columns (kilo-cfs): inflow(GC), infl(CJ)-Disch(GC), infl(W)-Disch(CJ), infl(RR)-Disch(W), infl(RI)-Disch(RR), infl(WNP)-Disch(RI), infl(PR)-Disch(WNP), infl(MN)-Disch(IH)-Disch(PR), infl(JD)-Disch(MN), infl(TD)-Disch(JD), infl(BNVL)-Disch(TD), inflow(BRW), infl(OX)-Disch(BRW), infl(HC)-Disch(OX), infl(LGRT)-Disch(HC), infl(LGS)-Disch(LGRT), infl(LMNT)-Disch(LG), infl(IH)-Disch(LMNT);
Plant.Data.Hydro.Inflow = Inflow;
%PowerGen = Generation (kW) every hour
Plant.Data.Hydro.PowerGen = PowerGen;
PowerGen(isnan(PowerGen)) = 0;
Plant.Data.Demand.E = sum(PowerGen,2);
%Still need 
%            storage data

%optimization options
Plant.optimoptions.Interval = 1;
Plant.optimoptions.Horizon = 24;
Plant.optimoptions.Resolution = 1;
Plant.optimoptions.Topt = 600;
Plant.optimoptions.Tmpc = 60;
Plant.optimoptions.nsSmooth = 0;
Plant.optimoptions.scaletime = 1;
Plant.optimoptions.fastsimulation = 1;
Plant.optimoptions.tspacing = 'constant';
Plant.optimoptions.sequential = 0;
Plant.optimoptions.excessHeat = 1;
Plant.optimoptions.thresholdSteps = 1;
Plant.optimoptions.Outputs = {'E'};
Plant.optimoptions.Buffer = 40;

Dam(1).Type = 'Hydro';
Dam(1).Name = {'Grand Coulee'};
Dam(1).Source = 'Water';
Dam(1).Output.E = 1;
Dam(1).Size = 9170;%This would be storage capacity in kilo-acre-feet
Dam(1).Enabled = 1;
Dam(1).VariableStruct.cityName = 'Grand Coulee, Wa';
Dam(1).VariableStruct.upriver = {};
Dam(1).VariableStruct.RM = 596.6;
Dam(1).VariableStruct.transTime = 87.37;
Dam(1).VariableStruct.MaxGenCapacity = 4414000; %power in kW
Dam(1).VariableStruct.RampUp = 189; %power change in kW/hr
Dam(1).VariableStruct.RampDown = 233;
Dam(1).VariableStruct.MaxGenFlow = 202.5;
Dam(1).VariableStruct.MaxSpillFlow = 0.2;
Dam(1).VariableStruct.MaxHead = 334.2;

Dam(2).Type = 'Hydro';
Dam(2).Name = {'Chief Joseph'}; 
Dam(2).Source = 'Water';
Dam(2).Output.E = 1;
Dam(2).Size = 588; %This would be storage capacity in kilo-acre-feet
Dam(2).Enabled = 1;
Dam(2).VariableStruct.cityName = {'Bridgeport, Wa'}; 
Dam(2).VariableStruct.upriver = {'Grand Coulee'};
Dam(2).VariableStruct.RM = 545.1; 
Dam(2).VariableStruct.transTime = 94.3; 
Dam(2).VariableStruct.MaxGenCapacity = 2180000; %power in kW
Dam(2).VariableStruct.RampUp =532; %power change in kW/hr
Dam(2).VariableStruct.RampDown =586; 
Dam(2).VariableStruct.MaxGenFlow = 179.3;
Dam(2).VariableStruct.MaxSpillFlow = 31.3; 
Dam(2).VariableStruct.MaxHead = 180.8; 

Dam(3).Type = 'Hydro';
Dam(3).Name = {'Wells'};
Dam(3).Source = 'Water';
Dam(3).Output.E = 1;
Dam(3).Size = 98; %This would be storage capacity in kilo-acre-feet
Dam(3).Enabled = 1;
Dam(3).VariableStruct.cityName = 'Azwell, Wa'; 
Dam(3).VariableStruct.upriver = {'Chief Joseph'}; 
Dam(3).VariableStruct.RM = 515.8; 
Dam(3).VariableStruct.transTime = 98.3; 
Dam(3).VariableStruct.MaxGenCapacity = 742000; %power in kW
Dam(3).VariableStruct.RampUp =533; %power change in kW/hr
Dam(3).VariableStruct.RampDown =608; 
Dam(3).VariableStruct.MaxGenFlow = 175.3;
Dam(3).VariableStruct.MaxSpillFlow = 57.96; 
Dam(3).VariableStruct.MaxHead = 75.19; 

Dam(4).Type = 'Hydro';
Dam(4).Name =  {'Rocky Reach'}; 
Dam(4).Source = 'Water';
Dam(4).Output.E = 1;
Dam(4).Size = 37; %This would be storage capacity in kilo-acre-feet
Dam(4).Enabled = 1;
Dam(4).VariableStruct.cityName = 'Wenatchee, Wa'; 
Dam(4).VariableStruct.upriver = {'Wells'};
Dam(4).VariableStruct.RM = 473.7; 
Dam(4).VariableStruct.transTime = 104.0; 
Dam(4).VariableStruct.MaxGenCapacity = 1012000; %power in kW
Dam(4).VariableStruct.RampUp =396; %power change in kW/hr
Dam(4).VariableStruct.RampDown =476; 
Dam(4).VariableStruct.MaxGenFlow = 159.8;
Dam(4).VariableStruct.MaxSpillFlow = 52.99;
Dam(4).VariableStruct.MaxHead =  95.19; 

Dam(5).Type = 'Hydro';
Dam(5).Name = {'Rock Island'}; 
Dam(5).Source = 'Water';
Dam(5).Output.E = 1;
Dam(5).Size = 12; %This would be storage capacity in kilo-acre-feet
Dam(5).Enabled = 1;
Dam(5).VariableStruct.cityName = 'Wenatchee, Wa'; 
Dam(5).VariableStruct.upriver = {'Rocky Reach'};
Dam(5).VariableStruct.RM = 453.4;
Dam(5).VariableStruct.transTime = 106.7; 
Dam(5).VariableStruct.MaxGenCapacity = 449000; %power in kW
Dam(5).VariableStruct.RampUp = 391; %power change in kW/hr
Dam(5).VariableStruct.RampDown = 478; 
Dam(5).VariableStruct.MaxGenFlow = 162.5;
Dam(5).VariableStruct.MaxSpillFlow = 73.4;
Dam(5).VariableStruct.MaxHead = 48.39; 

Dam(6).Type = 'Hydro';
Dam(6).Name = {'Wanapum'}; 
Dam(6).Source = 'Water';
Dam(6).Output.E = 1;
Dam(6).Size = 163; %This would be storage capacity in kilo-acre-feet
Dam(6).Enabled = 1;
Dam(6).VariableStruct.cityName = 'Beverly, Wa'; 
Dam(6).VariableStruct.upriver = {'Rock Island'}; 
Dam(6).VariableStruct.RM =  415.8; 
Dam(6).VariableStruct.transTime = 111.8; 
Dam(6).VariableStruct.MaxGenCapacity =  1517700; %power in kW
Dam(6).VariableStruct.RampUp =  924; %power change in kW/hr
Dam(6).VariableStruct.RampDown = 804; 
Dam(6).VariableStruct.MaxGenFlow = 173.7;
Dam(6).VariableStruct.MaxSpillFlow = 97.6;
Dam(6).VariableStruct.MaxHead = 85.03; 

Dam(7).Type = 'Hydro';
Dam(7).Name = {'Priest Rapids'};
Dam(7).Source = 'Water';
Dam(7).Output.E = 1;
Dam(7).Size = 44; %This would be storage capacity in kilo-acre-feet
Dam(7).Enabled = 1;
Dam(7).VariableStruct.cityName = 'Mattawa, Wa'; 
Dam(7).VariableStruct.upriver = {'Wanapum'}; 
Dam(7).VariableStruct.RM = 397.1; 
Dam(7).VariableStruct.transTime = 114.3; 
Dam(7).VariableStruct.MaxGenCapacity = 1553600; %power in kW
Dam(7).VariableStruct.RampUp = 913; %power change in kW/hr
Dam(7).VariableStruct.RampDown = 809; 
Dam(7).VariableStruct.MaxGenFlow = 162.8;
Dam(7).VariableStruct.MaxSpillFlow = 134.1;
Dam(7).VariableStruct.MaxHead = 86.14; 

Dam(8).Type = 'Hydro';
Dam(8).Name =  {'McNary'}; 
Dam(8).Source = 'Water';
Dam(8).Output.E = 1;
Dam(8).Size = 1345; %This would be storage capacity in kilo-acre-feet
Dam(8).Enabled = 1;
Dam(8).VariableStruct.cityName = 'Umatilla, Or'; 
Dam(8).VariableStruct.upriver = {'Priest Rapids', 'Ice Harbor'};
Dam(8).VariableStruct.RM = 292; 
Dam(8).VariableStruct.transTime = 128.5; 
Dam(8).VariableStruct.MaxGenCapacity = 1106540; %power in kW
Dam(8).VariableStruct.RampUp = 238; %power change in kW/hr
Dam(8).VariableStruct.RampDown = 378; 
Dam(8).VariableStruct.MaxGenFlow = 229.9;
Dam(8).VariableStruct.MaxSpillFlow = 176.8;
Dam(8).VariableStruct.MaxHead = 76.39; 

Dam(9).Type = 'Hydro';
Dam(9).Name = {'John Day'}; 
Dam(9).Source = 'Water';
Dam(9).Output.E = 1;
Dam(9).Size = 2530; %This would be storage capacity in kilo-acre-feet
Dam(9).Enabled = 1;
Dam(9).VariableStruct.cityName = 'Rufus, Or'; 
Dam(9).VariableStruct.upriver = {'McNary'}; 
Dam(9).VariableStruct.RM =  215.6; 
Dam(9).VariableStruct.transTime = 138.9; 
Dam(9).VariableStruct.MaxGenCapacity = 1932770; %power in kW
Dam(9).VariableStruct.RampUp = 221; %power change in kW/hr
Dam(9).VariableStruct.RampDown = 223; 
Dam(9).VariableStruct.MaxGenFlow = 260.5;
Dam(9).VariableStruct.MaxSpillFlow = 119;
Dam(9).VariableStruct.MaxHead = 106.42; 

Dam(10).Type = 'Hydro';
Dam(10).Name = {'The Dalles'}; 
Dam(10).Source = 'Water';
Dam(10).Output.E = 1;
Dam(10).Size = 330; %This would be storage capacity in kilo-acre-feet
Dam(10).Enabled = 1;
Dam(10).VariableStruct.cityName = 'The Dalles, Or'; 
Dam(10).VariableStruct.upriver = {'John Day'}; 
Dam(10).VariableStruct.RM = 191.5; 
Dam(10).VariableStruct.transTime = 142.1; 
Dam(10).VariableStruct.MaxGenCapacity = 2160000; %power in kW
Dam(10).VariableStruct.RampUp = 227; %power change in kW/hr
Dam(10).VariableStruct.RampDown = 228; 
Dam(10).VariableStruct.MaxGenFlow = 263.8;
Dam(10).VariableStruct.MaxSpillFlow = 128.1;
Dam(10).VariableStruct.MaxHead = 85.79; 

Dam(11).Type = 'Hydro';
Dam(11).Name = {'Bonneville'}; 
Dam(11).Source = 'Water';
Dam(11).Output.E = 1;
Dam(11).Size = 537; %This would be storage capacity in kilo-acre-feet
Dam(11).Enabled = 1;
Dam(11).VariableStruct.cityName = 'Bonneville, Or'; 
Dam(11).VariableStruct.upriver = {'The Dalles'}; 
Dam(11).VariableStruct.RM = 146.1; 
Dam(11).VariableStruct.transTime = 148.3; 
Dam(11).VariableStruct.MaxGenCapacity = 1227000; %power in kW
Dam(11).VariableStruct.RampUp = 0.319; %power change in kW/hr
Dam(11).VariableStruct.RampDown = 370;  
Dam(11).VariableStruct.MaxGenFlow = 268.2; 
Dam(11).VariableStruct.MaxSpillFlow = 164.5;
Dam(11).VariableStruct.MaxHead = 69.1; 

Dam(12).Type = 'Hydro';
Dam(12).Name = {'Brownlee'}; 
Dam(12).Source = 'Water';
Dam(12).Output.E = 1;
Dam(12).Size = [];%This would be storage capacity in kilo-acre-feet
Dam(12).Enabled = 1;
Dam(12).VariableStruct.cityName = 'Cambridge, Id'; 
Dam(12).VariableStruct.upriver = {};
Dam(12).VariableStruct.RM = 285; 
Dam(12).VariableStruct.transTime = 107.2; 
Dam(12).VariableStruct.MaxGenCapacity = []; %power in kW
Dam(12).VariableStruct.RampUp = []; %power change in kW/hr
Dam(12).VariableStruct.RampDown = [];
Dam(12).VariableStruct.MaxGenFlow = [];
Dam(12).VariableStruct.MaxSpillFlow = [];
Dam(12).VariableStruct.MaxHead = [];

Dam(13).Type = 'Hydro';
Dam(13).Name = {'Oxbow'}; 
Dam(13).Source = 'Water';
Dam(13).Output.E = 1;
Dam(13).Size = [];%This would be storage capacity in kilo-acre-feet
Dam(13).Enabled = 1;
Dam(13).VariableStruct.cityName = 'Oxbow, Or'; 
Dam(13).VariableStruct.upriver = {'Brownlee'}; 
Dam(13).VariableStruct.RM =273; 
Dam(13).VariableStruct.transTime = 108.8; 
Dam(13).VariableStruct.MaxGenCapacity = []; %power in kW
Dam(13).VariableStruct.RampUp = [];%power change in kW/hr
Dam(13).VariableStruct.RampDown = [];
Dam(13).VariableStruct.MaxGenFlow = [];
Dam(13).VariableStruct.MaxSpillFlow =[];
Dam(13).VariableStruct.MaxHead = [];

Dam(14).Type = 'Hydro';
Dam(14).Name = {'Hells Canyon'};
Dam(14).Source = 'Water';
Dam(14).Output.E = 1;
Dam(14).Size = [];%This would be storage capacity in kilo-acre-feet
Dam(14).Enabled = 1;
Dam(14).VariableStruct.cityName = 'Oxbow, Or'; 
Dam(14).VariableStruct.upriver = {'Oxbow'}; 
Dam(14).VariableStruct.RM =  247.6; 
Dam(14).VariableStruct.transTime = 112.2; 
Dam(14).VariableStruct.MaxGenCapacity = []; %power in kW
Dam(14).VariableStruct.RampUp = []; %power change in kW/hr
Dam(14).VariableStruct.RampDown = [];
Dam(14).VariableStruct.MaxGenFlow = [];
Dam(14).VariableStruct.MaxSpillFlow =[];
Dam(14).VariableStruct.MaxHead = [];

Dam(15).Type = 'Hydro';
Dam(15).Name =  {'Lower Granite'}; 
Dam(15).Source = 'Water';
Dam(15).Output.E = 1;
Dam(15).Size = 440; %This would be storage capacity in kilo-acre-feet
Dam(15).Enabled = 1;
Dam(15).VariableStruct.cityName = 'Almota, Wa'; 
Dam(15).VariableStruct.upriver = {'Hells Canyon'};
Dam(15).VariableStruct.RM = 107.5; 
Dam(15).VariableStruct.transTime = 131.2; 
Dam(15).VariableStruct.MaxGenCapacity = 627260; %power in kW
Dam(15).VariableStruct.RampUp = 215; %power change in kW/hr
Dam(15).VariableStruct.RampDown = 251; 
Dam(15).VariableStruct.MaxGenFlow = 89.6; 
Dam(15).VariableStruct.MaxSpillFlow = 46.5;
Dam(15).VariableStruct.MaxHead = 102.5; 

Dam(16).Type = 'Hydro';
Dam(16).Name = {'Little Goose'}; 
Dam(16).Source = 'Water';
Dam(16).Output.E = 1;
Dam(16).Size = 516; %This would be storage capacity in kilo-acre-feet
Dam(16).Enabled = 1;
Dam(16).VariableStruct.cityName = 'Starbuck, Wa'; 
Dam(16).VariableStruct.upriver = {'Lower Granite'}; 
Dam(16).VariableStruct.RM =  70.3; 
Dam(16).VariableStruct.transTime = 136.2; 
Dam(16).VariableStruct.MaxGenCapacity = 642820; %power in kW
Dam(16).VariableStruct.RampUp = 373; %power change in kW/hr
Dam(16).VariableStruct.RampDown = 427; 
Dam(16).VariableStruct.MaxGenFlow = 91; 
Dam(16).VariableStruct.MaxSpillFlow = 38.1;
Dam(16).VariableStruct.MaxHead = 98.96; 

Dam(17).Type = 'Hydro';
Dam(17).Name = {'Lower Monumental'};
Dam(17).Source = 'Water';
Dam(17).Output.E = 1;
Dam(17).Size = 432; %This would be storage capacity in kilo-acre-feet
Dam(17).Enabled = 1;
Dam(17).VariableStruct.cityName = 'Kahlotus, Wa'; 
Dam(17).VariableStruct.upriver = {'Little Goose'}; 
Dam(17).VariableStruct.RM = 41.6;
Dam(17).VariableStruct.transTime = 140.1;
Dam(17).VariableStruct.MaxGenCapacity = 699070; %power in kW
Dam(17).VariableStruct.RampUp = 274; %power change in kW/hr
Dam(17).VariableStruct.RampDown = 273; 
Dam(17).VariableStruct.MaxGenFlow = 95.3;
Dam(17).VariableStruct.MaxSpillFlow = 47.1;
Dam(17).VariableStruct.MaxHead = 441.41; 

Dam(18).Type = 'Hydro';
Dam(18).Name =  {'Ice Harbor'};
Dam(18).Source = 'Water';
Dam(18).Output.E = 1;
Dam(18).Size = 249; %This would be storage capacity in kilo-acre-feet
Dam(18).Enabled = 1;
Dam(18).VariableStruct.cityName = 'Pasco, Wa';
Dam(18).VariableStruct.upriver = {'Lower Monumental'};
Dam(18).VariableStruct.RM = 9.7;
Dam(18).VariableStruct.transTime =  144.4;
Dam(18).VariableStruct.MaxGenCapacity = 581060; %power in kW
Dam(18).VariableStruct.RampUp = 609; %power change in kW/hr
Dam(18).VariableStruct.RampDown = 510;
Dam(18).VariableStruct.MaxGenFlow = 82.9;
Dam(18).VariableStruct.MaxSpillFlow = 86.3;
Dam(18).VariableStruct.MaxHead = 102.32;


nH = length(Dam);%number of dams in network
Plant.Generator = [];
Plant.Network = [];
for dam = 1:1:nH
    Dam(dam).VariableStruct.diffTime = [];
    Dam(dam).VariableStruct.upriverNum = [];
    Plant.Network(dam).name= Dam(dam).VariableStruct.cityName; %names of each node are name of closest city
    Plant.Network(dam).Equipment = Dam(dam).Name; %equiptment Hydro.Name_of_Dam
    
    Plant.Network(dam).Electrical.connections = {};
    Plant.Network(dam).Electrical.Trans_Eff = [];
    Plant.Network(dam).Electrical.Trans_Limit = [];
    
    Plant.Network(dam).Hydro.Plant_Eff = 1; %Effeciency of plant is 100% for testing
    Plant.Network(dam).Hydro.connections = {};
    Plant.Network(dam).Hydro.InstreamFlow = [];
    
    for i = 1:1:nH
        if any(strcmp(Dam(dam).Name,Dam(i).VariableStruct.upriver))%find the downstream dam
            %% Currently assume an exact overlaying Electrical network (need to adjust this, assume 99% transmission)
            Plant.Network(dam).Electrical.connections(end+1) = Dam(i).Name;
            Plant.Network(dam).Electrical.Trans_Eff(end+1) = 1;
            Plant.Network(dam).Electrical.Trans_Limit = inf;
            
            Plant.Network(dam).Hydro.connections(end+1) = Dam(i).Name;
            Plant.Network(dam).Hydro.InstreamFlow(end+1) = 0;
        end
        for k = 1:1:length(Dam(dam).VariableStruct.upriver)
            diffTime = [];
            if strcmp(Dam(dam).VariableStruct.upriver{k},Dam(i).Name)
               Dam(dam).VariableStruct.diffTime(end+1)= Dam(dam).VariableStruct.transTime - Dam(i).VariableStruct.transTime;
               Dam(dam).VariableStruct.upriverNum(end+1) = i;
            end   
        end
    end
    name = Dam(dam).VariableStruct.upriver;
    for j = 1:1:length(name)
        Plant.Network(dam).Hydro.connections(end+1) = name(j);
        Plant.Network(dam).Hydro.Trans_Eff(end+1) = 0;
        Plant.Network(dam).Hydro.Trans_Limit = 0;
    end
    Plant.Generator(dam) = Dam;
end
Plant.Network(1).Electrical.Load = 1; %put all the load at the first node
calculateHistoricalFit %might need to update this to do extra hydro specific fitting of data