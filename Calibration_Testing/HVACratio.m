global Model_dir MapHandle
MapHandle =figure(1);

stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';}; 
buildTypeName = {'SmOff'; 'FFRest';'ware'; 'MRapt';'SDRest';'StMall';'Retail';'SmHotel';'MdOff'; 'Sch-pri';'OutP';'SMarket'; 'LgHotel';'Sch-sec'; 'LgOff';     'Hospital';     };
load('USbuildData.mat');
HVAC_ratio = zeros(50,1);
ratio = zeros(50,16);
h = waitbar(0,'Initializing');
for state = 1:1:50
    for j = 1:1:length(buildTypeName)
        name = strcat(char(buildTypeName(j)),'_',char(ClimateZone(stateAbrev(state),0 )),'_','New2010');
        load(fullfile(Model_dir,'System Library','Buildings',name))
        ratio(state,j) = sum(component.CoolingElectricalLoad)/sum(component.DemandE);
        waitbar((16*(state-1)+j)/800,h,strcat('Building # ',num2str(16*(state-1)+j)));
    end
    HVAC_ratio(state) = sum(buildData(state,:).*ratio(state,:))/sum(buildData(state,:));
end
HVAC_ratioBuild = zeros(16,1);
for i = 1:1:length(buildTypeName)
    HVAC_ratioBuild(i) = sum(buildData(:,i).*ratio(:,i))/sum(buildData(:,i));
end
close(h)