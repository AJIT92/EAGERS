function IC = manualInitialCondition
global Plant;
nG = length(Plant.Generator);
[~,n] = size(Plant.OneStep.organize);
nL = n-nG;
list = {};
sizes = {};
Index =[];
IC = zeros(1,nG+nL);
include = {'CHP Generator', 'Electric Generator', 'Chiller','Heater'};
for i = 1:1:length(Plant.Generator)
    if isfield(Plant.Generator(i).OpMatA,'Stor')
        list(end+1) = {strcat(Plant.Generator(i).Name,' --- % of max capacity')};
        sizes(end+1) = {'50'};
        Index(end+1) = i;
    elseif ismember(cellstr(Plant.Generator(i).Type),include)
        list(end+1) = {strcat(Plant.Generator(i).Name,' --- If greater than 0, IC must be between lower bound(',num2str(LB(i)),') and upper bound(',num2str(Plant.Generator(i).Size),').')};
        sizes(end+1) = {'0'};
        Index(end+1) = i;
    elseif ~isempty(strcmp(Plant.Generator(i).Source,'Renewable'))
        %renewable
    elseif ~isempty(strcmp(Plant.Generator(i).Type,'Utility'))
        %utility
     end
end
IC(Index) = str2double(inputdlg(list,'Specify Initial Condition (kW) or State of Charge (%)',1,sizes));
for i=1:1:nG
    if isfield(Plant.Generator(i).OpMatA,'Stor')
        IC(i) = IC(i)/100*Plant.Generator(i).OpMatA.Stor.UsableSize; % IC = halfway charged energy storage
    end
end
%% specify initial river flow and spillway flows: IC (nL)
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
for net = 1:1:length(networkNames)
    nLinet(net) = length(Plant.subNet.lineNames.(networkNames{net}));
end
nLcum = 0; %cumulative line #
for net = 1:1:length(networkNames)
    if ~strcmp(networkNames{net},'Hydro')
        nLcum = nLcum+nLinet(net); %all hydro lines have initial state
    else
        for i = 1:1:nLinet(net) 
            name = Plant.subNet.lineNames.Hydro{i};
            r = strfind(name,'_');

            if strcmp(name(r(1)+1:r(2)-1),'Hydro')
                name = name(1:r(1)-1);
                IC(nG+nLcum+i) = str2double(inputdlg(name,'Specify Initial Flow (1000 ft^3/s)',1,0));
            else
                name = name(1:r(1)-1);
                %need to find dam at this node
                IC(nG+nLcum+i) = str2double(inputdlg(name,'Specify Spillway Flow (1000 ft^3/s)',1,0));
            end
        end
    end
end