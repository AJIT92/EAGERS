function ScaleSystemComponents(DG,Cool,TESsize,BatSizeMin,Solar,ScaleSol,Wind,ScaleWind,Iterations)
global Project
MSG = msgbox('Re-sizing CHP System');
% Remove Battery, TES, and Renewables

if isfield(Project.System,'TES')
    TES = Project.System.TES;
    Project.System = rmfield(Project.System,'TES');
else TES = 1;
end
if isfield(Project.System,'Battery')
    Battery = Project.System.Battery;
    Project.System = rmfield(Project.System,'Battery');
else Battery = 1;
end
%determine if there are renewables to re-size
if ~isfield(Project,'Renewable')
    ScaleSol=0;
    ScaleWind=0;
elseif ~isfield(Project.Renewable,'Solar')
    ScaleSol=0;
elseif ~isfield(Project.Renewable,'Wind')
    ScaleWind=0;
end
% If renewables are demand based size them here
if ScaleSol==1%%%%when is ScaleSol initialized???
    SolarSizing(1,'demand','all',Solar)
end
if ScaleWind==1
    WindSizing(1,'demand','all',Wind)
end
% if chillers are present temporarily remove
if isfield(Project.System,'Chiller')
    Chillers = Project.System.Chiller;
    Project.System = rmfield(Project.System,'Chiller');
else Chillers = 1;
end
%Set new Stationary generator size
Project.System.CHP(1).OptSize = DG;
OptimalFCsizing(1,'OptSize');

% Set new chiller size
if isstruct(Chillers)
    Project.System.Chiller = Chillers;
    for chill =1:1:length(Project.System.Chiller)
        Project.System.Chiller(chill).OptSize = Cool;
    end
    OptimalChillerSizing('OptSize');
end
%Set new Stationary generator sizes
for chp =1:1:length(Project.System.CHP)
    Project.System.CHP(chp).OptSize = DG;
    OptimalFCsizing(chp,'OptSize');
end
% Add TES, Battery and Renewables back and set new sizes of each
if isstruct(TES)
    Project.System.TES = TES;
    for tes =1:1:length(Project.System.TES)
        Project.System.TES(tes).OptSize = TESsize;
    end
    OptimalTESsizing('OptSize');
end
if isstruct(Battery)
    Project.System.Battery = Battery;
    oldSizeMinutes = zeros(length(Battery),1);
    for bat =1:1:length(Project.System.Battery)
        oldSizeMinutes(bat) = Project.System.Battery(bat).SizeMin;
    end
    totalSizeMinutes= sum(oldSizeMinutes);
    for bat =1:1:length(Project.System.Battery)
        Project.System.Battery(bat).SizeMin = oldSizeMinutes(bat)*(BatSizeMin/totalSizeMinutes);
    end
    OptimalBatterySizing('editSizeMin')
end

%% Iterate sizes to reach solution
for repeat = 1:1:Iterations
    OptimalFCsizing(1,'OptSize');
    if ScaleSol==2
        SolarSizing(1,'generation','all',Solar)
    end
    if  ScaleWind==2
        WindSizing(1,'generation','all',Wind)
    end
    if isfield(Project.System,'Chiller')
        OptimalChillerSizing('OptSize');
    end
    if isfield(Project.System,'TES')
        OptimalTESsizing('OptSize');
    end
    if isfield(Project.System,'Battery')
        OptimalBatterySizing('editSizeMin')
    end
end
for chp =1:1:length(Project.System.CHP)
    OptimalFCsizing(chp,'OptSize');
end

CHPSize = 0;
for i = 1:1:length(Project.System.CHP)
    CHPSize = CHPSize+Project.System.CHP(i).SysSize(1);
end
Project.Economic.LifekWh = Project.Economic.LifeYrs*CHPSize*8760;

% if isdeployed
% else
%     discreteSizes();
% end

close(MSG)