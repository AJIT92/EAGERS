function updateGUIstatus(GenDisp,History)
global Plant Model_dir SYSINDEX
nG = length(Plant.Generator);
%% update status lights
for i = 1:1:nG
    x = [];
    if Plant.Generator(i).Enabled
        if GenDisp(2,i)>0 && GenDisp(1,i)==0  %just turned on
            [x,~] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
        end
        if (GenDisp(2,i)==0 && GenDisp(1,i)>0) || (isempty(History) && GenDisp(2,i)==0) %just turned off or 1st time running
        	[x,~] = imread(fullfile(Model_dir,'GUI','Graphics','yellow.png'));
        end
    end
    if ~isempty(x)
        set(Plant.GUIhandles.Switch,'Units','pixels');
        pos1 = get(Plant.GUIhandles.Switch,'Position');
        set(Plant.GUIhandles.Switch,'Units','characters');
        pos2 = get(Plant.GUIhandles.Switch,'Position');
        s = imresize(x,[3*pos1(3)/pos2(3) pos1(4)/pos2(4)]);
        set(Plant.GUIhandles.(strcat('GeneratorStat',num2str(i))),'cdata', s);
    end
end

%% Update status of selected generator in lower left box
if ~isempty(strfind(Plant.Generator(SYSINDEX).Type,'Storage'))
    power = (GenDisp(1,SYSINDEX)- GenDisp(2,SYSINDEX))/Plant.optimoptions.Resolution*Plant.Generator(SYSINDEX).OpMatA.Stor.DischEff;
    set(Plant.GUIhandles.GenStatus1,'string', num2str(power));
    set(Plant.GUIhandles.GenStatus2text,'string','State-Of-Charge (%)');
    SOC = GenDisp(2,SYSINDEX)/Plant.Generator(SYSINDEX).OpMatA.Stor.UsableSize*100;
    set(Plant.GUIhandles.GenStatus2,'string', num2str(SOC));
else
    set(Plant.GUIhandles.GenStatus1,'string', num2str(GenDisp(2,SYSINDEX)));
    set(Plant.GUIhandles.GenStatus2text,'string','Efficiency (%)');
    skip = false;
    if ~isempty(Plant.Generator(SYSINDEX).Output)
        cap = Plant.Generator(SYSINDEX).Output.Capacity*Plant.Generator(SYSINDEX).Size;
    end
    if strcmp(Plant.Generator(SYSINDEX).Type,'Electric Generator') || strcmp(Plant.Generator(SYSINDEX).Type,'CHP Generator')
        eff = Plant.Generator(SYSINDEX).Output.Electricity;
    elseif strcmp(Plant.Generator(SYSINDEX).Type,'Chiller') 
        eff = Plant.Generator(SYSINDEX).Output.Cooling;
    elseif strcmp(Plant.Generator(SYSINDEX).Type,'Heater')
        eff = Plant.Generator(SYSINDEX).Output.Heat;    
    elseif strcmp(Plant.Generator(SYSINDEX).Type,'Hydro')
        skip = true;
    else skip = true;
    end
    if ~skip
        Efficiency = interp1(cap,eff,GenDisp(2,SYSINDEX))*100;
        set(Plant.GUIhandles.GenStatus2,'string', num2str(Efficiency));
    end
end