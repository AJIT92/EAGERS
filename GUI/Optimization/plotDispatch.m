function plotDispatch(Timestamp,GenDisp,History)
%% plot dispatch into GUI, with historical operation
global Plant
nG = length(Plant.Generator);
stor = [];
Names = cell(nG,1);
for i = 1:1:nG
    Names(i) = {Plant.Generator(i).Name};
    if ~isempty(strfind(Plant.Generator(i).Type,'Storage')) && Plant.Generator(i).Enabled
        stor(end+1) = i;
    end
end
horizon = Plant.optimoptions.Horizon;

Time = (Timestamp-Timestamp(1))*24;
D = datevec(Timestamp(1));
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
nPlot = length(networkNames);

hoursF = D(4)+D(5)/60+D(6)/3600+Time;


%% Collate history and future data
if isempty(History)
    backSteps = 0;
    Data = GenDisp;
else
    backSteps = length(History(:,1))-1;
    Data = [History;GenDisp(2:end,:)];
end
dt = [Plant.optimoptions.Resolution*ones(backSteps,1); Time(2:end) - Time(1:end-1);];
if backSteps>0
    hoursB = ((hoursF(1) - backSteps*Plant.optimoptions.Resolution):Plant.optimoptions.Resolution:(hoursF(1)- Plant.optimoptions.Resolution))';
    hours = [hoursB;hoursF];
else hours = hoursF;
end

%% Make text strings to scroll across bottom axis
Dprev = datevec(Timestamp(1)-1);
Dnext = datevec(Timestamp(1)+1);
Dtimestamp = Timestamp(1)+hours/24;
months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Aug','Nov','Dec'};
dateText = ' ';
if backSteps<24/Plant.optimoptions.Resolution
    b = ceil(30*(24/Plant.optimoptions.Resolution-backSteps)/(24/Plant.optimoptions.Resolution));
    for i = 1:1:b
        dateText = strcat(dateText,{' '});
    end
end
if nnz(Dtimestamp<datenum([D(1),D(2),D(3)]))>0.3*length(hours) %more than 30% of window is previous day
    dateText = strcat(dateText,months(Dprev(2)),{' '},{num2str(Dprev(3))},{'  '},{num2str(Dprev(1))},{'                                '});
else  %append spaces
    c = ceil(30*nnz(Dtimestamp<datenum([D(1),D(2),D(3)]))/length(hours));
    for i = 1:1:c
        dateText = strcat(dateText,{' '});
    end
end
dateText = strcat(dateText,months(D(2)),{' '},{num2str(D(3))},{'  '},{num2str(D(1))});
if nnz(Dtimestamp>datenum([Dnext(1),Dnext(2),Dnext(3)]))>0.3*length(hours) %more than 30% of window is next day
    dateText = strcat(dateText,{'                         '},months(Dnext(2)),{' '},{num2str(Dnext(3))},{'  '},{num2str(Dnext(1))});
else  %append spaces
        dateText = strcat(dateText,{'                                      '});
end
dateText2 = strcat('  ',months(D(2)),{' '},{num2str(D(3))},{'  '},{num2str(D(1))});
if nnz(Dtimestamp>datenum([Dnext(1),Dnext(2),Dnext(3)]))>0.3*length(hours) %more than 30% of window is next day
    dateText2 = strcat(dateText2,{'                      '},months(Dnext(2)),{' '},{num2str(Dnext(3))},{'  '},{num2str(Dnext(1))});
else
    dateText2 = strcat(dateText2,{'                                '});
end

%% convert the saved SOC to power
StoragePower = 0*Data;
StorageState = 0*Data;
for i = 1:1:length(stor)
    StorageState(:,stor(i)) = Data(:,stor(i))+ ones(length(hours),1)*(Plant.Generator(stor(i)).OpMatA.Stor.Size - Plant.Generator(stor(i)).OpMatA.Stor.UsableSize); %add the unusable charge
    if strcmp(Plant.Generator(stor(i)).Type,'Hydro')
        %% Need the river segment and spill flow to calculate power
    else
        StoragePower(2:end,stor(i)) = (StorageState(1:end-1,stor(i)) - StorageState(2:end,stor(i)))./dt;  
    end
end
Data(:,stor) = StoragePower(:,stor);

%% Do actual Plotting
colorVec = Plant.Plotting.ColorMaps{1};
colorsPlot = interp1(linspace(0,1,length(colorVec)),colorVec,linspace(0,1,length(Names)));
if get(Plant.GUIhandles.StackedGraph,'Value')==1
    LINE = false;
else
    LINE = true;
end
axTick = (ceil(hours(1)):round((hours(end)-hours(1))/12):hours(end));
axIndex = mod(axTick,24);
axIndex([false,axIndex(2:end)==0]) = 24;
for q = 1:1:nPlot
    h = Plant.GUIhandles.(strcat('ResultPlot',num2str(q)));
    if q==1
        s1 = 1;
        tSize = 12;
    else
        s1 = length(Data(:,1))- length(hoursF) +1;
        tSize = 9;
    end
    if strcmp(get(h,'Visible'),'on')
        [name,name2,stor,plotBars,negBars,negBarIndex,storIndex] = sortForPlot(h,Data,s1,num2str(q));
        h = Plant.GUIhandles.(strcat('ResultPlot',num2str(q)));
        hold(h,'off');
        %% Plot    
        if ~isempty(plotBars)
            dataPlot(h,LINE,plotBars,negBars,hours,horizon,dt,s1,name,Names,colorsPlot,negBarIndex,D,tSize,axTick,axIndex)
            %% Storage
            if ~isempty(stor) && ~LINE
                addStorage2Plot(h,StorageState(s1:end,stor),hours,s1,colorsPlot,storIndex,name2,tSize,axTick,axIndex)
            end
            if q==1
                xlabel(h,dateText,'Color','k','FontSize',tSize)
            else xlabel(h,dateText2,'Color','k','FontSize',tSize)
            end
            ylabel(h,'Generation (kW)','Color','k','FontSize',tSize)  
        end
    end
end

function [name,name2,stor,plotBars,negBars,negBarIndex,storIndex] = sortForPlot(h,Data,s1,q)
global Plant
S = get(Plant.GUIhandles.(strcat('ResultName',q)),'String');
if strcmp(S,'Electrical')
    S = 'E';
elseif strcmp(S,'DistrictHeat')
    S = 'H';
elseif strcmp(S,'DistrictCool')
    S = 'C';
elseif strcmp(S,'Hydro')
    S = 'W';
elseif strcmp(S,'Steam')
    S = 'S';
end

name = {};
name2 ={};
stor =[];
plotBars =[];
negBars =[];
negBarIndex = [];
storIndex = [];
 %the following is necessary for Matlab 2015 and earlier due to issue with legend
pos = get(h,'Position');
delete(h);
h =  axes('Units','characters','Position', pos,'Tag', strcat('ResultPlot',q),'Parent', Plant.GUIhandles.uipanelMain1,'Visible','on');
Plant.GUIhandles.(strcat('ResultPlot',q)) = h;

for i = 1:1:length(Plant.Generator)
    include = false;
    dataI = Data(:,i);
    if strcmp(S,'W')
        %water flow is in lines, plot lines does this
    elseif isfield(Plant.Generator(i).OpMatA,'Stor') && isfield(Plant.Generator(i).OpMatA.output,S) && Plant.Generator(i).Enabled % energy storage
        stor(end+1) = i;
        include = true;
    elseif (strcmp(Plant.Generator(i).Type,'Utility')||~isempty(strfind(Plant.Generator(i).Type,'District'))) && isfield(Plant.Generator(i).OpMatA.output,S) %utilities
        include = true;
    elseif isfield(Plant.Generator(i).OpMatA.output,S) && ~(strcmp(Plant.Generator(i).Type,'Chiller') && strcmp(S,'E')) %generators
        include = true;
        if strcmp(S,'H') && isfield(Plant.Generator(i).OpMatA.output,'E')
            dataI = dataI*Plant.Generator(i).OpMatA.output.(S)(1); %Hratio for CHP generators
        end
    elseif strcmp(S,'E') && strcmp(Plant.Generator(i).Source,'Renewable')
        include = true;
    end
    if include
        name(end+1) = {Plant.Generator(i).Name};
        plotBars(:,end+1) = max(0,dataI(s1:end));
        if any(dataI<0) || isfield(Plant.Generator(i).OpMatA,'Stor')
            negBars(:,end+1) = -max(0,-dataI(s1:end));
            negBarIndex(end+1) = i;
        end
        if isfield(Plant.Generator(i).OpMatA,'Stor')% energy storage
            name2(end+1) ={Plant.Generator(i).Name};
            storIndex(end+1) = i;
        end
    end
end

function dataPlot(h,LINE,plotBars,negBars,hours,horizon,dt,s1,name,Names,colorsPlot,negBarIndex,D,tSize,axTick,axIndex)
hourNow = D(4) + D(5)/60 + D(6)/3600;
if strcmp(get(h,'Tag'),'ResultPlot1')
    s0 = nnz(hours<hourNow) +1;
else s0 = [];
end
plotTime = zeros(2*length(hours(s1:end))-2,1);
plotTime(1:2:2*length(hours(s1:end))-2) = [hours(s1);hours(s1+2:end)-.9999*dt(s1+1:end)];
plotTime(2:2:2*length(hours(s1:end))-2) = hours(s1+1:end);
posBars = zeros(length(plotBars(:,1))*2-2,length(plotBars(1,:)));
posBars(1:2:end,:) = plotBars(2:end,:);
posBars(2:2:end,:) = plotBars(2:end,:);
if LINE
    str2 = 'Color';
    h1 = plot(h,hours(s1:end),plotBars,'LineWidth',3);
else
    str2 = 'FaceColor';
    h1 = area(h,plotTime,posBars,'Linestyle','none');
end
for c = 1:1:length(h1)
    set(h1(c),str2,colorsPlot(strcmp(name(c),Names),:));
end
hold(h,'on');

if ~isempty(negBars)
    OoM = log10(max(sum(plotBars,2)+sum(negBars,2)));
else OoM = log10(max(sum(plotBars,2)));
end
if (OoM-floor(OoM))==0 %count in increments of 1, 10, 100 or 1000 etc
    Yspace = 10^(OoM-1);
    Ymax = 10^OoM;
elseif (OoM-floor(OoM))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
    Yspace = 10^floor(OoM);
    Ymax = 10^ceil(OoM);
elseif (OoM-floor(OoM))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
    Yspace = .5*10^floor(OoM);
    Ymax = .5*10^ceil(OoM);
else  %count in increments of 2, 20, 200 or 2000 etc
    Yspace = .2*10^floor(OoM);
    Ymax = .2*10^ceil(OoM);
end
negTicks = 0;
if ~isempty(negBars)
    negBarsPlot = zeros(length(negBars(:,1))*2-2,length(negBars(1,:)));
    negBarsPlot(1:2:end,:) = negBars(2:end,:);
    negBarsPlot(2:2:end,:) = negBars(2:end,:);
    if LINE
        h2 = plot(h,hours(s1:end),negBars,'LineWidth',3);
    else
        h2 = area(h,plotTime,negBarsPlot,'Linestyle','none');
    end
    for c = 1:1:length(h2)
        set(h2(c),str2,colorsPlot(negBarIndex(c),:));
    end
    negTicks = floor(min(min(negBars))/Yspace);
    if abs(negTicks)>3
        Yspace = 2*Yspace;
        negTicks = floor(min(min(negBars))/Yspace);
    end
end
Ymin = Yspace*negTicks;
if isempty(Ymin)
    Ymin = 0;
end
if ~isempty(s0)
    plot(h,[hours(s0),hours(s0)],[Ymin,Ymax],'c--')
    xlim(h,[hours(s0)-horizon, hours(end)])
else
    xlim(h,[hourNow, hours(end)])
end
ylim(h,[Ymin,Ymax])
set(h,'YTick',Ymin:Yspace:Ymax,'FontSize',tSize-2)
set(h,'XTick',axTick,'XTickLabel', {axIndex})

function addStorage2Plot(h,StorageState,hours,s1,colorsPlot,storIndex,name2,tSize,axTick,axIndex)
Xlim = get(h,'Xlim');
ticks = get(h,'YTick');
[AX, H1, H2] =plotyy(h,0,0,hours(s1:end),StorageState);
set(AX,{'ycolor'},{'k';'k'})
for c = 1:1:length(H2)
    set(H2(c),'Color',colorsPlot(storIndex(c),:),'LineStyle','-','LineWidth',2,'Marker','x','MarkerEdgeColor','k','MarkerSize',5)
end
xlim(AX(1),Xlim)
xlim(AX(2),Xlim)
set(AX(1),'XTick',axTick,'XTickLabel',{axIndex},'FontSize',tSize-1)
set(AX(2),'XTick',axTick,'XTickLabel',{axIndex},'FontSize',tSize-1)

ylim(AX(1),[min(ticks),max(ticks)])
set(AX(1),'YTick', ticks)

OoM = log10(max(max(StorageState)));
if (OoM-floor(OoM))==0 %count in increments of 1, 10, 100 or 1000 etc
    Ymax = 10^OoM;
elseif (OoM-floor(OoM))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
    Ymax = 10^ceil(OoM);
elseif (OoM-floor(OoM))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
    Ymax = .5*10^ceil(OoM);
else  %count in increments of 2, 20, 200 or 2000 etc
    Ymax = .2*10^ceil(OoM);
end
pTicks = nnz(get(h,'YTick')>0); % # of positive tick marks
Yspace = Ymax/pTicks;
negTicks = min(ticks)/(ticks(end)-ticks(end-1));
Ymin = Yspace*negTicks;  
set(AX(1), 'Box', 'off')%remove yticks from left side, so you can add correct ticks
ylim(AX(2),[Ymin,Ymax])
set(AX(2),'YTick', Ymin:Yspace:Ymax)
ylabel(AX(2),'State of Charge (kWh)','Color','k','FontSize',tSize)
if ~isempty(name2)
    legend(H2,name2,'Fontsize',tSize-1,'Orientation','Horizontal','Location','North','Box','off')
end