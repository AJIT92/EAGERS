function plotDispatch(GeneratorDispatch,Data,Time,D)
%% plot dispatch into GUI, with historical operation
global Plant
[m,n] = size(GeneratorDispatch);
nG = length(Plant.Generator);
stor = [];
dt = Time' - [0, Time(1:end-1)]';
MaxDODcapacity = zeros(1,nG);
for i = 1:1:nG
    if isfield(Plant.Generator(i).OpMatA,'Stor') && Plant.Generator(i).Enabled
        stor(end+1) = i;
        if isfield(Plant.Generator(i).VariableStruct,'MaxDOD')
            MaxDODcapacity(i) = Plant.Generator(i).Size*(1-Plant.Generator(i).VariableStruct.MaxDOD/100);
        end
    end
end
hoursF = D(4)+[0,Time];
[m2,n2] = size(Data);
backSteps = m2-m;
if backSteps>0
    hoursB = (hoursF(1) - backSteps*Plant.optimoptions.Resolution):Plant.optimoptions.Resolution:(hoursF(1)- Plant.optimoptions.Resolution);
    hours = [hoursB,hoursF];
else hours = hoursF;
end
StoragePower = zeros(m2,n2);
StorageState = zeros(m2,n2);
StorageState(:,stor) = Data(:,stor);
for i = 1:1:length(stor)
    StorageState(:,stor(i)) = StorageState(:,stor(i)) + MaxDODcapacity(i);
    StoragePower(backSteps+2:end,stor(i)) = (StorageState(backSteps+1:end-1,stor(i)) - StorageState(backSteps+2:end,stor(i)))./dt;  
end
StoragePower(2:backSteps+1,stor) = (StorageState(1:backSteps,stor) - StorageState(2:backSteps+1,stor))/Plant.optimoptions.Resolution; 
Data(:,stor) = StoragePower(:,stor);
nG = length(Plant.Generator);
Names = cell(nG,1);
for i = 1:1:nG
 Names(i) = {Plant.Generator(i).Name};
end
horizon = Plant.optimoptions.Horizon;

S = Plant.optimoptions.Outputs;
str = get(Plant.GUIhandles.MainTop,'String');
colorVec = Plant.Plotting.ColorMaps{1};
if get(Plant.GUIhandles.StackedGraph,'Value')==1
    LINE = false;
    str2 = 'FaceColor';
else
    LINE = true;
    str2 = 'Color';
end
for q = 1:1:length(S)
    if strcmp(S{q},'E')
        h = Plant.GUIhandles.ElectricGraph;
        primary = strcmp(str,'Electric Dispatch');
    elseif strcmp(S{q},'H')
        h = Plant.GUIhandles.HeatGraph;
        primary = strcmp(str,'Heating Dispatch');
    elseif strcmp(S{q},'C')
        h = Plant.GUIhandles.CoolGraph;
        primary = strcmp(str,'Cooling Dispatch');
    end
    if primary
        s1 = 1;
        axTick = hours(1:2:end);
        axIndex = mod(axTick,24);
        tSize = 12;
    else
        s1 = length(Data(:,1))- length(hoursF) +1;
        axTick = hours(1:4:end);
        axIndex = mod(axTick,24);
        tSize = 9;
    end
    if strcmp(get(h,'Visible'),'on')
        %the following is necessary for Matlab 2015 and earlier due to issue with legend
        pos = get(h,'Position');
        unit = get(h,'Units');
        parent = get(h,'Parent');
        delete(h);
        h = axes('Units',unit,'Position',pos,'Parent', parent);
        if strcmp(S{q},'E')
            Plant.GUIhandles.ElectricGraph = h;
        elseif strcmp(S{q},'H')
            Plant.GUIhandles.HeatGraph = h;
        elseif strcmp(S{q},'C')
            Plant.GUIhandles.CoolGraph = h;
        end
        %%%
        name = {};
        name2 ={};
        stor =[];
        plotBars =[];
        negBars =[];
        nCorrespond = [];
        for i = 1:1:length(Plant.Generator)
            include = false;
            if isfield(Plant.Generator(i).OpMatA,'Stor') && isfield(Plant.Generator(i).OpMatA.output,S{q}) && Plant.Generator(i).Enabled % energy storage
                stor(end+1) = i;
                include = true;
            elseif (strcmp(Plant.Generator(i).Type,'Utility')||~isempty(strfind(Plant.Generator(i).Type,'District'))) && isfield(Plant.Generator(i).OpMatA.output,S{q}) %utilities
                include = true;
            elseif isfield(Plant.Generator(i).OpMatA.output,S{q}) && ~(strcmp(Plant.Generator(i).Type,'Chiller') && strcmp(S{q},'E')) %generators
                include = true;
                if strcmp(S{q},'H') && isfield(Plant.Generator(i).OpMatA.output,'E')
                    Data(:,i) = Data(:,i)*Plant.Generator(i).OpMatA.output.(S{q})(1); %Hratio for CHP generators
                end
            elseif strcmp(S{q},'E') && strcmp(Plant.Generator(i).Source,'Renewable')
                include = true;
            end
            if include
                name(end+1) = {Plant.Generator(i).Name};
                plotBars(:,end+1) = max(0,Data(s1:end,i));
                if any(Data(:,i)<0)
                    negBars(:,end+1) = -max(0,-Data(s1:end,i));
                    nCorrespond(end+1) = i;
                end
                if isfield(Plant.Generator(i).OpMatA,'Stor') && isfield(Plant.Generator(i).OpMatA.output,S{q}) % energy storage
                    name2(end+1) ={Plant.Generator(i).Name};
                end
            end
        end

        %% Plot    
        if ~isempty(plotBars)
            colorsPlot = interp1(linspace(0,1,length(colorVec)),colorVec,linspace(0,1,length(Names)));
            plotTime = zeros(2*length(hours(s1:end))-2,1);
            plotTime(1:2:2*length(hours(s1:end))-2) = [hours(s1),hours(s1+2:end)-.9999]';
            plotTime(2:2:2*length(hours(s1:end))-2) = hours(s1+1:end);
            posBars = zeros(length(plotBars(:,1))*2-2,length(plotBars(1,:)));
            posBars(1:2:end,:) = plotBars(2:end,:);
            posBars(2:2:end,:) = plotBars(2:end,:);
            if LINE
                h1 = plot(h,hours(s1:end),plotBars,'LineWidth',3);
            else
                h1 = area(h,plotTime,posBars,'Linestyle','none');
            end
            for c = 1:1:length(h1)
                I = find(strcmp(name(c),Names));
                set(h1(c),str2,colorsPlot(I,:));
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
                    set(h2(c),str2,colorsPlot(nCorrespond(c),:));
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
            pTicks = Ymax/Yspace;
            
            if primary
                s0 = length(Data(:,1))- length(hoursF) +1;
                plot(h,[hours(s0),hours(s0)],[Ymin,Ymax],'c--')
                xlim(h,[D(4)-horizon, hours(end)])
            else
                xlim(h,[D(4), hours(end)])
            end
            ylim(h,[Ymin,Ymax])
            set(h,'YTick',Ymin:Yspace:Ymax,'FontSize',tSize-2)
            set(h,'XTick',axTick,'XTickLabel', {axIndex})
            %% Storage
            if ~isempty(stor) && ~LINE
                [AX, H1, H2] =plotyy(h,0,0,hours(s1:end),StorageState(s1:end,stor));
                set(AX,{'ycolor'},{'k';'k'})
                for c = 1:1:length(H2)
                    set(H2(c),'Color',colorsPlot(nCorrespond(c),:),'LineStyle','-','LineWidth',2,'Marker','x','MarkerEdgeColor','k','MarkerSize',5)
                end
                if primary
                    xlim(AX(1),[D(4)-horizon, hours(end)])
                    xlim(AX(2),[D(4)-horizon, hours(end)])
                else
                    xlim(AX(1),[D(4), hours(end)])
                    xlim(AX(2),[D(4), hours(end)])
                end
                set(AX(1),'XTick',axTick,'XTickLabel',{axIndex},'FontSize',tSize-1)
                set(AX(2),'XTick',axTick,'XTickLabel',{axIndex},'FontSize',tSize-1)
                ylim(AX(1),[Ymin,Ymax])
                set(AX(1),'YTick', Ymin:Yspace:Ymax)

                OoM = log10(max(max(StorageState(s1:end,stor))));
                if (OoM-floor(OoM))==0 %count in increments of 1, 10, 100 or 1000 etc
                    Ymax = 10^OoM;
                elseif (OoM-floor(OoM))> 0.6990 %count in increments of 1, 10, 100 or 1000 etc
                    Ymax = 10^ceil(OoM);
                elseif (OoM-floor(OoM))> 0.30103 %count in increments of 5, 50, 500 or 5000 etc
                    Ymax = .5*10^ceil(OoM);
                else  %count in increments of 2, 20, 200 or 2000 etc
                    Ymax = .2*10^ceil(OoM);
                end
                Yspace = Ymax/pTicks;
                Ymin = Yspace*negTicks;  
                set(AX(1), 'Box', 'off')%remove yticks from left side, so you can add correct ticks
                ylim(AX(2),[Ymin,Ymax])
                set(AX(2),'YTick', Ymin:Yspace:Ymax)
                ylabel(AX(2),'State of Charge (kWh)','Color','k','FontSize',tSize)
                if ~isempty(name2)
                    legend(H2,name2,'Fontsize',tSize-1,'Orientation','Horizontal','Location','North','Box','off')
                end
            end
            xlabel(h,'Time (hour)','Color','k','FontSize',tSize)
            ylabel(h,'Generation (kW)','Color','k','FontSize',tSize)
            if strcmp(S{q},'E')
                h = Plant.GUIhandles.CumulativeGraph;
                if strcmp(get(h,'Visible'),'on')
                    cla(h)
                    hold(h,'off')
                    s2 = length(Data(:,1))- length(hoursF) +1;
                    Tot = sum(plotBars(s2+1:end,:),1);
                    if ~isempty(negBars)
                        for n = 1:1:length(negBars(1,:))
                            Tot(nCorrespond(n)) = Tot(nCorrespond(n)) - sum(negBars(s2:end,n));
                        end
                    end
                    Tot = Tot/sum(Tot)*100;%convert to percent of net generation
                    h1 = area(h,[0,1],[Tot;Tot;],'Linestyle','none');
                    for c = 1:1:length(h1)
                        set(h1(c),'FaceColor',colorsPlot(c,:))
                    end
                    ylim(h,[0,100])
                    ylabel(h,'Percent of Generation (%)','Color','k','FontSize',tSize-1)
                    set(h,'YTick',0:20:100,'FontSize',tSize-2)
                end
            end
            
        end
    end
end