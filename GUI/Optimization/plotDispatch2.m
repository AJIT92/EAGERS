function plotDispatch2(GeneratorDispatch,Time)
%% plot dispatch
global Plant Si
%can't use Si here, because each step migh not be a whole number. Instead,
%round Si to the nearest whole number

dt = Time' - [0, Time(1:end-1)]';
Time = Time+(Si-2)*dt(1);
horizon = Plant.optimoptions.Horizon;
if horizon>24
    axStep = floor(horizon/10);
elseif horizon>12
    axStep=2;
else axStep=1;
end
axTick = ceil(Time(1)-dt(1)):axStep:Time(end);
axIndex = axTick;
while axIndex(end)>24
    axIndex(axIndex>24) = axIndex(axIndex>24)-24;
end

S = Plant.optimoptions.Outputs;
for q = 1:1:length(S)
    if strcmp(S{q},'E')
        h = Plant.GUIhandles.ElectricGraph;
    elseif strcmp(S{q},'H')
        h = Plant.GUIhandles.HeatGraph;
    elseif strcmp(S{q},'C')
        h = Plant.GUIhandles.CoolGraph;
    end
    if strcmp(get(h,'Visible'),'on')
        cla(h)
        hold(h,'off')
        if strcmp(S{q},'C')
            map = 'cool';
            tSize = 9;
        elseif strcmp(S{q},'H')
            map = 'autumn';
            tSize = 9;
        else
            map = 'parula';
            tSize = 12;
        end
        c = 1;
        colorVec = [];
        while isempty(colorVec) &&  c <=length(Plant.Plotting.ColorNames)
            if strcmp(map,Plant.Plotting.ColorNames{c});
                colorVec = Plant.Plotting.ColorMaps{c};
            else c = c+1;
            end
        end
        if isempty(colorVec) 
            colorVec = Plant.Plotting.ColorMaps{1};
        end

        name = {};
        stor =[];
        plotBars =[];
        negBars =[];
        StoragePower =[];
        StorageState = [];
        nCorrespond = [];
        for i = 1:1:length(Plant.Generator)
            if isfield(Plant.Generator(i).OpMatA,'Stor') && isfield(Plant.Generator(i).OpMatA.output,S{q}) % energy storage
                stor(end+1) = i;
                StoragePower(:,end+1) = (GeneratorDispatch(1:end-1,i)-GeneratorDispatch(2:end,i))./dt;
                plotBars(:,end+1) = [0;max(0,StoragePower(:,end));];
                negBars(:,end+1) = [0;-max(0,-StoragePower(:,end));];
                name(end+1) = cellstr(Plant.Generator(i).Name);
                nCorrespond(end+1) = length(name);
                if isfield(Plant.Generator(i).OpMatA.output,'E')
                    bat = Plant.Generator(i).VariableStruct;
                    DischEff = (bat.Voltage-bat.DischResist/10)/bat.Voltage;
                    MaxDODcapacity = Plant.Generator(i).Size*(1-bat.MaxDOD/100);
                    StorageState(:,end+1) = GeneratorDispatch(:,i)/DischEff+MaxDODcapacity;
                else
                    StorageState(:,end+1) = GeneratorDispatch(:,i)/Plant.Generator(i).VariableStruct.DischargeEff;
                end
            elseif (strcmp(Plant.Generator(i).Type,'Utility')||~isempty(strfind(Plant.Generator(i).Type,'District'))) && isfield(Plant.Generator(i).OpMatA.output,S{q}) %utilities
                plotBars(:,end+1) = max(0,GeneratorDispatch(:,i));
                name(end+1) = cellstr(Plant.Generator(i).Name);
                if min(GeneratorDispatch(:,i))<0
                    negBars(:,end+1)= -min(0,GeneratorDispatch(:,i));
                    nCorrespond(end+1) = length(name);
    %                 name(end+1) = cellstr(strcat(Plant.Generator(i).Name,'Sellback'));
                end
            elseif isfield(Plant.Generator(i).OpMatA.output,S{q}) && ~(strcmp(Plant.Generator(i).Type,'Chiller') && strcmp(S{q},'E')) %generators & renewables
                plotBars(:,end+1) = GeneratorDispatch(:,i)*Plant.Generator(i).OpMatA.output.(S{q})(1); %accounts for Hratio or such when gen has multiple outputs
                name(end+1) = cellstr(Plant.Generator(i).Name);
            end
        end

        %% Plot    
        if ~isempty(plotBars)
            colorsPlot = interp1(linspace(0,1,length(colorVec)),colorVec,linspace(0,1,length(name)));

            plotTime = zeros(2*length(Time),1);
            plotTime(1:2:end) = [Time(1)-dt(1),Time(2:end)-.9999]';
            plotTime(2:2:end) = Time;
            posBars = zeros(length(plotBars(:,1))*2-2,length(plotBars(1,:)));
            posBars(1:2:end,:) = plotBars(2:end,:);
            posBars(2:2:end,:) = plotBars(2:end,:);
            h1 = area(h,plotTime,posBars,'Linestyle','none');
            for c = 1:1:length(h1)
                set(h1(c),'FaceColor',colorsPlot(c,:))
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
                h2 = area(h,plotTime,negBarsPlot,'Linestyle','none');
                for c = 1:1:length(h2)
                    set(h2(c),'FaceColor',colorsPlot(nCorrespond(c),:));
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
            legend(h1,name,'Fontsize',tSize-1,'Orientation','Horizontal','Location','NorthOutside','Box','off')
            xlim(h,[Time(1)-dt(1), Time(end)])%[0,25]
            ylim(h,[Ymin,Ymax])
            set(h,'YTick',Ymin:Yspace:Ymax,'FontSize',tSize-2)
            set(h,'XTick',axTick,'XTickLabel', {axIndex})%{xticks(1:nS+1)})
            %% Storage
            if ~isempty(stor)
                [AX, H1, H2] =plotyy(h,0,0,[0,Time],StorageState');
                set(AX,{'ycolor'},{'k';'k'})
                set(H2(end),'Color','k','LineStyle','-','LineWidth',2)
                set(H2(1),'Color','r','LineStyle','-','LineWidth',2)
                if length(H2)>2
                    for c = 1:1:length(H2)
                        set(H2(i),'Color',colorsPlot(nCorrespond(c),:),'LineStyle','-','LineWidth',2)
                    end
                end
                xlim(AX(1),[Time(1)-dt(1),Time(end)])
                xlim(AX(2),[Time(1)-dt(1),Time(end)])
                set(AX(1),'XTick',axTick,'XTickLabel',{axIndex},'FontSize',tSize-1)
                set(AX(2),'XTick',axTick,'XTickLabel',{axIndex},'FontSize',tSize-1)
                ylim(AX(1),[Ymin,Ymax])
                set(AX(1),'YTick', Ymin:Yspace:Ymax)

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
                Yspace = Ymax/pTicks;
                Ymin = Yspace*negTicks;  
                set(AX(1), 'Box', 'off')%remove yticks from left side, so you can add correct ticks
                ylim(AX(2),[Ymin,Ymax])
                set(AX(2),'YTick', Ymin:Yspace:Ymax)
                ylabel(AX(2),'State of Charge (kWh)','Color','k','FontSize',tSize)
            end
            xlabel(h,'Time (hour)','Color','k','FontSize',tSize)
            ylabel(h,'Generation (kW)','Color','k','FontSize',tSize)
            if strcmp(S{q},'E')
                h = Plant.GUIhandles.CumulativeGraph;
                if strcmp(get(h,'Visible'),'on')
                    cla(h)
                    hold(h,'off')
                    Tot = sum(plotBars(2:end,:),1);
                    if ~isempty(negBars)
                        for n = 1:1;length(negBars(1,:))
                            Tot(nCorrespond(n)) = Tot(nCorrespond(n)) - sum(negBars(:,n));
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