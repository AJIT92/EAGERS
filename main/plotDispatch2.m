function plotDispatch2(GeneratorDispatch,Demand,Time,fig)
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

if isfield(Demand,'T')
    Demand = rmfield(Demand,'T');
end
S = fieldnames(Demand);
for q = 1:1:length(S)
    name = {};
    name2 = {};
    stor =[];
    plotBars =[];
    negBars =[];
    StoragePower =[];
    StorageState = [];
    for i = 1:1:length(Plant.Generator)
        if isfield(Plant.Generator(i).OpMatA,'Stor') && isfield(Plant.Generator(i).OpMatA.output,S{q}) % energy storage
            stor(end+1) = i;
            StoragePower(:,end+1) = (GeneratorDispatch(1:end-1,i)-GeneratorDispatch(2:end,i))./dt;
            plotBars(:,end+1) = [0;max(0,StoragePower(:,end));];
            negBars(:,end+1) = -[0;max(0,-StoragePower(:,end));];
            name(end+1) = {['Discharge ',Plant.Generator(i).Name]};
            name2(end+1) = {['Charge ',Plant.Generator(i).Name]};
            if isfield(Plant.Generator(i).OpMatA.output,'E')
                bat = Plant.Generator(i).VariableStruct;
                DischEff = (bat.Voltage-bat.DischResist/10)/bat.Voltage;
                MaxDODcapacity = Plant.Generator(i).Size*(1-bat.MaxDOD/100);
                StorageState(:,end+1) = GeneratorDispatch(:,i)/DischEff+MaxDODcapacity;
            else
                StorageState(:,end+1) = GeneratorDispatch(:,i)/Plant.Generator(i).VariableStruct.DischargeEff;
            end
        elseif strcmp(Plant.Generator(i).Type,'Utility') && isfield(Plant.Generator(i).OpMatA.output,S{q})
            plotBars(:,end+1) = max(0,GeneratorDispatch(:,i));
            name(end+1) = cellstr(Plant.Generator(i).Name);
            if min(GeneratorDispatch(:,i))<0
                negBars(:,end+1)= min(0,GeneratorDispatch(:,i));
                name2(end+1) = cellstr('Sellback');
            end
        elseif isfield(Plant.Generator(i).OpMatA.output,S{q}) && ~(strcmp(Plant.Generator(i).Type,'Chiller') && strcmp(S{q},'E')) %generators & renewables
            plotBars(:,end+1) = GeneratorDispatch(:,i)*Plant.Generator(i).OpMatA.output.(S{q})(1); %accounts for Hratio or such when gen has multiple outputs
            name(end+1) = cellstr(Plant.Generator(i).Name);
        end
    end
    
    %% Plot
    if ~isempty(plotBars)
%         clf(fig+q-1)
        figure(fig+q-1)
        hold off
        bar(Time,plotBars(2:end,:),'stacked','barwidth',1)
        hold on
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
            h = bar(Time,negBars(2:end,:),'stacked','barwidth',1);
            negTicks = floor(min(min(negBars))/Yspace);
            if abs(negTicks)>3
                Yspace = 2*Yspace;
                negTicks = floor(min(min(negBars))/Yspace);
            end
            set(h(length(h)),'FaceColor','c')
            set(h(1), 'FaceColor', 'm')
        end
        Ymin = Yspace*negTicks;
        pTicks = Ymax/Yspace;
        legend([name name2],'FontSize',16)
        xlim([Time(1)-dt(1), Time(end)])%[0,25]
        ylim([Ymin,Ymax])
        set(gca,'YTick',Ymin:Yspace:Ymax,'FontSize',14)
        set(gca, 'XTick',axTick,'XTickLabel', {axIndex})%{xticks(1:nS+1)})
        %% Storage
        if ~isempty(stor)
            [AX, H1, H2] =plotyy(0,0,Time,StorageState(2:end,:)');
            set(AX,{'ycolor'},{'k';'k'})
            set(H2(end),'Color','k','LineStyle','--','LineWidth',2)
            set(H2(1),'Color','r','LineStyle','--','LineWidth',2)
            xlim(AX(1),[Time(1)-dt(1),Time(end)])
            xlim(AX(2),[Time(1)-dt(1),Time(end)])
            set(AX(1),'XTick',axTick,'XTickLabel',{axIndex},'FontSize',16)
            set(AX(2),'XTick',axTick,'XTickLabel',{axIndex},'FontSize',16)
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
            ylim(AX(2),[Ymin,Ymax])
            set(AX(2),'YTick', Ymin:Yspace:Ymax)
            ylabel(AX(2),'State of Charge (kWh)','Color','k','FontSize',16)
        end
        xlabel('Time (hour)','Color','k','FontSize',16)
        ylabel('Generation (kW)','Color','k','FontSize',16)
        plotBars=[];
    end
end