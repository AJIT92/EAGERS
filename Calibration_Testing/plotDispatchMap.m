function plotDispatchMap(Generators,OptimalOutput,costMin,Range,Utility,GenType,ElecOrCool)
a = size(OptimalOutput);
netGeneration = zeros(a(1),1);
netSellBack = netGeneration;
name = {};
e = [];
i=1;
if strcmp(ElecOrCool,'Electric') || strcmp(ElecOrCool,'Both')
    while i <=length(GenType) %add up only electric generators
        if GenType(i)==2
            netGeneration = netGeneration+ OptimalOutput(:,i);
            name(end+1) = cellstr(Generators(i).Name);
            e(end+1) = i;
        elseif GenType(i) ==1
            name(end+1) = cellstr('Electric Grid');
            netGeneration = netGeneration+ OptimalOutput(:,i);
            e(end+1) = i;
            if max(Utility.SellBack)>0
                netSellBack = -OptimalOutput(:,i+1);
            end
            i = i+1; %handle both purchase & sellback states in one go
        end
        i=i+1;
    end
end
if strcmp(ElecOrCool,'Cooling') || strcmp(ElecOrCool, 'Both')
    i = 1;
    while i <=length(GenType) %add up only electric generators
        if GenType(i)==3
            netGeneration = netGeneration+ OptimalOutput(:,i);
            name(end+1) = cellstr(Generators(i).Name);
            e(end+1) = i;
        elseif GenType(i) ==9
            name(end+1) = cellstr('District Cooling');
            netGeneration = netGeneration+ max(0,OptimalOutput(:,i));
            e(end+1) = i;
            netSellBack = -max(0,-OptimalOutput(:,i));
            i = i+1; %handle both purchase & sellback states in one go
        end
        i=i+1;
    end
end
figure(1)
bar(Range,OptimalOutput(:,e),'stacked')
hold on
costkWh = costMin./netGeneration;
[AX, H1, H2] = plotyy(0,0,Range,costkWh);
set(AX,{'ycolor'},{'k';'k'})
set(H2,'Color','k','LineStyle','--','LineWidth',2)
set(H1,'Color','k','LineStyle','--','LineWidth',2)
name(end+1) = cellstr('Cost per kWh');
xlabel('Demand (kW)')
set(get(AX(1),'Ylabel'),'String','Net Generation (kW)')
set(get(AX(2),'Ylabel'),'String','Operating Cost ($/kWh of generation)')
ylim(AX(2),[0,.1])
set(AX(2),'YTick',0:.01:.1)
xlim(AX(1),[0,max(Range)])
xlim(AX(2),[0,max(Range)])
set(AX(1),'YTick',0:50:max(Range))
if max(abs(netSellBack))>0
    bar(AX(1),Range,netSellBack,'m')
    l = -100;%round(1.25*min(netSellBack)/10)*10;
    u = 400;%round(1.25*max(Range)/10)*10;
    ylim(AX(1),[l,u])
    set(AX(1),'YTick',l:50:u)
    l = round(-1.25*abs(min(netSellBack))/max(Range)*max(costkWh)*100)/100;
    u = .1;%round(1.25*max(costkWh)*100)/100;
    ylim(AX(2),[l,u])
    set(AX(2),'YTick',l:.01:u)
    if strcmp(ElecOrCool,'Electric')
        name(end+1) = cellstr('Grid Sellback');
    else name(end+1) = cellstr('DC Sellback');
    end
end
legend(name)