function EditSystem(handles)
%EditSystem(handles)

global Plant SYSINDEX

uiElmts = {'textEdit1';'textEdit1Units';'compText1';'textEdit2'; ...
    'textEdit2Units';'compText2';'textEdit3';'compText3';'textEdit4'; ...
    'compText4';'textFuel';'CompFuel';'textEff';'uitableEffCurve'; ...
    'EffCurve';'textResponse';'ResponseRate';};
if SYSINDEX > 0
    Gen = Plant.Generator(SYSINDEX);
else
    switch SYSINDEX
        case 0
            Gen = struct('Type', 'None', ...
                'Name', 'None');
        case -1
            Gen = struct('Type', 'Heating Demands', ...
                'Name', 'Heating Demands', ...
                'Demand', 100);
        case -2
            Gen = struct('Type', 'Hot Water Demands', ...
                'Name', 'Hot Water Demands', ...
                'Demand', 100);
        case -3
            Gen = struct('Type', 'Cooling Demands', ...
                'Name', 'Cooling Demands', ...
                'Demand', 100);
        case -4
            Gen = struct('Type', 'AC/DC Conversion', ...
                'Name', 'AC/DC Conversion', ...
                'Efficiency', 0.5);
    end
end
set(handles.CompName,'String',Gen.Name);
switch Gen.Type
    case {'CHP Generator';'ElectricGenerator'}
        show = {'textEdit1';'textEdit1Units';'compText1'; ...
            'uitableEffCurve';'EffCurve';'ResponseRate';'textResponse';};
        hide = setdiff(uiElmts,show); % everything else not in show
        for i = 1:1:length(show)
            set(handles.(show{i}),'Visible','on')
        end
        OutNames = {'Capacity';'Electricity';};
        Data(:,1) = Gen.Output.Capacity;
        Data(:,2) = Gen.Output.Electricity;
        if nnz(Gen.Output.Heat) > 0
            Data(:,3) = Gen.Output.Heat;
            OutNames = {'Capacity';'Electricity';'Heat';};
        end
        set(handles.compText1,'String',num2str(Gen.Size));
        set(handles.textEdit1,'String','Capacity ');
        set(handles.textEdit1Units,'String','(kW)');
        set(handles.uitableEffCurve,'Data',Data);
        set(handles.uitableEffCurve,'RowName',{});
        set(handles.uitableEffCurve,'ColumnName',OutNames);
        for i = 1:1:length(hide)
            set(handles.(hide{i}),'Visible','off')
        end
        plotGenEfficiency(Gen,handles)
        plotResponse(Gen,handles)
    case {'Utility'}
        if strcmp(Gen.Source,'Electricity')
            show = {'compText1';'textEdit1';'textEdit1Units'; ...
                'compText2';'textEdit2';'textEdit2Units';'compText3'; ...
                'textEdit3';'uitableEffCurve';};
            hide = setdiff(uiElmts,show); % everything else not in show
            for i = 1:1:length(show)
                set(handles.(show{i}),'Visible','on')
            end
            clearGenAxes(handles) % clear all figures
            % Populate text boxes
            set(handles.compText1,'String', ...
                strcat('[',num2str(Gen.VariableStruct.SumRates(:,1)'),']'));
            set(handles.textEdit1,'String','Summer Rates');
            set(handles.textEdit1Units,'String','[l m h]');
            set(handles.compText2,'String', ...
                strcat('[',num2str(Gen.VariableStruct.WinRates(:,1)'),']'));
            set(handles.textEdit2,'String','Winter Rates');
            set(handles.textEdit2Units,'String','[l m h]');
            set(handles.compText3,'String', ...
                num2str(Gen.VariableStruct.SellBackPerc));
            set(handles.textEdit3,'String','Sell Back %');
            set(handles.uitableEffCurve,'Data',Gen.VariableStruct.SumRateTable);
            set(handles.uitableEffCurve,'RowName', ...
                {'Sun';'Mon';'Tue';'Wed';'Thu';'Fri';'Sat';});
            hrs = cell(24,1);
            for i = 1:1:length(hrs)
                hrs{i} = num2str(i);
            end
            set(handles.uitableEffCurve,'ColumnName',hrs);
            for i = 1:1:length(hide)
                set(handles.(hide{i}),'Visible','off')
            end
        else %Gen.Source: NG (Natural Gas)
            show = {'compText1';'textEdit1';'textEdit1Units';};
            hide = setdiff(uiElmts,show); % everything else not in show
            for i = 1:1:length(show)
                set(handles.(show{i}),'Visible','on')
            end
            set(handles.compText1,'String',num2str(Gen.VariableStruct.Rate(1)));
            set(handles.textEdit1,'String','Fuel Rate');
            set(handles.textEdit1Units,'String','($/MMBTU)');
            for i = 1:1:length(hide)
                set(handles.(hide{i}),'Visible','off')
            end
        end
    case {'Solar'}
        show = {'compText1';'textEdit1';'textEdit1Units';};
        for i = 1:1:length(show)
            set(handles.(show{i}),'Visible','on')
        end
        set(handles.compText1,'String',num2str(Gen.VariableStruct.Size));
        set(handles.textEdit1,'String','Solar Capacity');
        set(handles.textEdit1Units,'String','(kW)');
        clearGenAxes(handles) % clear all figures
        hide = setdiff(uiElmts,show); % everything else not in show
        for i = 1:1:length(hide)
            set(handles.(hide{i}),'Visible','off')
        end
    case {'Electric Storage'}
        show = {'compText1';'textEdit1';'textEdit1Units';};
        for i = 1:1:length(show)
            set(handles.(show{i}),'Visible','on')
        end
        set(handles.compText1,'String',num2str(Gen.Size));
        set(handles.textEdit1,'String','Storage Capacity');
        set(handles.textEdit1Units,'String','(kW)');
        clearGenAxes(handles) % clear all figures
        hide = setdiff(uiElmts,show); % everything else not in show
        for i = 1:1:length(hide)
            set(handles.(hide{i}),'Visible','off')
        end
    case {'Heating Demands'}
        show = {'compText1';'textEdit1';'textEdit1Units'};
        for i = 1:1:length(show)
            set(handles.(show{i}),'Visible','on')
        end
        set(handles.compText1,'String',num2str(Gen.Demand));
        set(handles.textEdit1,'String','Demand');
        set(handles.textEdit1Units,'String','(kW)');
        clearGenAxes(handles) % clear all figures
        hide = setdiff(uiElmts,show); % everything else not in show
        for i = 1:1:length(hide)
            set(handles.(hide{i}),'Visible','off')
        end
    case {'Hot Water Demands'}
        show = {'compText1';'textEdit1';'textEdit1Units'};
        for i = 1:1:length(show)
            set(handles.(show{i}),'Visible','on')
        end
        set(handles.compText1,'String',num2str(Gen.Demand));
        set(handles.textEdit1,'String','Demand');
        set(handles.textEdit1Units,'String','(kW)');
        clearGenAxes(handles) % clear all figures
        hide = setdiff(uiElmts,show); % everything else not in show
        for i = 1:1:length(hide)
            set(handles.(hide{i}),'Visible','off')
        end
    case {'Cooling Demands'}
        show = {'compText1';'textEdit1';'textEdit1Units'};
        for i = 1:1:length(show)
            set(handles.(show{i}),'Visible','on')
        end
        set(handles.compText1,'String',num2str(Gen.Demand));
        set(handles.textEdit1,'String','Demand');
        set(handles.textEdit1Units,'String','(kW)');
        clearGenAxes(handles) % clear all figures
        hide = setdiff(uiElmts,show); % everything else not in show
        for i = 1:1:length(hide)
            set(handles.(hide{i}),'Visible','off')
        end
    case {'AC/DC Conversion'}
        show = {'compText1';'textEdit1';'textEdit1Units'};
        for i = 1:1:length(show)
            set(handles.(show{i}),'Visible','on')
        end
        set(handles.compText1,'String',num2str(Gen.Efficiency*100));
        set(handles.textEdit1,'String','Conversion Efficiency');
        set(handles.textEdit1Units,'String','(%)');
        clearGenAxes(handles) % clear all figures
        hide = setdiff(uiElmts,show); % everything else not in show
        for i = 1:1:length(hide)
            set(handles.(hide{i}),'Visible','off')
        end
    case {'None'} % for when everything has been removed from Plant
        clearGenAxes(handles) % clear all figures
        for i = 1:1:length(uiElmts)
            set(handles.(uiElmts{i}),'Visible','off')
        end
end

end

function plotGenEfficiency(Gen,handles)
axes(handles.EffCurve)
hold off
str = {};
if isfield(Gen.Output,'Electricity') && nnz(Gen.Output.Electricity) > 0
    c = Gen.Output.Capacity./Gen.Output.Electricity;
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Electricity, ...
        Gen.Output.Capacity(2:end),c(2:end));
    hold on
    set(H1,'Color','k','LineStyle','-','LineWidth',2,'Marker','o')
    str(end+1) = {'Electric'};
end
if isfield(Gen.Output,'Heat') && nnz(Gen.Output.Heat) > 0
    if strcmp(Gen.Type,'CHP Generator')
       plot(AX(1),Gen.Output.Capacity,Gen.Output.Heat,'r-o');
    else
        c = Gen.Output.Capacity./Gen.Output.Heat;
        [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Heat, ...
            Gen.Output.Capacity(2:end),c(2:end));
        set(H1,'Color','r','LineStyle','-','LineWidth',2,'Marker','o')
        hold on
    end
    str(end+1) = {'Heat'};
end
if isfield(Gen.Output,'Steam') && nnz(Gen.Output.Steam) > 0
    c = Gen.Output.Capacity./Gen.Output.Steam; 
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Steam, ...
        Gen.Output.Capacity(2:end),c(2:end));
    set(H1,'Color','m','LineStyle','-','LineWidth',2,'Marker','o')
    hold on
    str(end+1) = {'Steam'};
end
if isfield(Gen.Output,'Cooling') && nnz(Gen.Output.Cooling) > 0
    c = Gen.Output.Capacity./Gen.Output.Cooling; 
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Cooling, ...
        Gen.Output.Capacity(2:end),c(2:end));
    set(H1,'Color','b','LineStyle','-','LineWidth',2,'Marker','o')
    str(end+1) = {'Cooling'};
end
set(AX,{'ycolor'},{'k';'k'})
set(H2,'Color','g','LineStyle',':','LineWidth',3)
str(end+1) = {'Cost Curve'};

if isfield(Gen.Output,'Cooling') && max(Gen.Output.Cooling) > 1
    ylim([0 max(Gen.Output.Cooling)])
    a = round(max(Gen.Output.Cooling))+1;
    ylim(AX(1),[0,a])
    ylim(AX(2),[0,1])
else
    a = round(max(c))+1;
    ylim(AX(1),[0,1])
    ylim(AX(2),[0,a])
    set(AX(1),'YTick',0:.1:1)
    set(AX(2),'YTick',0:a/10:a)
end
set(get(AX(1),'Ylabel'),'String','Efficiency')
set(get(AX(2),'Ylabel'),'String','Cost Curve Shape')
xlabel('% of Capacity')
legend(str);
title('Efficiency / Cost')

set(handles.EffCurve,'UserData',AX)

end

function plotResponse(Gen,handles)
A = Gen.VariableStruct.StateSpace.A;
B = Gen.VariableStruct.StateSpace.B;
C = Gen.VariableStruct.StateSpace.C;
D = Gen.VariableStruct.StateSpace.D;
Dt = Gen.VariableStruct.StateSpace.Dt;
SS = ss(A,B,C,D,Dt);
x0 = [];
Names = {};
if isfield(Gen.Output,'Electricity') && nnz(Gen.Output.Electricity) > 0
    x0(end+1) = Gen.VariableStruct.Startup.Electricity(end);
    Names(end+1) = {'Electricity'};
end
if isfield(Gen.Output,'Cooling') && nnz(Gen.Output.Cooling) > 0

    x0(end+1) = Gen.VariableStruct.Startup.Cooling(end);
    Names(end+1) = {'Cooling'};
end
if isfield(Gen.Output,'Steam') && nnz(Gen.Output.Steam) > 0
    x0(end+1) = Gen.VariableStruct.Startup.Steam(end);
    Names(end+1) = {'Steam'};
end
if isfield(Gen.Output,'Heat') && nnz(Gen.Output.Heat) > 0
    x0(end+1) = Gen.VariableStruct.Startup.Heat(end);
    Names(end+1) = {'Heat'};
end
nS = round(3600/Dt)+1;
t = linspace(0, Dt*(nS-1),nS);
u = Gen.Size*linspace(1,1,nS);
[n,n2] = size(C);
X0 = zeros(n2,1);
for i = 1:1:n
    X0(find(C(i,:),1,'first'))=x0(i);
end
% r = size(SS(1,1).A);
% s = length(x0);
% x0 = [x0;zeros(r(1)-s,1);];
[y,t] = lsim(SS,u,t,X0);
RR = plot(handles.ResponseRate,t/60,y);
xlabel('Time (min)')
legend(Names)

set(handles.ResponseRate,'UserData',RR)

end

function clearGenAxes(handles)
AX_eff = get(handles.EffCurve,'UserData');
AX_res = get(handles.ResponseRate,'UserData');
if ~isempty(AX_eff)
    cla(AX_eff(1))
    cla(AX_eff(2))
    set(AX_eff(2),'Visible','off')
    legend(AX_eff(1),'hide')
end
if ~isempty(AX_res)
    set(AX_res,'Visible','off')
end

end
