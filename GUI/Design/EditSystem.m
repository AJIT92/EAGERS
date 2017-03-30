function EditSystem(handles)
global Plant SYSINDEX
Gen = Plant.Generator(SYSINDEX);
set(handles.CompName,'String',Gen.Name);
switch Gen.Type
    case {'CHP Generator';'ElectricGenerator'}
        show = {'MaxOutput';'textEdit1';'textEdit1Units';'uitableEffCurve';'EffCurve';'ResponseRate';'textResponse';};
        for i = 1:1:length(show)
            set(handles.(show{i}),'Visible','on')
        end
        OutNames = {'Capacity';'Electricity';};
        Data(:,1) = Gen.Output.Capacity;
        Data(:,2) = Gen.Output.Electricity;
        if nnz(Gen.Output.Heat)>0
            Data(:,3) = Gen.Output.Heat;
            OutNames = {'Capacity';'Electricity';'Heat';};
        end
        set(handles.MaxOutput,'String',num2str(Gen.Size));
        set(handles.textEdit1,'String','Capacity ');
        set(handles.textEdit1Units,'String','(kW)');
        set(handles.uitableEffCurve,'Data',Data);
        set(handles.uitableEffCurve,'ColumnName',OutNames);
        hide = {'MinOutput';'textEdit2';'textEdit2Units';'RampRate';'textEdit3';'MaxEff';'textEdit4';'textFuel';'Compfuel';};%'ResponseRate';'textResponse';};
        for i = 1:1:length(hide)
            set(handles.(hide{i}),'Visible','off')
        end
        plotGenEfficiency(Gen,handles)
        plotResponse(Gen,handles)
    case {'Utility'}
        if strcmp(Gen.Source,'Electricity')
            show = {'MaxOutput';'textEdit1';'textEdit1Units';'MinOutput';'textEdit2';'textEdit2Units';'RampRate';'textEdit3';'uitableEffCurve';'EffCurve';'ResponseRate';'textResponse';};
            for i = 1:1:length(show)
                set(handles.(show{i}),'Visible','on')
            end
            set(handles.MaxOutput,'String',strcat('[',num2str(Gen.VariableStruct.SumRates(:,1)'),']'));
            set(handles.textEdit1,'String','Summer Rates');
            set(handles.textEdit1Units,'String','[l m h]');
            set(handles.MinOutput,'String',strcat('[',num2str(Gen.VariableStruct.WinRates(:,1)'),']'));
            set(handles.textEdit2,'String','Winter Rates');
            set(handles.textEdit2Units,'String','[l m h]');
            set(handles.RampRate,'String',num2str(Gen.VariableStruct.SellBackPerc));
            set(handles.textEdit3,'String','Sell Back %');
            set(handles.uitableEffCurve,'Data',Gen.VariableStruct.SumRateTable);
            set(handles.uitableEffCurve,'ColumnName',{'Sun';'Mon';'Tue';'Wed';'Thu';'Fri';'Sta';});
            hide = {'MaxEff';'textEdit4';'textFuel';'Compfuel';'EffCurve';'ResponseRate';'textResponse';};
            for i = 1:1:length(hide)
                set(handles.(hide{i}),'Visible','off')
            end
        else %fuel
            show = {'MaxOutput';'textEdit1';'textEdit1Units';};
            for i = 1:1:length(show)
                set(handles.(show{i}),'Visible','on')
            end
            set(handles.MaxOutput,'String',num2str(Gen.VariableStruct.Rate(1)));
            set(handles.textEdit1,'String','Fuel Rate');
            set(handles.textEdit1Units,'String','$/MMBTU');
            hide = {'MinOutput';'textEdit2';'textEdit2Units';'RampRate';'textEdit3';'uitableEffCurve';'EffCurve';'ResponseRate';'textResponse';'MaxEff';'textEdit4';'textFuel';'Compfuel';'EffCurve';'ResponseRate';'textResponse';};
            for i = 1:1:length(hide)
                set(handles.(hide{i}),'Visible','off')
            end
        end
    case {'Solar'}
        show = {'MaxOutput';'textEdit1';'textEdit1Units';};
        for i = 1:1:length(show)
            set(handles.(show{i}),'Visible','on')
        end
        set(handles.MaxOutput,'String',num2str(Gen.VariableStruct.Size));
        set(handles.textEdit1,'String','Solar Capacity');
        set(handles.textEdit1Units,'String','kW');
        hide = {'MinOutput';'textEdit2';'textEdit2Units';'RampRate';'textEdit3';'uitableEffCurve';'EffCurve';'ResponseRate';'textResponse';'MaxEff';'textEdit4';'textFuel';'Compfuel';'EffCurve';'ResponseRate';'textResponse';};
        for i = 1:1:length(hide)
            set(handles.(hide{i}),'Visible','off')
        end
end

function plotGenEfficiency(Gen,handles)
axes(handles.EffCurve)
hold off
str = {};
if isfield(Gen.Output,'Electricity') && nnz(Gen.Output.Electricity)>0
    c = Gen.Output.Capacity./Gen.Output.Electricity;
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Electricity,Gen.Output.Capacity(2:end),c(2:end));
    hold on
    set(H1,'Color','k','LineStyle','-','LineWidth',2,'Marker','o')
    str(end+1) = {'Electric'};
end
if isfield(Gen.Output,'Heat') && nnz(Gen.Output.Heat)>0
    if strcmp(Gen.Type,'CHP Generator')
       plot(AX(1),Gen.Output.Capacity,Gen.Output.Heat,'r-o');
    else
        c = Gen.Output.Capacity./Gen.Output.Heat;
        [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Heat,Gen.Output.Capacity(2:end),c(2:end));
        set(H1,'Color','r','LineStyle','-','LineWidth',2,'Marker','o')
        hold on
    end
    str(end+1) = {'Heat'};
end
if isfield(Gen.Output,'Steam') && nnz(Gen.Output.Steam)>0
    c = Gen.Output.Capacity./Gen.Output.Steam; 
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Steam,Gen.Output.Capacity(2:end),c(2:end));
    set(H1,'Color','m','LineStyle','-','LineWidth',2,'Marker','o')
    hold on
    str(end+1) = {'Steam'};
end
if isfield(Gen.Output,'Cooling') && nnz(Gen.Output.Cooling)>0
    c = Gen.Output.Capacity./Gen.Output.Cooling; 
    [AX, H1, H2] = plotyy(Gen.Output.Capacity,Gen.Output.Cooling,Gen.Output.Capacity(2:end),c(2:end));
    set(H1,'Color','b','LineStyle','-','LineWidth',2,'Marker','o')
    str(end+1) = {'Cooling'};
end
set(AX,{'ycolor'},{'k';'k'})
set(H2,'Color','g','LineStyle',':','LineWidth',3)
str(end+1) = {'Cost Curve'};

if isfield(Gen.Output,'Cooling') && max(Gen.Output.Cooling)>1
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
legend(str)
title('Efficiency / Cost')


function plotResponse(Gen,handles)
A = Gen.VariableStruct.StateSpace.A;
B = Gen.VariableStruct.StateSpace.B;
C = Gen.VariableStruct.StateSpace.C;
D = Gen.VariableStruct.StateSpace.D;
Dt = Gen.VariableStruct.StateSpace.Dt;
SS = ss(A,B,C,D,Dt);
x0 = [];
Names = {};
if isfield(Gen.Output,'Electricity') && nnz(Gen.Output.Electricity)>0
    x0(end+1) = Gen.VariableStruct.Startup.Electricity(end);
    Names(end+1) = {'Electricity'};
end
if isfield(Gen.Output,'Cooling') && nnz(Gen.Output.Cooling)>0

    x0(end+1) = Gen.VariableStruct.Startup.Cooling(end);
    Names(end+1) = {'Cooling'};
end
if isfield(Gen.Output,'Steam') && nnz(Gen.Output.Steam)>0
    x0(end+1) = Gen.VariableStruct.Startup.Steam(end);
    Names(end+1) = {'Steam'};
end
if isfield(Gen.Output,'Heat') && nnz(Gen.Output.Heat)>0
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
plot(handles.ResponseRate,t/60,y);
xlabel('Time (min)')
legend(Names)