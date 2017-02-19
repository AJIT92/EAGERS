function CTmap(Comp,Turb,n)
%% CTmap, draws the compressor and turbine maps in a single figure
global TagInf
Line(1,1:3) = 'k- ';
Line(2,1:3) = 'k-.';
Line(3,1:3) = 'k: '; 
Line(4,1:3) = 'k--';
Line(5,1:3) = 'k- ';
Line(6,1:3) = 'k-.';
Line(7,1:3) = 'k: '; 
Line(8,1:3) = 'k--';
Line(9,1:3) = 'k- ';
Line(10,1:3) = 'k-.';

if isfield(TagInf,'Time') && ~isempty(TagInf.Time)
    point1 = 1;
    while TagInf.Time(point1)<1
        point1 = point1+1;
    end
    point2 =length(TagInf.Time);
else
    point1 = [];
end

figure(n)
Ax(1) = axes('position',[0.1300    0.1100    0.765    0.8150],'Visible','on');
hold on
%% Begin Comp Map%%%%
if ~isempty(Comp)
    size = length(Comp.NflowGMap);
    labels = {};
    for j=1:1:length(Comp.RPM)
        plot(Ax(1),Comp.NflowGMap(j,1:size),(Comp.Pdesign-1)*(Comp.PressRatio(j,1:size)-1)+1,Line(j,1:3),'linewidth',2)
        labels(end+1) = {strcat('nRPM = ',num2str(Comp.RPM(j)))};
    end
    %%%%%-- Make the surge line----%%%
    X=linspace(.3,1.2);
    Y= Comp.surgefit(1)*X.^4+Comp.surgefit(2)*X.^3+Comp.surgefit(3)*X.^2+Comp.surgefit(4)*X+Comp.surgefit(5);
    plot(Ax(1),X,(Comp.Pdesign-1)*(Y-1)+1,'r-','linewidth',3)
    labels(end+1) = {'Surge Line'};
    %%%%%-----------%%%%%%%%%
    contour(Comp.NflowGMap,(Comp.Pdesign-1)*(Comp.PressRatio-1)+1,Comp.Efficiency*Comp.PeakEfficiency,logspace(log10(Comp.PeakEfficiency+.05),log10(.8*Comp.PeakEfficiency),12))
    labels(end+1) = {'Efficiency'};
    if ~isempty(point1) %plot the actual operation including surge margin
        for k=1:1:(length(TagInf.Time(point1:point2)))%% The following Finds The Surge Margin--%%
            SurgePR= Comp.surgefit(1)*TagInf.(Comp.name).Nflow(k-1+point1)^4+Comp.surgefit(2)*TagInf.(Comp.name).Nflow(k-1+point1)^3+Comp.surgefit(3)*TagInf.(Comp.name).Nflow(k-1+point1)^2+Comp.surgefit(4)*TagInf.(Comp.name).Nflow(k-1+point1)+Comp.surgefit(5);
            SurgeGap(k) = SurgePR-TagInf.(Comp.name).PR(k-1+point1)/Comp.Pdesign;
            SurgeMargin(k) =SurgeGap(k)/SurgePR*100;
        end
        %%%Plot the operation & final point%%
        plot(Ax(1),TagInf.(Comp.name).Nflow(point1:point2),TagInf.(Comp.name).PR(point1:point2),'b--','linewidth',2);
        plot(Ax(1),TagInf.(Comp.name).Nflow(point2),TagInf.(Comp.name).PR(point2),'go','linewidth',4);
        labels(end+1) = {'Operating Point'};
        labels(end+1) = {'Final Point'};
    end
end
%% Begin Turb Map
if ~isempty(Turb)
    if isfield(Turb,'PdesignLP') %2-spool turbine
        size = length(Turb.NflowGMap);
        for j=1:1:length(Turb.RPM)
            plot(Turb.NflowSpeed(j,1:size)+.5,(Turb.Pdesign*Turb.PdesignLP-1)*(Turb.PressRatio(j,1:size)-1)+1,Line(j,1:3))
        end
        contour(Turb.NflowSpeed+.5,(Turb.Pdesign*Turb.PdesignLP-1)*(Turb.PressRatio-1)+1,Turb.Efficiency*Turb.PeakEfficiency,linspace(Turb.PeakEfficiency,.8*Turb.PeakEfficiency,8))%,[.925 .92 .915 .90 .88 .8])
        %%Plot Turb operation
        plot(TagInf.(Turb.name).Nflow(point1:point2).*TagInf.(Turb.name).NRPM(point1:point2)+.5,TagInf.(Turb.name).PR(point1:point2).*TagInf.(Turb.name).PR(point1:point2),'b--','linewidth',2);
        plot(TagInf.(Turb.name).Nflow(length(TagInf.(Turb.name).Nflow))*TagInf.(Turb.name).NRPM(length(TagInf.(Turb.name).NRPM))+.5,TagInf.(Turb.name).PR(length(TagInf.(Turb.name).PR))*TagInf.(Turb.name).PR(length(TagInf.(Turb.name).PR)),'go','linewidth',4);
    %% Begin Turb2 Map
        size = length(Turb.NflowGMap_LP);
        for j=1:1:length(Turb.RPM_LP)
            plot(Turb.NflowSpeed_LP(j,1:size)+.8,(Turb.PdesignLP-1)*(Turb.PressRatio_LP(j,1:size)-1)+1,Line(j,1:3))
        end
        contour(Turb.NflowSpeed_LP+.8,(Turb.PdesignLP-1)*(Turb.PressRatio_LP-1)+1,TurbEfficiency_LP*Turb.PeakEfficiency_LP,linspace(Turb.PeakEfficiency_LP,.8*Turb.PeakEfficiency_LP,8))%,[.925 .92 .915 .90 .88 .8])
        %%Plot Turb operation
        plot(Turb.Nflow_LP(point1:point2).*Turb.NRPM_LP(point1:point2)+.8,Turb.PR_LP(point1:point2),'b--','linewidth',2);
        plot(Turb.Nflow_LP(length(TurbNflow_LP))*Turb.NRPM_LP(length(Turb.NRPM_LP))+.8,Turb.PR_LP(length(Turb.PR_LP)),'go','linewidth',4);
        axis([0 1.8 0 Turb.Pdesign*Turb.Pdesign_LP])
    else
        size = length(Turb.NflowGMap);
        for j=1:1:length(Turb.RPM)
            plot(Turb.NflowSpeed(j,1:size)+.5,(Turb.Pdesign-1)*(Turb.PressRatio(j,1:size)-1)+1,Line(j,1:3))
        end
        contour(Turb.NflowSpeed+.5,(Turb.Pdesign-1)*(Turb.PressRatio-1)+1,Turb.Efficiency*Turb.PeakEfficiency,linspace(Turb.PeakEfficiency,.8*Turb.PeakEfficiency,8))%,[.925 .92 .915 .90 .88 .8])
        if ~isempty(point1) %plot the actual operation including surge margin
            plot(TagInf.(Turb.name).Nflow(point1:point2).*TagInf.(Turb.name).NRPM(point1:point2)+.5,TagInf.(Turb.name).PR(point1:point2),'b--','linewidth',2);
        end
            plot(TagInf.(Turb.name).Nflow(length(TagInf.(Turb.name).Nflow))*TagInf.(Turb.name).NRPM(length(TagInf.(Turb.name).NRPM))+.5,TagInf.(Turb.name).PR(length(TagInf.(Turb.name).PR)),'go','linewidth',4);
    end
    if isempty(Comp)
        labels = {};
        for j=1:1:length(Turb.RPM)
            labels(end+1) = {strcat('nRPM = ',num2str(Turb.RPM(j)))};
        end
        if ~isempty(point1) %plot the actual operation including surge margin
            labels(end+1) = {'Operating Point'};
            labels(end+1) = {'Final Point'};
        end
    end
end
set(gca,'Fontsize',12), legend(labels,'Location','NorthWest')%
ylabel('Pressure Ratio','fontsize',16,'fontweight','l')
xlabel('Normalized Flow','fontsize',16,'fontweight','l')
set(Ax(1),'XColor','k','YColor','k','Fontsize',12)
set(Ax(1),'XTick',0:.2:1.5) 
% %Changing fonts for axis labels and titles
% xlabel('Power [MW]','FontName','AvantGarde','fontsize',16,'fontweight','l'); 
% ylabel('Compressor Pressure [kPa]','FontName','AvantGarde','fontsize',16,'fontweight','l');             
% set(gca,'FontName','Helvetica');
% set(gcf,'PaperPositionMode','auto')