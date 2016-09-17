function PlotData(TimeVec,Y,forecast,currentstr,Data,HistDisp,RunDisp,BuildProf)
cla reset
if ~isempty(TimeVec)
   X = TimeVec;
elseif ~isempty(HistDisp)
    X = HistDisp.X;
elseif ~isempty(RunDisp)
    X = RunDisp.X;
elseif ~isempty(BuildProf)
    X = BuildProf.X;
end

D = datevec(X(1));
days = max(1,round(X(end)-X(1)));
Steps = round(1/(X(2)-X(1)));%points in 1 day of data
dNum = X(1);
monthvec = [0 31 59 90 120 151 181 212 243 273 304 334 365];
leapyear = mod((D(1)-2004),4)==0;%if it is a leap year it will have a remainder of 0
if leapyear
    monthvec(3:end) = monthvec(3:end)+1;
end
if days>366 %if you are plotting multiple years make the ticks in years
    E = datevec(X(end));
    plotAxis = dNum+linspace(0,round(X(end)-X(1)),E(1)-D(1)+1);
elseif days>=365 %if you are plotting a year, make the ticks in months
    plotAxis = dNum+[monthvec(D(2):end),monthvec(12)+monthvec(1:D(2))]-monthvec(D(2))-D(3);%make sure to add the right number of days given the month and day of the year.
elseif days > 31
    [y,m1,~] = datevec(X(1));
    [~,m2,~] = datevec(X(end));
    M=0;
    for i=1:1:(m2-m1)
        d = datenum([y,m1+i,1])-datenum([y,m1+i-1,1]);%days in month
        M(end+1) = d;
    end
    plotAxis = dNum+M;
elseif days>1
    plotAxis = dNum+linspace(0,days,days+1);
else
    h = floor(24*(X(end)-X(1)));
    plotAxis = dNum+linspace(0,h/24,h+1);  
end
%% Plot
AX(1) = gca;
if strcmp(currentstr,'SOC')
    if isempty(RunDisp)
        [AX, H1, H2] =plotyy(X,Y(:,1),X,Y(:,2));
    else [AX, H1, H2] =plotyy(RunDisp.X,RunDisp.Y(:,1),RunDisp.X,RunDisp.Y(:,2));
    end
    set(AX,{'ycolor'},{'k';'k'})
    if ~isempty(HistDisp)
        hold on
        plot(AX(1),HistDisp.X,HistDisp.Y(:,1))
        plot(AX(2),HistDisp.X,HistDisp.Y(:,2))
    end
    ylabel(AX(1),'Output (kW)','Color','k')
    ylabel(AX(2),'State of Charge (kWh)','Color','k')
    hold on
elseif strcmp(currentstr,'Disp') 
    [m,n] = size(Y);
    color = {'g','r','c','m','k'};
    for i = 1:1:n
        plot(X,Y(:,i),color{i});
        hold on
    end        
elseif forecast
    color = {'k';'g';'r';'b'};
    Y(2:4,:) = CreateForecast(currentstr,X(1),[(X-X(1)).*24]',(length(X)-1)/Steps,[],Data,'HiLow');
    for i = 1:1:2 %actual, prediction, high, low predictions
        plot(X,Y(i,:),char(color(i)))
        if i ==1
            hold on
        end
    end
    legend('Actual', 'Prediction')%,'Pred High', 'Pred Low')         
elseif ~isempty(Y)
    plot(X,Y);
end
if ~isempty(HistDisp)
    plot(HistDisp.X,HistDisp.Y,'b');
end
if ~isempty(RunDisp) && ~strcmp(currentstr, 'SOC') %don't repeat the plot if you are showing SOC
    plot(RunDisp.X,RunDisp.Y,'k');
end 
if ~isempty(BuildProf)
    plot(BuildProf.X,BuildProf.Y,'b');
end 
z=1;
if strcmp(currentstr,'T') 
    ylabel('Temperature (C)')
elseif strcmp(currentstr,'SOC') 
    z=2;
elseif strcmp(currentstr,'Disp') 
    ylabel('Output (kW)')
else%if ~strcmp(currentstr, 'SOC')
    ylabel('Demand (kW)')
end
         

for i = 1:1:z
    if i==1
        ax = AX(1);
    else ax = AX(2);
    end
    xlim(ax,[X(1) X(end)])
    set(ax,'xtick',plotAxis)
    monthLabel = ['January  '; 'February '; 'March    '; 'April    '; 'May      '; 'June     '; 'July     '; 'August   '; 'September'; 'October  ';'November ' ;'December ';];
    if days>366
        datetick(ax,'x','yyyy','keeplimits')
        xlabel(ax,'Year') 
    elseif days==365|| days==366
        datetick(ax,'x','mmmdd','keeplimits')
        xlabel(ax,num2str(D(1)))
    elseif days>=28 && days<=31
        datetick(ax,'x','dd','keeplimits')
        xlabel(ax,strcat(monthLabel(D(2),:),'  ',num2str(D(1))))
    elseif days==7
        datetick(ax,'x','dd','keeplimits')
        xlabel(ax,strcat(monthLabel(D(2),:),'  ',num2str(D(1))))
    elseif days ==1
        datetick(ax,'x','HH','keeplimits','keepticks')
        xlabel(ax,strcat(['Hours of ', monthLabel(D(2),:), num2str(D(3))]))
    end
end