function PlotData(h,TimeVec,Y,Ylab,color)
D = datevec(TimeVec(1));
days = max(1,round(TimeVec(end)-TimeVec(1)));
dNum = TimeVec(1);
monthvec = [0 31 59 90 120 151 181 212 243 273 304 334 365];
leapyear = mod((D(1)-2004),4)==0;%if it is a leap year it will have a remainder of 0
if leapyear
    monthvec(3:end) = monthvec(3:end)+1;
end
if days>sum(monthvec) %if you are plotting multiple years make the ticks in years
    E = datevec(TimeVec(end));
    plotAxis = dNum+linspace(0,round(TimeVec(end)-TimeVec(1)),E(1)-D(1)+1);
elseif days==sum(monthvec) %if you are plotting a year, make the ticks in months
    plotAxis = dNum+[monthvec(D(2):end),monthvec(12)+monthvec(1:D(2))]-monthvec(D(2))-D(3);%make sure to add the right number of days given the month and day of the year.
elseif days > 31
    [y,m1,~] = datevec(TimeVec(1));
    [~,m2,~] = datevec(TimeVec(end));
    M=0;
    for i=1:1:(m2-m1)
        d = datenum([y,m1+i,1])-datenum([y,m1+i-1,1]);%days in month
        M(end+1) = d;
    end
    plotAxis = dNum+M;
elseif days>1
    plotAxis = dNum+linspace(0,days,days+1);
else
    hours = floor(24*(TimeVec(end)-TimeVec(1)));
    plotAxis = dNum+linspace(0,hours/24,floor(hours/3)+1);  
end

[m,n] = size(Y);
for i = 1:1:n
    plot(h,TimeVec,Y(:,i),color{i});
    hold on
end
ylabel(h,Ylab)

xlim(h,[TimeVec(1) TimeVec(end)])
set(h,'xtick',plotAxis)
monthLabel = ['January  '; 'February '; 'March    '; 'April    '; 'May      '; 'June     '; 'July     '; 'August   '; 'September'; 'October  ';'November ' ;'December ';];
if days>366
    datetick(h,'x','yyyy','keeplimits')
    xlabel(h,'Year') 
elseif days==365|| days==366
    datetick(h,'x','mmmdd','keeplimits')
    xlabel(h,num2str(D(1)))
elseif days>=28 && days<=31
    datetick(h,'x','dd','keeplimits')
    xlabel(h,strcat(monthLabel(D(2),:),'  ',num2str(D(1))))
elseif days==7
    datetick(h,'x','dd','keeplimits')
    xlabel(h,strcat(monthLabel(D(2),:),'  ',num2str(D(1))))
elseif days ==1
    datetick(h,'x','HH','keeplimits','keepticks')
    xlabel(h,strcat(['Hours of ', monthLabel(D(2),:), num2str(D(3))]))
end