function calculateHistoricalFit
global Plant
h=waitbar(0,'Recalculating surface fit');
Steps = round(1/(Plant.Data.Timestamp(2)-Plant.Data.Timestamp(1)));%points in 1 day of data
Resolution =24/Steps;
[~,m,~,hour,minutes] = datevec(Plant.Data.Timestamp(2:end-1));%avoid checking month 1st and last point if starting or ending at 00
m = [m(1) m m(end)];%make m the same length as data
hour = [max(0,floor(hour(1)-Resolution)) hour min(24,hour(end)+floor((minutes(end)+60*Resolution)/60))]';%make hour the same lenght as data
minutes = [max(0,minutes(1)-60*Resolution) minutes mod(minutes(end)+60*Resolution,60)]';%make hour the same lenght as data
months = unique(m);
Xs = 1+round(Steps*(ceil(Plant.Data.Timestamp(1))-Plant.Data.Timestamp(1)));

%% fit temperature data
Z = ones(12,24);
for i = 1:1:length(months)
    waitbar(i/length(months),h,'Calculating Temperature fit');
    D = datevec(Plant.Data.Timestamp(Xs));
    days = floor(nnz(m==D(2))/Steps); %data points in month/points per day
    Total = zeros(24,1);
    points = zeros(24,1);
    for d = 1:1:days
        for j = 1:1:24
            XFs = Xs+round(Steps/24)-1;
            Total(j) = Total(j)+sum(Plant.Data.Temperature(Xs:XFs));
            points(j) = points(j) +(XFs-Xs+1);
            Xs = XFs+1;
        end
    end
    Z(D(2),:) = (Total./points)';
end
if length(months)<12
    if nnz(months==1)==1 && nnz(months==12)==1
        j =1;
        for i= 2:1:11
            if nnz(months==i)==1
                j=j+1;
            else Z(i,:) = Z(j,:);
            end
        end
    else
        for i= 1:1:min(months)-1
            Z(i,:) = Z(min(months),:);
        end
        for i= max(months)+1:1:12
            Z(i,:) = Z(max(months),:);
        end
    end
end
Plant.Data.HistProf.Temperature = Z;

monthNames = cellstr({'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';});
dateDay = round(Plant.Data.Timestamp);
BusDay = zeros(length(dateDay),1);
for i = 1:1:length(dateDay)
    BusDay(i) = isbusday(dateDay(i),Plant.Data.Holidays,[1 0 0 0 0 0 1]);
end
%% make a surface fit for Demand.E, Demand.H, and Demand.C
Outs = fieldnames(Plant.Data.Demand);
for s = 1:1:length(Outs)
    if length(months)<=4 %make 1 surface fit
        for k = 1:1:length(Plant.Data.Demand.(Outs{s})(1,:))
            Z = Plant.Data.Demand.(Outs{s})(:,k);
            Y = Plant.Data.Temperature;
            X = hour+minutes/60;
            if Steps>100 % too many data points for a surface fit, average each hour
                r = ceil(Steps/24);
                n2 = floor(length(X)/r);
                X2 = zeros(n2,1);
                Y2 = zeros(n2,1);
                Z2 = zeros(n2,1);
                for i=1:1:n2%turn while loop into for loop to make loop faster
                    X2(i)=sum(X((i-1)*r+1:i*r))/r;
                    Y2(i)=sum(Y((i-1)*r+1:i*r))/r;
                    Z2(i)=sum(Z((i-1)*r+1:i*r))/r;
                end
                X = X2;
                Y = Y2;
                Z = Z2;
            end
            Plant.Data.HistProf.(Outs{s})(k).Annual = fit([X Y],Z,'lowess');
        end
    else %split surface fit by month and weekday/weekend
        for i = 1:1:length(months)
            waitbar((i-1)/length(months),h,strcat('Calculating surface fit for ',Outs{s},' for the month of ',monthNames{months(i)}));
            for k = 1:1:length(Plant.Data.Demand.(Outs{s})(1,:))
                Z = Plant.Data.Demand.(Outs{s})(:,k);
                Y = Plant.Data.Temperature;
                X = hour+minutes/60;
                if Steps>150 % too many data points for a surface fit, average each 10 min
                    r = ceil(Steps/(24*6));
                    n2 = floor(length(X)/r);
                    X2=zeros(n2,1);
                    Y2=zeros(n2,1);
                    Z2=zeros(n2,1);
                    for j=1:1:n2%turn the while loop into a for loop to make loop faster
                        X2(j)=sum(X((j-1)*r+1:j*r))/r;
                        Y2(j)=sum(Y((j-1)*r+1:j*r))/r;
                        Z2(j)=sum(Z((j-1)*r+1:j*r))/r;
                    end
                     X = X2;
                    Y = Y2;
                    Z = Z2;
                end
                % business days
                H1 = BusDay((find(m==months(i))));%find which days are businessdays
                H1non0=find(H1);%find the indexes of all the nonzeros (all the businessdays are nonzero)
                Z1 = Z(H1non0);%remove the days that are non business days
                Y1 = Y(H1non0);
                X1 = X(H1non0);%cannot use nonzeros function because there are some elements in X1 that are 0
                % weekends/holidays
                H2 = 1-H1;
                H2non0=find(H2);
                Z2 = Z(H2non0);
                Y2 = Y(H2non0);
                X2 = X(H2non0);
                Plant.Data.HistProf.(Outs{s})(k).(strcat(monthNames{months(i)},'WeekDay')) = fit([X1, Y1],Z1,'lowess');%'poly23');%'loess');
                Plant.Data.HistProf.(Outs{s})(k).(strcat(monthNames{months(i)},'WeekEnd')) = fit([X2, Y2],Z2,'lowess');%'poly23');%'loess');
            end
        end
    end
end
close(h)