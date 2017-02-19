function [Rgen,renew] = RenewableOutput(Date,Time,str)
%Date is the start date
%Time is the Time vector in hours
%str chooses between forecast and actual data
global Plant
Rgen =[];
nG = length(Plant.Generator);
include = {'Solar', 'Wind'};
renew = [];
for i = 1:1:nG 
    if ismember(cellstr(Plant.Generator(i).Type),include) && Plant.Generator(i).Enabled == 1 
        renew(end+1) = i;%only fill GeneratorProfile if there are renewable generators
        Gen = Plant.Generator(i).VariableStruct;
        A = datevec(Date);
        B = datenum(A(1),1,1);
        expectedT = Date:1/96:Date+Time(end)/24;%the data for the renewable power is recorded in 15 min intervals starting at Jan 1st
        %online loops might be in 10 minute intervals or shorter, so add a
        %time point at the end
        if ~(expectedT>=(Date+Time(end)/24))
            expectedT(end+1) = expectedT(end)+1/96;
        end
        if length(expectedT)==1 && Time(end)~=0
            expectedT = [Date,Date+1/96];
        end
        Xts = round((Date-B)*96)+1;
        XFts = round((expectedT(end)-B)*96)+1;
        if strcmp(str,'Actual')%% Actual Renewable Generation        
            if strcmp(Gen.Type,'Solar') 
                if strcmp(Gen.Tracking,'fixed')
                    Power = Gen.Sizem2*(Gen.Irrad/1000).*cosd(Gen.SunZen-Gen.Tilt).*(cosd(Gen.SunAz-Gen.Azimuth))*Gen.Eff;
                elseif strcmp(Gen.Tracking,'1axis')
                    Power = Gen.Sizem2*(Gen.Irrad/1000).*cosd(Gen.SunZen-Gen.Tilt)*Gen.Eff;
                else Power = Gen.Sizem2*(Gen.Irrad/1000)*Gen.Eff; %Dual axis
                end
                Power = Power(Xts:XFts);
            elseif strcmp(Gen.Type,'Wind')
                %% need to add wind
            end
            if length(Time)==1
                Rgen = Power;
            else
                Rgen(:,end+1) = interp1(expectedT,Power,Date+Time/24)';
            end
        elseif strcmp(str,'Predict')
            %% Predict Renewable Generation (average 4 days before and 4 days after)
            %% in this case a "day" is one horizon
            if strcmp(Gen.Type,'Solar') 
                if strcmp(Gen.Tracking,'fixed')
                    Power = Gen.Sizem2*(Gen.Irrad/1000).*cosd(Gen.SunZen-Gen.Tilt).*(cosd(Gen.SunAz-Gen.Azimuth))*Gen.Eff;
                elseif strcmp(Gen.Tracking,'1axis')
                    Power = Gen.Sizem2*(Gen.Irrad/1000).*cosd(Gen.SunZen-Gen.Tilt)*Gen.Eff;
                else Power = Gen.Sizem2*(Gen.Irrad/1000)*Gen.Eff; %Dual axis
                end
                n = length(Power);
                if n<XFts+3*96
                    Power = [Power;Power];
                end
                day_4 = Power(max(1,Xts-4*96):max((XFts-Xts+1),XFts-4*96));
                day_3 = Power(max(1,Xts-3*96):max((XFts-Xts+1),XFts-3*96));
                day_2 = Power(max(1,Xts-2*96):max((XFts-Xts+1),XFts-2*96));
                day_1 = Power(max(1,Xts-1*96):max((XFts-Xts+1),XFts-1*96));
                day1 = Power(min(n-(XFts-Xts),Xts):min(n,XFts));
                day2 = Power(min(n-(XFts-Xts),Xts+1*96):min(n,XFts+1*96));
                day3 = Power(min(n-(XFts-Xts),Xts+2*96):min(n,XFts+2*96));
                day4 = Power(min(n-(XFts-Xts),Xts+3*96):min(n,XFts+3*96));
                AveragedDays = (day_4+day_3+day_2+day_1+day1+day2+day3+day4)/8;
            elseif strcmp(Gen.Type,'Wind')
                %% need to add wind
            end
            Rgen(:,end+1) = interp1(expectedT, AveragedDays, Date+Time/24)';
        end
    end
end