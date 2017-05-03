function Ren= RenewableOutput(Gen,Date,Time,str)
%Date is the start date
%Time is the Time vector in hours
%str chooses between forecast and actual data
A = datevec(Date);
B = datenum(A(1),1,1);
RenTimestamp = ((Date-B):1/96:(Date-B+Time(end)/24+1/96))';%the data for the renewable power is recorded in 15 min intervals starting at Jan 1st
if strcmp(str,'Actual')%% Actual Renewable Generation        
    if strcmp(Gen.Type,'Solar') 
        if strcmp(Gen.Tracking,'fixed')
            Power = Gen.Sizem2*(Gen.Irrad/1000).*cosd(Gen.SunZen-Gen.Tilt).*(cosd(Gen.SunAz-Gen.Azimuth))*Gen.Eff;
        elseif strcmp(Gen.Tracking,'1axis')
            Power = Gen.Sizem2*(Gen.Irrad/1000).*cosd(Gen.SunZen-Gen.Tilt)*Gen.Eff;
        else Power = Gen.Sizem2*(Gen.Irrad/1000)*Gen.Eff; %Dual axis
        end
        Ren = interp1(RenTimestamp+B,Power(round(RenTimestamp*96+1)),Date+Time/24);
    elseif strcmp(Gen.Type,'Wind')
        %% need to add wind
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
        if n<(RenTimestamp(end)*96+1+3*96)
            Power = [Power;Power(1:10*96)];%append with another 10 days if it crosses from Dec 31 to Jan1
        end
        Xts = round(RenTimestamp(1)*96)+1;
        XFts = round(RenTimestamp(end)*96)+1;
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
    Ren = interp1(RenTimestamp+B, AveragedDays, Date+Time/24);
end