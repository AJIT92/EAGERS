function [CO2avg, NOxavg, SO2avg] = AvgEmissionProfile(State,varargin)
global Model_dir
Profile = 'avgDayCombust';
State = char(State);
load(fullfile(Model_dir,'Emissions','EmissionProfileByState', strcat(State, '.mat')));
CO2=eval([State '.CO2' Profile]);
NOx=eval([State '.NOx' Profile]);
SO2=eval([State '.SO2' Profile]);
monthEnd = [31 59 90 120 151 181 212 243 273 304 334 365];
month = 1;
for i=1:1:365
    if i > monthEnd(month)
        month = month + 1;
    end
   CO2avg(1+24*(i-1):24*i) = CO2(month,:);
   NOxavg(1+24*(i-1):24*i) = NOx(month,:);
   SO2avg(1+24*(i-1):24*i) = SO2(month,:);
end
