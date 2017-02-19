function [CO2, NOx, SO2] = EmissionProfile(State,varargin)
global Model_dir
if isempty(varargin)
    ZipCode =0;
    Profile = 'profileCombust';
elseif length(varargin)==1
    ZipCode = varargin{1};
    Profile = 'profileCombust';
else ZipCode = varargin{1};
    GridMix = varargin{2};
    if GridMix ==1
        Profile = 'profile';
    else Profile = 'profileCombust';
    end
end
if ZipCode==0 % use state profiles
    State = char(State);
    load(fullfile(Model_dir,'Emissions','EmissionProfileByState', strcat(State, '.mat')));
    CO2=eval([State '.CO2' Profile]);
    NOx=eval([State '.NOx' Profile]);
    SO2=eval([State '.SO2' Profile]);
else %% use egrid region database
    EmissionRegion = EmissionZone( char(State),ZipCode );
    load(fullfile(Model_dir,'Emissions','EmissionProfileByRegion',char(EmissionRegion)));
    CO2 = eval([char(EmissionRegion) '.CO2' Profile]);
    NOx = eval([char(EmissionRegion) '.NOx' Profile]);
    SO2 = eval([char(EmissionRegion) '.SO2' Profile]);
end