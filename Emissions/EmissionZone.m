function [ eGridZone ] = EmissionZone( State,ZipCode )
%EmissionZone Returns an eGrid Emission Region code based on State (two letter abbreviation string) and Zip Code (5 digit number or string) inputs.
%   The function first looks for State variable, and if it's not in the
%   list of states and territories in the EmissionZonesByState folder, it
%   checks for a valid zip code.  If neither are valid, a emission region 
%   of NWPP is returned.
global Model_dir
if ischar(ZipCode)
    ZipCode=str2double(ZipCode);
end
if ~ischar(State)
    State=char(State);
end

if exist(fullfile(Model_dir,'Emissions','EmissionZonesByState',strcat(State,'.mat')),'file')==2
    State=upper(State);
    load(fullfile(Model_dir,'Emissions','EmissionZonesByState',strcat(State,'.mat')))
    
    zipInd=eval([State '.ZipInfo.ZipCode'])==ZipCode;
    County=eval([State '.ZipInfo.County(zipInd)']);
    eGridZone=cell2mat(eval([State '.Major']));
    
    if ~isempty(County)
        excInd=strcmpi(eval([State '.Exceptions(1:end,1)']),County);
        if any(excInd)
           eGridZone =eval([State '.Exceptions{excInd,2}']);
        end
    end
else
    eGridZone='NWPP';
end
