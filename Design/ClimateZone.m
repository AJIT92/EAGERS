function ASHRAEZone = ClimateZone(State,ZipCode)
%ClimateZone Returns an ASHRAE Climate Zone based on State (two letter abbreviation string) and Zip Code (5 digit number or string) inputs.
%   The function first looks for State variable, and if it's not in the
%   list of states and territories in the ClimateZonesByState folder, it
%   checks for a valid zip code.  If neither are valid, a climate zone of
%   5A is returned.
global Model_dir
if ischar(ZipCode)
    ZipCode=str2double(ZipCode);
end
if iscell(State)
    State = char(State{1});
end
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};         
if length(State)>2
    Num = find(strcmp(State,stateName));
    State = char(stateAbrev(Num));
end
if exist(fullfile(Model_dir,'Design','ClimateZonesByState',strcat(State,'.mat')),'file')==2
    State=upper(State);
    goodInput=1;
else
    
    load(fullfile(Model_dir,'data','ZipState.mat'))
    zipIndx=find(ZipCodeInfo.ZipCode==ZipCode);
    if ~isempty(zipIndx)
        State=ZipCodeInfo.State{zipIndx};
        goodInput=1;
    else
        goodInput=0;
    end
end

if goodInput
    load(fullfile(Model_dir,'Design','ClimateZonesByState',strcat(State,'.mat')))
    
    zipInd=eval([State '.ZipInfo.ZipCode'])==ZipCode;
    County=eval([State '.ZipInfo.County(zipInd)']);
    ASHRAEZone=cell2mat(eval([State '.Major']));
    
    if ~isempty(County)
        excInd=strcmpi(eval([State '.Exceptions(1:end,1)']),County);
        if any(excInd)
            ASHRAEZone=eval([State '.Exceptions{excInd,2}']);
        end
    end
else
    ASHRAEZone='5A';
end

switch ASHRAEZone
    case '1A'
        ASHRAEZone='1A';
    case '2A'
        ASHRAEZone='2A';
    case '2B'
        ASHRAEZone='2B';
    case '3A'
        ASHRAEZone='3A';
    case '3B'
        if strcmpi(State,'CA')
            ASHRAEZone='3B-Coast';
        else
            ASHRAEZone='3B';
        end
    case '3C'
        ASHRAEZone='3C';
    case '4A'
        ASHRAEZone='4A';
    case '4B'
        ASHRAEZone='4B';
    case '4C'
        ASHRAEZone='4C';
    case '5A'
        ASHRAEZone='5A';
    case '5B'
        ASHRAEZone='5B';
    case '5C'
        ASHRAEZone='4C';
    case '6A'
        ASHRAEZone='6A';
    case '6B'
        ASHRAEZone='6B';
    case '7'
        ASHRAEZone='7';
    case '8'
        ASHRAEZone='8';
    otherwise
        ASHRAEZone='5A';
        
end