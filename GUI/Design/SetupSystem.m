function SetupSystem(hObject,handles)
%SetupSystem(hObject,handles)
%Edits the selected component in System Spec tab

global Plant SYSINDEX

% Put Plant.Generator information into lists
nG = length(Plant.Generator);
isType = zeros(nG,1);
isFC = zeros(nG,1);
type = cell(nG,1);
name = cell(nG,1);
source = cell(nG,1);
for i = 1:nG
    type(i) = {Plant.Generator(i).Type};
    name(i) = {Plant.Generator(i).Name};
    source(i) = {Plant.Generator(i).Source};
    if strcmp(type(i),'CHP Generator') || ...
            strcmp(type(i),'Electric Generator')
        if Plant.Generator(i).VariableStruct.isFuelCell
            isFC(i) = 1;
        end
    end
end

% Find indices of relevant components
str = get(hObject,'String');
switch str
    case {'Fuel Cell';'ICE / mGT'}
        if strcmp(str,'Fuel Cell')
            isType = strcmp('CHP Generator',type).*isFC + ...
                strcmp('Electric Generator',type).*isFC;
        else
            isType = strcmp('CHP Generator',type).*(~isFC) + ...
                strcmp('Electric Generator',type).*(~isFC);
        end
    case {'Electric Utility'}
        isType = strcmp('Utility',type);% .* strcmp('Electricity',source);
    case {'Solar PV'}
        isType = strcmp('Solar',type);
    case {'Battery'}
        isType = strcmp('Electric Storage',type);
    case {'Heating Demands'}
        SYSINDEX = -1;
    case {'Hot Water Demands'}
        SYSINDEX = -2;
    case {'Cooling Demands'}
        SYSINDEX = -3;
    case {'AC / DC Conversion'}
        SYSINDEX = -4;
end
isType = nonzeros(linspace(1,nG,nG)'.*isType);

% Decide whether to show the user a list of the chosen type
if length(isType) > 1
    % s - selection index; OK - whether an option was selected
    [s,OK] = listdlg('PromptString','Select desired component', ...
        'ListSize',[160 100], ...
        'SelectionMode','single', ...
        'ListString',name(isType));
    if OK
        SYSINDEX = isType(s);
        EditSystem(handles)
    end
elseif length(isType) == 1
    SYSINDEX = isType;
    EditSystem(handles)
else
    EditSystem(handles)
end

