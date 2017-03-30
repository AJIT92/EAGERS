function SetupSystem(hObject,handles)
%this function edits the selected component
global Plant SYSINDEX
%Edits the selected component in System Spec tab
% Get Generator.Type corresponding to button clicked
nG = length(Plant.Generator);
isType = zeros(nG,1);
isFC = zeros(nG,1);
type = cell(nG,1);
name = cell(nG,1);
for i = 1:nG
    type(i) = {Plant.Generator(i).Type};
    name(i) = {Plant.Generator(i).Name};
    if strcmp(Plant.Generator(i).Type,'CHP Generator') || strcmp(Plant.Generator(i).Type,'Electric Generator')
        if Plant.Generator(i).VariableStruct.isFuelCell
            isFC(i) = 1;
        end
    end
end
cat = get(hObject,'String');
switch cat
    case {'Fuel Cell';'ICE / mGT'}
        if strcmp(cat,'Fuel Cell')
            isType = strcmp('CHP Generator',type).*isFC + strcmp('Electric Generator',type).*isFC;
        else
            isType = strcmp('CHP Generator',type).*(~isFC) + strcmp('Electric Generator',type).*(~isFC);
        end
    case {'Solar PV'}
        isType = strcmp('Solar',type);
end
isType = nonzeros(linspace(1,nG,nG)'.*isType);
% Decide whether to show the user a list of the chosen type
if length(isType)>1
    % s - index of selection list; v - whether OK was clicked
    [s,v] = listdlg('PromptString','Select the desired component...',...
        'SelectionMode','single',...
        'ListString',name(isType));
    SYSINDEX = isType(s);
elseif ~isempty(isType)
    SYSINDEX = isType;
end
EditSystem(handles)

        