function AddSystem(hObject,handles)
%this function adds a system to the plant
global Plant Model_dir SYSINDEX
nG = length(Plant.Generator);
SYSINDEX = nG+1;
files=dir(fullfile(Model_dir, 'System Library',get(hObject,'Tag'),'*.mat'));
list=strrep({files.name},'.mat','');
[s,OK] = listdlg('PromptString','Select Model', 'SelectionMode','single','ListString',list);
if OK~=1
    disp('Invalid selection. Exiting...')
else
    componentName = list{s};
    load(fullfile(Model_dir, 'System Library',get(hObject,'Tag'),strcat(componentName,'.mat')));
    Plant.Generator(SYSINDEX) = component;
    Plant.Network(1).Equipment(end+1) = {strcat(Plant.Generator(SYSINDEX).Type,'.',Plant.Generator(SYSINDEX).Name)};
end
EditSystem(handles)
updateSystemRep(hObject,[], handles)