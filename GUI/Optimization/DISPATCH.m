function varargout = DISPATCH(varargin)
% DISPATCH MATLAB code for DISPATCH.fig
%      DISPATCH, by itself, creates a new DISPATCH or raises the existing
%      singleton*.
%
%      H = DISPATCH returns the handle to a new DISPATCH or the handle to
%      the existing singleton*.
%
%      DISPATCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPATCH.M with the given input arguments.
%
%      DISPATCH('Property','Value',...) creates a new DISPATCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DISPATCH_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DISPATCH_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
 
% Edit the above text to modify the response to help DISPATCH
 
% Last Modified by GUIDE v2.5 28-Apr-2017 13:28:12
 
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DISPATCH_OpeningFcn, ...
                   'gui_OutputFcn',  @DISPATCH_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
 
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
 
 
% --- Executes just before DISPATCH is made visible.
function DISPATCH_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for EAGERS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Plant Model_dir
set(gcf,'Name','DISPATCH')
movegui(gcf,'center');
if ~isfield(Plant.optimoptions,'method') || strcmp(Plant.optimoptions.method,'Planning')
    Plant.optimoptions.method = 'Dispatch';
end
if ~isfield(Plant.optimoptions,'MixedInteger')
    Plant.optimoptions.MixedInteger = true;
end
if ~isfield(Plant.optimoptions,'SpinReserve')
    Plant.optimoptions.SpinReserve = false;
    Plant.optimoptions.SpinReservePerc = 0;
end

%Set handles for communications
set(handles.editCommandOnOff,'string','--')
set(handles.editCommandSet,'string','--')
set(handles.editMeasureOnOff,'string','--')
set(handles.editMeasureInput,'string','--')
set(handles.editMeasurePrimary,'string','--')
set(handles.editMeasureSecondary,'string','--')

setupTabs(hObject, handles);

%make gen list for tabs 1 and 5
list={};
for i=1:length(Plant.Generator)
    if Plant.Generator(i).Enabled ==1
        list(end+1) = {Plant.Generator(i).Name};
    elseif Plant.Generator(i).Enabled ==0
        list(end+1) = {Plant.Generator(i).Name};
    end 
end
set(handles.uipanelMain1,'UserData',list)

Plant.Plotting.ColorNames = {'parula';'autumn';'cool';'spring';'summer';'winter';};
for i = 1:1:length(Plant.Plotting.ColorNames)
    colormap(handles.ElectricGraph,Plant.Plotting.ColorNames{i});
    Plant.Plotting.ColorMaps{i} = colormap(handles.ElectricGraph);
end

%make swap button
set(handles.SwitchChart,'Units','pixels')
pos = get(handles.SwitchChart,'Position');
[x,map] = imread(fullfile(Model_dir,'GUI','Graphics','swap.png'));
s = imresize(x,[pos(3) pos(4)]);
set(handles.SwitchChart,'cdata',s)

handles = GenList_Make(handles);
set(handles.uipanelMain1,'Visible','off')
set(handles.uipanelMain5,'Visible','on')
handles = GenList_Make(handles);
set(handles.uipanelMain5,'Visible','off')
set(handles.uipanelMain1,'Visible','on')

Plant.GUIhandles = handles;

set(handles.constant, 'value', strcmp(Plant.optimoptions.tspacing,'constant'));
set(handles.linear, 'value', strcmp(Plant.optimoptions.tspacing, 'linear'));
set(handles.logarithm, 'value', strcmp(Plant.optimoptions.tspacing, 'logarithm'));
set(handles.manual, 'value', strcmp(Plant.optimoptions.tspacing, 'manual'));
set(handles.Interval,'string',Plant.optimoptions.Interval);
set(handles.Horizon, 'string', Plant.optimoptions.Horizon);
set(handles.Resolution, 'string', Plant.optimoptions.Resolution);
set(handles.scaletime, 'string', Plant.optimoptions.scaletime);

set(handles.sequential, 'value', Plant.optimoptions.sequential);
set(handles.simultaneous, 'value', ~Plant.optimoptions.sequential);
set(handles.excessHeat, 'value', Plant.optimoptions.excessHeat);
set(handles.nsSmooth, 'string', Plant.optimoptions.nsSmooth);

set(handles.fastsimulation, 'value', Plant.optimoptions.fastsimulation);
set(handles.slowsimulation, 'value', ~Plant.optimoptions.fastsimulation);

set(handles.NoMixedInteger, 'value', ~Plant.optimoptions.MixedInteger);
set(handles.MixedInteger, 'value', Plant.optimoptions.MixedInteger);

set(handles.noSpinReserve, 'value', ~Plant.optimoptions.SpinReserve);
set(handles.SpinReserve, 'value', Plant.optimoptions.SpinReserve);
set(handles.SpinReservePerc, 'string', Plant.optimoptions.SpinReservePerc);
set(handles.editBuffer, 'string', Plant.optimoptions.Buffer);

set(handles.Control, 'value', strcmp(Plant.optimoptions.method,'Control'));
set(handles.Dispatch, 'value', strcmp(Plant.optimoptions.method,'Dispatch'));
set(handles.Topt, 'string', Plant.optimoptions.Topt);
set(handles.Tmpc, 'string', Plant.optimoptions.Tmpc);

MainWindow_Setup(hObject, eventdata, handles);

Plant.Plotting.ColorNames = {'parula';'autumn';'cool';'spring';'summer';'winter';};
for i = 1:1:length(Plant.Plotting.ColorNames)
    colormap(handles.ElectricGraph,Plant.Plotting.ColorNames{i});
    Plant.Plotting.ColorMaps{i} = colormap(handles.ElectricGraph);
end

%Update forecast handles
days = round(Plant.Data.Timestamp(end)-Plant.Data.Timestamp(1));
set(handles.sliderZoom,'Min',1,'Max',4,'Value',1,'SliderStep',[1/3,1/3])
set(handles.sliderDate,'Min',1,'Max',2,'Value',1,'SliderStep',[1/(days-1),1/(days-1)])
set(handles.sliderDate,'Max',2)

%%put something into axes
% InitialDispatch

% --- Outputs from this function are returned to the command line.
function varargout = DISPATCH_OutputFcn(hObject, eventdata, handles) 
%ask user to save?

% Contains setup for tabs
function setupTabs(hObject, handles)
%% Main Tabs
% Assumptions:
% 1. Tags of main tab static text boxes are of form, 'MainTab1',
% 'MainTab2', etc.
% 2. Tags of main tab panels are of form, 'uipanelMain1', 'uipanelMain2',
% etc.
TabText = {'Main Window';'Market Services';'Historian/Forecast';'Control Options';'Communication'};
set(hObject,'UserData',TabText);

for i = 1:length(TabText)
    j = num2str(i);
    % panel management
    set(handles.(strcat('MainTab',j)),'Units','characters','String',TabText{i});
    set(handles.(strcat('uipanelMain',j)),'Units','characters','BorderType','none')
    if i ==1
        pan1pos = get(handles.uipanelMain1,'Position');
        pos = get(handles.MainTab1','Position');
        set(handles.MainTab1,'Position',[pos(1),pos(2),pos(3),pos(4)+.5])
    else
        pos2 = get(handles.(strcat('uipanelMain',j)),'Position');
        set(handles.(strcat('uipanelMain',j)),'Position',[pan1pos(1), pan1pos(2)+(pan1pos(4)-pos2(4)), pos2(3), pos2(4)])
        set(handles.(strcat('uipanelMain',j)),'Visible','off')
    end
end

% Main tabs callback
function mainTab_Callback(hObject, eventdata, handles)
global Plant
n = get(hObject,'Tag');
n = n(end);
m = [];
i = 1;
% Find out which tab is currently selected
while isempty(m) && isfield(handles,strcat('uipanelMain',num2str(i)))
    if strcmp(get(handles.(strcat('uipanelMain',num2str(i))),'Visible'),'on')
        m = i;
    else i = i+1;
    end
end
m = num2str(m);

% CRUCIAL IN NEXT 3 STEPS: m, then n.
% Change color
bColor = get(handles.(strcat('MainTab',m)),'BackgroundColor');
set(handles.(strcat('MainTab',m)),'BackgroundColor',max(0,bColor-.1))
bColor = get(handles.(strcat('MainTab',n)),'BackgroundColor');
set(handles.(strcat('MainTab',n)),'BackgroundColor',min(1,bColor+.1))

% Change dimensions
pos = get(handles.(strcat('MainTab',m)),'Position');
set(handles.(strcat('MainTab',m)),'Position',[pos(1),pos(2),pos(3),pos(4)-.5])
pos = get(handles.(strcat('MainTab',n)),'Position');
set(handles.(strcat('MainTab',n)),'Position',[pos(1),pos(2),pos(3),pos(4)+.5])

% Change visibility
set(handles.(strcat('uipanelMain',m)),'Visible','off')
set(handles.(strcat('uipanelMain',n)),'Visible','on')

%% actions specific to one tab or another
if strcmp(n,'3')
    ForecastPlot(handles)
elseif strcmp(n,'1')
    MainWindow_Setup(hObject, eventdata, handles)
elseif strcmp(n,'4')
    set(handles.Interval,'string',Plant.optimoptions.Interval);
    set(handles.Horizon, 'string', Plant.optimoptions.Horizon);
    set(handles.Resolution, 'string', Plant.optimoptions.Resolution);
    set(handles.Topt, 'string', Plant.optimoptions.Topt);
    set(handles.Tmpc, 'string', Plant.optimoptions.Tmpc);
    set(handles.nsSmooth, 'string', Plant.optimoptions.nsSmooth);
    set(handles.scaletime, 'string', Plant.optimoptions.scaletime);
    set(handles.fastsimulation, 'value', Plant.optimoptions.fastsimulation);
    set(handles.slowsimulation, 'value', ~Plant.optimoptions.fastsimulation);
    set(handles.constant, 'value', strcmp(Plant.optimoptions.tspacing,'constant'));
    set(handles.linear, 'value', strcmp(Plant.optimoptions.tspacing, 'linear'));
    set(handles.logarithm, 'value', strcmp(Plant.optimoptions.tspacing, 'logarithm'));
    set(handles.manual, 'value', strcmp(Plant.optimoptions.tspacing, 'manual'));
    set(handles.sequential, 'value', Plant.optimoptions.sequential);
    set(handles.simultaneous, 'value', ~Plant.optimoptions.sequential);
    set(handles.excessHeat, 'value', Plant.optimoptions.excessHeat);
end

% --- Populates the Main Window Tab 
function MainWindow_Setup(hObject, eventdata, handles)
global Plant
Plant.optimoptions.method = 'Dispatch';

if ismember('E',Plant.optimoptions.Outputs)
    set(handles.checkboxElectric,'value',1)
    set(handles.ElecTitle,'Visible','on')
    set(Plant.GUIhandles.ElectricGraph,'Visible','on')
    set(handles.Forecast,'Visible','on')
    MakeMain('Electric Dispatch',handles)%Always make Electric the default if it exists
else
    set(handles.checkboxElectric,'value',0)
    set(handles.ElecTitle,'Visible','off')
    set(Plant.GUIhandles.ElectricGraph,'Visible','off')
    set(handles.Forecast,'Visible','off')
    %% move and re-size heat or cooling dispatch axes & title
end
if ismember('H',Plant.optimoptions.Outputs)
    set(handles.checkboxHeat,'Value',1)
    set(handles.HeatFore,'Visible','on')
    set(handles.HeatTitle,'Visible','on')
    set(Plant.GUIhandles.HeatGraph,'Visible','on')
    if length(Plant.optimoptions.Outputs) == 1 || (length(Plant.optimoptions.Outputs) == 2 && ~ismember('E',Plant.optimoptions.Outputs))
        MakeMain('Heating Dispatch',handles)
    else
        set(handles.HeatTitle,'Position',[150,42,40,2])
        set(Plant.GUIhandles.HeatGraph,'Position',[150,24,59,17.5])
    end
else
    set(handles.checkboxHeat,'Value',0)
    set(handles.HeatTitle,'Visible','off')
    set(Plant.GUIhandles.HeatGraph,'Visible','off')
    set(handles.HeatFore,'Visible','off')
end
if ismember('C',Plant.optimoptions.Outputs)
    set(handles.checkboxCooling,'Value',1)
    set(handles.CoolFore,'Visible','on')
    set(handles.CoolTitle,'Visible','on')
    set(Plant.GUIhandles.CoolGraph,'Visible','on')    
    if length(Plant.optimoptions.Outputs) == 1
        MakeMain('Cooling Dispatch',handles)
    elseif length(Plant.optimoptions.Outputs) == 2
        set(handles.CoolTitle,'Position',[150,42,40,2])
        set(Plant.GUIhandles.CoolGraph,'Position',[150,24,59,17.5])
    else
        set(handles.CoolTitle,'Position',[150,20,40,2])
        set(Plant.GUIhandles.CoolGraph,'Position',[150,2,59,17.5])
    end
else
    set(handles.checkboxCooling,'Value',0)
    set(handles.CoolTitle,'Visible','off')
    set(Plant.GUIhandles.CoolGraph,'Visible','off')
    set(handles.CoolFore,'Visible','off')
end
if ismember('S',Plant.optimoptions.Outputs)
    set(handles.checkboxSteam,'Value',1) 
else
    set(handles.checkboxSteam,'Value',0)
end


% --- Executes on button press in SwitchChart.
function SwitchChart_Callback(hObject, eventdata, handles)
global Plant
order = Plant.optimoptions.Outputs;
if ismember('E',order)
    e = get(Plant.GUIhandles.ElectricGraph,'Position');
    et = get(Plant.GUIhandles.ElecTitle,'Position');
end
if ismember('H',order)
    h = get(Plant.GUIhandles.HeatGraph,'Position');
    ht = get(Plant.GUIhandles.HeatTitle,'Position');
end
if ismember('C',order)
    c = get(Plant.GUIhandles.CoolGraph,'Position');
    ct = get(Plant.GUIhandles.CoolTitle,'Position');
end
if length(order) == 3
        set(Plant.GUIhandles.ElecTitle,'Position',ct);
        set(Plant.GUIhandles.HeatTitle,'Position',et);
        set(Plant.GUIhandles.CoolTitle,'Position',ht);
        set(Plant.GUIhandles.ElectricGraph,'Position',c);
        set(Plant.GUIhandles.HeatGraph,'Position',e);
        set(Plant.GUIhandles.CoolGraph,'Position',h);
elseif length(order) ==2
    if ismember('E',order) && ismember('H',order)
        set(Plant.GUIhandles.ElecTitle,'Position',ht);
        set(Plant.GUIhandles.HeatTitle,'Position',et);
        set(Plant.GUIhandles.ElectricGraph,'Position',h);
        set(Plant.GUIhandles.HeatGraph,'Position',e);
    elseif ismember('E',order) && ismember('C',order)
        set(Plant.GUIhandles.ElecTitle,'Position',ct);
        set(Plant.GUIhandles.CoolTitle,'Position',et);
        set(Plant.GUIhandles.CoolGraph,'Position',c);
        set(Plant.GUIhandles.ElectricGraph,'Position',e);
    elseif ismember('H',order) && ismember('C',order)
        set(Plant.GUIhandles.HeatTitle,'Position',ct);
        set(Plant.GUIhandles.CoolTitle,'Position',ht);
        set(Plant.GUIhandles.CoolGraph,'Position',c);
        set(Plant.GUIhandles.HeatGraph,'Position',h);
    end
end

%If only 1 output, make it main graph
function MakeMain(hObject, handles)
global Plant
maint = [42,45,50,2];%Main Title Position
main = [42,21,86,23];
if strcmp(hObject,'Electric Dispatch')
    set(Plant.GUIhandles.ElectricGraph,'Position',main);
    set(handles.ElecTitle,'Position',maint);
elseif strcmp(hObject,'Heating Dispatch')
    set(Plant.GUIhandles.HeatGraph,'Position',main);
    set(handles.HeatTitle,'Position',maint);
else
    set(Plant.GUIhandles.CoolGraph,'Position',main);
    set(handles.CoolTitle,'Position',maint);
end

% --- Executes on button press in ShowCumulative.
function ShowCumulative_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.CumulativeGraph,'Visible','on')
    set(handles.HeatGraph,'Visible','off')
    set(handles.CoolGraph,'Visible','off')
else
    set(handles.CumulativeGraph,'Visible','off')
    if ismember('H',Plant.optimoptions.Outputs)
        set(handles.HeatGraph,'Visible','on')
    end
    if ismember('H',Plant.optimoptions.Outputs)
        set(handles.CoolGraph,'Visible','on')
    end
end


function handles = GenList_Make(handles)
global Plant Model_dir
list = get(handles.uipanelMain1,'UserData');
if strcmp(get(handles.uipanelMain1,'Visible'),'on')
    r = 'Main';
    p = 1;
elseif strcmp(get(handles.uipanelMain5,'Visible'),'on')
    r = 'Comm';
    p = 5;
end
if length(list)>10
    set(handles.NextGen1,'Visible','on')
    set(handles.NextGen5,'Visible','on')
end
%Makes buttons for GUI && status buttons next to corresponding generator
nG = length(list);
colorVec = Plant.Plotting.ColorMaps{1};
colorsPlot = interp1(linspace(0,1,length(colorVec)),colorVec,linspace(0,1,nG));
for i=1:1:nG 
    num = num2str(i);
    curtab = floor((i-1)/10)+1;
    prev = 10*(curtab-1);
    if curtab==1
        vis = 'on';
    else vis = 'off';
    end
    if p ==1
        pos = 45 - 2*(i-prev);
    elseif p ==5
        pos = 26 - 2*(i-prev);
    end
    callback = strcat('@(hObject,eventdata)DISPATCH(Gen',r,'_Callback,hObject,eventdata,guidata(hObject))');
    btn = uicontrol('Style', 'pushbutton', 'String', list{i},...
    'Units','characters',...
    'Position', [1 pos 25 1.8],...
    'Tag', strcat('Generator',num),...
    'FontSize', 10,...
    'Parent', handles.(strcat('uipanelMain',num2str(p))),...
    'Callback',eval(callback),...
    'Visible',vis,...
    'UserData',i);
    handles.(strcat('Generator',num)) = btn;
    set(handles.(strcat('Generator',num)),'BackgroundColor',colorsPlot(i,:))
    if p ==1 %Only make Status buttons on Main Window
        pos = 45.5 - 2*(i-prev);
        callback = strcat('@(hObject,eventdata)DISPATCH(Status_Callback,hObject,eventdata,guidata(hObject))');
        if Plant.Generator(i).Enabled
            [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
        else
            [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
        end
        pSize = pixelSize(handles);
        s = imresize(x,[3*pSize(1) pSize(2)]);
        if Plant.Generator(i).Enabled
            enableGen  = 'bold';
        else enableGen  = 'normal';
        end
        btn = uicontrol('Style', 'pushbutton', 'String', '',...
        'Units','characters',...
        'Position', [27 pos 3 1],...
        'Tag', strcat('GeneratorStat',num),...
        'cdata', s,...
        'FontWeight',enableGen,...
        'Parent', handles.uipanelMain1,...
        'Callback',eval(callback),...
        'Visible',vis,...
        'UserData',i);
        handles.(strcat('GeneratorStat',num)) = btn;
    end
end
% --- Executes on button press in PrevGen1.
function PrevGen_Callback(hObject, eventdata, handles)
if strcmp(get(handles.uipanelMain1,'Visible'),'on')
    panel = 'uipanelMain1';
    button1 = 'PrevGen1';
    button2 = 'NextGen1';
    mult = 2;%status and button made on main window
elseif strcmp(get(handles.uipanelMain5,'Visible'),'on')
    panel = 'uipanelMain5';
    button1 = 'PrevGen5';
    button2 = 'NextGen5';
    mult = 1;
end
list = get(handles.uipanelMain1,'UserData');
childHandles = get(handles.(panel),'Children');
vis = [];
i = 0;
while isempty(vis) && i<length(list)*mult%Find out last visible generator
    i = i+1;
    if strcmp(get(childHandles(i),'Visible'),'on')
        vis = get(childHandles(i),'UserData');
    end
end
if vis == length(list) 
    r = floor(vis/10)*10;%round down to nearest 10
    if (vis-r)<10
        start = mult*(vis-r);%only turn off the visible, if less then 10 are showing
    end
else
    start = i + (10*mult-1);
end
new = start + (10*mult);
for j = i:start%current visible buttons off
    set(childHandles(j),'Visible','off')
end

for j = start+1:new%current visible buttons off
    set(childHandles(j),'Visible','on')
end
if new == mult*length(list)
    set(handles.(button1),'Visible','off')
    set(handles.(button2),'Visible','on')
else
    set(handles.(button1),'Visible','on')
    set(handles.(button2),'Visible','on')
end

% --- Executes on button press in NextGen1.
function NextGen_Callback(hObject, eventdata, handles)
if strcmp(get(handles.uipanelMain1,'Visible'),'on')
    panel = 'uipanelMain1';
    button1 = 'PrevGen1';
    button2 = 'NextGen1';
    mult = 2;%status and button made on main window
elseif strcmp(get(handles.uipanelMain5,'Visible'),'on')
    panel = 'uipanelMain5';
    button1 = 'PrevGen5';
    button2 = 'NextGen5';
    mult = 1;
end
list = get(handles.uipanelMain1,'UserData');
childHandles = get(handles.(panel),'Children');
vis = [];
i = 0;
while isempty(vis) && i<length(list)*mult%Find out last visible generator
    i = i+1;
    if strcmp(get(childHandles(i),'Visible'),'on')
        vis = get(childHandles(i),'UserData');
    end
end
dif = length(list)-vis;
for j = i:(mult*length(list))%current visible buttons off
    set(childHandles(j),'Visible','off')
end
if dif<10 || (dif==10 && length(list)==10)%For less than (or exactly 10 left) visible components
    set(handles.(button2),'Visible','off')
    set(handles.(button1),'Visible','on')
    for j = 1:(mult*dif)%show rest of buttons 
        set(childHandles(j),'Visible','on')
    end
else
    set(handles.(button1),'Visible','on')
    set(handles.(button2),'Visible','on')
    n=dif-10;%For plants with more components
    for j = mult*n+1:i-1%show next 10 buttons 
        set(childHandles(j),'Visible','on')
    end
end

%When Component Buttons on the main tab are clicked
function GenMain_Callback(hObject, eventdata, handles)
global Plant
gen = get(hObject,'String');
i = get(hObject,'UserData');
size = num2str(Plant.Generator(i).Size);
set(handles.SelGen,'Title',Plant.Generator(i).Name,'UserData',i)
if ~isempty(strfind(gen,'Utility'))
    set(handles.GenSpec1,'String','Inf')
else
    set(handles.GenSpec1,'String',size)
end
if Plant.Generator(i).Enabled == 1
    set(handles.GenEnable,'Value',1)
    set(handles.GenDisable,'Value',0)
elseif Plant.Generator(i).Enabled == 0
    set(handles.GenEnable,'Value',0)
    set(handles.GenDisable,'Value',1)
end
 

%When Status colors are clicked
function Status_Callback(hObject, eventdata,handles)
global Plant Model_dir
i = get(hObject,'UserData');
size = num2str(Plant.Generator(i).Size);
if ~isempty(strfind(Plant.Generator(i).Name,'Utility'))
    set(handles.GenSpec1,'String','Inf')
else
    set(handles.GenSpec1,'String',size)
end
set(handles.SelGen,'Title',Plant.Generator(i).Name,'UserData',i)
if strcmp(get(hObject,'FontWeight'),'bold')
    Plant.Generator(i).Enabled = 0;
    [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
else
    Plant.Generator(i).Enabled = 1;
    [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
end
set(handles.GenEnable,'value',Plant.Generator(i).Enabled)
set(handles.GenDisable,'Value',~Plant.Generator(i).Enabled)
pSize = pixelSize(handles);
s = imresize(x,[3*pSize(1) pSize(2)]);
set(hObject,'FontWeight','bold','cdata',s)
    
function pSize = pixelSize(handles)
set(handles.Switch,'Units','pixels');
pos1 = get(handles.Switch,'Position');
set(handles.Switch,'Units','characters');
pos2 = get(handles.Switch,'Position');
pSize(1) = pos1(3)/pos2(3);
pSize(2) = pos1(4)/pos2(4);

% --- Executes on button press in Start.
function Start_Callback(hObject, eventdata, handles)
global Virtual RealTime DispatchWaitbar 
Virtual = 1;
RealTime = 0;
set(handles.Start,'Value',1);%reset start button
set(handles.Stop,'Value',0);%reset stop button
DispatchWaitbar=waitbar(0,'Running Dispatch','Visible','off');
RunOptimization
waitfor(DispatchWaitbar)
set(handles.Stop,'Value',1);%reset stop button
set(handles.Start,'Value',0);%reset start button

% --- Executes on button press in Stop.
function Stop_Callback(hObject, eventdata, handles)
%stops dispatch
 global DispatchWaitbar Virtual
 Virtual = 0;
 close(DispatchWaitbar)
 DispatchWaitbar=[];

 
% --- Executes on button press in Switch.
function Switch_Callback(hObject, eventdata, handles)
%%send user back to EPT, pass along the plant generators 
close
MainScreen1

% --- Executes on button press in GenEnable.
function GenEnable_Callback(hObject, eventdata, handles)
global Plant Model_dir
i = get(handles.SelGen,'UserData');
Plant.Generator(i).Enabled = 1;
[x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
pSize = pixelSize(handles);
s = imresize(x,[3*pSize(1) pSize(2)]);
num = 2*length(Plant.Generator) - (2*i - 1);
set(handles.uipanelMain1.Children(num),'FontWeight','normal','cdata',s)

% --- Executes on button press in GenDisable.
function GenDisable_Callback(hObject, eventdata, handles)
global Plant Model_dir
i = get(handles.SelGen,'UserData');
Plant.Generator(i).Enabled = 0;
[x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
pSize = pixelSize(handles);
s = imresize(x,[3*pSize(1) pSize(2)]);
num = 2*length(Plant.Generator) - (2*i - 1);
set(handles.uipanelMain1.Children(num),'FontWeight','normal','cdata',s)
        

function GenStatus1_Callback(hObject, eventdata, handles)
function GenStatus1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function GenStatus2_Callback(hObject, eventdata, handles)
function GenStatus2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function GenSpec1_Callback(hObject, eventdata, handles)
global Plant
type = get(handles.SelGen,'Title');
size = str2double(get(hObject,'string'));
set(handles.GenSpec1,'Value',str2double(get(hObject,'string')));
for i = 1:length(Plant.Generator)
    if strcmp(Plant.Generator(i).Name,type) && ~strcmp(Plant.Generator(i).Type,'Utility')%Cant update a utility size
        Plant.Generator(i).Size = size;
    end
end

function GenSpec1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function GenSpec2_Callback(hObject, eventdata, handles)
function GenSpec2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in LineGraph.
function LineGraph_Callback(hObject, eventdata, handles)
if get(handles.LineGraph,'Value')==1
    a = get(handles.StackedGraph,'BackgroundColor');
    b = get(handles.LineGraph,'BackgroundColor');
    c = get(handles.StackedGraph,'ForegroundColor');
    d = get(handles.LineGraph,'ForegroundColor');
    set(handles.LineGraph,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.StackedGraph,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.LineGraph,'Value',1); %was already pressed
end

% --- Executes on button press in StackedGraph.
function StackedGraph_Callback(hObject, eventdata, handles)
if get(handles.StackedGraph,'Value')==1
    a = get(handles.StackedGraph,'BackgroundColor');
    b = get(handles.LineGraph,'BackgroundColor');
    c = get(handles.StackedGraph,'ForegroundColor');
    d = get(handles.LineGraph,'ForegroundColor');
    set(handles.LineGraph,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.StackedGraph,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.StackedGraph,'Value',1); %was already pressed
end

% --- Executes on button press in AutoControl.
function AutoControl_Callback(hObject, eventdata, handles)
if get(handles.AutoControl,'Value')==1
    a = get(handles.ManualControl,'BackgroundColor');
    b = get(handles.AutoControl,'BackgroundColor');
    c = get(handles.ManualControl,'ForegroundColor');
    d = get(handles.AutoControl,'ForegroundColor');
    set(handles.AutoControl,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.ManualControl,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.AutoControl,'Value',1); %was already pressed
end

% --- Executes on button press in ManualControl.
function ManualControl_Callback(hObject, eventdata, handles)
if get(handles.ManualControl,'Value')==1
    a = get(handles.ManualControl,'BackgroundColor');
    b = get(handles.AutoControl,'BackgroundColor');
    c = get(handles.ManualControl,'ForegroundColor');
    d = get(handles.AutoControl,'ForegroundColor');
    set(handles.AutoControl,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.ManualControl,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.ManualControl,'Value',1); %was already pressed
end

% --- Executes on button press in VirtualMode.
function VirtualMode_Callback(hObject, eventdata, handles)
if get(handles.VirtualMode,'Value')==1
    a = get(handles.VirtualMode,'BackgroundColor');
    c = get(handles.VirtualMode,'ForegroundColor');
    if get(handles.ObserverMode,'Value')==1
        b = get(handles.ObserverMode,'BackgroundColor');
        d = get(handles.ObserverMode,'ForegroundColor');
        set(handles.ObserverMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    elseif get(handles.ControllerMode,'Value')==1
        b = get(handles.ControllerMode,'BackgroundColor');
        d = get(handles.ControllerMode,'ForegroundColor');
        set(handles.ControllerMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    end
    set(handles.VirtualMode,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.VirtualMode,'Value',1); %was already pressed
end

% --- Executes on button press in ObserverMode.
function ObserverMode_Callback(hObject, eventdata, handles)
if get(handles.ObserverMode,'Value')==1
    a = get(handles.ObserverMode,'BackgroundColor');
    c = get(handles.ObserverMode,'ForegroundColor');
    if get(handles.VirtualMode,'Value')==1
        b = get(handles.VirtualMode,'BackgroundColor');
        d = get(handles.VirtualMode,'ForegroundColor');
        set(handles.VirtualMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    elseif get(handles.ControllerMode,'Value')==1
        b = get(handles.ControllerMode,'BackgroundColor');
        d = get(handles.ControllerMode,'ForegroundColor');
        set(handles.ControllerMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    end
    set(handles.ObserverMode,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.ObserverMode,'Value',1); %was already pressed
end

% --- Executes on button press in ControllerMode.
function ControllerMode_Callback(hObject, eventdata, handles)
if get(handles.ControllerMode,'Value')==1
    a = get(handles.ControllerMode,'BackgroundColor');
    c = get(handles.ControllerMode,'ForegroundColor');
    if get(handles.ObserverMode,'Value')==1
        b = get(handles.ObserverMode,'BackgroundColor');
        d = get(handles.ObserverMode,'ForegroundColor');
        set(handles.ObserverMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    else
        b = get(handles.VirtualMode,'BackgroundColor');
        d = get(handles.VirtualMode,'ForegroundColor');
        set(handles.VirtualMode,'Value',0,'BackgroundColor',a,'ForegroundColor',c);
    end
    set(handles.ControllerMode,'Value',1,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.ControllerMode,'Value',1); %was already pressed
end

%% Historian/Forecast Tab
function ForecastPlot(handles)
global Plant
S = fieldnames(Plant.Data.Demand);
Ylab = 'Demand (kW)';
PlotIndex = PlotWindow(handles.sliderZoom,handles.sliderDate,Plant.Data.Timestamp);
TimeVec = Plant.Data.Timestamp(PlotIndex(1):PlotIndex(2))';
Steps = round(1/(TimeVec(2)-TimeVec(1)));%points in 1 day of data
if Plant.Data.Timestamp(1)<= TimeVec(1)-1
    PrevDay.Timestamp = Plant.Data.Timestamp(PlotIndex(1)-Steps:PlotIndex(1)-1);
    PrevDay.T = Plant.Data.Temperature(PlotIndex(1)-Steps:PlotIndex(1)-1);
else
    PrevDay.Timestamp = Plant.Data.Timestamp(PlotIndex(1):PlotIndex(2));
    PrevDay.T = Plant.Data.Temperature(PlotIndex(1):PlotIndex(2));
end
%temperature forecast?
%     Y = Plant.Data.Temperature(PlotIndex(1):PlotIndex(2))';
%     PlotData(h,TimeVec,Y,0,'E',Plant.Data,Plant.Data.HistProf,[],[])
for i = 1:1:length(S)
    if strcmp(S{i},'E')
        h = handles.Forecast;
    elseif strcmp(S{i},'H')
        h = handles.HeatFore;
    elseif strcmp(S{i},'C')
        h = handles.CoolFore;
    end
    Y = Plant.Data.Demand.(S{i})(PlotIndex(1):PlotIndex(2))';
    if Plant.Data.Timestamp(1)<= TimeVec(1)-1
        PrevDay.(S{i}) = Plant.Data.Demand.(S{i})(PlotIndex(1)-Steps:PlotIndex(1)-1);
    else
        PrevDay.(S{i}) = Plant.Data.Demand.(S{i})(PlotIndex(1):PlotIndex(2));
    end
    cla(h)
    color = {'k';'g';'r';'b';'c';'m';'y';};
    %% need to go back and understand/clean up create forecast
%     Y(:,2:4) = CreateForecast(S{i},TimeVec(1),((TimeVec-TimeVec(1)).*24)',(length(TimeVec)-1)/Steps,PrevDay,Plant.Data,'HiLow');
    PlotData(h,TimeVec,Y,Ylab,color)
end

function PlotIndex = PlotWindow(zoom,date,timestamp)
a = datevec(timestamp(1));
b = datevec(timestamp(end));
D1 = datenum([a(1) a(2) a(3)]);
months = max(1,12*(b(1)-a(1))+b(2)-a(2));
years = max(1,b(1)-a(1));
days = round(timestamp(end)-timestamp(1));
weeks = ceil(days/7);
Z = max(1,round(get(zoom,'Value')));
if Z == 1
    day = 1+round((get(date,'Value')-1)*(days-1));
    endday = day;
elseif Z ==2
    day = 1+7*round((get(date,'Value')-1)*(weeks-1));
    endday = day+6;
elseif Z ==3
    month = 1+ round((get(date,'Value')-1)*(months-1));
    day = datenum(a(1),month,1)-datenum(a(1),1,1)+1;
    endday = day+(datenum(a(1),month+1,1)-datenum(a(1),month,1)-1);
elseif Z ==4
    year = a(1)+ round((get(date,'Value')-1)*(years-1));
    day = 1+round(datenum(year,1,1)-round(timestamp(1)));
    endday = round(datenum(year+1,1,1)-round(timestamp(1)));
end
PlotIndex(1) = max(1,nnz(timestamp<=(D1+(day-1))));
PlotIndex(2) = nnz(timestamp<=(D1+endday));
if PlotIndex(2)>length(timestamp)
    set(zoom,'Value',Z-1)
    PlotWindow(zoom,date,timestamp)
end   

function sliderDate_Callback(hObject, eventdata, handles)
ForecastPlot(handles)
function sliderDate_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function sliderZoom_Callback(hObject, eventdata, handles)
ForecastPlot(handles)
function sliderZoom_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%% Control Options Tab
function checkboxElectric_Callback(hObject, eventdata, handles)
checkOutputs
function checkboxHeat_Callback(hObject, eventdata, handles)
checkOutputs
function checkboxCooling_Callback(hObject, eventdata, handles)
checkOutputs
function checkboxSteam_Callback(hObject, eventdata, handles)
checkOutputs

function checkOutputs
global Plant
Plant.Optimoptions.Outputs = {};
if get(handles.checkboxElectric,'Value')
    Plant.Optimoptions.Outputs(end+1) = 'E';
end
if get(handles.checkboxHeat,'Value')
    Plant.Optimoptions.Outputs(end+1) = 'H';
end
if get(handles.checkboxCooling,'Value')
    Plant.Optimoptions.Outputs(end+1) = 'C';
end
if get(handles.checkboxSteam,'Value')
    Plant.Optimoptions.Outputs(end+1) = 'S';
end

function changingtimesteps_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'constant'
        Plant.optimoptions.tspacing = 'constant';
    case 'manual'
        Plant.optimoptions.tspacing = 'manual';
        prompt = {'Specify a vector of times out of one horizon for each timestep: (ex: .05,.1,.25,.75,1 would result in one timestep at 5hr, 10hr, 25hr, 75hr, and 100hr for a 100 hour horizon)'};
        dlg_title = 'Manual Timesteps';
        num_lines = 1;
        def_ans = {'0.0035,0.0070,0.0105,0.0140,0.0175,0.0210,0.0417,0.0625,0.0833,0.1042,0.125,0.1667,0.2083,0.25,0.3333,0.4167,0.5,0.5833,0.6667,0.75,0.8333,0.9167,1'};
        a = inputdlg(prompt,dlg_title,num_lines,def_ans);
        Plant.optimoptions.manualT = str2double(strsplit(a{:}, ','));
    case 'linear'
        Plant.optimoptions.tspacing = 'constant';
    case 'logarithm'
        Plant.optimoptions.tspacing = 'logarithm';
end
function Interval_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Interval = str2double(get(handles.Interval, 'String'));
function Interval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Resolution_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Resolution = str2double(get(handles.Resolution, 'String'));
function Resolution_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Horizon_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Horizon = str2double(get(handles.Horizon, 'String'));
function Horizon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function scaletime_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.scaletime = str2double(get(handles.scaletime, 'String'));
function scaletime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in ChillingHeating.
function ChillingHeating_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'sequential'
        Plant.optimoptions.sequential = 1;
    case 'simultaneous'
        Plant.optimoptions.sequential = 0;
end
function excessHeat_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.excessHeat = get(hObject, 'Value');

function nsSmooth_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.nSSmooth = str2double(get(handles.nsSmooth, 'String'));
function nsSmooth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in simulationspeed.
function simulationspeed_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'fastsimulation'
        Plant.optimoptions.fastsimulation = 1;
    case 'slowsimulation'
        Plant.optimoptions.fastsimulation = 0;
end

% --- Executes when selected object is changed in uipanelMixedInteger.
function uipanelMixedInteger_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'NoMixedInteger'
        Plant.optimoptions.MixedInteger = 0;
    case 'MixedInteger'
        Plant.optimoptions.MixedInteger = 1;
end

% --- Executes when selected object is changed in SpinningReserve.
function SpinningReserve_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'noSpinReserve'
        Plant.optimoptions.SpinReserve = false;
    case 'SpinReserve'
        Plant.optimoptions.SpinReserve = true;
        Plant.optimoptions.SpinReservePerc = str2double(get(handles.SpinReservePerc, 'String'));
end
function SpinReservePerc_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.SpinReservePerc = str2double(get(handles.SpinReservePerc, 'String'));
function SpinReservePerc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editBuffer_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Buffer = str2double(get(handles.editBuffer, 'String'));

% --- Executes during object creation, after setting all properties.
function editBuffer_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in ControlDispatch.
function ControlDispatch_SelectionChangeFcn(hObject, eventdata, handles)
global Plant
switch get(eventdata.NewValue,'Tag')
    case 'Control'
        Plant.optimoptions.method = 'Control';
    case 'Dispatch'
        Plant.optimoptions.method = 'Dispatch';
end

function Topt_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Topt = str2double(get(handles.Topt, 'String'));
function Topt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Tmpc_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Tmpc = str2double(get(handles.Tmpc, 'String'));
function Tmpc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NoMixedInteger_Callback(hObject, eventdata, handles)
function MixedInteger_Callback(hObject, eventdata, handles)
function noSpinReserve_Callback(hObject, eventdata, handles)
function SpinReserve_Callback(hObject, eventdata, handles)
function Control_Callback(hObject, eventdata, handles)
function Dispatch_Callback(hObject, eventdata, handles)

%% Communication Tab
%When Component Buttons on the Communication tab are clicked
function GenComm_Callback(hObject, eventdata, handles)
global Plant
genHandles = get(handles.uipanelMain5,'Children');
j = 1;
while ~get(genHandles(j),'Value')
    j = j+1;
end
name = get(genHandles(j),'String');
i = 1;
while ~strcmp(name,Plant.Generator(i).Name)
    i = i+1;
end
set(handles.CurrentComm,'String',name)
if isfield(Plant.Generator(i).VariableStruct,'Comm')
    Comm = Plant.Generator(i).VariableStruct.Comm;
else
    Comm.OnOff = 0;
    Comm.Set = 0;
end
if isfield(Plant.Generator(i).VariableStruct,'Measure')
    Measure = Plant.Generator(i).VariableStruct.Measure;
else
    Measure.OnOff = 0;
    Measure.Input = 0;
    Measure.Electric = 0;
    Measure.Thermal = 0;
end
set(handles.editCommandOnOff,'string',num2str(Comm.OnOff))
set(handles.editCommandSet,'string',num2str(Comm.Set))
set(handles.editMeasureOnOff,'string',num2str(Measure.OnOff))
set(handles.editMeasureInput,'string',num2str(Measure.Input))
set(handles.editMeasurePrimary,'string',num2str(Measure.Electric))
set(handles.editMeasureSecondary,'string',num2str(Measure.Thermal))


function editCommandOnOff_Callback(hObject, eventdata, handles)
recordComm(handles)
function editCommandOnOff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editCommandSet_Callback(hObject, eventdata, handles)
recordComm(handles)
function editCommandSet_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureSecondary_Callback(hObject, eventdata, handles)
recordComm(handles)
function editMeasureSecondary_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasurePrimary_Callback(hObject, eventdata, handles)
recordComm(handles)
function editMeasurePrimary_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureInput_Callback(hObject, eventdata, handles)
recordComm(handles)
function editMeasureInput_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editMeasureOnOff_Callback(hObject, eventdata, handles)
recordComm(handles)

function editMeasureOnOff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function recordComm(handles)
global Plant
nG = length(Plant.Generator);
click = get(handles.CurrentComm,'String');

Comm.OnOff = str2double(get(handles.editCommandOnOff,'String'));
Com.Set = str2double(get(handles.editCommandSet,'String'));

Measure.OnOff = str2double(get(handles.editMeasureOnOff,'String'));
Measure.Input = str2double(get(handles.editMeasureInput,'String'));
Measure.Primary = str2double(get(handles.editMeasurePrimary,'String'));
Measure.Secondary = str2double(get(handles.editMeasureSecondary,'String'));
for j = 1:nG
    if strcmp(click, Plant.Generator(j).Name)%May need to change if multiple generators have the same name
        Plant.Generator(j).VariableStruct.Comm = Comm;
        Plant.Generator(j).VariableStruct.Measure = Measure;
    end
end
