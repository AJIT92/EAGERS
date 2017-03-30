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
 
% Last Modified by GUIDE v2.5 27-Feb-2017 21:43:12
 
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
Plant.optimptions.method = 'Dispatch';
%Set up Control Spec values
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

if ismember('E',Plant.optimoptions.Outputs)
    set(handles.checkboxElectric,'value',1)
else
    set(handles.checkboxElectric,'value',0)
    set(handles.MainTop,'Visible','off')
    set(handles.ElectricGraph,'Visible','off')
    set(handles.Forecast,'Visible','off')
    %% move and re-size heat or cooling dispatch axes & title
end
if ismember('H',Plant.optimoptions.Outputs)
    set(handles.checkboxHeat,'Value',1)
else
    set(handles.checkboxHeat,'Value',0)
    set(handles.SideTop,'Visible','off')
    set(handles.HeatGraph,'Visible','off')
    set(handles.HeatFore,'Visible','off')
end
if ismember('C',Plant.optimoptions.Outputs)
    set(handles.checkboxCooling,'Value',1)
else
    set(handles.checkboxCooling,'Value',0)
    set(handles.SideBottom,'Visible','off')
    set(handles.CoolGraph,'Visible','off')
    set(handles.CoolFore,'Visible','off')
end
if ismember('S',Plant.optimoptions.Outputs)
    set(handles.checkboxSteam,'Value',1) 
else
    set(handles.checkboxSteam,'Value',0)
end

%Set handles for communications
set(handles.editCommandOnOff,'string','--')
set(handles.editCommandSet,'string','--')
set(handles.editMeasureOnOff,'string','--')
set(handles.editMeasureInput,'string','--')
set(handles.editMeasureElectric,'string','--')
set(handles.editMeasureThermal,'string','--')

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
[x,map] = imread(fullfile(Model_dir,'GUI','Graphics','swap.png'));
s = imresize(x,[48 48]);
set(handles.SwitchChart,'cdata',s)

handles = GenList_Make(handles);
set(handles.uipanelMain1,'Visible','off')
set(handles.uipanelMain5,'Visible','on')
handles = GenList_Make(handles);
set(handles.uipanelMain5,'Visible','off')
set(handles.uipanelMain1,'Visible','on')

Plant.GUIhandles = handles;


%Update forecast handles
days = round(Plant.Data.Timestamp(end)-Plant.Data.Timestamp(1));
set(handles.sliderZoom,'Min',1,'Max',4,'Value',1,'SliderStep',[1/3,1/3])
set(handles.sliderDate,'Min',1,'Max',2,'Value',1,'SliderStep',[1/(days-1),1/(days-1)])
set(handles.sliderDate,'Max',2)

%update optimization method
P = path;
if strfind(P,fullfile(Model_dir,'Optimization','ComplementaryQP'))
    N = {'cQP';'NN';'Network';};
elseif strfind(P,fullfile(Model_dir,'Optimization','NeuralNetwork'))
    N = {'NN';'cQP';'Network';};
elseif strfind(P,fullfile(Model_dir,'Optimization','NetworkQP'))
    N = {'Network';'NN';'cQP';};
else
    addpath(fullfile(Model_dir,'Optimization','NetworkQP'))
    N = {'Network';'NN';'cQP';};
end
if get(handles.(N{1}),'Value')==1
    %do nothing
else
    a = get(handles.(N{1}),'BackgroundColor');
    c = get(handles.(N{1}),'ForegroundColor');
    if get(handles.(N{2}),'Value')==1
        b = get(handles.(N{2}),'BackgroundColor');
        d = get(handles.(N{2}),'ForegroundColor');
        set(handles.(N{2}),'Value',0,'BackgroundColor',a,'ForegroundColor',c);
        set(handles.(N{1}),'Value',0,'BackgroundColor',b,'ForegroundColor',d);
    elseif get(handles.(N{3}),'Value')==1
        b = get(handles.(N{3}),'BackgroundColor');
        d = get(handles.(N{3}),'ForegroundColor');
        set(handles.(N{3}),'Value',0,'BackgroundColor',a,'ForegroundColor',c);
        set(handles.(N{1}),'Value',0,'BackgroundColor',b,'ForegroundColor',d);
    end
end

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
TabText = {'Main Window';'Running Costs';'Forecasting';'Control Spec';'Communication'};
set(hObject,'UserData',TabText);

for i = 1:length(TabText)
    j = num2str(i);
    % panel management
    set(handles.(strcat('MainTab',j)),'Units','normalized','String',TabText{i});
    set(handles.(strcat('uipanelMain',j)),'Units','normalized','BorderType','none')
    if i ==1
        pan1pos = get(handles.uipanelMain1,'Position');
    else
        set(handles.(strcat('uipanelMain',j)),'Position',pan1pos)
        set(handles.(strcat('uipanelMain',j)),'Visible','off')
    end
end

% Main tabs callback
function mainTab_Callback(hObject, eventdata, handles)
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
set(handles.(strcat('MainTab',m)),'Position',[pos(1),pos(2),pos(3),pos(4)-.003])
pos = get(handles.(strcat('MainTab',n)),'Position');
set(handles.(strcat('MainTab',n)),'Position',[pos(1),pos(2),pos(3),pos(4)+.003])

% Change visibility
set(handles.(strcat('uipanelMain',m)),'Visible','off')
set(handles.(strcat('uipanelMain',n)),'Visible','on')

%% actions specific to one tab or another
if strcmp(n,'3')
    ForecastPlot(handles)
end

% --- Executes on button press in SwitchChart.
function SwitchChart_Callback(hObject, eventdata, handles)
global Plant
order = Plant.optimoptions.Outputs;
n = get(handles.MainTop,'String');
if strcmp(n,'Electric Dispatch')
    m = get(handles.ElectricGraph,'Position');
elseif strcmp(n,'Heating Dispatch')
    m = get(handles.HeatGraph,'Position');
elseif strcmp(n,'Cooling Dispatch')
    m = get(handles.CoolGraph,'Position');
end

if m(4) > 23%keeps the same size graphs if simulation hasn't been run yet
    main = [35.625,19.7619,90,23.9524];
    sidetop = [143,24.095,60.00,17.3];
    sidebot = [143,2.571,60.00,17.3];
else
    main = [35.625,19.9,90,20];
    sidetop = [143,24,60.00,14.2];
    sidebot = [143,3,60.00,14.2];
end

if length(order) == 3
    if strcmp(get(handles.MainTop,'String'),'Electric Dispatch')
        set(handles.MainTop,'String','Heating Dispatch');
        set(handles.SideTop,'String','Cooling Dispatch');
        set(handles.SideBottom,'String','Electric Dispatch');
        set(handles.ElectricGraph,'Position',sidebot);
        set(handles.HeatGraph,'Position',main);
        set(handles.CoolGraph,'Position',sidetop);
    elseif strcmp(get(handles.MainTop,'String'),'Heating Dispatch')
        set(handles.MainTop,'String','Cooling Dispatch');
        set(handles.SideTop,'String','Electric Dispatch');
        set(handles.SideBottom,'String','Heating Dispatch');
        set(handles.ElectricGraph,'Position',sidetop);
        set(handles.HeatGraph,'Position',sidebot);
        set(handles.CoolGraph,'Position',main);
    elseif strcmp(get(handles.MainTop,'String'),'Cooling Dispatch')
        set(handles.MainTop,'String','Electric Dispatch');
        set(handles.SideTop,'String','Heating Dispatch');
        set(handles.SideBottom,'String','Cooling Dispatch');
        set(handles.ElectricGraph,'Position',main);
        set(handles.HeatGraph,'Position',sidetop);
        set(handles.CoolGraph,'Position',sidebot);
    end
elseif length(order) ==2
    if strcmp(get(handles.MainTop,'String'),'Electric Dispatch')
        if strcmp(get(handles.SideTop,'String'),'Heating Dispatch')
            set(handles.MainTop,'String','Heating Dispatch');
            set(handles.SideTop,'String','Electric Dispatch');
            set(handles.ElectricGraph,'Position',sidetop);
            set(handles.HeatGraph,'Position',main);
        else
            set(handles.MainTop,'String','Cooling Dispatch');
            set(handles.SideTop,'String','Electric Dispatch');
            set(handles.CoolGraph,'Position',main);
            set(handles.ElectricGraph,'Position',sidetop);
        end
    elseif strcmp(get(handles.MainTop,'String'),'Heating Dispatch')
        if strcmp(get(handles.SideTop,'String'),'Electric Dispatch')
            set(handles.MainTop,'String','Electric Dispatch');
            set(handles.SideTop,'String','Heating Dispatch');
            set(handles.ElectricGraph,'Position',main);
            set(handles.HeatGraph,'Position',sidetop);
        else
            set(handles.MainTop,'String','Cooling Dispatch');
            set(handles.SideTop,'String','Heating Dispatch');
            set(handles.CoolGraph,'Position',main);
            set(handles.HeatGraph,'Position',sidetop);
        end
    elseif strcmp(get(handles.MainTop,'String'),'Cooling Dispatch')
        if strcmp(get(handles.SideTop,'String'),'Heating Dispatch')
            set(handles.MainTop,'String','Heating Dispatch');
            set(handles.SideTop,'String','Cooling Dispatch');
            set(handles.HeatGraph,'Position',main);
            set(handles.CoolGraph,'Position',sidetop);
        else
            set(handles.MainTop,'String','Electric Dispatch');
            set(handles.SideTop,'String','Cooling Dispatch');
            set(handles.ElectricGraph,'Position',main);
            set(handles.CoolGraph,'Position',sidetop);
        end
    end
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

%Makes buttons for GUI && status buttons next to corresponding generator
nG = length(list);
colorVec = Plant.Plotting.ColorMaps{1};
colorsPlot = interp1(linspace(0,1,length(colorVec)),colorVec,linspace(0,1,nG));
for i=1:min(10,nG) %Set first 10 generators visible
    num = num2str(i);
    pos = 567 - 30*(i-1);
    name = strcat(char(39),'Gen',r,'_Callback',char(39));
    callback = strcat('@(hObject,eventdata)DISPATCH(',name,',hObject,eventdata,guidata(hObject))');
    btn = uicontrol('Style', 'pushbutton', 'String', list{i},...
    'Position', [3 pos 100 21],...
    'Tag', strcat('Generator',num),...
    'FontSize', 10,...
    'Parent', handles.(strcat('uipanelMain',num2str(p))),...
    'Callback',eval(callback),...
    'Visible','on',...
    'UserData',i);
    handles.(strcat('Generator',num)) = btn;
    set(handles.(strcat('Generator',num)),'BackgroundColor',colorsPlot(i,:))
    if p ==1 %Only make Status buttons on Main Window
        pos = 569 - (30*(i-1));
        name = strcat(char(39),'Status_Callback',char(39));
        callback = strcat('@(hObject,eventdata)DISPATCH(',name,',hObject,eventdata,guidata(hObject))');
        if Plant.Generator(i).Enabled
            [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
        else
            [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
        end
        s = imresize(x,[18 18]);
        if Plant.Generator(i).Enabled
            enableGen  = 'bold';
        else enableGen  = 'normal';
        end
        btn = uicontrol('Style', 'pushbutton', 'String', '',...
        'Position', [116 pos 18 18],...
        'Tag', strcat('GeneratorStat',num),...
        'cdata', s,...
        'FontWeight',enableGen,...
        'Parent', handles.uipanelMain1,...
        'Callback',eval(callback),...
        'Visible','on',...
        'UserData',i);
        handles.(strcat('GeneratorStat',num)) = btn;
    end
end
if length(list)>10
    set(handles.NextGen1,'Visible','on')
    set(handles.NextGen5,'Visible','on')
    tab = ceil(length(list)/10);%Determine how many sets of 10 generators need to be displayed
    curtab = 2;
    prev = 10;
    %%%%% Make buttons for remaining generators, but have visibility off
    for j = curtab:tab
        for i = (prev+1):min(curtab*10,length(list))%For every set of 10, creates the buttons back to the top
            num = num2str(i);
            pos = 567 - 30*(i-(prev+1));
            name = strcat(char(39),'Gen',r,'_Callback',char(39));
            callback = strcat('@(hObject,eventdata)DISPATCH(',name,',hObject,eventdata,guidata(hObject))');
            btn = uicontrol('Style', 'pushbutton', 'String', list{i},...
            'Position', [23 pos 100 21],...
            'Tag', strcat('Generator',num),...
            'FontSize', 10,...
            'Parent', handles.(strcat('uipanelMain',num2str(p))),...
            'Callback',eval(callback),...
            'Visible','off',...
            'UserData',i);
            handles.(strcat('Generator',num)) = btn;
            if p ==1 %Only make Status buttons on Main Window
                pos = 569 - 30*(i-(prev+1));
                name = strcat(char(39),'Status_Callback',char(39));
                callback = strcat('@(hObject,eventdata)DISPATCH(',name,',hObject,eventdata,guidata(hObject))');
                if Plant.Generator(i).Enabled
                    [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
                else
                    [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
                end
                s = imresize(x,[18 18]);
                if Plant.Generator(i).Enabled
                    enableGen  = 'bold';
                else enableGen  = 'normal';
                end
                btn = uicontrol('Style', 'pushbutton', 'String', '',...
                'Position', [116 pos 18 18],...
                'Tag', strcat('GeneratorStat',num),...
                'cdata', s,...
                'FontWeight',enableGen,...
                'Parent', handles.uipanelMain1,...
                'Callback',eval(callback),...
                'Visible','off',...
                'UserData',i);
                handles.(strcat('GeneratorStat',num)) = btn;
            end
        end
        prev = i;
        curtab = curtab+1;
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
for i = 1:(length(list)*mult)%Find out last visible generator
    if strcmp(get(handles.(panel).Children(i),'Visible'),'on')
        vis = get(handles.(panel).Children(i),'UserData');
        last = i;%last visible child
        break
    end
end
if vis == length(list) 
    r = floor(vis/10)*10;%round down to nearest 10
    if (vis-r)<10
        start = mult*(vis-r);%only turn off the visible, if less then 10 are showing
    end
else
    start = last + (10*mult-1);
end
new = start + (10*mult);
for i = last:start%current visible buttons off
    set(handles.(panel).Children(i),'Visible','off')
end

for i = start+1:new%current visible buttons off
    set(handles.(panel).Children(i),'Visible','on')
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
for i = 1:(length(list)*2)%Find out last visible generator
    if strcmp(get(handles.(panel).Children(i),'Visible'),'on')
        vis = get(handles.(panel).Children(i),'UserData');
        last = i;%last visible child
        break
    end
end
dif = length(list)-vis;
for i = last:(mult*length(list))%current visible buttons off
    set(handles.(panel).Children(i),'Visible','off')
end
if dif<10 || (dif==10 && length(list)==10)%For less than (or exactly 10 left) visible components
    set(handles.(button2),'Visible','off')
    set(handles.(button1),'Visible','on')
    for j = 1:(mult*dif)%show rest of buttons 
        set(handles.(panel).Children(j),'Visible','on')
    end
else
    set(handles.(button1),'Visible','on')
    set(handles.(button2),'Visible','on')
    n=dif-10;%For plants with more components
    for j = mult*n+1:last-1%show next 10 buttons 
        set(handles.(panel).Children(j),'Visible','on')
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
    set(handles.GenEnable,'Value',0)
    set(handles.GenDisable,'Value',1)
    Plant.Generator(i).Enabled = 0;
    [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
    s = imresize(x,[18 18]);
    set(hObject,'FontWeight','normal','cdata',s)
else
    set(handles.GenEnable,'value',1)
    set(handles.GenDisable,'Value',0)
    Plant.Generator(i).Enabled = 1;
    [x,map] = imread(fullfile(Model_dir,'GUI','Graphics','green.png'));
    s = imresize(x,[18 18]);
    set(hObject,'FontWeight','bold','cdata',s)
end

% --- Executes on button press in Start.
function Start_Callback(hObject, eventdata, handles)
global Virtual RealTime Model_dir DispatchWaitbar 
Virtual = 1;
RealTime = 0;
P = path;
if strfind(P,fullfile(Model_dir,'Optimization','ComplementaryQP'))
    if get(handles.NN,'Value')==1 || get(handles.Network,'Value')==1
        rmpath(fullfile(Model_dir,'Optimization','ComplementaryQP'));
    end
end
if strfind(P,fullfile(Model_dir,'Optimization','NeuralNetwork'))
    if get(handles.NN,'Value')==1 || get(handles.cQP,'Value')==1
        rmpath(fullfile(Model_dir,'Optimization','NeuralNetwork'));
    end
end
if strfind(P,fullfile(Model_dir,'Optimization','NeuralNetwork'))
    if get(handles.cQP,'Value')==1 || get(handles.Network,'Value')==1
        rmpath(fullfile(Model_dir,'Optimization','NeuralNetwork'));
    end
end
if get(handles.cQP,'Value')==1
    addpath(fullfile(Model_dir,'Optimization','ComplementaryQP'));
elseif get(handles.NN,'Value')==1
    addpath(fullfile(Model_dir,'Optimization','NeuralNetwork'));
elseif get(handles.Network,'Value')==1
    addpath(fullfile(Model_dir,'Optimization','NetworkQP'));   
end
DispatchWaitbar=waitbar(0,'Running Dispatch','Visible','off');
RunOptimization
waitfor(DispatchWaitbar)


% --- Executes on button press in Stop.
function Stop_Callback(hObject, eventdata, handles)
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
s = imresize(x,[18 18]);
num = 2*length(Plant.Generator) - (2*i - 1);
set(handles.uipanelMain1.Children(num),'FontWeight','normal','cdata',s)

% --- Executes on button press in GenDisable.
function GenDisable_Callback(hObject, eventdata, handles)
global Plant Model_dir
i = get(handles.SelGen,'UserData');
Plant.Generator(i).Enabled = 0;
[x,map] = imread(fullfile(Model_dir,'GUI','Graphics','red.png'));
s = imresize(x,[18 18]);
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

%% Forecast Tab

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


%% Optimization Options/ Control Spec
function initialize_gui(fig_handle, handles, isreset)
function scaletime_Callback(hObject, eventdata, handles)
function scaletime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function nsSmooth_Callback(hObject, eventdata, handles)
function nsSmooth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function checkboxElectric_Callback(hObject, eventdata, handles)
function checkboxHeat_Callback(hObject, eventdata, handles)
function checkboxCooling_Callback(hObject, eventdata, handles)
function checkboxSteam_Callback(hObject, eventdata, handles)


% --- Executes on button press in cQP.
function cQP_Callback(hObject, eventdata, handles)
a = [0,0.8,0.4];
b = [1,1,1];
c = [0,0,0];
d = [0.4,0.4,0.4];
if get(handles.cQP,'Value')==1
    set(handles.cQP,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.NN,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
    set(handles.Network,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.cQP,'Value',1); %was already pressed
end

% --- Executes on button press in NN.
function NN_Callback(hObject, eventdata, handles)
a = [0,0.8,0.4];
b = [1,1,1];
c = [0,0,0];
d = [0.4,0.4,0.4];
if get(handles.NN,'Value')==1
    set(handles.cQP,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
    set(handles.NN,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
    set(handles.Network,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
else set(handles.NN,'Value',1); %was already pressed
end

% --- Executes on button press in Network.
function Network_Callback(hObject, eventdata, handles)
a = [0,0.8,0.4];
b = [1,1,1];
c = [0,0,0];
d = [0.4,0.4,0.4];
if get(handles.Network,'Value')==1
    set(handles.cQP,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
    set(handles.NN,'Value',0,'BackgroundColor',b,'ForegroundColor',d);
    set(handles.Network,'Value',1,'BackgroundColor',a,'ForegroundColor',c);
else set(handles.Network,'Value',1); %was already pressed
end


% --- Executes on button press in SAVE.
function SAVE_Callback(hObject, eventdata, handles)
global Plant
Plant.optimoptions.Interval = str2double(get(handles.Interval, 'String'));
Plant.optimoptions.Topt = str2double(get(handles.Topt, 'String'));
Plant.optimoptions.Tmpc = str2double(get(handles.Tmpc, 'String'));
Plant.optimoptions.nsSmooth = str2double(get(handles.nsSmooth, 'String'));
Plant.optimoptions.scaletime = str2double(get(handles.scaletime, 'String'));
Plant.optimoptions.thresholdSteps = str2double(get(handles.thresholdSteps, 'String'));
if get(handles.excessHeat, 'Value') ==1
    Plant.optimoptions.excessHeat = 1;%excessHeat is 1 if you can produce escess heat
else Plant.optimoptions.excessHeat = 0;%it is zero if all heat produced must be used
end
if get(handles.fastsimulation,'Value')==1
    Plant.optimoptions.fastsimulation = 1;
else Plant.optimoptions.fastsimulation = 0;
end
if get(handles.sequential,'Value')==1
    %% Only valid question if a) there are chillers and electric generators & b) there is cold thermal storage so that the chiller dispatch affects the electric dispatch
    %% Combined has the drawback of only single value for chiller efficiency (like hratio) rather than a quadratic fit
    Plant.optimoptions.sequential = 1;
else Plant.optimoptions.sequential = 0;
end
%% also build the list of outputs to consider (i.e. E, H, S, H2, C for electricity, heat steam hydrogen, cooling..)
Plant.optimoptions.Outputs = {};
if get(handles.checkboxElectric,'Value')==1
    Plant.optimoptions.Outputs(end+1) ={'E'}; 
end
if get(handles.checkboxHeat,'Value')==1
    Plant.optimoptions.Outputs(end+1) ={'H'}; 
end
if get(handles.checkboxCooling,'Value')==1
    Plant.optimoptions.Outputs(end+1) ={'C'}; 
end
if get(handles.checkboxSteam,'Value')==1
    Plant.optimoptions.Outputs(end+1) ={'S'}; 
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

function excessHeat_Callback(hObject, eventdata, handles)
function Interval_Callback(hObject, eventdata, handles)
function Interval_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Resolution_Callback(hObject, eventdata, handles)
Plant.optimoptions.Resolution = str2double(get(handles.Resolution, 'String'));
function Resolution_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Horizon_Callback(hObject, eventdata, handles)
Plant.optimoptions.Horizon = str2double(get(handles.Horizon, 'String'));

function Horizon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Topt_Callback(hObject, eventdata, handles)
function Topt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Tmpc_Callback(hObject, eventdata, handles)
function Tmpc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function thresholdSteps_Callback(hObject, eventdata, handles)
function thresholdSteps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Communication Tab

%When Component Buttons on the Communication tab are clicked
function GenComm_Callback(hObject, eventdata, handles)
global Plant Generator Comm Measure EditCommHandle
nG = length(Plant.Generator);
genHandles = get(handles.uipanelMain5,'Children');
EditCommHandle = handles.uipanelMain5;
for j = 1:nG
    click = get(genHandles(j),'Value');
    if click == 1
        for i = 1:length(Plant.Generator)
            gen = get(genHandles(j),'String');
            if strcmp(gen,Plant.Generator(i).Name)
                set(handles.CurrentComm,'String',gen)
                if isfield(Plant.Generator(i).VariableStruct,'Comm')
                    Comm = Plant.Generator(i).VariableStruct.Comm;
                else Comm.OnOff = 0;
                    Comm.Set = 0;
                end
                if isfield(Plant.Generator(i).VariableStruct,'Measure')
                    Measure = Plant.Generator(i).VariableStruct.Measure;
                else Measure.OnOff = 0;
                    Measure.Input = 0;
                    Measure.Electric = 0;
                    Measure.Thermal = 0;
                end
            set(handles.editCommandOnOff,'string',num2str(Comm.OnOff))
            set(handles.editCommandSet,'string',num2str(Comm.Set))
            set(handles.editMeasureOnOff,'string',num2str(Measure.OnOff))
            set(handles.editMeasureInput,'string',num2str(Measure.Input))
            set(handles.editMeasureElectric,'string',num2str(Measure.Electric))
            set(handles.editMeasureThermal,'string',num2str(Measure.Thermal))
            end
        end
    end
end
function editCommandOnOff_Callback(hObject, eventdata, handles)
global Comm
Comm.OnOff = str2double(get(hObject,'String'));
function editCommandOnOff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editCommandSet_Callback(hObject, eventdata, handles)
global Comm
Comm.Set = str2double(get(hObject,'String'));
function editCommandSet_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editMeasureThermal_Callback(hObject, eventdata, handles)
global Measure
Measure.Thermal = str2double(get(hObject,'String'));
function editMeasureThermal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editMeasureElectric_Callback(hObject, eventdata, handles)
global Measure
Measure.Electric = str2double(get(hObject,'String'));
function editMeasureElectric_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editMeasureInput_Callback(hObject, eventdata, handles)
global Measure
Measure.Input = str2double(get(hObject,'String'));
function editMeasureInput_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function editMeasureOnOff_Callback(hObject, eventdata, handles)
global Measure
Measure.OnOff = str2double(get(hObject,'String'));
function editMeasureOnOff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pushbuttonFinish_Callback(hObject, eventdata, handles)
global Plant Comm Measure
nG = length(Plant.Generator);
genHandles = get(handles.uipanelMain5,'Children');
click = get(handles.CurrentComm,'String');
for j = 1:nG
    if strcmp(click, Plant.Generator(j).Name)%May need to change if multiple generators have the same name
        Plant.Generator(j).VariableStruct.Comm = Comm;
        Plant.Generator(j).VariableStruct.Measure = Measure;
    end
end
