function varargout = SensitivityResult(varargin)
% SENSITIVITYRESULT M-file for SensitivityResult.fig
%      SENSITIVITYRESULT, by itself, creates a new SENSITIVITYRESULT or raises the existing
%      singleton*.
%
%      H = SENSITIVITYRESULT returns the handle to a new SENSITIVITYRESULT or the handle to
%      the existing singleton*.
%
%      SENSITIVITYRESULT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SENSITIVITYRESULT.M with the given input arguments.
%
%      SENSITIVITYRESULT('Property','Value',...) creates a new SENSITIVITYRESULT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SensitivityResult_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SensitivityResult_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SensitivityResult

% Last Modified by GUIDE v2.5 01-Apr-2013 11:07:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SensitivityResult_OpeningFcn, ...
                   'gui_OutputFcn',  @SensitivityResult_OutputFcn, ...
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

% --- Executes just before SensitivityResult is made visible.
function SensitivityResult_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SensitivityResult (see VARARGIN)

% Choose default command line output for SensitivityResult
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
global Model_dir;

resultsdir=fullfile(Model_dir, 'results');% strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\results');
% filename=fullfile(resultsdir,'Project1');
x=dir(fullfile(resultsdir,'*.mat'));
listTemp=strrep({x.name},'.mat','');
i = 0;
for j = 1:1:length(listTemp')
    if  strncmp(listTemp(j),'sensitivity',10)
        i = i+1;
        list(i) = strrep(listTemp(j),'sensitivity','');
    end
end
if ~isempty(find(strcmp('Project1',list)));
    currentProject = find(strcmp('Project1',list));
else currentProject=1;
end
set(handles.popupmenuResultsFile,'string',list,'value',currentProject)

popupmenuResultsFile_Callback(handles.popupmenuResultsFile,eventdata,handles);


% --- Executes on selection change in popupmenuResultsFile.
function popupmenuResultsFile_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuResultsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuResultsFile contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuResultsFile
global result Model_dir
str=get(hObject,'string');
val=get(hObject,'value'); 
resultsdir=fullfile(Model_dir, 'results'); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\results');
load(fullfile(resultsdir, strcat('sensitivity',str{val})))
result = sensitivity;
if result.StudyType == 1 || result.StudyType == 2
varList = {'NPV $/kW demand';
            'NPV Total';
            'NPV savings %';
            'Payback Period'; 
            'Proportion of Self Generation';
            'FC/CHP Capacity Factor';
            'FC/CHP size'; 
            'Chiller size'; 
            'TES size';
            'NPV $ / kWhr';
            };
else varList = {'Total NPV $';
            'Payback Period'; 
            'NPV savings %';};
end

set(handles.popupmenu1,'string',varList,'value',1)
set(handles.popupmenu2,'string',varList,'value',1)

plotAxes(1,handles)
plotAxes(2,handles)



% --- Executes during object creation, after setting all properties.
function popupmenuResultsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuResultsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Outputs from this function are returned to the command line.
function varargout = SensitivityResult_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAxes(1,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotAxes(2,handles)

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonUpdate.
function pushbuttonUpdate_Callback(hObject, eventdata, handles)
global vintage climate buildType filter result
filter = filterData(result.BuildName, vintage, climate, buildType);
plotAxes(1,handles)
plotAxes(2,handles)

% --- Executes on button press in pushbuttonClose.
function pushbuttonClose_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)


function plotAxes(fig,handles)
global result 

if fig==1
    Fig.Handle =handles.axes1;
    listVal=get(handles.popupmenu1,'value');
else Fig.Handle = handles.axes2;
    listVal=get(handles.popupmenu2,'value');
end

if result.StudyType == 1 || result.StudyType == 2
varList = {'NPV $/kW demand';
            'NPV Total';
            'NPV savings %';
            'Payback Period'; 
            'Proportion of Self Generation';
            'FC/CHP Capacity Factor';
            'FC/CHP size'; 
            'Chiller size'; 
            'TES size';
            'NPV $ / kWhr';
            };
else varList = {'Total NPV $';
            'Payback Period'; 
            'NPV savings %';};
end

label = char(varList(listVal));

switch label
    case 'NPV $/kW demand'
        Fig.Y = [result.NPVdispatch(:,1)./result.AvgDemand' result.NPVdispatch(:,2)./result.AvgDemand' result.NPVdispatch(:,3)./result.AvgDemand' result.NPVdispatch(:,4)./result.AvgDemand' result.NPVdispatch(:,5)./result.AvgDemand'];
        Fig.ylabel=('Fuel, Grid, and O&M/Finance Costs ($/kW)');
    case 'NPV Total'
        Fig.Y = [result.NPVdispatch(:,1)  result.NPVdispatch(:,2)  result.NPVdispatch(:,3)  result.NPVdispatch(:,4)  result.NPVdispatch(:,5) ];
        Fig.ylabel=('Fuel, Grid, O&M, and Finance Costs ($)');
    case 'Total NPV $'
        Fig.Y = sum(result.NPVdispatch,2);
        Fig.ylabel=('Fuel, Grid, and O&M/Finance Costs ($)');
    case 'NPV savings %'
        Fig.Y = (1-sum(result.NPVdispatch,2)./sum(result.NPVbaseline,2))*100;
        Fig.ylabel=('% Reduction in NPV from baseline');
    case 'Payback Period'
        Fig.Y = result.Payback;
        Fig.ylabel=('Payback Period (Yrs)');
    case 'Proportion of Self Generation'
        Fig.Y = result.SelfGen;
        Fig.ylabel=('Proportion of Self Generation');
    case 'FC/CHP Capacity Factor'
        Fig.Y = result.CapacityFactor;
        Fig.ylabel=('FC/CHP Capacity Factor');
    case 'FC/CHP size'
        Fig.Y = result.CapacitykW;
        Fig.ylabel=('FC/CHP size (kW)');
    case 'Chiller size'
        Fig.Y = result.ChillerSize;
        Fig.ylabel=('Chiller size (kW)');
    case 'TES size'
        Fig.Y = result.TESsize;
        Fig.ylabel=('TES size (kWh)');
    case 'NPV $ / kWhr'
        Fig.Y = sum(result.NPVdispatch,2)./result.AvgDemand'/175200;
        Fig.ylabel=('NPV $ / kWhr');
end
axes(Fig.Handle);
legpos{1}=[0.180 0.741 0.157 0.148];
legpos{2}=[0.694 0.745 0.152 0.141];
if result.StudyType == 1
    bar(Fig.Y,'stacked')
    set(gca,'XTickLabel',result.xTicks, 'xtick', 1:17)
    rotateticklabel(gca, 90, 10);
    if strcmp(label,'NPV $/kW demand')
        l=legend('Demand Charges', 'Grid Energy Charges', 'Fuel Cost', 'O&M', ...
            'Financing');
        set(l, 'position', legpos{fig} , 'fontsize', 8);
    end
else 
    if strcmp(char(varList(listVal)),'NPV $/kW demand') || strcmp(char(varList(listVal)),'NPV Total')
        Fig.Y = sum(Fig.Y,2);
    end
    plot(result.X,Fig.Y)
%     if strcmp(char(varList(listVal)),'NPV $/kW demand')
%         legend('Demand Charges', 'Grid Energy Charges', 'Fuel Cost', 'O&M', 'Financing')
%     end
end
xlabel(Fig.Handle,result.xlabel)
ylabel(Fig.Handle,Fig.ylabel)

function th=rotateticklabel(h,rot,fontSize)
%set the default rotation if user doesn't specify
if nargin==1
    rot=90;
end
%make sure the rotation is in the range 0:360 (brute force method)
while rot>360
    rot=rot-360;
end
while rot<0
    rot=rot+360;
end
%get current tick labels
a=get(h,'XTickLabel');
%erase current tick labels from figure
set(h,'XTickLabel',[]);
%get tick label positions
b=get(h,'XTick');
c=get(h,'YTick');
yLimits=get(gca,'ylim');
%make new tick labels
if ~exist('fontSize', 'var')
    fontSize=20;
end
if rot<180
    %th=text(b,repmat(c(1)-.1*(c(2)-c(1)),length(b),1),a,'HorizontalAlignment','right','rotation',rot,...
     %   'fontsize', fontSize);
     th=text(b,repmat(yLimits(1)-0.15*(diff(yLimits)),length(b),1),a,'HorizontalAlignment','right','rotation',rot,...
        'fontsize', fontSize);
else
    th=text(b,repmat(c(1)-.1*(c(2)-c(1)),length(b),1),a,'HorizontalAlignment','left','rotation',rot,...
        'fontsize', fontSize);
end





