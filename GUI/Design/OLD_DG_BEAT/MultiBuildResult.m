function varargout = MultiBuildResult(varargin)
% MULTIBUILDRESULT M-file for MultiBuildResult.fig
%      MULTIBUILDRESULT, by itself, creates a new MULTIBUILDRESULT or raises the existing
%      singleton*.
%
%      H = MULTIBUILDRESULT returns the handle to a new MULTIBUILDRESULT or the handle to
%      the existing singleton*.
%
%      MULTIBUILDRESULT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTIBUILDRESULT.M with the given input arguments.
%
%      MULTIBUILDRESULT('Property','Value',...) creates a new MULTIBUILDRESULT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MultiBuildResult_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MultiBuildResult_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MultiBuildResult

% Last Modified by GUIDE v2.5 18-Mar-2014 10:00:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MultiBuildResult_OpeningFcn, ...
                   'gui_OutputFcn',  @MultiBuildResult_OutputFcn, ...
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

% --- Executes just before MultiBuildResult is made visible.
function MultiBuildResult_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MultiBuildResult (see VARARGIN)

% Choose default command line output for MultiBuildResult
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

global Model_dir

resultdir=fullfile(Model_dir, 'results'); %strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\result');
x=dir(fullfile(resultdir,'*.mat'));
listTemp=strrep({x.name},'.mat','');
i = 0;
for j = 1:1:length(listTemp')
    if  strncmp(listTemp(j),'multiBuild',10)
        i = i+1;
        list(i) = strrep(listTemp(j),'multiBuild','');
        date(i) = x(j).datenum;
    end
end
[Y, currentProject] = max(date);

set(handles.popupmenuResultsFile,'string',list,'value',currentProject)

set(handles.sliderZoomX1,'Min',0,'Max',1,'Value',1)
set(handles.sliderZoomY1,'Min',0,'Max',1,'Value',1)
set(handles.sliderCenterX1,'Min',0,'Max',1,'Value',.5)
set(handles.sliderCenterY1,'Min',0,'Max',1,'Value',.5)

varList = {'NPV savings %';
            'Payback Period'; 
            'Demand Charge Savings $/kW';
            'Energy Charge Savings $/kW';
            'Added Fuel Costs $/kW'
            'Proportion of Self Generation';
            'FC/CHP Capacity Factor';
            'FC/CHP size'; 
            'Chiller size'; 
            'TES size'; 
            'Electric demand load factor'; 
            'Electric demand std dev'; 
            'Electric to cooling demand ratio';
            'Electric to heating demand ratio'; 
            'Cooling to heating demand ratio';
            'Average Building Load (kW)'};

set(handles.popupmenu1X,'string',varList,'value',6)
set(handles.popupmenu1Y,'string',varList,'value',1)

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

resultdir=fullfile(Model_dir, 'results');  % strrep(which('NREL_FCModel.m'),'\main\NREL_FCModel.m','\result');
load(fullfile(resultdir, strcat('multiBuild',str{val})))
result = multiBuild;

stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
set(handles.popupmenuState, 'string', stateName(:,1), 'value',5)

buildlist =result.BuildType(1);
count = 1;
for i = 2:1:length(result.BuildType)
    if max(strcmp(result.BuildType(i),buildlist))==0
        count = count+1;
        buildlist(count) = result.BuildType(i);
    end
end
if length(buildlist)>1
    buildlist(2:end+1) = buildlist;
    buildlist(1) = cellstr('All');
end
set(handles.popupmenuBuilding, 'string', buildlist, 'value',1)

climlist =result.ClimateType(1);
count = 1;
for i = 2:1:length(result.ClimateType)
    if max(strcmp(result.ClimateType(i),climlist))==0
        count = count+1;
        climlist(count) = result.ClimateType(i);
    end
end
if length(climlist)>1
    climlist(2:end+1) = climlist;
    climlist(1) = cellstr('All');
end
set(handles.popupmenuClimate, 'string', climlist, 'value',1)

statelist =result.State(1);
count = 1;
for i = 2:1:length(result.State)
    if max(strcmp(result.State(i),statelist))==0
        count = count+1;
        statelist(count) = result.State(i);
    end
end
if length(statelist)>1
statelist = [cellstr('All'); statelist'];
end
set(handles.popupmenuState, 'string', statelist, 'value',1)

popupmenuAxes1_Callback(hObject, eventdata, handles)
 


% --- Executes during object creation, after setting all properties.
function popupmenuResultsFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuResultsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Outputs from this function are returned to the command line.
function varargout = MultiBuildResult_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1X.
function popupmenu1X_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenu1X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu1Y.
function popupmenu1Y_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function popupmenu1Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function sliderZoomX1_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZoomX1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function sliderZoomX1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZoomX1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderCenterX1_Callback(hObject, eventdata, handles)
% hObject    handle to sliderCenterX1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function sliderCenterX1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderCenterX1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderZoomY1_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZoomY1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function sliderZoomY1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZoomY1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderCenterY1_Callback(hObject, eventdata, handles)
% hObject    handle to sliderCenterY1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function sliderCenterY1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderCenterY1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbuttonClose.
function pushbuttonClose_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(gcf)


function popupmenuAxes1_Callback(hObject, eventdata, handles)
global result
Building = get(handles.popupmenuBuilding,'Value');
buildTypes = get(handles.popupmenuBuilding,'String');
if strcmp(buildTypes(Building),'All')
    buildTypes = buildTypes(2:end);
else buildTypes =buildTypes(Building);
end
Climate = get(handles.popupmenuClimate,'Value');
climateTypes = get(handles.popupmenuClimate,'String');
if strcmp(climateTypes(Climate),'All')
    climateTypes = climateTypes(2:end);
else climateTypes = climateTypes(Climate);
end
State = get(handles.popupmenuState,'Value');
stateTypes = get(handles.popupmenuState,'String');
if strcmp(stateTypes(Climate),'All')
    stateTypes = stateTypes(2:end);
else stateTypes = stateTypes(State);
end

listValX=get(handles.popupmenu1X,'value');
listValY=get(handles.popupmenu1Y,'value');

[resultX, Fig.xlabel] = result2(listValX,handles);
[resultY, Fig.ylabel] = result2(listValY,handles);
zoom = [get(handles.sliderZoomX1,'Value') get(handles.sliderZoomY1,'Value') get(handles.sliderCenterX1,'Value') get(handles.sliderCenterY1,'Value')];
Fig.Handle =handles.axes1;
% ColorBuild = {'bo-';'go-';'ro-';'co-';'mo-';'bx:';'gx:';'rx:';'cx:';'mx:';'b+-.';'g+-.';'r+-.';'c+-.';'m+-.';'b*--'};
count = 0;
for i = 1:1:length(resultX)
    if max(strcmp(result.BuildType(i),buildTypes))==1 && max(strcmp(result.ClimateType(i),climateTypes))==1 && max(strcmp(result.State(i),stateTypes))==1
        count = count +1;
        Fig.X(count) = resultX(i);
        Fig.Y(count) = resultY(i);
        Fig.color(count,:) = mod(i,16)/15;
        TableBuild(count) = result.BuildType(i);
        TableClimate(count) = result.ClimateType(i);
        TableState(count) = result.State(i);
    end
end
if length(Fig.X) ==   1
    if max(Fig.X)>0
        Fig.xLim = [0 2*max(Fig.X)];
    else Fig.xLim = [2*max(Fig.X) 0];
    end
    if max(Fig.Y)>0
        Fig.yLim = [0 2*max(Fig.Y)];
    else Fig.yLim = [2*max(Fig.Y) 0];
    end
else
    R = (max(Fig.X)-min(Fig.X));
    M = (max(Fig.X)+min(Fig.X))/2+(R-R*zoom(1))*(zoom(3)-.5);
    Fig.xLim = [M-.5*R*zoom(1),M+.5*R*zoom(1)] ;
    R = (max(Fig.Y)-min(Fig.Y));
    M = (max(Fig.Y)+min(Fig.Y))/2+(R-R*zoom(2))*(zoom(4)-.5);
    Fig.yLim = [M-.5*R*zoom(2),M+.5*R*zoom(2)] ;
end
if Fig.xLim(1) == Fig.xLim(2)
    Fig.xLim(1) = Fig.xLim(1)-1;
    Fig.xLim(2) = Fig.xLim(2)+1;
end
if Fig.yLim(1) == Fig.yLim(2)
    Fig.yLim(1) = Fig.yLim(1)-1;
    Fig.yLim(2) = Fig.yLim(2)+1;
end
axes(Fig.Handle);
scatter(Fig.X,Fig.Y,5,Fig.color);
xlabel(Fig.Handle,Fig.xlabel)
ylabel(Fig.Handle,Fig.ylabel)
axis([Fig.xLim Fig.yLim])

Data = cell(count+1,4);
j = 1;
if length(buildTypes)>1
    Data(1,j) =  cellstr('Building');
    Data(2:end,j) = TableBuild';
    j = j+1;
end
if length(climateTypes)>1
    Data(1,j) =  cellstr('Climate');
    Data(2:end,j) = TableClimate';
    j = j+1;
end
if length(stateTypes)>1
    Data(1,j) =  cellstr('State');
    Data(2:end,j) = TableState';
    j = j+1;
end
VarList = get(handles.popupmenu1X,'string');
Data(1,j) =  VarList(listValX);
Data(2:end,j) = cellstr(num2str(Fig.X'));
j = j+1;
Data(1,j) =  VarList(listValY);
Data(2:end,j) = cellstr(num2str(Fig.Y'));
set(handles.uitable1,'Data',Data)



function [Result, label] = result2(listVal,handles)
global result
varList = get(handles.popupmenu1X,'string');

label = varList(listVal);
n = length(result.BuildType);
Result = zeros(n,1);
for i = 1:1:n
    if listVal == 1
        NPVbaseline = sum([result.costOut(i).NPVbaselineDemCharges result.costOut(i).NPVbaselineUseCharges result.costOut(i).NPVbaselineFuelCost result.costOut(i).NPVbaselineOandM result.costOut(i).NPVbaselineFinance]);
        NPVdispatch = sum([result.costOut(i).NPVnewDemCharges result.costOut(i).NPVnewUseCharges result.costOut(i).NPVnewFuelCost result.costOut(i).NPVnewOandM result.costOut(i).NPVnewFinance]);
        Result(i) = (NPVbaseline-NPVdispatch)/NPVbaseline*100;
    elseif listVal ==2
        Result(i) = result.costOut(i).Payback;
    elseif listVal ==3
        Result(i) = (result.costOut(i).Year1baseCharges(1) - result.costOut(i).Year1dispatchCharges(1))/result.eOut(i).SysSize;
    elseif listVal ==4
        Result(i) = (result.costOut(i).Year1baseCharges(2) - result.costOut(i).Year1dispatchCharges(2))/result.eOut(i).SysSize;
    elseif listVal ==5
        Result(i) = (result.costOut(i).Year1dispatchCharges(3) - result.costOut(i).Year1baseCharges(3))/result.eOut(i).SysSize;
    elseif listVal ==6
        Result(i) = result.eOut(i).SelfGen;
    elseif listVal ==7
        Result(i) = result.Dispatch(i).ElecTotProd/(8760*result.eOut(i).SysSize);
    elseif listVal ==8
        Result(i) = result.eOut(i).SysSize;
    elseif listVal ==9
        Result(i) = sum(result.eOut(i).ChillerSize);
    elseif listVal ==10
        Result(i) = result.eOut(i).TESsize;
    elseif listVal ==11
        Result(i) = result.LoadFactor(i);
    elseif listVal ==12
        Result(i) = result.LoadScatter(i);
    elseif listVal ==13
        Result(i) = result.LoadEtoC(i);
    elseif listVal ==14
        Result(i) = result.LoadEtoH(i);
    elseif listVal ==15
        Result(i) = result.LoadCtoH(i);
    elseif listVal ==16
        Result(i) = sum(result.Baseline(i).Elec)/8760;
    end
end

% --- Executes on selection change in popupmenuBuilding.
function popupmenuBuilding_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuBuilding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes1_Callback(hObject, eventdata, handles)
 

% --- Executes during object creation, after setting all properties.
function popupmenuBuilding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuBuilding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuClimate.
function popupmenuClimate_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuClimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes1_Callback(hObject, eventdata, handles)
 

% --- Executes during object creation, after setting all properties.
function popupmenuClimate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuClimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuState.
function popupmenuState_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popupmenuAxes1_Callback(hObject, eventdata, handles)
 

% --- Executes during object creation, after setting all properties.
function popupmenuState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonExport.
function pushbuttonExport_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global result 
