function varargout = Economics(varargin)
% ECONOMICS MATLAB code for Economics.fig
%      ECONOMICS, by itself, creates a new ECONOMICS or raises the existing
%      singleton*.
%
%      H = ECONOMICS returns the handle to a new ECONOMICS or the handle to
%      the existing singleton*.
%
%      ECONOMICS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ECONOMICS.M with the given input arguments.
%
%      ECONOMICS('Property','Value',...) creates a new ECONOMICS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Economics_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Economics_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Economics

% Last Modified by GUIDE v2.5 17-Mar-2014 15:46:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Economics_OpeningFcn, ...
                   'gui_OutputFcn',  @Economics_OutputFcn, ...
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


% --- Executes just before Economics is made visible.
function Economics_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Economics (see VARARGIN)

% Choose default command line output for Economics
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global Project COMPONENT 
COMPONENT=Project.Economic;

CHP=Project.System.CHP;
COMPONENT.CHPSize = 0;
for i = 1:1:length(CHP)
    COMPONENT.CHPSize = COMPONENT.CHPSize + CHP(i).SysSize(1);
end

COMPONENT.Payback = Project.Result.costOut.Payback;
if isempty(COMPONENT.LifekWh)
    COMPONENT.LifekWh = COMPONENT.CHPSize*8760*COMPONENT.LifeYrs;
end
set(handles.editName,'string',COMPONENT.Name)
set(handles.editInstall,'string',round(COMPONENT.InstallCost))
set(handles.editCHP_OM,'string',round(COMPONENT.CHP_OM))
set(handles.editElecChill,'string',round(COMPONENT.ElecChill))
set(handles.editChillOM,'string',COMPONENT.ChillOM)
set(handles.editColdStorage,'string',COMPONENT.ColdStore)
set(handles.editColdOM,'string',COMPONENT.ColdOM)
set(handles.editBattery,'string',COMPONENT.Battery)
set(handles.editBatteryOM,'string',COMPONENT.BatteryOM)
set(handles.editAbsorbChill,'string',COMPONENT.AbsorbChill)
set(handles.checkboxAbsorpChill,'value',COMPONENT.CurveFit)
set(handles.editAbChillOM,'string',COMPONENT.AbChillOM)
set(handles.editSolar,'string',COMPONENT.Solar)
set(handles.editSolarOM,'string',COMPONENT.SolarOM)
set(handles.editWind,'string',COMPONENT.Wind)
set(handles.editWindOM,'string',COMPONENT.WindOM)
set(handles.editFinance,'string',COMPONENT.FinanceYrs)
set(handles.editDepreciate,'string',COMPONENT.DepreciateYrs)
set(handles.editInflation,'string',COMPONENT.Inflation)
set(handles.editInterest,'string',COMPONENT.Interest)
set(handles.editPayback,'string',num2str(COMPONENT.Payback,2))
set(handles.editLifeYrs,'string',COMPONENT.LifeYrs)
set(handles.editLifekWh,'string',num2str(COMPONENT.LifekWh,2))
set(handles.editReplace,'string',COMPONENT.StackReplaceCost)
set(handles.editIncentive,'string',COMPONENT.Incentive)
set(handles.editSellBack,'string',Project.Utilities.Grid.SellBackRate)
set(handles.uitableProduction,'Data',COMPONENT.ProdCosts);
axes(findobj(gcf,'tag','axes1'))
plot(COMPONENT.ProdCosts(:,1),COMPONENT.ProdCosts(:,2),'b-o');
set(gca,'tag','axes1')
xlabel('Quantity of FC Produced (kW)')
ylabel('Production Cost [$/kW]')
set(handles.editSellBack,'string',num2str(Project.Utilities.Grid.SellBackRate))
if ~isfield(Project.Utilities.Grid,'SellBackRatePerc');
    Project.Utilities.Grid.SellBackPerc =100;
end
set(handles.editSellbackPerc,'string',num2str(Project.Utilities.Grid.SellBackPerc));
if Project.Utilities.Grid.SellBackRate == 0
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonNoSellBack)
elseif Project.Utilities.Grid.SellBackRate == -1
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonReverseMeter)
else
    set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonFixedSellBack)
end


% --- Outputs from this function are returned to the command line.
function varargout = Economics_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.Name=get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonSaveOnly.
function pushbuttonSaveOnly_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project COMPONENT 
Project.Economic=COMPONENT;
clear global COMPONENT SYSINDEX
close(gcf)

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global COMPONENT SYSINDEX
close(gcf)

% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project COMPONENT Model_dir
[f,p]=uiputfile('*.mat','Save As Economics Component',fullfile(Model_dir,'System Library','Economic'));
if f==0;return;end
component=COMPONENT;
save([p f],'component')
Project.Economic=COMPONENT;
clear global COMPONENT
close(gcf)

function editInstall_Callback(hObject, eventdata, handles)
% hObject    handle to editInstall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.InstallCost=str2double(get(hObject,'string'));     
COMPONENT.Payback = calcPayback();
set(handles.editPayback,'string',COMPONENT.Payback)

% --- Executes during object creation, after setting all properties.
function editInstall_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editInstall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editCHP_OM_Callback(hObject, eventdata, handles)
% hObject    handle to editCHP_OM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.CHP_OM=str2double(get(hObject,'string'));
COMPONENT.Payback = calcPayback();
set(handles.editPayback,'string',COMPONENT.Payback)

% --- Executes during object creation, after setting all properties.
function editCHP_OM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCHP_OM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editElecChill_Callback(hObject, eventdata, handles)
% hObject    handle to editElecChill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.ElecChill = str2double(get(hObject,'String'));
set(handles.editElecChill,'string',COMPONENT.ElecChill)

% --- Executes during object creation, after setting all properties.
function editElecChill_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editElecChill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editChillOM_Callback(hObject, eventdata, handles)
% hObject    handle to editChillOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.ChillOM = str2double(get(hObject,'String'));
set(handles.editChillOM,'string',COMPONENT.ChillOM)

% --- Executes during object creation, after setting all properties.
function editChillOM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editChillOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editColdStorage_Callback(hObject, eventdata, handles)
% hObject    handle to editColdStorage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.ColdStore = str2double(get(hObject,'String'));
set(handles.editColdStorage,'string',COMPONENT.ColdStore)

% --- Executes during object creation, after setting all properties.
function editColdStorage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editColdStorage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editColdOM_Callback(hObject, eventdata, handles)
% hObject    handle to editColdOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.ColdOM = str2double(get(hObject,'String'));
set(handles.editColdOM,'string',COMPONENT.ColdOM)

% --- Executes during object creation, after setting all properties.
function editColdOM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editColdOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editBattery_Callback(hObject, eventdata, handles)
% hObject    handle to editBattery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.Battery = str2double(get(hObject,'String'));
set(handles.editBattery,'string',COMPONENT.Battery)

% --- Executes during object creation, after setting all properties.
function editBattery_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBattery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editBatteryOM_Callback(hObject, eventdata, handles)
% hObject    handle to editBatteryOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.BatteryOM = str2double(get(hObject,'String'));
set(handles.editBatteryOM,'string',COMPONENT.BatteryOM)

% --- Executes during object creation, after setting all properties.
function editBatteryOM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBatteryOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editAbsorbChill_Callback(hObject, eventdata, handles)
% hObject    handle to editAbsorbChill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.AbsorbChill = str2double(get(hObject,'String'));
set(handles.editAbsorbChill,'string',COMPONENT.AbsorbChill)

% --- Executes during object creation, after setting all properties.
function editAbsorbChill_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAbsorbChill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editAbChillOM_Callback(hObject, eventdata, handles)
% hObject    handle to editAbChillOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.AbChillOM = str2double(get(hObject,'String'));
set(handles.editAbChillOM,'string',COMPONENT.AbChillOM)

% --- Executes during object creation, after setting all properties.
function editAbChillOM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAbChillOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkboxAbsorpChill.
function checkboxAbsorpChill_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAbsorpChill (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.CurveFit = str2double(get(hObject,'String'));

function editSolar_Callback(hObject, eventdata, handles)
% hObject    handle to editSolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.Solar = str2double(get(hObject,'String'));
set(handles.editSolar,'string',COMPONENT.Solar)

% --- Executes during object creation, after setting all properties.
function editSolar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSolar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSolarOM_Callback(hObject, eventdata, handles)
% hObject    handle to editSolarOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.SolarOM = str2double(get(hObject,'String'));
set(handles.editSolarOM,'string',COMPONENT.SolarOM)

% --- Executes during object creation, after setting all properties.
function editSolarOM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSolarOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editWind_Callback(hObject, eventdata, handles)
% hObject    handle to editWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.Wind = str2double(get(hObject,'String'));
set(handles.editWind,'string',COMPONENT.Wind)

% --- Executes during object creation, after setting all properties.
function editWind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editWindOM_Callback(hObject, eventdata, handles)
% hObject    handle to editWindOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.WindOM = str2double(get(hObject,'String'));
set(handles.editWindOM,'string',COMPONENT.WindOM)

% --- Executes during object creation, after setting all properties.
function editWindOM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWindOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editFinance_Callback(hObject, eventdata, handles)
% hObject    handle to editFinance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.FinanceYrs=str2double(get(hObject,'string'));
COMPONENT.Payback = calcPayback();
set(handles.editPayback,'string',COMPONENT.Payback)

% --- Executes during object creation, after setting all properties.
function editFinance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFinance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editDepreciate_Callback(hObject, eventdata, handles)
% hObject    handle to editDepreciate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.DepreciateYrs=str2double(get(hObject,'string'));
COMPONENT.Payback = calcPayback();
set(handles.editPayback,'string',COMPONENT.Payback)

% --- Executes during object creation, after setting all properties.
function editDepreciate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDepreciate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editInflation_Callback(hObject, eventdata, handles)
% hObject    handle to editInflation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.Inflation=str2double(get(hObject,'string'));
COMPONENT.Payback = calcPayback();
set(handles.editPayback,'string',COMPONENT.Payback)

% --- Executes during object creation, after setting all properties.
function editInflation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editInflation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editInterest_Callback(hObject, eventdata, handles)
% hObject    handle to editInterest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.Interest=str2double(get(hObject,'string'));
COMPONENT.Payback = calcPayback();
set(handles.editPayback,'string',COMPONENT.Payback)

% --- Executes during object creation, after setting all properties.
function editInterest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editInterest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editPayback_Callback(hObject, eventdata, handles)
% hObject    handle to editPayback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT Project
COMPONENT.Payback=str2double(get(hObject,'string'));
OriginCost = Project.Economic.InstallCost;
tempPayback = calcPayback();
error = COMPONENT.Payback-tempPayback;
while abs(error)>.1
    if tempPayback ==20 && Project.Economic.InstallCost<50
        error = 0;
        COMPONENT.Payback = 20;
        set(handles.editPayback,'string',COMPONENT.Payback)
        Project.Economic.InstallCost = OriginCost;
    else
        COMPONENT.InstallCost = COMPONENT.InstallCost*(1+error/tempPayback);
        tempPayback = calcPayback();
        error = COMPONENT.Payback-tempPayback;
    end
end
set(handles.editInstall,'string',COMPONENT.InstallCost)

% --- Executes during object creation, after setting all properties.
function editPayback_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPayback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function payBack = calcPayback()
global Project COMPONENT
Project.Economic = COMPONENT;
Project.Result = runAnalyses(Project);
payBack = Project.Result.costOut.Payback;

function editLifeYrs_Callback(hObject, eventdata, handles)
% hObject    handle to editLifeYrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT 
COMPONENT.LifeYrs = str2double(get(hObject,'String'));
COMPONENT.LifekWh = COMPONENT.LifeYrs*COMPONENT.CHPSize*8760;
set(handles.editLifekWh,'string',COMPONENT.LifekWh)

% --- Executes during object creation, after setting all properties.
function editLifeYrs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLifeYrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editLifekWh_Callback(hObject, eventdata, handles)
% hObject    handle to editLifekWh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT 
COMPONENT.LifekWh = str2double(get(hObject,'String'));
COMPONENT.LifeYrs = COMPONENT.LifekWh/(COMPONENT.CHPSize*8760);
set(handles.editLifeYrs,'string',COMPONENT.LifeYrs)

% --- Executes during object creation, after setting all properties.
function editLifekWh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLifekWh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editReplace_Callback(hObject, eventdata, handles)
% hObject    handle to editReplace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.StackReplaceCost=str2double(get(hObject,'string'));      
COMPONENT.Payback = calcPayback();
set(handles.editPayback,'string',COMPONENT.Payback)

% --- Executes during object creation, after setting all properties.
function editReplace_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editReplace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editIncentive_Callback(hObject, eventdata, handles)
% hObject    handle to editIncentive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.Incentive=str2double(get(hObject,'string'));
COMPONENT.Payback = calcPayback();
set(handles.editPayback,'string',COMPONENT.Payback)

% --- Executes during object creation, after setting all properties.
function editIncentive_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editIncentive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes when entered data in editable cell(s) in uitableProduction.
function uitableProduction_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitableProduction (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
COMPONENT.ProdCosts = get(handles.uitableProduction,'Data');
axes(findobj(gcf,'tag','axes1'))
plot(COMPONENT.ProdCosts(:,1),COMPONENT.ProdCosts(:,2),'b-o');
set(gca,'tag','axes1')
xlabel('Quantity of FC Produced (kW)')
ylabel('Production Cost [$/kW]')


function editSellBack_Callback(hObject, eventdata, handles)
% hObject    handle to editSellBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Utilities.Grid.SellBackRate = str2double(get(hObject,'String'));
set(handles.uipanelSellBack,'SelectedObject',handles.radiobuttonFixedSellBack)

% --- Executes during object creation, after setting all properties.
function editSellBack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSellBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editSellbackPerc_Callback(hObject, eventdata, handles)
% hObject    handle to editSellbackPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
Project.Utilities.Grid.SellBackPerc =str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function editSellbackPerc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSellbackPerc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in uipanelChillType.
function uipanelSellBack_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelChillType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global Project 
switch get(eventdata.NewValue,'Tag')
    case 'radiobuttonNoSellBack'
        Project.Utilities.Grid.SellBackRate = 0;
        set(handles.editSellBack,'string','0.00')
    case 'radiobuttonFixedSellBack'
        if Project.Utilities.Grid.SellBackRate < 0;
            Project.Utilities.Grid.SellBackRate = 0;
        else
            Project.Utilities.Grid.SellBackRate = str2double(get(handles.editSellBack,'String'));
        end
        set(handles.editSellBack,'string',num2str(Project.Utilities.Grid.SellBackRate))
    case 'radiobuttonReverseMeter'
        Project.Utilities.Grid.SellBackRate = -1;
        set(handles.editSellBack,'string','-1')
end
