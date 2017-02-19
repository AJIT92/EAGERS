function varargout = NatGasSetup(varargin)
% NATGASSETUP MATLAB code for NatGasSetup.fig
%      NATGASSETUP, by itself, creates a new NATGASSETUP or raises the existing
%      singleton*.
%
%      H = NATGASSETUP returns the handle to a new NATGASSETUP or the handle to
%      the existing singleton*.
%
%      NATGASSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NATGASSETUP.M with the given input arguments.
%
%      NATGASSETUP('Property','Value',...) creates a new NATGASSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NatGasSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NatGasSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NatGasSetup

% Last Modified by GUIDE v2.5 07-Apr-2014 11:51:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NatGasSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @NatGasSetup_OutputFcn, ...
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


% --- Executes just before NatGasSetup is made visible.
function NatGasSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NatGasSetup (see VARARGIN)

% Choose default command line output for NatGasSetup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NatGasSetup wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global Project COMPONENT

COMPONENT=Project.Utilities.NatGas;

set(handles.editName,'string',COMPONENT.Name)
MonthNames = cellstr(['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';]);
MonthNames(13:length(COMPONENT.fuelPriceYears)) = cellstr('-');
if ~isfield(COMPONENT,'monthPrice')
    COMPONENT.monthPrice(1:12,1) = COMPONENT.fuelPrice(1);
end
monthprice = num2cell([COMPONENT.monthPrice(1:12,1)', linspace(0,0,length(COMPONENT.fuelPriceYears)-12)]);
data =[num2cell(COMPONENT.fuelPriceYears), num2cell(COMPONENT.fuelPrice),MonthNames,monthprice'];
set(handles.uitablePrice,'Data',data);
listUser = {'Residential'; 'Commercial';'Industrial';'Electric Utility';};
set(handles.popupmenuUser,'string',listUser,'value',2)
stateName = {'Alabama';'Alaska';'Arizona';'Arkansas';'California';'Colorado';'Connecticut';'Delaware';'Florida';'Georgia';
             'Hawaii';'Idaho';'Illinois';'Indiana';'Iowa';'Kansas';'Kentucky';'Louisiana';'Maine';'Maryland';
             'Massachusetts';'Michigan';'Minnesota';'Mississippi';'Missouri';'Montana';'Nebraska';'Nevada';'NewHampshire';'NewJersey';
             'NewMexico';'NewYork';'NorthCarolina';'NorthDakota';'Ohio';'Oklahoma';'Oregon';'Pennsylvania';'RhodeIsland';'SouthCarolina';
             'SouthDakota';'Tennessee';'Texas';'Utah';'Vermont';'Virginia';'Washington';'WestVirginia';'Wisconsin';'Wyoming';};
val = find(strcmp(Project.State,stateName),1);
if isempty(val)
    val = 1;
end
set(handles.popupmenuState,'string',stateName,'value',val)

plotGasRate(handles);

% --- Outputs from this function are returned to the command line.
function varargout = NatGasSetup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1} = handles.output;


% --- Executes when entered data in editable cell(s) in uitablePrice.
function uitablePrice_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitablePrice (see GCBO)
% eventdata  structure with the following fields (see UITABLEPRICE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT
Data=get(handles.uitablePrice,'Data');
COMPONENT.fuelPriceYears=cell2mat(Data(:,1));
COMPONENT.fuelPrice=cell2mat(Data(:,2));
COMPONENT.monthPrice = cell2mat(Data(1:12,4));
axes(findobj(gcf,'tag','axes1'))
MonthNames = cellstr(['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';]);
MonthNames(13:length(COMPONENT.fuelPriceYears)) = cellstr('-');
if ~isfield(COMPONENT,'monthPrice')
    COMPONENT.monthPrice(1:12,1) = COMPONENT.fuelPrice(1);
end
monthprice = num2cell([COMPONENT.monthPrice(1:12,1)', linspace(0,0,length(COMPONENT.fuelPriceYears)-12)]);
data =[num2cell(COMPONENT.fuelPriceYears), num2cell(COMPONENT.fuelPrice),MonthNames,monthprice'];
set(handles.uitablePrice,'Data',data);
axes(findobj(gcf,'tag','axes1'));
plot(COMPONENT.fuelPriceYears,COMPONENT.fuelPrice,'b-o');
set(gca,'tag','axes1')
xlabel('Year')
ylabel('Price')
title('Natural Gas Price')


% --- Executes when selected cell(s) is changed in uitablePrice.
function uitablePrice_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitablePrice (see GCBO)
% eventdata  structure with the following fields (see UITABLEPRICE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonSaveOnlyToProject.
function pushbuttonSaveOnlyToProject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveOnlyToProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project COMPONENT 
Project.Utilities.NatGas=COMPONENT;
clear global COMPONENT
close(gcf)
uiresume

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global COMPONENT
close(gcf)


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

% --- Executes on button press in pushbuttonSaveAs.
function pushbuttonSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project COMPONENT Model_dir
[f,p]=uiputfile('*.mat','Save As Gas Utility Structure',fullfile(Model_dir,'System Library','NatGas'));
if f==0;return;end
component=COMPONENT;
save([p f],'component')
Project.Utilities.NatGas=COMPONENT;
clear global COMPONENT 
close(gcf)
uiresume


% --- Executes on selection change in popupmenuState.
function popupmenuState_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Project
list = get(handles.popupmenuState,'string');
val = get(handles.popupmenuState,'value');
Project.State = list(val);
plotGasRate(handles);

% --- Executes during object creation, after setting all properties.
function popupmenuState_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuUser.
function popupmenuUser_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuUser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotGasRate(handles);

% --- Executes during object creation, after setting all properties.
function popupmenuUser_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuUser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonSetProjection.
function pushbuttonSetProjection_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSetProjection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COMPONENT Project
str = get(handles.popupmenuUser,'string');
val = get(handles.popupmenuUser,'Value');
user = str{val};
StateNames = get(handles.popupmenuState,'string');
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};
I = find(strcmp(Project.State,StateNames),1);
state = char(stateAbrev(I));
Projection = ProjectUtilityCost('gas',user,21,state);
for i = 1:1:floor(length(Projection)/12)
    FuelPrice(i) = mean(Projection(12*(i-1)+1:12*i));
end
COMPONENT.fuelPrice = FuelPrice';
COMPONENT.fuelPriceYears = linspace(2014,2034,21)';
COMPONENT.monthPrice = Projection(1:12)';

MonthNames = cellstr(['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';]);
MonthNames(13:length(COMPONENT.fuelPriceYears)) = cellstr('-');
monthprice = num2cell([COMPONENT.monthPrice(1:12,1)', linspace(0,0,length(COMPONENT.fuelPriceYears)-12)]);
data =[num2cell(COMPONENT.fuelPriceYears), num2cell(COMPONENT.fuelPrice),MonthNames,monthprice'];
set(handles.uitablePrice,'Data',data);
plotGasRate(handles);

%% 
function plotGasRate(handles)
global Project Model_dir COMPONENT
load(fullfile(Model_dir,'System Library', 'NatGas','RateData','GasRate'))
str = get(handles.popupmenuUser,'string');
val = get(handles.popupmenuUser,'Value');
user = str{val};

StateNames = get(handles.popupmenuState,'string');
stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';};
I = find(strcmp(Project.State,StateNames),1);
state = char(stateAbrev(I));
Projection = ProjectUtilityCost('gas',user,10,state);
Date1 = GasRate.Date(:,1);
Date2 = GasRate.Date(1:nnz(GasRate.Date(:,2)),2);
Date3 = GasRate.Date(1:nnz(GasRate.Date(:,3)),3);
CityGate = GasRate.(state).CityGate;
Residential = GasRate.(state).Residential;
Commercial = GasRate.(state).Commercial;
Industrial = GasRate.(state).Industrial;
ForElectricGeneration = GasRate.(state).ForElectricGeneration;
axes(handles.axes1)
cla reset
set(gca,'tag','axes1')
% plot(COMPONENT.fuelPriceYears,COMPONENT.fuelPrice,'b-o');
plot(Date1,CityGate,'r')
hold on
plot(Date1,Residential,'g')
plot(Date1,Commercial,'b')
plot(Date2,Industrial,'m')
plot(Date3,ForElectricGeneration,'c')

if strcmp(user,'Residential')
    colPlot = 'g';
elseif strcmp(user,'Commercial')
    colPlot = 'b';
elseif strcmp(user,'Industrial')
    colPlot = 'm';
elseif strcmp(user,'Electric Utility')
    colPlot = 'c';
end
CurrentProjection1 = interp1(COMPONENT.fuelPriceYears,COMPONENT.fuelPrice,2014:2023);
for year = 1:10
    for month = 1:1:12
        j = 12*(year-1)+month;
        Date4(j) =  datenum(2014, j,01);
        CurrentProjection(j)=CurrentProjection1(year);%+COMPONENT.monthPrice(month)-mean(COMPONENT.monthPrice);
    end
end
plot(Date4,CurrentProjection,'k--','Linewidth',4)
plot(Date4,Projection,colPlot)
DatePlot = [Date1(1:24:end)' Date4(1:24:end)];
set(gca,'xtick',DatePlot)
datetick('x','yy','keepticks')
xlabel('Year')
ylabel('Natural Gas Cost [$/mmBTU]')
legend('CityGate','Residential' ,'Commercial','Industrial','Elecric Utility','Current Projection', 'New Projection', 'Location','NorthWest')
title(char(Project.State))
