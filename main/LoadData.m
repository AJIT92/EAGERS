function LoadData(varargin)
% LOADDATA MATLAB code for LoadData.fig
%      LOADDATA, by itself, creates a new LOADDATA or raises the existing
%      singleton*.
%
%      H = LOADDATA returns the handle to a new LOADDATA or the handle to
%      the existing singleton*.
%
%      LOADDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOADDATA.M with the given input arguments.
%
%      LOADDATA('Property','Value',...) creates a new LOADDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LoadData_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LoadData_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LoadData

% Last Modified by GUIDE v2.5 09-Apr-2015 11:37:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LoadData_OpeningFcn, ...
                   'gui_OutputFcn',  @LoadData_OutputFcn, ...
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


% --- Executes just before LoadData is made visible.
function LoadData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LoadData (see VARARGIN)

% Choose default command line output for LoadData
% Update handles structure
guidata(hObject, handles);
handles.output = hObject;
global Plant DataOriginal figHandle2
figHandle2 = handles.figure1;
DataOriginal = Plant.Data;
S = fieldnames(Plant.Data.Demand);
str = {'Date'};
str(end+1) = {'Temperature'};
for i = 1:1:length(S)
    str(end+1) = cellstr(strcat('Demand.',char(S(i))));
end
set(handles.lbExistingData,'string',str,'value',2);
set(handles.popupmenuColumn,'string',str(2:end),'value',1)
fits = {};
S = fieldnames(Plant.Data.HistProf);
if nnz(strcmp(S,'Temperature'))
    fits(end+1) = {'Temperature'};
end
if nnz(strcmp(S,'E'))
    fits(end+1) = {'Electricity'};
end
if nnz(strcmp(S,'H'))
    fits(end+1) = {'Heating'};
end
if nnz(strcmp(S,'C'))
    fits(end+1) = {'Cooling'};
end

set(handles.sliderZoom,'Min',1,'Max',4,'Value',1)
set(handles.sliderZoom,'SliderStep',[1/3,1/3])
set(handles.sliderDate,'Min',1,'Max',2,'Value',1)
lbExistingData_Callback(hObject, eventdata, handles)
sliderZoom_Callback(handles.sliderZoom, eventdata, handles)
readMATfiles(handles)
set(handles.lbXLSfiles, 'string', {'empty list'},'value',1)


% --- Outputs from this function are returned to the command line.
function varargout = LoadData_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in lbExistingData.
function lbExistingData_Callback(hObject, eventdata, handles)
% hObject    handle to lbExistingData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%build list for Setup Listbox
global plot2
list=get(handles.lbExistingData,'string');
val=get(handles.lbExistingData,'value');
compSel=char(list{val});
set(handles.popupmenuColumn,'string',list,'value',val)

[plot2.X,plot2.Y] = pullData(compSel,'all');
plot2.scale=1;
plot2chart(handles,2)


function [time, data] = pullData(sel,window)
global Plant
if ischar(window) && strcmp(window,'all')
    Xs = 1;
    XFs = length(Plant.Data.Timestamp);
else Xs = window(1);
    XFs = window(2);
end

time = Plant.Data.Timestamp(Xs:XFs);
if strcmp(sel,'Date')
    data = Plant.Data.Timestamp(Xs:XFs);
elseif strncmp(sel,'Demand',6)
    data = Plant.Data.Demand.(sel(end))(Xs:XFs);
elseif strcmp(sel,'Temperature')
    data = Plant.Data.Temperature(Xs:XFs);
end


function putData(sel,data,window)
global Plant
if ischar(window) && strcmp(window,'all')
    Xs = 1;
    XFs = length(Plant.Data.Timestamp);
else
    if isempty(Plant.Data.Timestamp)
        Xs = 1;
        points = 1+96*(window(2)-window(1));
        Plant.Data.Timestamp(1:points) = window(1)+1/96*(1:1:points);
        Plant.Data.Temperature(1:points) = ones(1,points);
    elseif window(1)>=Plant.Data.Timestamp(1) && window(2)<=Plant.Data.Timestamp(end)
        Xs = nnz(Plant.Data.Timestamp<=window(1));
    elseif window(1)<Plant.Data.Timestamp(1) %% shift all existing data later
        preAdd = 96*(Plant.Data.Timestamp(1)-window(1));
        Plant.Data.Timestamp(preAdd+1:end+preAdd) = Plant.Data.Timestamp(1:end);
        Plant.Data.Timestamp(1:preAdd) = linspace(window(1),Plant.Data.Timestamp(1)-1/96,preAdd);
        Plant.Data.Temperature(preAdd+1:end+preAdd) = Plant.Data.Temperature(1:end);
        Plant.Data.Temperature(1:preAdd) = mean(Plant.Data.Temperature(preAdd+1:end+preAdd));
        S = fieldnames(Plant.Data.Demand);
        for i = 1:1:length(S)
            Plant.Data.Demand.(char(S(i)))(preAdd+1:end+preAdd) = Plant.Data.Demand.(char(S(i)))(1:end);
            Plant.Data.Demand.(char(S(i)))(1:preAdd) = mean(Plant.Data.Demand.(char(S(i)))(preAdd+1:end+preAdd));
        end
        Xs = 1;
    elseif window(2)>Plant.Data.Timestamp(end) %% extend all data variables & start with average values
        postAdd = 96*(window(2) - Plant.Data.Timestamp(end));
        Plant.Data.Timestamp(end+1:end+postAdd) = Plant.Data.Timestamp(end)+1/96*(1:1:postAdd);
        Plant.Data.Temperature(end+1:end+postAdd) = mean(Plant.Data.Temperature)*ones(1,postAdd);
        S = fieldnames(Plant.Data.Demand);
        for i = 1:1:length(S)
            Plant.Data.Demand.(char(S(i)))(end+1:end+postAdd) = mean(Plant.Data.Demand.(char(S(i))))*ones(1,postAdd);
        end
        Xs = nnz(Plant.Data.Timestamp<=window(1));
    end
    XFs = Xs+round(96*(window(2)-window(1)));
end
if strcmp(sel(1:6),'Demand')
    Plant.Data.Demand.(sel(end))(Xs:XFs) = data;
elseif strcmp(sel,'Temperature')
    Plant.Data.Temperature(Xs:XFs) = data;
end


% --- Executes during object creation, after setting all properties.
function lbExistingData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbExistingData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function readMATfiles(handles)
% find names for .mat structures
global Model_dir
MATfiles=dir(fullfile(Model_dir,'data','*.mat'));
name = {};
for i = 1:1:length(MATfiles)
    A = open(fullfile(Model_dir,'data',MATfiles(i).name));
    A = A.(strrep(MATfiles(i).name,'.mat',''));
    if isstruct(A)
        B = fieldnames(A);
%             h=waitbar(0,'Loading .mat files');
        for j = 1:1:length(B)
%             waitbar(j/length(B),h,'Loading .mat files');
            if isstruct(A.(char(B(j))))
                C = fieldnames(A.(char(B(j))));
                for k = 1:1:length(C)
                    if ~(strcmp(char(C(k)),'Timestamp')||strcmp(char(C(k)),'Holidays'))
                        name(end+1) = cellstr(strcat(strrep(MATfiles(i).name,'.mat',''),'.',char(B(j)),'.',char(C(k))));
                    end
                end
            else
                if ~(strcmp(char(B(j)),'Timestamp')||strcmp(char(B(j)),'Holidays'))
                    name(end+1) = cellstr(strcat(strrep(MATfiles(i).name,'.mat',''),'.',char(B(j))));
                end
            end
        end
%         close(h)
    else name(end+1) = cellstr(strrep(MATfiles(i).name,'.mat',''));
    end
end
MATfiles = {};
for i = 1:1:length(name)
    if ~(strcmp(char(name(i)),'Timestamp')||strcmp(char(name(i)),'Holidays'))
        MATfiles(end+1) = name(i);
    end
end
set(handles.lbMATfiles, 'string', MATfiles,'value',1)
if ~isempty(name)
    loadFromFile(char(name(1)),'.mat',0,handles)
end

% --- Executes on button press in pushbuttonLoadXLS.
function pushbuttonLoadXLS_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadXLS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
readXLSfiles(handles)


function readXLSfiles(handles)
% find names for .xls sheets
global Model_dir
XLSfiles=dir(fullfile(Model_dir,'data','*.xls'));
name = {};
h=waitbar(0,'Loading .xls files');
for i = 1:1:length(XLSfiles)
    waitbar(i/length(XLSfiles),h,'Loading .xls files');
    [status,sheets] = xlsfinfo(XLSfiles(i).name);
    if length(sheets)==12 && strcmp(sheets(1),'Jan')
        name(end+1) = cellstr(strcat(strrep(XLSfiles(i).name,'.xls',''),'.','_12monthsheets'));
    else
        for j = 1:1:length(sheets)
            if ~isempty(status(j))
                name(end+1) = cellstr(strcat(strrep(XLSfiles(i).name,'.mat',''),'.',sheets(j)));
            end
        end 
    end
end
close(h)
XLSfiles = {};
for i = 1:1:length(name)
    if ~(strcmp(char(name(i)),'Timestamp')||strcmp(char(name(i)),'Holidays'))
        XLSfiles(end+1) = name(i);
    end
end
set(handles.lbXLSfiles, 'string', XLSfiles,'value',1)

%%% load data from selected file
function loadFromFile(name,type,plot,handles)
global plot1 plot2 Model_dir
%load 1st variable into plot 1 (type = .mat or .xls)
ind= strfind(name,'.');

if strcmp(type,'.xls')
    filename = strcat(name(1:ind(1)-1),'.xls');
    sheetname = name(ind(1)+1:end);
    if strcmp(sheetname,'_12monthsheets')
        months = cellstr({'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';});
        [A, labels] = xlsread(filename,char(months(1)));     
    else [A, labels] = xlsread(filename,sheetname);
    end
    if A(1,1)>=datenum(1900,1,1) && A(1,1)<=datenum(2100,1,1) && A(end,1)>=datenum(1900,1,1) && A(end,1)<=datenum(2100,1,1)
        B = datevec(A(1,1));
    else B= datevec(plot2.X(1));
    end
    if strcmp(sheetname,'_12monthsheets')
        plot1.X = linspace(datenum(B(1),1,1,.25,0,0),datenum(B(1)+1,1,1),96*(datenum(B(1)+1,1,1)-datenum(B(1),1,1)));
        for i = 2:1:12
            waitbar(i/12,handles.h,'Loading .xls sheets');
            [C] = xlsread(filename,char(months(i)));
            steps = length(C(:,1));
            A(end+1:end+steps,:) = C(1:steps,1:length(C(1,:)));
        end
    else plot1.X = linspace(datenum(B),datenum(B)+max(size(A))/96,max(size(A)));
    end
    if A(1,1)>=datenum(1900,1,1) && A(1,1)<=datenum(2100,1,1) && A(end,1)>=datenum(1900,1,1) && A(end,1)<=datenum(2100,1,1)
        labels = labels(1,2:end);
        A = A(:,2:end);
    end
    yearpoints15 = 96*(datenum(B(1)+1,1,1)-datenum(B(1),1,1));
    if length(A(:,1))< yearpoints15 %repeat hour or 2 hour values for 4 or 8 15 min data points
        n= yearpoints15/length(A(:,1));
        B = zeros(yearpoints15,length(A(1,:)));
        for i=1:1:n
            B(i:n:yearpoints15,:)=A;
        end
    elseif length(A(:,1))> yearpoints15 %average intermediate point to get 15 min values
        n= length(A(:,1))/yearpoints15;
        B = zeros(yearpoints15,length(A(1,:)));
        C = B;
        for i=1:1:n
            C=C+A(i:n:end,:);
        end
        B(1:yearpoints15,:) = C/n;
    else B = A;
    end
    plot1.Y = B';
    
elseif strcmp(type,'.mat')
    labels = {};
    if ~isempty(ind)
        A = open(fullfile(Model_dir,'data',strcat(name(1:ind(1)-1),type)));
        if length(ind)==1
            plot1.Y = A.(name(1:ind(1)-1)).(name(ind(1)+1:end));
        else plot1.Y = A.(name(1:ind(1)-1)).(name(ind(1)+1:ind(2)-1)).(name(ind(2)+1:end));
        end
    else plot1.Y = open(fullfile(Model_dir,'data',strcat(name,type)));
    end
    if ~isempty(ind) && isfield(A.(name(1:ind(1)-1)),'Timestamp')
        plot1.X = A.(name(1:ind(1)-1)).Timestamp;
    else plot1.X = linspace(plot2.X(1),plot2.X(1)+length(plot1.Y),length(plot1.Y)/96);
    end
end
plot1.scale = 1;
str = {'all'};
A = size(plot1.Y);
alphabet = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'X';'Y';'Z';};
for i = 1:1:min(A)
    if ~isempty(labels)
        str(end+1) = labels(1,i);
    else str(end+1) = cellstr(strcat('column',alphabet(i)));
    end
end
if plot==1
    if strcmp(type,'.xls')
        set(handles.popupmenuColumn,'string',str,'value',2)
    else set(handles.popupmenuColumn,'string',str,'value',1)
    end
    plot2chart(handles,1)    
end

%%% Plot to figure
function plot2chart(handles,new)
global plot1 plot2 plotC
if new ==1
    plotC = plot1;
    plotC.plot = 1;
elseif new ==2 
    plotC = plot2;
    plotC.plot = 2;
end
axes(handles.axes1)
Zoom = max(1,round(get(handles.sliderZoom,'Value')));
column = get(handles.popupmenuColumn,'value');
columnList = get(handles.popupmenuColumn,'string');
if ischar(columnList)
    columnList = cellstr(columnList);
end
columnName = columnList{column};
days = round(plotC.X(end)-plotC.X(1));
weeks = ceil(days/7);
if Zoom == 1
    Xs = max(1,nnz(plotC.X<=(ceil(plotC.X(1))+ round((get(handles.sliderDate,'Value')-1)*(days-1)))));
    XFs = nnz(plotC.X<=(plotC.X(Xs)+1));
elseif Zoom ==2
    Xs = max(1,nnz(plotC.X<=(ceil(plotC.X(1))+ 7*round((get(handles.sliderDate,'Value')-1)*(weeks-1)))));
    XFs = nnz(plotC.X<=(plotC.X(Xs)+7));
elseif Zoom ==3
    month = 1+ round((get(handles.sliderDate,'Value')-1)*11);
    D = datevec(plotC.X(1));
    if month ==1
        days = 0;
    else 
        days = datenum([D(1) D(2)+(month-1) 1 0 0 0])-plotC.X(1);%days until start of selected month
    end
    days2 = datenum([D(1) D(2)+(month) 1 0 0 0])-plotC.X(1);%days until end of selected month
    Xs = max(1,nnz(plotC.X<=(plotC.X(1)+days)));
    XFs = max(1,nnz(plotC.X<=(plotC.X(1)+days2)));
else
    Xs = 1;
    XFs = length(plotC.X);
end

if strcmp(columnName,'all')
    B = plotC.Y(:,Xs:XFs);
elseif strncmp(columnName,'column',6)
    B = plotC.Y(column-1,Xs:XFs);
else B = plotC.Y(Xs:XFs);
end 

% apply scaling factor for unit conversion
if isnumeric(plotC.scale)
    Y = B*plotC.scale;
elseif strcmp(plotC.scale,'CF')
    Y = 9/5*B+32;
elseif strcmp(plotC.scale,'FC')
    Y = (B-32)*5/9;
end

X = plotC.X(Xs:XFs);
hold off
plot(X,Y)
xlim([X(1) X(end)])

D = datevec(X(1));
days = round(X(end)-X(1));
dNum = X(1);
if days>366
    E = datevec(X(end));
    plotAxis = dNum+linspace(0,round(X(end)-X(1)),E(1)-D(1)+1);
elseif days == 365
    plotAxis = dNum+[0 31 59 90 120 151 181 212 243 273 304 334 365];
elseif days==366
    plotAxis = dNum+[0 31 60 91 121 152 182 213 244 274 305 335 366];
elseif days>=28 && days<=31
    plotAxis = dNum+linspace(0,days,days+1);
elseif days ==7
    plotAxis = dNum+linspace(0,7,8);
elseif days ==1
    plotAxis = dNum+linspace(0,1,25);  
end
set(gca,'xtick',plotAxis)
monthLabel = ['January  '; 'February '; 'March    '; 'April    '; 'May      '; 'June     '; 'July     '; 'August   '; 'September'; 'October  ';'November ' ;'December ';];
if days>366
    datetick('x','yyyy','keeplimits')
    xlabel('Year') 
elseif days==365|| days==366
    datetick('x','mmmdd','keeplimits')
    xlabel(num2str(D(1))) 
elseif days>=28 && days<=31
    datetick('x','dd','keeplimits')
    xlabel(strcat(monthLabel(D(2),:),'  ',num2str(D(1))))
elseif days==7
    datetick('x','dd','keeplimits')
    xlabel(strcat(monthLabel(D(2),:),'  ',num2str(D(1))))
elseif days ==1
    datetick('x','HH','keeplimits','keepticks')
    xlabel(strcat(['Hours of ', monthLabel(D(2),:), num2str(D(3))]))
end

% --- Executes on slider movement.
function sliderZoom_Callback(hObject, eventdata, handles)
% hObject    handle to sliderZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant
Zoom = max(1,round(get(handles.sliderZoom,'Value')));
a = datevec(Plant.Data.Timestamp(1));
b = datevec(Plant.Data.Timestamp(end));
months = max(1,12*(b(1)-a(1))+b(2)-a(2));
years = max(1,b(1)-a(1));
days = Plant.Data.Timestamp(end)-Plant.Data.Timestamp(1)+1;

if Zoom ==1
    set(handles.sliderDate,'SliderStep',[1/days,10/days])
elseif Zoom == 2
    set(handles.sliderDate,'SliderStep',[7/days,70/days])
elseif Zoom == 3
    set(handles.sliderDate,'SliderStep',[1/months,3/months])
else set(handles.sliderDate,'SliderStep',[1/years,1/years])
end
plot2chart(handles,0)


% --- Executes during object creation, after setting all properties.
function sliderZoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function sliderDate_Callback(hObject, eventdata, handles)
% hObject    handle to sliderDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot2chart(handles,0)

% --- Executes during object creation, after setting all properties.
function sliderDate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in lbMATfiles.
function lbMATfiles_Callback(hObject, eventdata, handles)
% hObject    handle to lbMATfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list=get(handles.lbMATfiles,'string');
val=get(handles.lbMATfiles,'value');
compSel=list{val};
name = char(compSel);
loadFromFile(name,'.mat',1,handles)

% --- Executes during object creation, after setting all properties.
function lbMATfiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbMATfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in lbXLSfiles.
function lbXLSfiles_Callback(hObject, eventdata, handles)
% hObject    handle to lbXLSfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list=get(handles.lbXLSfiles,'string');
val=get(handles.lbXLSfiles,'value');
compSel=list{val};
name = char(compSel);
handles.h=waitbar(0,'Loading .xls file: please wait');
loadFromFile(name,'.xls',1,handles)
close(handles.h)

% --- Executes during object creation, after setting all properties.
function lbXLSfiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lbXLSfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonDone.
function pushbuttonDone_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global figHandle2
close(figHandle2)

% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant
saveFile(Plant)

function saveFile(Plant) %don't save it as global
global Model_dir
[f,p]=uiputfile(fullfile(Model_dir,'Plant','PlantNew.mat'),'Save Plant As...');
if f==0; return; end
Plant.Name=strrep(f,'.mat','');
save([p,f],'Plant')


% --- Executes on button press in pushbuttonRevert.
function pushbuttonRevert_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRevert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant DataOriginal
Plant.Data = DataOriginal;
LoadData_OpeningFcn(hObject, eventdata, handles)


% --- Executes on button press in pushbuttonClear.
function pushbuttonClear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Plant
Plant.Data.Timestamp = [];
Plant.Data.Temperature = [];
Plant.Data.Demand = [];

% --- Executes on selection change in popupmenuColumn.
function popupmenuColumn_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuColumn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global plotC plot2
list=get(handles.lbExistingData,'string');
val=get(handles.lbExistingData,'value');
compSel=char(list{val});
list2=get(handles.popupmenuColumn,'string');
val2=get(handles.popupmenuColumn,'value');
if plotC.plot ==2
    [plot2.X,plot2.Y] = pullData(compSel,'all');
    plotC = plot2;
    plotC.scale = 1;
    plotC.plot =2;
end
plot2chart(handles,0)

% --- Executes during object creation, after setting all properties.
function popupmenuColumn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuColumn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonAdd.
function pushbuttonAdd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global plotC
val = get(handles.lbExistingData,'value');
VarList = get(handles.lbExistingData,'string');
sel = VarList{val};
ColNames = get(handles.popupmenuColumn,'string');
c = get(handles.popupmenuColumn,'value');
A = size(plotC.Y);
if c==1
    columns = A(1);
else columns = 1;
end
D = menu('Specify Duration to copy','All','Year','Month','Week','Day', 'Cancel');

A = datevec(plotC.X(1));
A(4) = A(4)+round(A(5)/15)/4;
A = A(1:4)';
if D ==1
    Start = plotC.X(1);
    EndDate = plotC.X(end);
else
    prompt = {'Enter year','Enter month:','Enter day:','Enter hour (in .25 increments):'};
    dlg_title = 'Specify the start of the period to copy';
    num_lines = 1;
    A = str2double(inputdlg(prompt, dlg_title, num_lines,num2str(A)));
    Start = datenum(A(1),A(2),A(3),A(4),0,0);
    if D == 2
        EndDate = datenum(A(1)+1,A(2),A(3),A(4),0,0)-1/96;
    elseif D == 3
        EndDate = datenum(A(1),A(2)+1,A(3),A(4),0,0)-1/96;
    elseif D == 4
        EndDate = datenum(A(1),A(2),A(3)+7,A(4),0,0)-1/96;
    elseif D == 5
        EndDate = datenum(A(1),A(2),A(3)+1,A(4),0,0)-1/96;
    end
    if (EndDate)>plotC.X(end)
        %error message: insufficient data (fill with average)
        plotC.X(end+1:end+96*((EndDate)-plotC.X(end)))=linspace(plotC.X(end+1/96),EndDate,((EndDate)-plotC.X(end))*96);
        plotC.Y(end+1:end+96*((EndDate)-plotC.X(end)))=mean(plotC.Y);
    end
end
Xs = nnz(plotC.X<=Start);
XFs = nnz(plotC.X<=EndDate);
data = plotC.Y(Xs:XFs);
% apply scaling factor for unit conversion
if isnumeric(plotC.scale)
    data = data*plotC.scale;
elseif strcmp(plotC.scale,'CF')
    data = 9/5*data+32;
elseif strcmp(plotC.scale,'FC')
    data = (data-32)*5/9;
end

list = get(handles.lbExistingData,'string');
answer = char(list(listdlg('PromptString','Select variable to copy to','SelectionMode','single','ListString',list)));
dlg_title = 'Specify when to copy to';
prompt = {'Enter year','Enter month:','Enter day:','Enter hour (in .25 increments):'};
sDate = inputdlg(prompt, dlg_title, 1,cellstr(num2str(A)));
year2 = str2double(sDate{1});
month2 = str2double(sDate{2});
day2 = str2double(sDate{3});
hour2 = str2double(sDate{4});
Start2 = datenum(year2,month2,day2,hour2,0,0);
EndDate2 = Start2+(EndDate-Start);
putData(answer,data,[Start2,EndDate2])



% --- Executes on button press in pushbuttonScale.
function pushbuttonScale_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global plotC
CF = menu('Conversion:','C --> F','F-->C', '*1000','*0.001','Other');
if CF == 1
    plotC.scale = 'CF';
elseif CF == 2
    plotC.scale = 'FC';
elseif CF == 3
    plotC.scale = plotC.scale*1000;
elseif CF ==4
    plotC.scale = plotC.scale*.001;
elseif CF ==5
    s = str2double(char(inputdlg('Scaling Factor', 'Specify the scaling factor in decimal form (i.e. .001 for 1/1000th)', 1,{'1000'})));
    plotC.scale = plotC.scale*s;
end
plot2chart(handles,0)

