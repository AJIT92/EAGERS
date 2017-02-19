function ExportResult
global Model_dir Plant Dispatch DataLog
if ~isempty(DataLog)
    [f,p]=uiputfile(fullfile(Model_dir,'results','LoggedData','Test1.mat'),'Save logged test data...');
    save([p,f],'DataLog')
end
K = menu('Save Dialog','Export Dispatch to Excel','Do not Export');
if K ==1
    [f,p]=uiputfile(fullfile(Model_dir,'results',strcat(Plant.Name,'.xls')),'Save Plant As...');
    if f==0; return; end
    filename = fullfile(p,f);
    MSG = msgbox('Exporting...');

    ColumnLabelSheet1 = {'Name';'Type';'Size';};
    xlswrite(filename,ColumnLabelSheet1,'System','A1:A3');
    Row1Sheet1 = {Plant.Name;};
    Row2Sheet1 = {'Fuel  -->'};
    Row3Sheet1 = {'Size (kW) -->'};
    Row4Sheet1 = {'Day'};

    D = datevec(Dispatch.Dispatch.Timestamp');
    Data = D(:,3);
    Row4Sheet1(end+1) = {'Hour'};
    Data(:,end+1) = D(:,4) + round(4*D(:,5)/60)/4;
    Row4Sheet1(end+1) = {'Temperature'};
    Data(:,end+1) = Dispatch.Dispatch.Temperature;
    if isfield(Dispatch.Dispatch.Demand,'E')
        Row4Sheet1(end+1) = {'Electric Demand (kW)'};
        Data(:,end+1) = Dispatch.Dispatch.Demand.E';
    end
    if isfield(Dispatch.Dispatch.Demand,'H')
        Row4Sheet1(end+1) = {'Heating Demand (kW)'};
        Data(:,end+1) = Dispatch.Dispatch.Demand.H';
    end
    if isfield(Dispatch.Dispatch.Demand,'C')
        Row4Sheet1(end+1) = {'Cooling Demand (kW)'};
        Data(:,end+1) = Dispatch.Dispatch.Demand.C';
    end
    s = size(Data);
    j = s(2);
    for i = 1:1:length(Dispatch.ActiveGenerators.Name)
        j = j+1;
        Row1Sheet1(j) = Dispatch.ActiveGenerators.Name(i);
        Row2Sheet1(j) = Dispatch.ActiveGenerators.Source(i);
        Row3Sheet1(j) = Dispatch.ActiveGenerators.Size(i);
        Data(:,j) = Dispatch.Dispatch.GeneratorState(:,i);
        if nnz(Dispatch.GenType(i)==[1 9 10])
            Row4Sheet1(j) = {'Energy Purchase (kW)'};
        elseif nnz(Dispatch.GenType(i)==-1)
            Row4Sheet1(j) = {'Sellback (kW)'};
        elseif nnz(Dispatch.GenType(i)==[6 7 8])
            Row4Sheet1(j) = {'Usable Energy Stored (kWh)'};
        elseif nnz(Dispatch.GenType(i)==[2 3 4])
            Row4Sheet1(j) = {'Ouput (kW)'};
            j = j+1;
            Row4Sheet1(j) = {'Input (kW)'};
            Data(:,j) = Dispatch.Dispatch.GeneratorInput(:,i);
        end
    end
    xlswrite(filename,Row1Sheet1,'System','B1');
    xlswrite(filename,Row2Sheet1,'System','B2');
    xlswrite(filename,Row3Sheet1,'System','B3');
    xlswrite(filename,Row4Sheet1,'System','B4');
    r = size(Data);
    Data2 = {};
    for i = 1:1:r(2)
        Data2(:,i) = cellstr(num2str(Data(:,i),3));
    end
    xlswrite(filename,Data2,'System','B5');
    close(MSG)
end