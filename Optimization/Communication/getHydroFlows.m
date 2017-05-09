function Flow = getHydroFlows(Date,line)
%Date is a vector of timstamps that you request the flow rate at. If we are
%at hour 5, and you need the previous 8 hours of flow data, 4 points will
%be from previous dispatchs, and 4 from the historical record. The last
%value in Date should be 1 timestep ago.
% line is the index of hydro line names
global Plant
nG = length(Plant.Generator);
NodeNames = cell(length(Plant.Network),1);
for i = 1:1:length(Plant.Network)
    NodeNames(i) = Plant.Network(i).name;
end
name = Plant.subNet.lineNames.Hydro{line};
r = strfind(name,'_');
equip = Plant.Network(strcmp(name(1:r(1)-1),NodeNames)).Equipment;%all equipment at the node specified in the first part of the line name
for k = 1:1:length(equip)
    I = find(strcmp(equip(k),Plant.Data.Hydro.Equipment));%column index of this dam in the stored matrices of Data.Hydro
    if ~isempty(I)
        if ~isfield(Plant,'Dispatch')
            tKnown = Date(end);
        else tKnown = Plant.Dispatch.Timestamp(1);
        end
        Time = [];
        Flow = [];
        if ~isfield(Plant,'Dispatch') || any(Date<Plant.Dispatch.Timestamp(1)) %historical data
            x1 = max(1,nnz(Plant.Data.Hydro.Timestamp<Date(1))); 
            n = nnz(Plant.Data.Hydro.Timestamp(x1:end)<tKnown);%take care of any initial conditions before data exists
            if n == 0
                Time = [0 Date(end)];  
                if strcmp(name(r(1)+1:r(2)-1),'Hydro')
                    Flow = [Plant.Data.Hydro.Outflow(1,I), Plant.Data.Hydro.Outflow(1,I)];  
                else Flow = [Plant.Data.Hydro.Spillflow(1,I), Plant.Data.Hydro.Spillflow(1,I)]; 
                end
            else
                Time = Plant.Data.Hydro.Timestamp(x1:x1+n); 
                if strcmp(name(r(1)+1:r(2)-1),'Hydro')
                    Flow = Plant.Data.Hydro.Outflow(x1:x1+n,I);
                else Flow = Plant.Data.Hydro.Spillflow(x1:x1+n,I);
                end
            end
        end
        if isfield(Plant,'Dispatch')
            j = 0;
            A = nnz(Plant.Dispatch.Timestamp>0);
            x1 = max(1,nnz(Plant.Dispatch.Timestamp(1:A)<tKnown));
            while (x1+j)<=A && Plant.Dispatch.Timestamp(x1+j-1)<Date(end)
                Time(end+1) = Plant.Dispatch.Timestamp(x1+j);%first time is 1 step before first date, last step is after last date
                Flow(end+1) = Plant.Dispatch.GeneratorState(x1+j,nG+nLcum+line);
                j = j+1;
            end
        end
%         x1 = max(1,nnz(Plant.Data.Hydro.Timestamp<Date(1))); 
%         if length(Date) == 1
%             r = (Date - Plant.Data.Hydro.Timestamp(x1))/(Plant.Data.Hydro.Timestamp(x1+1) - Plant.Data.Hydro.Timestamp(x1));
%             if strcmp(name(r(1)+1:r(2)-1),'Hydro')
%                 Flow = (1-r)*Plant.Data.Hydro.Outflow(x1,I) + r*Plant.Data.Hydro.Outflow(x1+1,I);
%             else
%                 Flow = (1-r)*Plant.Data.Hydro.Spillflow(x1,I) + r*Plant.Data.Hydro.Spillflow(x1+1,I);
%             end
%         else
%             x2 = nnz(Plant.Data.Hydro.Timestamp<Date(end))+1; 
%             n = nnz(Date<=Plant.Data.Hydro.Timestamp(x1));%take care of any initial conditions before data exists
% 
%             if strcmp(name(r(1)+1:r(2)-1),'Hydro')
%                 Flow(1:n) = Plant.Data.Hydro.Outflow(x1);
%                 Flow(n+1:end) = interp1(Plant.Data.Hydro.Timestamp(x1:x2),Plant.Data.Hydro.Outflow(x1:x2+1,I),Date);
%             else
%                 Flow(1:n) = Plant.Data.Hydro.Spillflow(x1);
%                 Flow(n+1:end) = interp1(Plant.Data.Hydro.Timestamp(x1:x2),Plant.Data.Hydro.Spillflow(x1:x2+1,I),Date);
%             end
%         end
        Flow = interp1(Time,Flow,Date);
    end
end