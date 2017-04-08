function connectPorts
global modelParam Outlet Inlet
CompNames = fieldnames(modelParam.Components);
controls = fieldnames(modelParam.Controls);
list = [CompNames;controls;];
for k = 1:1:length(list)
    block = list{k};

    list2 = modelParam.(block).InletPorts;
    for i = 1:1:length(list2)
        port = list2{i};
        if isfield(modelParam.(block).(port),'IC')
                Inlet.(block).(port) = modelParam.(block).(port).IC; %use initial condition, update later if connected to an outlet
            else Inlet.(block).(port) = []; %don't have an IC for this inlet port, will get updated later
        end
        if length(modelParam.(block).connections)<i || isempty(modelParam.(block).connections{i})
            modelParam.(block).(port).connected={};
        else
            if ischar(modelParam.(block).connections{i})
                modelParam.(block).(port).connected = modelParam.(block).connections(i);
            else
                modelParam.(block).(port).IC = modelParam.(block).connections{i};
                modelParam.(block).(port).connected={};
            end
        end
    end
    list3 = modelParam.(block).OutletPorts;
    for i = 1:1:length(list3)
        port = list3{i};
        Outlet.(block).(port) = modelParam.(block).(port).IC;
    end   
end