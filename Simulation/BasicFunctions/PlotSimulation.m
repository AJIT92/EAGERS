function PlotSimulation(T,Y,plotFC,animateFC,PlotTurbos)
global TagInf  modelParam

time = TagInf.Time;
if isfield(modelParam,'Plot')
    n = length(modelParam.Plot);
    for i = 1:1:n
        tagName = modelParam.Plot{i};
        r = strfind(tagName,'.');
        block = tagName(1:r-1);
        tag = tagName(r+1:end);
        figure(i)
        if isfield(TagInf,block) && isfield(TagInf.(block),tag)
            plot(time,TagInf.(block).(tag))
            hold on
            ylabel(tagName);
        end
    end
end
if plotFC
    CompNames = fieldnames(modelParam.Components);
    for i = 1:1:length(CompNames)
        if strcmp(modelParam.Components.(CompNames{i}).type,'FuelCell') || strcmp(modelParam.Components.(CompNames{i}).type,'Electrolyzer')
            block = modelParam.(CompNames{i});
            CellMap(Y(end,:),block,n+1);
            n = n+1;
        end
    end 
end
if animateFC
    CompNames = fieldnames(modelParam.Components);
    for i = 1:1:length(CompNames)
        if strcmp(modelParam.Components.(CompNames{i}).type,'FuelCell') || strcmp(modelParam.Components.(CompNames{i}).type,'Electrolyzer')
            block = modelParam.(CompNames{i});
            Animate(T,Y,block,n+1)
            n = n+1;
        end
    end
end
if PlotTurbos
    CompNames = fieldnames(modelParam.Components);
    for i = 1:1:length(CompNames)
        if strcmp(modelParam.Components.(CompNames{i}).type,'Blower') 
            CTmap(modelParam.(CompNames{i}),[],n+1)
            n = n+1;
        end
        if strcmp(modelParam.Components.(CompNames{i}).type,'Compressor') 
            comp = modelParam.(CompNames{i});
            CTmap(comp,[],n+1)
            n = n+1;
        end
        if strcmp(modelParam.Components.(CompNames{i}).type,'Turbine') 
            turb = modelParam.(CompNames{i});
            CTmap([],turb,n+1)
            n = n+1;
        end
    end
end