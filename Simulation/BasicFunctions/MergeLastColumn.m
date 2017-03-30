function Out = MergeLastColumn(Flow, Dir, Cells)
spec = fieldnames(Flow);
for i = 1:1:length(spec)
    if strcmp(spec{i},'T')
        Out.T  = mean(Flow.T(Dir(:,end))); %temperature 
    else
        Out.(spec{i}) = max(0,sum(Flow.(spec{i})(Dir(:,end)))*Cells);%avoid sending negative outflows to other blocks
    end
end