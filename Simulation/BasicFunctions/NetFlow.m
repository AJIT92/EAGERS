function netflow = NetFlow(Flow)
speciesName = fieldnames(Flow);
netflow=0;
for i = 1:1:length(speciesName)
    if ~strcmp(speciesName{i},'T')
        netflow = netflow+ Flow.(speciesName{i});
    end
end