function build_subNet
%identify generators and lines and their position in the network
%organize into sub-networks
% Group nodes between which transmission losses don't occur
% only create lines where transmission losses do occur
global Plant 
nG = length(Plant.Generator);
nodes = length(Plant.Network);
genNames = cell(nG,1);
nodeNames = cell(nodes,1);
networkNames = fieldnames(Plant.Network);
networkNames = networkNames(~strcmp('name',networkNames));
networkNames = networkNames(~strcmp('Equipment',networkNames));
subNet.lineNames = {};
subNet.lineEff = [];   %added space for line effeciencies
subNet.lineLimit = [];
for i = 1:1:nG
    genNames(i,1) = {Plant.Generator(i).Name};
end
for i = 1:1:nodes
    nodeNames(i) = {Plant.Network(i).name};
end
for net = 1:1:length(networkNames)
    subNet.(networkNames{net}) = [];
    n = 0;
    for i = 1:1:nodes
        if ~isempty(Plant.Network(i).(networkNames{net}))
            %first check and see if this node is already part of a subNet node
            j = 0;
            I = false;
            while ~I && j<n
                j = j+1;
                I = any(strcmp(nodeNames(i),subNet.(networkNames{net})(j).nodes));
            end
            if I
                %node was already pulled into a previus sub-net node because transmission factor was = 1
            else
                %add a new subnet node, and lump in any nodes with a perfect transmission factor
                j = 0;
                load = [];
                nName = nodeNames(i); %name of current node
                subNet.(networkNames{net})(n+1).nodes = nName;
                subNet.(networkNames{net})(n+1).connections = {};
                if isfield(Plant.Network(i).(networkNames{net}),'Load') && ~isempty(Plant.Network(i).(networkNames{net}).Load)%%note if there is a demand at this node
                    load(end+1) = Plant.Network(i).(networkNames{net}).Load;
                end
                connect = Plant.Network(i).(networkNames{net}).connections;
                TransEff = Plant.Network(i).(networkNames{net}).Trans_Eff;
                TransLimit = Plant.Network(i).(networkNames{net}).Trans_Limit;
                while j<length(connect)
                    j = j+1;
                    if TransEff(j) ==1 
                        %perfect transmission, don't create a line, and pull in new connections from this adjacent node
                        I = find(strcmp(connect{j},nodeNames),1,'first');
                        c = 0;
                        while c<length(Plant.Network(I).(networkNames{net}).connections)
                            c = c+1;
                            %make sure the connection doesn't already exist, and that the connection is not to a node already in this subNet node
                            if nnz(strcmp(Plant.Network(I).(networkNames{net}).connections(c),{subNet.(networkNames{net})(n+1).connections;subNet.(networkNames{net})(n+1).nodes}))==0
                                connect(end+1) = Plant.Network(I).(networkNames{net}).connections(c);
                                TransEff(end+1) = Plant.Network(I).(networkNames{net}).Trans_Eff(c);
                                TransLimit(end+1) = Plant.Network(I).(networkNames{net}).Trans_Limit(c);
                            end
                        end
                        subNet.(networkNames{net})(n+1).nodes(end+1) = nodeNames(I,1); %add this node to the current subNet node
                        if isfield(Plant.Network(I).(networkNames{net}),'Load') && ~isempty(Plant.Network(I).(networkNames{net}).Load)%%note if there is a demand at this node
                            load(end+1) = Plant.Network(I).(networkNames{net}).Load;
                        end
                    else
                        %imperfect transmission, need a line, make sure one does not already exist
                        I = false;
                        k = 0;
                        while ~I && k<n
                            k = k+1;
                            I = any(strcmp(connect{j},subNet.(networkNames{net})(k).nodes));
                        end
                        if ~I %connected node has not yet been collected into any subNet node
                            I = true;
                            pconnected = connect{j};
                        else
                            pconnected = subNet.(networkNames{net})(k).nodes{1}; %name of subNet Node connected to (may not be same as connect if nodes were agregated)
                            reverseName = strcat(pconnected,'_',networkNames{net},'_',nName{1});
                            if ismember(reverseName,subNet.lineNames)
                                I = false; %don't create new line (it already exists)
                            end
                        end
                        if I
                            subNet.lineNames(end+1,1) = {strcat(nName{1},'_',networkNames{net},'_',pconnected)};
                            subNet.lineEff(end+1,1) = TransEff(j);
                            subNet.lineLimit(end+1,1) = TransLimit(j);
                        end
                        subNet.(networkNames{net})(n+1).connections(end+1) = {pconnected};
                    end
                end
                subNet.(networkNames{net})(n+1).Load = load;
                n = n+1;
            end
        end
    end
    %identify equipment at each subNet node
    for m = 1:1:n
        for k = 1:1:length(subNet.(networkNames{net})(m).nodes)
            i = find(strcmp(subNet.(networkNames{net})(m).nodes(k),nodeNames),1,'first');
            gen = [];
            equip = Plant.Network(i).Equipment;
            for j = 1:1:length(equip)
                s = strfind(equip{j},'.');
                I = find(strcmp(equip{j}(s+1:end),genNames),1,'first');
                if ~isempty(I)
                    if ~isempty(Plant.Generator(I).OpMatA.states) || strcmp(Plant.Generator(I).Source,'Renewable') %avoid things like gas utility with no states
                        gen(end+1) = I;
                    end
                else disp(strcat('error, generator is not in library',equip{j}))
                end
            end
        end
        subNet.(networkNames{net})(m).Equipment = gen;
    end
end
Plant.subNet = subNet;