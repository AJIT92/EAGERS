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

for i = 1:1:nG
    genNames(i,1) = {Plant.Generator(i).Name};
end
for i = 1:1:nodes
    nodeNames(i) = {Plant.Network(i).name};
end
for net = 1:1:length(networkNames)
    subNet.(networkNames{net}) = [];
    if strcmp(networkNames{net},'Hydro')
        subNet.lineNames.(networkNames{net}) = {};
        subNet.lineMinimum.(networkNames{net}) = []; 
    else
        subNet.lineNames.(networkNames{net}) = {};
        subNet.lineEff.(networkNames{net}) = [];   %added space for line effeciencies
        subNet.lineLimit.(networkNames{net}) = [];
    end
    n = 0;
    for i = 1:1:nodes
        if ~isempty(Plant.Network(i).(networkNames{net}))
            %first check and see if this node is already part of a subNet node
            %nodes with perfect transmission are agregated into the first node in the nameList that they have perfect 2-way connection with
            [I,aNodes,connect] = agregatedNode(nodeNames{i},networkNames{net});
            if I == i%add a new subnet node
                subNet.(networkNames{net})(n+1).nodes = aNodes;
                subNet.(networkNames{net})(n+1).connections = {};
                subNet.(networkNames{net})(n+1).Load = [];
                for j = 1:1:length(aNodes)
                    I = find(strcmp(aNodes{j},nodeNames),1,'first');
                    if isfield(Plant.Network(I).(networkNames{net}),'Load') && ~isempty(Plant.Network(I).(networkNames{net}).Load)%%note if there is a demand at this node
                        subNet.(networkNames{net})(n+1).Load(end+1) = Plant.Network(I).(networkNames{net}).Load;
                    end
                end
                for j=1:1:length(connect(:,1))
                    if ~any(strcmp(connect{j,2},aNodes))%imperfect transmission, need a line
                        [J, cNodes,~] = agregatedNode(connect{j,2},networkNames{net});
                        pconnected = nodeNames{J};%name of node that the connected node will be agregated into if it is perfectly connected to any others
                        subNet.(networkNames{net})(n+1).connections(end+1) = {pconnected};
                        if J>i %new line connection, otherwise this was handled previously in the reverse direction
                            if strcmp(networkNames{net},'Hydro')
                                %% Make sure it is upriver node to downriver node
                                subNet.lineNames.(networkNames{net})(end+1,1) = (strcat(nodeNames(i),'_',networkNames{net},'_',pconnected));
                                %% this needs to be corrected
                                subNet.lineMinimum.(networkNames{net})(end+1,1) = Plant.Network(i).Hydro.InstreamFlow;
                            else
                                [eff, limit,dir] = lineProp(subNet.(networkNames{net})(n+1).nodes,cNodes,networkNames{net});%find forward & reverse transmission efficiency & limit
                                if strcmp(dir,'none') %no transmission (zero efficiency)
                                    %do nothing
                                else
                                    if strcmp(dir,'reverse')
                                        subNet.lineNames(end+1,1) = (strcat(pconnected,'_',networkNames{net},'_',nodeNames(i)));
                                    else
                                        subNet.lineNames.(networkNames{net})(end+1,1) = (strcat(nodeNames(i),'_',networkNames{net},'_',pconnected));
                                    end
                                    if strcmp(dir,'dual')
                                        subNet.lineEff.(networkNames{net})(end+1,1:2) = eff;
                                        subNet.lineLimit.(networkNames{net})(end+1,1:2) = limit;
                                    else
                                        subNet.lineEff.(networkNames{net})(end+1,1) = eff;
                                        subNet.lineLimit.(networkNames{net})(end+1,1) = limit;
                                    end
                                end
                            end
                        end
                    end
                end
                n = n+1;
            end
        end
    end
    %identify equipment at each subNet node
    if strcmp(networkNames{net},'Electrical')
        out = 'E';
    elseif strcmp(networkNames{net},'DistrictHeat')
        out = 'H';
    elseif strcmp(networkNames{net},'DistrictCool')
        out = 'C';
    end
    for m = 1:1:n
        for k = 1:1:length(subNet.(networkNames{net})(m).nodes)
            i = find(strcmp(subNet.(networkNames{net})(m).nodes{k},nodeNames),1,'first');
            gen = [];
            equip = Plant.Network(i).Equipment;
            for j = 1:1:length(equip)
                s = strfind(equip{j},'.');
                I = find(strcmp(equip{j}(s+1:end),genNames),1,'first');
                if ~isempty(I)
                    if isfield(Plant.Generator(I).OpMatA.output,out) && (~isempty(Plant.Generator(I).OpMatA.states) || strcmp(Plant.Generator(I).Source,'Renewable')) %avoid things like gas utility with no states
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

function [TransEff,TransLimit,dir] = lineProp(node1,node2,net)
%find the transmission efficiency and limit between 2 connected nodes
%if one of the nodes is perfectly connected to another node there may be
%more than one pathway connecting them, so agregate the lines
global Plant
nodes = length(Plant.Network);
nodeNames = cell(nodes,1);
for i = 1:1:nodes
    nodeNames(i) = {Plant.Network(i).name};
end
TransEff = zeros(1,2);
TransLimit = zeros(1,2);
node1 = unique(node1);
node2 = unique(node2);
for j = 1:1:length(node1)
    I = find(strcmp(node1{j},nodeNames),1,'first');
    for k = 1:1:length(node2)
        J = find(strcmp(node2{k},nodeNames),1,'first');
        
        %forward direction efficieny & limit
        c = find(strcmp(node2{k},Plant.Network(I).(net).connections));
        if ~isempty(c)
            weight = Plant.Network(I).(net).Trans_Limit(c)/(TransLimit(1,1) + Plant.Network(I).(net).Trans_Limit(c));
            TransEff(1,1) = (1-weight)*TransEff(1,1) + weight*Plant.Network(I).(net).Trans_Eff(c);%weighted efficiency of 2 concurrent lines
            TransLimit(1,1) = TransLimit(1,1) + Plant.Network(I).(net).Trans_Limit(c);
        end
        %reverse direction efficiency and limit
        c = find(strcmp(node1{j},Plant.Network(J).(net).connections));
        if ~isempty(c)
            weight = Plant.Network(J).(net).Trans_Limit(c)/(TransLimit(1,2) + Plant.Network(J).(net).Trans_Limit(c));
            TransEff(1,2) = (1-weight)*TransEff(1,2) + weight*Plant.Network(J).(net).Trans_Eff(c);
            TransLimit(1,2) = TransLimit(1,2) + Plant.Network(J).(net).Trans_Limit(c);
        end
    end
end 
if TransEff(1,1)==0 && TransEff(1,2)>0
    dir = 'reverse';
    TransEff = TransEff(1,2);
    TransLimit = TransLimit(1,2);
elseif TransEff(1,1)>0 && TransEff(1,2)==0
    dir = 'forward';
    TransEff = TransEff(1,1);
    TransLimit = TransLimit(1,1);
elseif TransEff(1,1)==0 && TransEff(1,2)==0
    dir = 'none';
    TransEff = [];
    TransLimit = [];
else
    dir = 'dual';
end

function [I,aNodes,connect] = agregatedNode(node,net)
%Any connected nodes with perfect bi-directional transfer are agregated into the node earliest in the list of names
%This function finds which node in the list that is
global Plant
nodes = length(Plant.Network);
nodeNames = cell(nodes,1);
for i = 1:1:nodes
    nodeNames(i) = {Plant.Network(i).name};
end
I = find(strcmp(node,nodeNames),1,'first');
aNodes = {node};
connect = cell(length(Plant.Network(I).(net).connections),2);
connect(:,2) = Plant.Network(I).(net).connections;
connect(:,1) = {node};
[m,~] = size(connect);
k = 0;
while k<m
    k = k+1;
    [TransEff, ~,~] = lineProp(connect(k,1),connect(k,2),net);
    if ~isempty(TransEff) && length(TransEff) == 2 && min(TransEff)==1 && ~strcmp(net,'Hydro')%perfect bi-directional energy transfer, hydro lines are river segments, can't agregate
        J = find(strcmp(connect{k,2},nodeNames),1,'first');
        aNodes(end+1) = nodeNames(J);%add to list of agregated nodes
        %%add additional connections to check
        c = length(Plant.Network(J).(net).connections);
        connect(end+1:end+c,1) = connect(k,2);
        connect(end+1:end+c,2) = Plant.Network(J).(net).connections;
        I = min(J,I); %keep lowest number (index in list of node names
    end
    [m,~] = size(connect);
end