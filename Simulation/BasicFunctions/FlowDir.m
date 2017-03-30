function block = FlowDir(block,flows)
% Script which orients the nodes relative to the flow directions
% direction = 1: co-flow, direction = 2: counter-flow, direction = 3: cross-flow
nodes = block.nodes;
if isfield(block,'rows')
    rows = block.rows;
    columns = block.columns;
else
    columns = nodes;
    rows = 1;
end
for j = 1:1:columns
    block.Flow1Dir(:,j) = (j:columns:nodes)';
end
if flows>=2
    switch block.direction
        case 'coflow'
            for j = 1:1:columns
                block.Flow2Dir(:,j) = (j:columns:nodes)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
            end
        case 'counterflow'
            for j = 1:1:columns
                block.Flow2Dir(:,j) = (columns-j+1:columns:nodes)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
            end
        case 'crossflow'
            for j = 1:1:rows
                block.Flow2Dir(:,j) = (1+columns*(j-1):j*columns)'; % matrix where first column is all nodes that recieve fresh air, second column is all nodes that those feed into....
            end
    end
end
if flows>=3 %reformer/methanator
    for j = 1:1:columns
        block.Flow3Dir(:,j) = (columns-j+1:columns:nodes)';%assumed opposite of flow 1
    end
end
block.HTadjacent = zeros(nodes,4);
for i = 1:1:nodes
    block.HTadjacent(i,1) = i-1;%previous node
    block.HTadjacent(i,2) = i+1;%next node
    block.HTadjacent(i,3) = i-columns;%node to left
    block.HTadjacent(i,4) = i+columns;%node to right
end
block.HTadjacent(1:columns:end,1)=linspace(1,nodes-columns+1,rows);%first node in each row has nothing before it
block.HTadjacent(columns:columns:end,2)=linspace(columns,nodes,rows)';%last node in each row has nothing after it
block.HTadjacent(1:columns,3)=linspace(1,columns,columns);%first row has nothing to left
block.HTadjacent(end-columns+1:end,4)=linspace(nodes-columns+1,nodes,columns);%last row has nothing to right