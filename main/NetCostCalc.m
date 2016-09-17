function Cost = NetCostCalc(Dispatch,Time,Input)
global Plant
nG = length(Plant.Generator);
dt = (Time - [0, Time(1:end-1)])';
scaleCost = updateGeneratorCost(Time).*(dt*ones(1,nG));%% All costs were assumed to be 1 when building matrices
stor =[];
for i = 1:1:nG
    if isfield(Plant.Generator(i).OpMatB,'Stor')
       stor(end+1) = i; 
    end
end
scaleCost(:,stor)=0; %remove cost of storage
Cost = sum(Input(2:end,:).*scaleCost,2);

%% ___ %%% alternate: uses the cost in the optimization (FitB)
% dt = Time(2:end)-Time(1:end-1);
% Cost = zeros(length(Time)-1,1);
% StorPower =[];
% for i = 1:1:length(Plant.Generator)
%     if nnz(Dispatch(:,i))>0
%         s = Plant.Generator(i).OpMatB.states;
%         if isfield(Plant.Generator(i).OpMatB,'Stor')
%             StorPower(:, end+1) = Dispatch(1:end-1,i) - Dispatch(2:end,i); %no cost
%         elseif isfield(Plant.Generator(i).OpMatB,'constCost') %all of these cost terms need to be scaled later on
%             I = Plant.Generator(i).OpMatB.(s{1}).ub;
%             b1 = Plant.Generator(i).OpMatB.(s{1}).f;
%             b2 = Plant.Generator(i).OpMatB.(s{2}).f;
%             b3 = Plant.Generator(i).OpMatB.(s{2}).H;
%             b4 = Plant.Generator(i).OpMatB.constCost;
%             Cost = Cost + (min(Dispatch(2:end,i),I)*b1 + max(Dispatch(2:end,i)-I,0)*b2 + max(Dispatch(2:end,i)-I,0).^2*b3+b4).*dt.*scaleCost(:,i);
%         elseif~isempty(s) %single state generators (linear cost term) & utilities
%             b1 = Plant.Generator(i).OpMatB.(s{1}).f;
%             Cost = Cost + Dispatch(2:end,i)*b1.*scaleCost(:,i).*dt;
%             %need to handle case of sell-back
%         elseif ~isempty(strcmp(Plant.Generator(i).Source,'Renewable'))
%             %renewable
%         end
%     end
% end
