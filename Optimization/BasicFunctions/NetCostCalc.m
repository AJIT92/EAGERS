function Cost = NetCostCalc(Var1,Timestamp,method)
global Plant
if strcmp(method,'Dispatch')
    Dispatch = Var1;
    Input = 0*Dispatch;
    for i = 1:1:length(Plant.Generator)
        chiller = 0;
        if ~isempty(Plant.Generator(i).Output)
            cap = Plant.Generator(i).Output.Capacity*Plant.Generator(i).Size;
        end
        eff = [];
        if strcmp(Plant.Generator(i).Type,'Electric Generator') || strcmp(Plant.Generator(i).Type,'CHP Generator')
            eff = Plant.Generator(i).Output.Electricity;
        elseif strcmp(Plant.Generator(i).Type,'Chiller') && ~isfield(Plant.Data.Demand,'E')%don't include cost if it shows up in generator demand
            eff = Plant.Generator(i).Output.Cooling;
            if ~Plant.optimoptions.sequential
                chiller = 1;
            end
        elseif strcmp(Plant.Generator(i).Type,'Heater')
            eff = Plant.Generator(i).Output.Heat;    
        end
        if ~isempty(eff) && ~chiller %dont add the cost of a chiller if you ran E and C simultaneously, or you will double count the chiller demand
            Input(:,i) = Dispatch(:,i)./interp1(cap,eff,Dispatch(:,i));
        elseif strcmp(Plant.Generator(i).Type,'Utility')
            Input(:,i) = Dispatch(:,i);
        end
    end
    %% ___ %%% alternate: uses the cost in the optimization (FitB)
    % dt = Timestamp(2:end)-Timestamp(1:end-1);
    % Cost = zeros(length(Timestamp)-1,1);
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
elseif strcmp(method,'Input')
    Input = Var1;
    Dispatch = 0*Input;
    Dispatch(Input>1) = 1;
end
run = nnz(Timestamp);
nG = length(Plant.Generator);
dt = (Timestamp(2:run) - Timestamp(1:run-1))*24;
scaleCost = updateGeneratorCost(Timestamp(1:run-1)).*(dt*ones(1,nG));%% All costs were assumed to be 1 when building matrices
stor =[];
startupcost = zeros(run-1,nG);
for i = 1:1:nG
    if isfield(Plant.Generator(i).VariableStruct, 'StartCost')
        nStartups = (Dispatch(2:run,i)>0).*(Dispatch(1:run-1,i)==0);
        startupcost(:,i) = nStartups*Plant.Generator(i).VariableStruct.StartCost;
    end
    if isfield(Plant.Generator(i).OpMatB,'Stor')
       stor(end+1) = i; 
    end
end
scaleCost(:,stor)=0; %remove cost of storage
Cost = sum(Input(2:run,1:nG).*scaleCost + startupcost,2);