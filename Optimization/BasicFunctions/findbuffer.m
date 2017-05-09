%% calculate the upper and lower buffer thresholds for each storage system
global Plant
nG = length(Plant.Generator);
BuffPerc = Plant.optimoptions.Buffer; % percentage for buffer on storage
for j = 1:1:nG
    if isfield(Plant.Generator(j).OpMatA,'Stor')
        Out = fieldnames(Plant.Generator(j).OpMatA.output);
        if any(strcmp('W',Out))
            %hydro
            dischargeCapacity = (Plant.Generator(j).VariableStruct.MaxGenFlow + Plant.Generator(j).VariableStruct.MaxSpillFlow)/12.1; %flow rate in 1000 ft^3 converted to 1000 acre ft (1000 acre-ft = 12.1 x 1000 ft^3/s * 1 hr)
            dischargeCapacity = dischargeCapacity*Plant.optimoptions.Horizon/10; %amount the reservio can discharge in 10% of the dispatch horizon.
            Buffer = min((BuffPerc/100)*Plant.Generator(j).OpMatA.Stor.UsableSize,dischargeCapacity);
        else
            chargeCapacity = 0; 
            if any(strcmp('E',Out))
                include = {'Electric Generator','CHP Generator'};
            elseif isfield(Plant.Generator(j).OpMatA.output, 'C')
                include = {'Chiller'};
            elseif isfield(Plant.Generator(j).OpMatA.output, 'H')
                include = {'Heater';'CHP Generator';};
            end
            for i = 1:1:nG
                if ismember(Plant.Generator(i).Type,include)
                    if strcmp(Out(1),'H')
                        chargeCapacity = chargeCapacity+Plant.Generator(i).Size*Plant.Generator(i).OpMatA.output.(Out{1}); %second part is to multiply CHP generators by the heat ratio
                    else
                        chargeCapacity = chargeCapacity+Plant.Generator(i).Size;
                    end
                end
            end
            chargeCapacity = chargeCapacity*Plant.optimoptions.Horizon/10; %amount the plant can charge the storage in 10% of the dispatch horizon.
            Buffer = min((BuffPerc/100)*Plant.Generator(j).OpMatA.Stor.UsableSize,chargeCapacity);
        end
        Plant.Generator(j).OpMatA.link.bineq(end-1) = -Buffer; %lower buffer ineq :  -SOC - W <= -Buffer becomes W>= buffer -SOC
        Plant.Generator(j).OpMatA.link.bineq(end) = Plant.Generator(j).OpMatA.Stor.UsableSize-Buffer; %upper buffer ineq :  SOC - Z <= (UB-Buffer)
        Plant.Generator(j).OpMatA.Z.ub = Buffer;
        Plant.Generator(j).OpMatA.W.ub = Buffer;
        Plant.Generator(j).OpMatB = Plant.Generator(j).OpMatA;
    end
end