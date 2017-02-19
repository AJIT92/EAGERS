%storage needs to be loaded last so that the buffers can be calculated
global Plant UB dischEff
nG = length(Plant.Generator);
BuffPerc = Plant.optimoptions.Buffer; % percentage for buffer on storage
for j = 1:1:nG
    if isfield(Plant.Generator(j).OpMatA,'Stor')
        chargeCapacity = 0;
        if isfield(Plant.Generator(j).OpMatA.output, 'E')
            include = {'Electric Generator','CHP Generator'};
            for i = 1:1:nG
                if ismember(Plant.Generator(i).Type,include)
                    chargeCapacity = chargeCapacity+Plant.Generator(i).Size;
                end
            end
        elseif isfield(Plant.Generator(j).OpMatA.output, 'C')
            include = {'Chiller'};
            for i = 1:1:nG
                if ismember(Plant.Generator(i).Type,include)
                    chargeCapacity = chargeCapacity+Plant.Generator(i).Size;
                end
            end
        elseif isfield(Plant.Generator(j).OpMatA.output, 'H')
            include1 = {'Heater'};
            include2 = {'CHP Generator'};
            for i = 1:1:nG
                if ismember(Plant.Generator(i).Type,include1)
                    chargeCapacity = chargeCapacity+Plant.Generator(i).Size;
                elseif ismember(Plant.Generator(i).Type,include2)
                    chargeCapacity = chargeCapacity+Plant.Generator(i).Size.*Plant.Generator(i).OpMatA.output.H;%OpMatA.output.H is the heat ratio
                end
            end
        end
        Buffer = min((BuffPerc/100)*UB(j), (chargeCapacity))*dischEff(j);
        Plant.Generator(j).OpMatA.link.bineq(end-1) = -Buffer; %lower buffer ineq :  -SOC - W <= -Buffer becomes W>= buffer -SOC
        Plant.Generator(j).OpMatA.link.bineq(end) = UB(j)*dischEff(j)-Buffer; %upper buffer ineq :  SOC - Z <= (UB-Buffer)
        Plant.Generator(j).OpMatA.Z.ub = Buffer;
        Plant.Generator(j).OpMatA.W.ub = Buffer;
        Plant.Generator(j).OpMatB = Plant.Generator(j).OpMatA;
    end
end