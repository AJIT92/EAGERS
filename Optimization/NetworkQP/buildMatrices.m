function QPall = buildMatrices(fit,Time)
global Model_dir Plant Time
% Model_dir=strrep(which('Line_Losses_Method1.m'),fullfile('main','Line_Losses_Method1.m'),'');
% AddPaths()
% load(fullfile(Model_dir,'Plant','NetworkTest'))

% Time = [1,2,3]; %just a placeholder for vector of Time so I can control nS
dt = Time - [0, Time(1:end-1)]; %change in time per timestep
nS = length(Time); %total number of time steps
T_max = 200; %Maximum Transmission across each line


generators = length(Plant.Generator);
for i = 1:1:generators
    genNames(i,1) = {Plant.Generator(i).Name};
end

lineNames = {};
lineEff = [];   %added space for line effeciencies
gS = 0; %generator States
rC = 0; %ramping Constraints; 2 rC per ramping generator
ic = 0; %Number of Initial conditions
st = 0; %Number of storage states
nodes = length(Plant.Network);
    
for i = 1:1:nodes
    nodeNames(i,1) = {Plant.Network(i).name};
    connect  = Plant.Network(i).connections;
    for j = 1:1:length(connect) %identify all of the individual lines connecting the nodes
        if ~ismember(strcat(nodeNames{i},'_',connect{j}),lineNames) &&  ~ismember(strcat(connect{j},'_',nodeNames{i}),lineNames) 
            lineNames(end+1,1) = strcat(nodeNames(i),'_',connect{j});
            lineEff(end+1,1) = Plant.Network(i).line_eff(j);
        end
    end
    gen = Plant.Network(i).gen;
    for j = 1:1:length(gen)
        s = strfind(gen{j},'.');
        genType(i,j) = {gen{j}(1:s-1)};
        genName(i,j) = {gen{j}(s+1:end)};
        I = find(strcmp(genName{i,j},genNames),1,'first');
        genStates(i,j) = {gS+1:gS+length(Plant.Generator(I).OpMatB.states)}; %note when setting this up in buildmatrices, this may have to revert to OpMatA
        gS = gS + length(Plant.Generator(I).OpMatB.states);
        if isfield(Plant.Generator(I).OpMatB,'Ramp')
            rC = rC + 2; %add 2 ramping constraints
            ic = ic+1; %adding 1 initial condition for every ramp state; storage has a ramp rate;
        end
        if isfield(Plant.Generator(I).OpMatA,'Stor')
            st = st+1; %adding 1 for every storage state 
        end 
    end
end

%Creating global variables that are used in updateMatrices, updateMatrices1Step, & buildMatrices1Step
Plant.nodalVar.ic = ic;
Plant.nodalVar.nodes = nodes;
Plant.nodalVar.genNames = genNames;
Plant.nodalVar.genStates = genStates;

%%% This is only for Electricity; ie no heat, cooling, etc
Op = strcat('OpMatB');
Outs = Plant.optimoptions.Outputs;
Organize.States={};
Organize.Equalities = {};
Organize.Inequalities = {};
nG = generators;
Organize.IC = zeros(ic,1);
seq = 1;
while seq<= 1 %Currently just working with 'E'
    xL = 0;
    req = 0; % row index of the Aeq matrix and beq vector
    r = 0; % row index of the A matrix & b vector
    if strcmp('E',Outs{seq})>0
        Organize.E.Stor = zeros(1,nG);
        Organize.E.Demand = {'beq',(1:nS)};
        include = {'CHP Generator', 'Electric Generator'};
    end
    
    thisSeq = zeros(nG,1);
    stor = zeros(1,nG);
    utility = zeros(1,nG);
    Hratio = zeros(1,nG);
    renew = zeros(1,nG);
    
    QP.organize = cell(nS+1,nG); 
    ineqbuildcell = [];
    buildcell = [];
        
    for t = 1:1:nS
        for i = 1:1:nodes
        gen = Plant.Network(i).gen;
            for j = 1:1:length(gen)
                strf = strfind(gen{j},'.');
                genName = gen{j}(strf+1:end);
                I = find(strcmp(genName,genNames),1,'first');
                statesIndex = genStates{i,j};
                Gen = Plant.Generator(I).(Op);
                %% identify system type
                if ismember(Plant.Generator(I).Type,include)      
                    thisSeq(I) = I;
                    if strcmp(Plant.Generator(I).Type,'CHP Generator')
                        Hratio(I) = Plant.Generator(I).OpMatB.output.H;
                    end
                elseif strcmp(Plant.Generator(I).Type,'Utility')
                    if isfield(Plant.Generator(I).OpMatB.output,Outs{seq}) 
                        utility(I) =I;
                    end
                elseif isfield(Plant.Generator(I).OpMatB,'Stor')
                    if isfield(Plant.Generator(I).OpMatB.output,Outs{seq})
                        stor(I) = I; %storage in this seq
                    end
                elseif strcmp('E',Outs{seq}) && strcmp(Plant.Generator(I).Source,'Renewable')
                    renew(I) = I;
                end

                if thisSeq(I) || stor(I) || utility(I)
                    if t == 1
                        if isfield(Gen,'Ramp') 
                            Organize.IC(I) = req+1;
                            QP.organize{1,I} = xL; %output state organized into matrix of time vs. generator (IC)   
                            if stor(I)
                                QP.organize{1,I} = xL-1;
                            end 
                            neq = eval(Gen.req);
                            req = req + neq;
                        end 
                    end
                    s = length(Gen.states);%generator with multiple states
                    xLend = xL + s +ic;
                    if stor(I) 
                        s = 1; % storage only outputs the 1st state
                    end
                    QP.organize{t+1,I} = xL+ic+1:xL+(s-1)+ic+1; %output state organized into matrix of time vs. generator
                    if stor(I)
                        s = length(Gen.states);
                    end
                   if s >= 1
                        if t == 1;
                            ineqbuildcell(I,1) = req;
                            ineqbuildcell(I,end+1:end+s) = req+ic+1:req+(s-1)+ic+1;
                        else
                            ineqbuildcell(I,end+1:end+s) = xL+req+1:xL+req+(s-1)+1;
                        end  
                   end 
                    Organize.Equalities(I) = {[]}; %leaving blank; mutli-gens at each node means no associate per gen per line Aeq
                    ineq = rC*nS+1+(3*(t-1));
                    if ~utility(I)
                          if stor(I)
                              buildcell(I,end+1:end+2) =(r+1:r+2);
                              buildcell(I,end+1:end+3) = (ineq:ineq+2);
                          else
                              buildcell(I,end+1:end+2) =(r+1:r+2);
                          end 
                        r = r + 2;
                    end
                    if t == nS
                        Organize.Inequalities(1,I) = {sort((nonzeros(buildcell(I,:)))')};
                        Organize.States(1,I) = {(nonzeros(ineqbuildcell(I,:)))'};
                    end
                    xL = xLend-ic; %subtract ic to prevent it from being added to every generator loop at every time step
                else %if this is not a member of this sequence, but in case it is the last generator
                    Organize.States(I) = {[]};
                    Organize.Equalities(I) = {[]};
                    Organize.Inequalities(I) = {[]}; 
                end    
            end      
        end  
        xL = (gS + 3*length(lineNames))*t;
    end  
 
   
    thisSeq = nonzeros(thisSeq)';
    utility = nonzeros(utility)';
    stor = nonzeros(stor)';
    renew = nonzeros(renew)';
    CHPindex = nonzeros((1:nG).*(Hratio>0))';
    Organize.(Outs{seq}).thisSeq = thisSeq;
    Organize.(Outs{seq}).stor = stor;
    Organize.(Outs{seq}).utility = utility;
    Organize.(Outs{seq}).renew = renew;
    Organize.(Outs{seq}).CHPindex = CHPindex;
    Organize.(Outs{seq}).Hratio = nonzeros(Hratio)';
    allStor = [stor]; 
    allGen = [thisSeq];
    allUtility = [utility];
   

    xL = gS + 3*length(lineNames); %gen states + line states 
    req = length(nodeNames);
    r = 4*length(lineNames);

    QP.f = zeros(xL,1);
    QP.H = zeros(xL);      
    QP.Aeq = zeros(req,xL);
    QP.beq = zeros(req*nS,1);
    A = zeros(r,xL);
    b = zeros(r*nS,1);
    QP.lb = zeros(xL,1);
    QP.ub = Inf*ones(xL,1);

    S = zeros(3*st,xL);           
    b_stor = zeros(3*st,1);   

    Aeq_ic = zeros(ic,ic+(xL*nS));
    beq_ic = zeros(ic,1);
    ub_ic = Inf*ones(ic,1);


    z = 0;
    g = 1;
    %Sets up Aeq,lb,ub,f,S,b and H matrices
    for i = 1:1:nodes
        gen = Plant.Network(i).gen;
        for j = 1:1:length(gen)
            s = strfind(gen{j},'.');
            genName = gen{j}(s+1:end);
            I = find(strcmp(genName,genNames),1,'first');
            states = Plant.Generator(I).OpMatB.states;
            statesIndex = genStates{i,j};
            for k = 1:1:length(states)
                QP.H(statesIndex(k),statesIndex(k)) = Plant.Generator(I).OpMatB.(states{k}).H;
                QP.f(statesIndex(k),1) = Plant.Generator(I).OpMatB.(states{k}).f; 
                QP.ub(statesIndex(k),1) = Plant.Generator(I).OpMatB.(states{k}).ub;
                QP.lb(statesIndex(k),1) = Plant.Generator(I).OpMatB.(states{k}).lb;
                if ~isfield(Plant.Generator(I).OpMatB,'link')
                    QP.Aeq(i,statesIndex(k)) = Plant.Generator(I).OpMatB.output.E;
                else 
                    if k == 1
                        %No equalities
                        if isfield(Plant.Generator(I).OpMatB.link,'ineq')
                            [m,p] = size(Plant.Generator(I).OpMatB.link.ineq);
                            S(1+z:m+z,statesIndex(k):statesIndex(k)+p-1) = Plant.Generator(I).OpMatB.link.ineq(1:m,1:p);   
                            b_stor(1+z:m+z) = Plant.Generator(I).OpMatB.link.bineq(1:m);
                        end
                        QP.Aeq(i,statesIndex(k)) = 1/dt(t); %If member is 'X' then add negative SOC(t) state
                        if ismember('Y',states)
                            QP.Aeq(i,statesIndex(k+1)) = -1/dt(t); %charging state (negative load)
                        end
                        z = z+m;
                    end 
                end  
            end
            if isfield(Plant.Generator(I).OpMatB,'Ramp') %creating IC upper bounds to be the maximum size of generator
                Aeq_ic(g,g) = 1;
                beq_ic(g,1) = 0.5*Plant.Generator(I).Size;
                ub_ic(g,1) = Plant.Generator(I).Size;
                g = g+1;
            end 
        end 
        connect = Plant.Network(i).connections;
        for j=1:1:length(connect) %identify all of the individual lines connecting the nodes
            I = find(strcmp(strcat(nodeNames(i),'_',connect{j}),lineNames),1,'first');
            if ~isempty(I)
                QP.Aeq(i,gS+3*I-2) = -1; %negative value for state of power down line from current node to connected node
                QP.Aeq(i,gS+3*I-1) = -1; %negative value for state of power penalty down line from current node to connected node
            else I = find(strcmp(strcat(connect{j},'_',nodeNames(i)),lineNames),1,'first');
                QP.Aeq(i,gS+3*I-2) = 1; %positive value for state of power from connected node to current node
                QP.Aeq(i,gS+3*I) = -1; %negative value for state of power penalty from connected node to current node
            end
        end
    end

    %Sets up A matrices with ub and lb restrictions
    for i = 1:1:length(lineNames)
        I = find(strcmp(lineNames(i),lineNames),1,'first');
        A(4*i-3,gS+3*I-2) = -(1-lineEff(i)); % efficiency for state of power down line from connected node to current node
        A(4*i-3,gS+3*I-1) = -1; % efficiency for state of power penalty down line from  connected node to current node
        A(4*i-1,gS+3*I) = -1; % inequality that state of power penalty down line from current node to connected node 
        A(4*i-1,gS+3*I-2) = (1-lineEff(i)); % efficiency for state of power  from current node to connected node 
        A(4*i-2,gS+3*I-1) = -1; % efficiency for state of power penalty from connected node to current node
        A(4*i,gS+3*I) = -1; % inequality that state of power penalty from current node to connected node 
        QP.lb(gS+3*(I-1)+1:gS+3*I) = -T_max;
        QP.ub(gS+3*(I-1)+1:gS+3*I) = T_max;
    end



    QP.f = kron(ones(nS,1),QP.f);
    QP.H = kron(eye(nS),QP.H);
    QP.Aeq = kron(eye(nS),QP.Aeq);
    A = kron(eye(nS),A);
    QP.lb = kron(ones(nS,1),QP.lb);
    QP.ub = kron(ones(nS,1),QP.ub);    
    S = kron(eye(nS),S);
    b_stor = kron(ones(nS,1),b_stor);    

    R = zeros(rC*nS,xL*nS);
    b_ramp = zeros(rC,1);

    %Adding intial condition states to front of vectors/matrices
    QP.f = [zeros(ic,1);QP.f];
    [m,n] = size(QP.H);
    QP.H = [zeros(ic), zeros(ic,n);zeros(m,ic),QP.H];
    QP.ub = [ub_ic; QP.ub];
    QP.lb = [zeros(ic,1); QP.lb];
    QP.Aeq = [zeros(req*nS,ic), QP.Aeq];
    A = [zeros(r*nS,ic), A];
    S = [zeros(3*st*nS,ic),S];
    R = [zeros(rC*nS,ic),R];


    %ramping constraints                               
    h = 0;
    p = 1;
    for t = 1:1:nS
        if t == 1
            c = 1+ic;
        else
            c = ((xL*(t-2))+1+ic);
        end
        z = 0;
        for i = 1:1:nodes
            gen = Plant.Network(i).gen;
            for j = 1:1:length(gen)
                s = strfind(gen{j},'.');
                genName = gen{j}(s+1:end);
                I = find(strcmp(genName,genNames),1,'first');
                states = Plant.Generator(I).OpMatB.states;
                statesIndex = genStates{i,j};
                k = length(states);
                if statesIndex >= 1
                    if isfield(Plant.Generator(I).OpMatB,'Ramp')
                        if t ==1
                            ir = 1+z;
                            b_ramp(ir:ir+1,1) = Plant.Generator(I).OpMatB.Ramp.b(1:2,1);
                            o = [1;-1];                 %Ramping constraints for single state
                            R(ir:ir+1,p) = -1*o;
                            if ~isfield(Plant.Generator(I).OpMatB,'link') %Prevents four states being added for Storage
                                o = kron(ones(1,k),o);
                                R(ir:ir+1,c:c+(k-1)) = o;
                            else
                                R(ir:ir+1,c) = o;
                            end
                            c = c + k;
                            p = p + 1;
                            if length(gen) >= 1
                                z = z+2;
                            end
                        else
                            ir = 1+(rC*(t-1))+z;
                            o = [1;-1];                 %Ramping constraints for single state
                            if ~isfield(Plant.Generator(I).OpMatB,'link') %Prevents four states being added for Storage
                                o = kron(ones(1,k),o);       %Adds a matrix of [o] that is length(k) columns 
                                R(ir:ir+1,c:c+(k-1)) = -1*o; %Adds matrix values for previous timestep
                                c = c+xL;
                                R(ir:ir+1,c:c+(k-1)) = o;   %Adds matrix for current timestep
                                c = c+k-xL;
                            else
                                R(ir:ir+1,c) = -1*o; %Adds matrix values for previous timestep
                                c = c+xL;
                                R(ir:ir+1,c) = o;   %Adds matrix for current timestep
                                c = c+k-xL;
                            end 
                            if length(gen) >= 1 
                                z = z+2;
                            end
                        end
                    else          
                        c = c+k;
                    end
                end  
            end
        end 
    end 

    b_ramp = kron(ones(nS,1),b_ramp);

    %Storage initial Conditions
    k = 1;
    r = 0;
    z = 1;
    for i = 1:1:nodes
        gen = Plant.Network(i).gen;
        for j = 1:1:length(gen)
            s = strfind(gen{j},'.');
            genName = gen{j}(s+1:end);
            I = find(strcmp(genName,genNames),1,'first');
            statesIndex = genStates{i,j};
            if isfield(Plant.Generator(I).OpMatB,'Ramp')
                if isfield(Plant.Generator(I).OpMatB,'link')
                    [m,p] = size(Plant.Generator(I).OpMatB.link.ineq);
                    for t = 1:1:nS
                        if t == 1
                            S(r+1,z) = -1/dt(t); %This is IC for A matrix for storage
                            QP.Aeq(i,z) = -1/dt(t); %This is IC for AEQ matrix
                        else
                            S(r+3*(t-1)+1,statesIndex(k)+(xL*(t-2))+ic) = -1/dt(t); %State of Charge for Previous time step; SOC(t-1), for first charging penalty ineq.
                            QP.Aeq(req*(t-1)+I,statesIndex(k)+(xL*(t-2))+ic) = -1/dt(t);%SOC(t-1)
                        end  
                    end  
                    r = r+m;
                    k = k+p;
                else
                z = z+1;
                end
            end 
        end 
    end 

    QP.A = [ R; S; A];
    QP.b = [b_ramp; b_stor; b];
    QP.Aeq = [Aeq_ic; QP.Aeq];
    QP.beq = [beq_ic; QP.beq];
    
    QPall.QP.(Outs{seq})= QP;% quadratic programing matrices (H, f, Aeq, beq, A, b, ub, lb) -- parts may need to be updated later 
    seq = seq+1;
end 
QPall.Organize = Organize; %indices (rows and columns) associated with each generator, allowing generators to be removed later
QPall.Timestamp = Time(1);


 %end
%  
%     
%     %adding small cost to storage at last time step resolves
%     %issues with QuadProg
%     for i = 1:1:nodes
%         gen = Plant.Network(i).gen;
%         for j = 1:1:length(gen)
%             s = strfind(gen{j},'.');
%             genName = gen{j}(s+1:end);
%             I = find(strcmp(genName,genNames),1,'first');
%             states = Plant.Generator(I).OpMatB.states;
%             statesIndex = genStates{i,j};
%             if isfield(Plant.Generator(I).OpMatB,'link')
%                 k = 1;
%                 QP.f(((nS-1)*xL)+ic+statesIndex(k),1) = 5.1159; %Checkin to see if adding marginal cost to gen5 to see if it helps
%             end
%         end
%     end 
%     
% %This has been put into Update Matrices to Update beq for each time step%%%%
%  for t = 1:1:nS 
%     for i = 1:1:nodes
%          if ~isempty(Plant.Network(i).demand)
%              QP.beq((t-1)*req + i + ic) = Plant.Data.Demand(Plant.Network(i).demand).E(t);
% %             else
% %               LineLosses_updateMatrices %function to update IC, f, H, and beq matrices
%          end
%     end
%  end