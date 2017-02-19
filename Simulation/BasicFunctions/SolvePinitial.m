function SolvePinitial(list)
%% Solve for initial Pressures
global modelParam Outlet
r =0;
Pplus = {};
Pminus = {};
Pdifference = [];
OutPort = {};
for k = 1:1:length(list) %look at all blocks that have pressure inlets and outlets
    block = list{k};
    [m,n] = size(modelParam.(block).P_Difference);%some blocks like heat exchangers can have 2 Pin and 2 Pout, some have none
    for i = 1:1:m
        r = r+1;
        Port1 = modelParam.(block).P_Difference{i,1}; %The higher pressure term, usually Pin
        Port2 = modelParam.(block).P_Difference{i,2}; %The lower pressure term, usually Pout
        if strcmp(modelParam.(block).(Port1).type,'out')
            OutPort{r} = {'Port1',block,Port1};
        elseif strcmp(modelParam.(block).(Port2).type,'out')
            OutPort{r} = {'Port2',block,Port2};
        end
        Pdifference(r) = modelParam.(block).(Port1).IC - modelParam.(block).(Port2).IC; %the difference in pressure
        if isfield(modelParam.(block),'dMdP')
            dMdP(r,1:2) = modelParam.(block).dMdP(i,1:2); %if it is a known function of flow rate (compressor or turbine, scale by dM/dP
            mFlow(r) = modelParam.(block).mFlow(i); % the nominal flow rate being targeted
        else dMdP(r,1:2) = [0,0];
        end
        %find the model state associated with the pressure of port 1 & port 2
        Pplus{r} = findPstate(block,Port1);
        Pminus{r} = findPstate(block,Port2);
    end
end
A = zeros(r,r);
b = zeros(r,1);
stateP = zeros(r,1);
OutPort2 = {};
for i = 1:1:r
    if dMdP(i,1)~=0 %&& ~ischar(Pplus{i}) && ~ischar(Pminus{i}) %both Pin and Pout must be states to solve for mass flow
        a = dMdP(i,:); % equality constraint to get a mass flow
        b(i) = mFlow(i);
    else
        a = [1, 1];
        b(i) = Pdifference(i); %equality constraint to get a pressure drop
    end
    if ischar(Pplus{i})
        b(i) = b(i)-str2double(Pplus{i})*a(1);%move a constant pressure term to the right side
    else
        if nnz(stateP==Pplus{i})>0 %state was already in list, put value in correct column
            k = find(stateP==Pplus{i},1);
        else%first time using this state, put in next column
            k = find(stateP==0,1);
            stateP(k) = Pplus{i};%record so we know what column of A is associated with what state of model.IC
        end
        A(i,k) = a(1);
        if ~isempty(OutPort{i}) && strcmp(OutPort{i}{1},'Port1')
            OutPort2{k}(1:2) = OutPort{i}(2:3); %keep track of what state in Ax = b corresponds to a certain outlet initial condition
        end
    end
    if ischar(Pminus{i})
        b(i) = b(i) + str2double(Pminus{i})*a(2);
    else
        if nnz(stateP==Pminus{i})>0 %state was already in list, put value in correct column
            k = find(stateP==Pminus{i},1);
        else %first time using this state, put in next column
            k = find(stateP==0,1);
            stateP(k) = Pminus{i};%record so we know what column of A is associated with what state of model.IC
        end
        A(i,k) = -a(2);
        if ~isempty(OutPort{i}) && strcmp(OutPort{i}{1},'Port2')
            OutPort2{k}(1:2) = OutPort{i}(2:3); %keep track of what state in Ax = b corresponds to a certain outlet initial condition
        end
    end
end
lA = find(stateP==0,1)-1;
if ~isempty(lA)
    A = A(:,1:lA); %remove columns from a becasue some pressure states are constant
else lA = length(stateP);
end
Pguess = lsqnonneg(A,b); %solve nonnegative least squares to guess initial pressures.

%% update outlets 
for i = 1:1:lA
    modelParam.Scale(stateP(i)) = Pguess(i);% update model initial condition
    if ~isempty(OutPort2{i})
        block = OutPort2{i}{1};
        port = OutPort2{i}{2};
        Outlet.(block).(port) = Pguess(i); %update outlets based on calculated pressure
        modelParam.(block).(port).IC = Pguess(i);
    end
end
modelParam.Pstates = stateP(1:lA);
modelParam.Poutlets = OutPort2(1:lA);
for k = 1:1:length(list) %push pressure states back into sub-block scaling factors and port initial conditions
    block = list{k};
    modelParam.(block).Scale = modelParam.Scale(modelParam.(block).States);
end

function StateNum = findPstate(block,Port)
%find the model state associated with the pressure of port 1
global modelParam
if ~isempty(modelParam.(block).(Port).Pstate)
    StateNum = modelParam.(block).States(modelParam.(block).(Port).Pstate); %the state # in the IC vector corresponding to the pressure state
elseif isempty(modelParam.(block).(Port).connected)%unconected inlet pressure port (stays at initial condition)
    StateNum = num2str(modelParam.(block).(Port).IC);%put into a string so that is moved to the b side of Ax = b
else %go to level #2
    BlockPort = char(modelParam.(block).(Port).connected);
    ni = strfind(BlockPort,'.');
    if isempty(ni) %connected to a function or lookup (not a component)
        StateNum = num2str(feval(BlockPort,0));%finds value of look-up function at time = 0 & converts to string
    else
        connectedBlock = BlockPort(1:ni-1);
        connectedPort = BlockPort(ni+1:end);
        if ~isempty(modelParam.(connectedBlock).(connectedPort).Pstate)
            StateNum = modelParam.(connectedBlock).States(modelParam.(connectedBlock).(connectedPort).Pstate);%the state # in the IC vector corresponding to the pressure state
        elseif isempty(modelParam.(connectedBlock).(connectedPort).connected)%unconected inlet pressure port (stays at initial condition)
            StateNum = num2str(modelParam.(connectedBlock).(connectedPort).IC);%put into a string so that is moved to the b side of Ax = b
        else %go to level #3
            BlockPort = char(modelParam.(connectedBlock).(connectedPort).connected);
            ni = strfind(BlockPort,'.');
            if isempty(ni)
                StateNum = num2str(feval(BlockPort,0));%finds value of look-up function at time = 0 & converts to string
            else
                connectedBlock = BlockPort(1:ni-1);
                connectedPort = BlockPort(ni+1:end);
                if ~isempty(modelParam.(connectedBlock).(connectedPort).Pstate)
                    StateNum = modelParam.(connectedBlock).States(modelParam.(connectedBlock).(connectedPort).Pstate);%the state # in the IC vector corresponding to the pressure state
                elseif isempty(modelParam.(connectedBlock).(connectedPort).connected)%unconected inlet pressure port (stays at initial condition)
                    StateNum = num2str(modelParam.(connectedBlock).(connectedPort).IC);%put into a string so that is moved to the b side of Ax = b
                else%go to level #4
                    BlockPort = char(modelParam.(connectedBlock).(connectedPort).connected);
                    ni = strfind(BlockPort,'.');
                    if isempty(ni)
                        StateNum = num2str(feval(BlockPort,0));%finds value of look-up function at time = 0 & converts to string
                    else
                        connectedBlock = BlockPort(1:ni-1);
                        connectedPort = BlockPort(ni+1:end);
                        if ~isempty(modelParam.(connectedBlock).(connectedPort).Pstate)
                            StateNum = modelParam.(connectedBlock).States(modelParam.(connectedBlock).(connectedPort).Pstate);%the state # in the IC vector corresponding to the pressure state
                        else isempty(modelParam.(connectedBlock).(connectedPort).connected)%unconected inlet pressure port (stays at initial condition)
                            StateNum = num2str(modelParam.(connectedBlock).(connectedPort).IC);%put into a string so that is moved to the b side of Ax = b
                        end
                    end
                end
            end
        end
    end
end