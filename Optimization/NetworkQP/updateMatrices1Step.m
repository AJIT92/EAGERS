function QP = updateMatrices1Step(QP,scaleCost,marginal,IC,dt,Organize)
global Plant UB 
nG = length(UB);
allStor = [Organize.stor];
allGen = [Organize.thisSeq, Organize.utility];

nodes = Plant.nodalVar.nodes; 
genNames = Plant.nodalVar.genNames;
genStates = Plant.nodalVar.genStates; 

%% update costs
H = diag(QP.H);%convert to vector
for i = 1:1:nodes
    gen = Plant.Network(i).gen;
    for j = 1:1:length(gen)
        s = strfind(gen{j},'.');
        genName = gen{j}(s+1:end);
        I = find(strcmp(genName,genNames),1,'first');
        states = Plant.Generator(I).OpMatB.states;
        statesIndex = genStates{i,j};
        k = Organize.States{I}; %states of this generator
        if ismember(I,allGen)
            H(k) = H(k)*scaleCost(I);
            QP.f(k) = QP.f(k)*scaleCost(I);
        elseif ismember(I,allStor)
            Output = fieldnames(Plant.Generator(I).OpMatB.output);
            if isempty(IC) && ismember(I,Organize.stor) % first initialization give arbitrarily high cost to storage (but not thermal storage if in conjunction with electric dispatch)
                QP.f(k) = marginal.(Output{1});
                H(k) = 1e8*marginal.(Output{1});
            elseif isempty(IC)%thermal storage 1st time
                QP.f(k) = marginal.(Output{1});
                H(k) = 0;
            else
                MaxCharge =min((UB(I)-IC(I))/dt,Plant.Generator(I).OpMatB.Ramp.b(1)); %minimum of peak Charge, and SOC/time (completely charging storage in next step)
                QP.f(k) = marginal.(Output{1});
                if MaxCharge == 0
                    H(k) = 0;
                else
                    H(k) = 2*QP.f(k)/MaxCharge;  %factor of 2 because its solving C = 0.5*x'*H*x + f'*x
                end
            end
        end
    end 
end 
QP.H = diag(H);%convert back to matrix
end