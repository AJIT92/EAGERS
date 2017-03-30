function [GenDisp, Cost, Feasible] = FindFeasible(QPmain,Locked)
%remove ramping constraints and add them back in until it fails to isolate
%what constraint is causing it to be infeasible
%%currently unsure of what the fix is.

noRamp = QPmain.Organize.Dispatchable;
%first check it is feasible without ramping
QP = removeRamping(QPmain,noRamp);
[GenDisp, Cost, Feasible] = DispatchQP(QP,Locked);%this is the dispatch with fit B

% %now try fixing problem
% for i = 1:1:nnz(noRamp)
%     QP = removeRamping(QPmain,noRamp);
%     [GenDisp, Cost, Feasible] = DispatchQP(QP,Locked);%this is the dispatch with fit B
% end

function QP = removeRamping(QP,noRamp)
nG = length(QP.Organize.Dispatchable);
[r,~] = size(QP.A);
rkeep = linspace(1,r,r)';
for i = 1:1:nG
    if QP.Organize.Dispatchable(i) ==1 && noRamp(i) ==1 %remove ramping from this generator
        rows = Organize.Ramping{i};
        rkeep(rows) = 0;
    end
end
rkeep = nonzeros(rkeep);
if r>0 && nnz(rkeep)>0
    QP.A = QP.A(rkeep,xkeep);
    QP.b = QP.b(rkeep);
else QP.A = [];
    QP.b = [];
end