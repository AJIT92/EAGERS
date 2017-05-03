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
ramplimit = zeros(r,1);
for i = 1:1:nG
    if QP.Organize.Dispatchable(i) ==1 && noRamp(i) ==1 %remove ramping from this generator
        rows = QP.Organize.Ramping{i};
        ramplimit(rows(1):rows(end)) = rows(1):rows(end);
    end
end
ramplimit = nonzeros(ramplimit);
QP.b(ramplimit) = inf;
