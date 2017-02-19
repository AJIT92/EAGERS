%% Electric Chillers (kW/ton)

figure(10)
load('2E_2C_1TES.mat')
i = 3;
x = Generators(i).Output.Capacity(2:end)*Generators(i).Size;

LB = .1*x(1);
UB = x(end);
COP =  Generators(i).Output.Cooling(2:end);
%ADD theta (startup cost)
SU = 0*Generators(i).VariableStruct.StartCost/LB;
%%%%%

c = x./COP*.06; %kWe cost converted to $ by $.06/kWe
[~, I] = max(COP);
P = x(I); % point of maximum efficiency


alpha = max(0,x-P);
beta = max(0,P-x);
options2 = optimoptions(@lsqlin,'Algorithm','active-set','MaxIter',50,'Display','none');
C = [x alpha beta];
LT  = lsqlin(C,c,[],[],[],[],[],[],[],options2);
linFit = LT(1)*x + LT(2)*alpha + LT(3)*beta; 
plot(x,linFit,'m')
hold on

c1 = LT(1)*LB + LT(2)*0 + LT(3)*(P-LB); 
x = [0;LB;x];
alpha = max(0,x-P);
beta = max(0,x-LB);
gamma = min(0,x-LB);
c2 = [0;c1+SU*LB;c+SU*LB];
C = [.5*alpha.^2 gamma.*beta x alpha beta];
QT = lsqlin(C,c2,[],[],[0 -1 -10 0 0; 0 0 1 0 0;],[0;c2(2)/x(2)],[],[],[],options2);
plot(x,[0; c1;c],'g')
% plot(x,c2,'b')
% x = linspace(0,x(end))';
% alpha = max(0,x-P);
% beta = max(0,x-LB);
% gamma = min(0,x-LB);
theta = min(x,LB);
newFit = .5*QT(1)*alpha.^2 + QT(2)*gamma.*beta + QT(3)*x + QT(4)*alpha + QT(5)*beta -SU*theta;
plot(x,newFit,'r')



H = [0 0 0 0 0; 0 QT(1) 0 0 0; 0 0 0 QT(2) 0; 0 0 QT(2) 0 0; 0 0 0 0 0;];
f = [QT(3) QT(4) QT(5) 0 -max(1e-8,SU)];
A = [1 -1 0 0 0; -1 0 0 1 0; 0 0 0 -1 1;];
b = [P -LB LB];
Aeq = [1 0 -1 -1 0; 1 0 0 0 0;]; 
x2 = [x(1)+1e-8; x(2:end-1); x(end)-1e-8 ];
beq = [0*x+LB x2];
lb = [0;0;0;-LB;0;];
ub = [UB;UB-P;UB-LB;0;LB;];

options = optimoptions(@quadprog,'Algorithm','active-set','MaxIter',50);%,'Display','none');
Cost = zeros(length(x),1);
alpha = zeros(length(x),1);
beta = zeros(length(x),1);
gamma = zeros(length(x),1);
theta = zeros(length(x),1);
for j = 1:1:length(x)
[Test, Cost(j)] = quadprog(H,f,A,b,Aeq,beq(j,:),lb,ub,[],options);
alpha(j) = Test(2);
beta(j) = Test(3);
gamma(j) = Test(4);
theta(j) = Test(5);
end
plot(x,Cost,'k--')
eig(H)
error = sum(abs(newFit-Cost));


figure(11)
linFit2 = LT(1)*x + LT(2)*alpha; 
c3 = c2-linFit2;
plot(x,[0;c1;c]-linFit2,'g')
hold on
% plot(x,c3,'m')

alpha = max(0,x-P);
beta = max(0,x-LB);
gamma = min(0,x-LB);
C = [.5*alpha.^2 gamma.*beta x alpha beta];
QT2 = lsqlin(C,c3,[],[],[0 -1 -10 0 0; 0 0 1 0 0;],[0;c2(2)/x(2)],[],[],[],options2);
% x = linspace(0,x(end))';
% alpha = max(0,x-P);
% beta = max(0,x-LB);
% gamma = min(0,x-LB);
theta = min(x,LB);
newFit = .5*QT2(1)*alpha.^2 + QT2(2)*gamma.*beta + QT2(3)*x + QT2(4)*alpha + QT2(5)*beta -SU*theta;
plot(x,newFit,'r')


H = [0 0 0 0 0; 0 QT2(1) 0 0 0; 0 0 0 QT2(2) 0; 0 0 QT2(2) 0 0; 0 0 0 0 0;];
f = [QT2(3) QT2(4) QT2(5) 0 -max(1e-8,SU)];
A = [1 -1 0 0 0; -1 0 0 1 0; 0 0 0 -1 1;];
b = [P -LB LB];
Aeq = [1 0 -1 -1 0; 1 0 0 0 0;]; 
x2 = [x(1)+1e-8; x(2:end-1); x(end)-1e-8 ];
beq = [0*x+LB x2];
lb = [0;0;0;-LB;0;];
ub = [UB;UB-P;UB-LB;0;LB;];
for j = 1:1:length(x)
[Test, Cost(j)] = quadprog(H,f,A,b,Aeq,beq(j,:),lb,ub,[],options);
alpha(j) = Test(2);
beta(j) = Test(3);
gamma(j) = Test(4);
theta(j) = Test(5);
end
plot(x,Cost,'k--')