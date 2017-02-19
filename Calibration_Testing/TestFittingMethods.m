options2 = optimoptions(@lsqlin,'Algorithm','active-set','MaxIter',50,'Display','none');
options3 = optimoptions(@linprog,'Algorithm','active-set','MaxIter',50,'Display','none');
options4 = optimoptions(@fmincon,'Algorithm','active-set','MaxIter',50,'Display','none');
figure(1)
load('3E.mat')
i = 1;
x = Generators(i).Output.Capacity(2:end)*Generators(i).Size;

LB = 0.1*x(1);
UB = x(end);
n =  Generators(i).Output.Electricity(2:end);
%ADD theta (startup cost)
SU = 2*Generators(i).VariableStruct.StartCost/LB;
%%%%%

c = x./n*4.15/293.15; %kWe cost with $4.15/1000 cu ft gas
costTerms.Convex(i,:) = polyfit(x,c,2); %%Quadratic fit to cost
x = [0;LB;x];
c1 = costTerms.Convex(i,1)*LB^2+costTerms.Convex(i,2)*LB+costTerms.Convex(i,3);
c2 = [0;c1+SU*LB;c+SU*LB];
beta = max(0,x-LB);
gamma = min(0,x-LB);
C = [.5*beta.^2  x beta ];
QT = lsqlin(C,c2,[],[],[0 1 0;],[c2(2)/x(2)],[],[],[],options2);
plot(x,[0;c1;c],'g')
hold on
% plot(x,c2,'m')
% x = [linspace(0,1,10) linspace(1,x(end))]';
beta = max(0,x-LB);
gamma = min(0,x-LB);
theta = min(LB,x);
newFit = .5*QT(1)*beta.^2 + -10*QT(2)*gamma.*beta + QT(2)*x + QT(3)*beta -SU*theta;
plot(x,newFit,'r')
QToriginal = QT;


H = [0 0 0 0; 0 QT(1) -10*QT(2) 0; 0 -10*QT(2) 0 0; 0 0 0 0;];
f = [QT(2) QT(3) 0 -max(1e-8,SU)];
A = [-1 0 1 0; 0 0 -1 1;];
b = [-LB LB];
Aeq = [1 -1 -1 0; 1 0 0 0;]; 
x2 = [x(1)+1e-8; x(2:end-1); x(end)-1e-8 ];
beq = [0*x+LB x2];
lb = [0;0;-LB;0;];
ub = [UB;UB-LB;0;LB;];

options = optimoptions(@quadprog,'Algorithm','active-set','MaxIter',50);%,'Display','none');
Cost = zeros(length(x),1);
beta = zeros(length(x),1);
gamma = zeros(length(x),1);
theta = zeros(length(x),1);
for j = 1:1:length(x)
[Test, Cost(j)] = quadprog(H,f,A,b,Aeq,beq(j,:),lb,ub,[],options);
beta(j) = Test(2);
gamma(j) = Test(3);
theta(j) = Test(4);
end
plot(x,Cost,'k--')
eig(H)
error = sum(abs(newFit-Cost));

% figure(15)
% H = [0 0 0 0; 0 QT(1) -10*QT(2) 0; 0 -10*QT(2) 0 0; 0 0 0 0;];
% f = [QT(2) QT(3) 0 -max(1e-8,SU)];
% A = [-1 0 1 0; 0 0 -1 1;-1 1 1 0; -1 0 0 0;];
% b = [0*x-LB 0*x+LB 0*x-LB -x2];
% % Aeq = [1 -1 -1 0; 1 0 0 0;]; 
% % beq = [0*x+LB x2];
% lb = [0;0;-LB;0;];
% ub = [UB;UB-LB;0;LB;];
% 
% % options = optimoptions(@quadprog,'Algorithm','trust-region-reflective','MaxIter',50);%,'Display','none');
% Cost = zeros(length(x),1);
% beta = zeros(length(x),1);
% gamma = zeros(length(x),1);
% theta = zeros(length(x),1);
% for j = 1:1:length(x)
% [Test, Cost(j)] = quadprog(H,f,A,b(j,:),[],[],lb,ub,[],options);
% beta(j) = Test(2);
% gamma(j) = Test(3);
% theta(j) = Test(4);
% end
% plot(x,Cost,'k--')
% eig(H)
% error = sum(abs(newFit-Cost));


% xsquare1 = zeros(length(x2));
% beta = max(0,x2-LB);
% gamma = min(0,x2-LB);
% theta = min(LB,x2);
% newFit = .5*QT(1)*beta.^2 + -10*QT(2)*gamma.*beta + QT(2)*x2 + QT(3)*beta -SU*theta;
% xsquare1(1,:) = newFit;
% for i = 2:1:length(x2)
%     xsquare1(i,:) = xsquare1(1,:)+1.1*newFit(i);
% end
% figure(27)
% surf(x2,x2,xsquare1);


figure(19) %reduced necessary constraints to 1.
QT = lsqlin(C,[0;c1;c],[],[],[0 1 0;],c1/x(2),[],[],[],options2);
H = [0 0 0; 0 QT(1) -10*QT(2); 0 -10*QT(2) 0;];
f = [QT(2) QT(3) -1e-5];
Aeq = [1 -1 -1; 1 0 0;]; 
beq = [0*x+LB x2];
lb = [0;0;-LB;];
ub = [UB;inf;0;];

Cost = zeros(length(x),1);
beta = zeros(length(x),1);
gamma = zeros(length(x),1);
for j = 1:1:length(x)
[Test, Cost(j)] = quadprog(H,f,[],[],Aeq,beq(j,:),lb,ub,[],options);
beta(j) = Test(2);
gamma(j) = Test(3);
end
plot(x,Cost,'k--')

% figure(20) %make convex?.
% beta = max(0,x-LB);
% gamma = min(0,x-LB);
% zeta = (-gamma)+beta;
% C = [.5*beta.^2  x beta ];
% QT = lsqlin(C,[0;c1;c],[],[],[0 1 0;],c1/x(2),[],[],[],options2);
% newFit = .5*QT(1)*beta.^2 + QT(2)*x + QT(3)*beta + -QT(3)*zeta -QT(3)*zeta;
% plot(x,[0;c1;c],'g')
% hold on
% plot(x,newFit,'r')
% H = [0 0 0 0; 0 QT(1) 0 0; 0 0 0 0;];
% f = [QT(2) QT(3) -QT(3)-1e-4 -QT(3)];
% Aeq = [1 -1 -1; 1 0 0;]; 
% beq = [0*x+LB x2];
% lb = [0;0;-LB;];
% ub = [UB;inf;0;];
% 
% Cost = zeros(length(x),1);
% beta = zeros(length(x),1);
% gamma = zeros(length(x),1);
% for j = 1:1:length(x)
% [Test, Cost(j)] = quadprog(H,f,[],[],Aeq,beq(j,:),lb,ub,[],options);
% beta(j) = Test(2);
% gamma(j) = Test(3);
% end
% plot(x,Cost,'k--')

% figure(3) %try with only 2 states (both quadratic)
% plot(x,[0;c1;c],'b')
% zeta = UB-x;
% C = [.5*zeta.^2 .5*beta.^2  zeta 0*x 0*x+1];
% QT = lsqlin(C,[0;c1;c],[],[],[],[],[],[],[],options2);
% x = [linspace(0,1,10) linspace(1,x(end))]';
% zeta = UB-x;
% beta = max(0,x-LB);
% newFit2 = .5*QT(1)*zeta.^2 + .5*QT(2)*beta.^2 + QT(3)*zeta + QT(4)*beta +QT(5);
% hold on
% plot(x,newFit2,'r')
% 
% H = [QT(1) 0; 0 QT(2);];
% f = [QT(3) QT(4)];
% A = [-1 -1; 1 0;]; 
% b = [0*x+LB-UB UB-x];
% lb = [-inf;0;];
% ub = [UB;inf;];
% Cost = zeros(length(x),1);
% beta = zeros(length(x),1);
% zeta = zeros(length(x),1);
% for j = 1:1:length(x)
% [Test, Cost(j)] = quadprog(H,f,A,b(j,:),[],[],lb,ub,[],options);
% beta(j) = Test(2);
% zeta(j) = Test(1);
% end
% plot(x,Cost+QT(5),'k--')

% %% try optimizing 2 generators with 3 states per gen
% figure(7)
% beta = max(0,x-LB);
% gamma = min(0,x-LB);
% C = [.5*beta.^2  x beta ];
% QT = lsqlin(C,[0;c1;c],[],[],[0 1 0;],c1/x(2),[],[],[],options2);
% QT2 = lsqlin(C,[0;1.25*c1;1.25*c],[],[],[0 1 0;],c1*1.25/x(2),[],[],[],options2);
% H = [0 0 0 0 0 0; 0 QT(1) -10*QT(2) 0 0 0; 0 -10*QT(2) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 QT2(1) -10*QT(2); 0 0 0 0 -10*QT(2) 0;];
% f = [QT(2) QT(3) -1e-5 QT2(2) QT2(3) -1e-5];
% Aeq = [1 -1 -1 0 0 0; 0 0 0 1 -1 -1; 1 0 0 1 0 0;]; 
% beq = [0*x+LB 0*x+LB 2*x2];
% lb = [0;0;-LB;0;0;-LB;];
% ub = [UB;inf;0;UB;inf;0;];
% x = [linspace(0,1,10) linspace(1,x(end))]';
% x2 = [x(1)+1e-8; x(2:end-1); x(end)-1e-8 ];
% beta = max(0,x-LB);
% gamma = min(0,x-LB);
% beq = [0*x+LB 0*x+LB 2*x2];
% newFit = .5*QT(1)*beta.^2 + -10*QT(2)*gamma.*beta + QT(2)*x + QT(3)*beta;
% C2gens = zeros(length(x));
% C2gens(1,:) = newFit;
% for i = 2:1:length(x)
%     C2gens(i,:) = C2gens(1,:)+1.25*newFit(i);
% end
% surf(x,x,C2gens);
% hold on
% 
% % QT2 = lsqlin(C,[0;1.25*c1;1.25*c],[],[],[],[],[],[],[],options2);%% try optimizing 2 generators with 2 states per gen
% % H = [QT(1) 0 0 0; 0 QT(2) 0 0; 0 0 QT2(1) 0; 0 0 0 QT2(2);];
% % f = [QT(3) QT(4) QT2(3) QT2(4)];
% % A = [-1 -1 0 0; 0 0 -1 -1; 1 0 1 0;]; 
% % b = [0*x+LB-UB 0*x+LB-UB 2*UB-2*x];
% % lb = [0;0;0;0;];
% % ub = [UB;UB-LB;UB;UB-LB;];
% Cost = zeros(length(x),1);
% Cost2 = zeros(length(x),1);
% xsoln1 =[];
% xsoln2 =[];
% for j = 1:1:length(x)
% % x0 = linprog(f,A,b(j,:),[],[],lb,ub,[],options3);
% % [Test, Cost(j)] = quadprog(H,f,A,b(j,:),[],[],lb,ub,x0,options);
% % xsoln1(j) = UB-Test(1);
% % xsoln2(j) = UB-Test(3);
% x0 = linprog(f.*[1 0 1 1 1 1],[],[],Aeq,beq(j,:),lb,ub,[],options3);
% xsoln1a(j) = x0(1);
% xsoln2a(j) = x0(4);
% [Test, Cost(j)] = quadprog(H,f,[],[],Aeq,beq(j,:),lb,ub,x0,options);
% xsoln1(j) = Test(1);
% xsoln2(j) = Test(4);
% x0 = linprog(f.*[1 1 1 1 0 1],[],[],Aeq,beq(j,:),lb,ub,[],options3);
% xsoln1b(j) = x0(1);
% xsoln2b(j) = x0(4);
% [Test, Cost2(j)] = quadprog(H,f,[],[],Aeq,beq(j,:),lb,ub,x0,options);
% xsoln3(j) = Test(1);
% xsoln4(j) = Test(4);
% % QP = @(x) .5*x'*H*x+f*x;
% % MsQP =createOptimProblem('fmincon','objective',QP,'Aeq',Aeq,'beq',beq(j,:),'lb',lb,'ub',ub,'x0',x0,'options',options4);
% % ms = MultiStart('StartPointsToRun','bounds-ineqs','UseParallel','never','Display','off');
% % [Test, Cost(j)] = run(ms,MsQP,10);
% if Cost2(j)<Cost(j)
%     xsoln1(j) = xsoln3(j);
%     xsoln2(j) = xsoln4(j);
% end
% end
% % plot3(xsoln1,xsoln2,Cost+QT(5)+QT2(5)+1,'k--','LineWidth',3)
% Cost3 = min(Cost,Cost2);
% plot3(xsoln1a,xsoln2a,Cost+1,'r--','LineWidth',3)
% plot3(xsoln1b,xsoln2b,Cost2+1,'b--','LineWidth',3)
% plot3(xsoln1,xsoln2,Cost3+1,'k--','LineWidth',3)
% 
% % figure(11) %view local minimum in problem
% % output = 230;
% % r = linspace(0,output);
% % cr = interp1([x(1); x(11:end)],[newFit(1); newFit(11:end)],r)+interp1([x(1); x(11:end)],1.25*[newFit(1); newFit(11:end)],output-r);
% % plot(r,cr)

%% try optimizing 3 generators at a demand of 400
% figure(12)
% beta = max(0,x-LB);
% gamma = min(0,x-LB);
% C = [.5*beta.^2  x beta ];
% QT = lsqlin(C,[0;c1;c],[],[],[0 1 0;],c1/x(2),[],[],[],options2);
% QT2 = lsqlin(C,[0;1.25*c1;1.25*c],[],[],[0 1 0;],c1*1.25/x(2),[],[],[],options2);
% QT3 = lsqlin(C,[0;1.4*c1;1.4*c],[],[],[0 1 0;],c1*1.4/x(2),[],[],[],options2);
% H = [0 0 0 0 0 0 0 0 0; 0 QT(1) -10*QT(2) 0 0 0 0 0 0; 0 -10*QT(2) 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 QT2(1) -10*QT(2) 0 0 0; 0 0 0 0 -10*QT(2) 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 QT3(1) -10*QT3(2); 0 0 0 0 0 0 0 -10*QT3(2) 0;];
% f = [QT(2) QT(3) -1e-5 QT2(2) QT2(3) -1e-5 QT3(2) QT3(3) -1e-5];
% Aeq = [1 -1 -1 0 0 0 0 0 0; 0 0 0 1 -1 -1 0 0 0; 0 0 0 0 0 0 1 -1 -1; 1 0 0 1 0 0 1 0 0;]; 
% beq = [LB LB LB 400];
% lb = [0;0;-LB;0;0;-LB;0;0;-LB;];
% ub = [UB;inf;0;UB;inf;0;UB;inf;0;];
% x = [linspace(0,1,10) linspace(1,x(end))]';
% x2 = [x(1)+1e-8; x(2:end-1); x(end)-1e-8 ];
% beta = max(0,x-LB);
% gamma = min(0,x-LB);
% newFit = .5*QT(1)*beta.^2 + -10*QT(2)*gamma.*beta + QT(2)*x + QT(3)*beta;
% C3gens = zeros(length(x));
% for i = 1:1:length(x)
%     for j = 1:1:length(x)
%         z = 400-x(i) - x(j);
%         if z>UB
%             C3gens(i,j) = NaN;
%         elseif z>LB
%             C3gens(i,j) = newFit(i)+1.25*newFit(j)+.5*QT3(1)*(z-LB).^2 + QT3(2)*z + QT3(3)*(z-LB);
%         else C3gens(i,j) =  newFit(i)+1.25*newFit(j)+QT3(2)*z;
%         end
%     end
% end
% surf(x,x,C3gens');
% hold on
% x0 = linprog(f,[],[],Aeq,beq,lb,ub,[],options3);
% x0 = [0 0 -LB 200 200-LB 0 200 200-LB 0]';
% [Test, Cost] = quadprog(H,f,[],[],Aeq,beq,lb,ub,x0,options);
% plot3(Test(1),Test(4),Cost+.1,'ko','MarkerSize',10)
% plot3(x0(1),x0(4),Cost+1,'go','MarkerSize',6)
% 
% QP = @(x) .5*x'*H*x+f*x;
% MsQP =createOptimProblem('fmincon','objective',QP,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'x0',x0,'options',options4);
% ms = MultiStart('StartPointsToRun','bounds-ineqs','UseParallel','never','Display','off');
% [Test, Cost] = run(ms,MsQP,3);

% 
% 
% figure(77)
% H = [0 0 0 0 0 0 0 0 0; 0 QT(1) 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 QT2(1) 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 QT3(1) 0; 0 0 0 0 0 0 0 0 0;];
% f = [QT(2) QT(3) -1e-5 QT2(2) QT2(3) -1e-5 QT3(2) QT3(3) -1e-5];
% Aeq = [1 -1 -1 0 0 0 0 0 0; 0 0 0 1 -1 -1 0 0 0; 0 0 0 0 0 0 1 -1 -1; 1 0 0 1 0 0 1 0 0;]; 
% beq = [LB LB LB 400];
% lb = [0;0;-LB;0;0;-LB;0;0;-LB;];
% ub = [UB;inf;0;UB;inf;0;UB;inf;0;];
% % QP = @(x) .5*x'*H*x+f*x;
% nleq = @(x) deal(0, -x(3:3:end)+x(2:3:end-1)-abs(x(1:3:end-2)-[LB LB LB]'));
% [Test, Cost] = fmincon(QP,x0,[],[],Aeq,beq,lb,ub,nleq,options4);
