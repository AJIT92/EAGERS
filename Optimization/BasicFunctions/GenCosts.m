function costTerms = GenCosts(Gen,LB,UB,out)
capacity = Gen.Output.Capacity*UB;
efficiency = Gen.Output.(out);
operationRange = find(capacity>=LB);
x = capacity(operationRange);
options = optimset('Algorithm','active-set','LargeScale','off','MaxIter',50,'Display','none');
c = x./efficiency(operationRange); %cost of generator in terms of input
Index = find(efficiency(operationRange)==max(efficiency(operationRange)),1,'last');

%% FIT A: piecewise convex curve with zero y-intercept (may be same as linear)
costTerms.P = x(Index);
costTerms.Convex(1) = c(Index)/x(Index); %linear segment
if Index < length(c)%if index is equal or greater than c, then there is no convex portion
    alpha = x(Index:end)-x(Index);%this is (x-P). In the paper alpha is gamma
    C=[alpha .5*alpha.^2];
    c_alpha = c(Index:end)-c(Index);
    costTerms.Convex(2:3) = lsqlin(C,c_alpha,[-1 0;0 -1],[-costTerms.Convex(1);0],[C(end,1), C(end,2)],c_alpha(end),[],[],[],options);%%Quadratic fit to cost (made for alpha)
else %if there is no convex portion (Index>= length(c))
    costTerms.Convex(2) = 0;% no quadratic term
    costTerms.Convex(3) = 0; %no quadratic term
end
costTerms.Convex(3) = max(0,costTerms.Convex(3));%ensure a positive H value

%% Fit B: piecewise convex with non-zero y-intercept, first find point I beyond which cost curve is convex
dc_dxi = (c(2:end)-c(1:end-1))./(x(2:end)-x(1:end-1));
convex = and((dc_dxi(2:end)>1.001*(dc_dxi(1:end-1))),(dc_dxi(1:end-1)>0));
k = length(convex);
P_i = find(convex==1,1,'first');
if ~isempty(P_i)%possibly a convex section
    MostlyConvex = 0*convex+1;
    for j = 1:1:k-4
        MostlyConvex(j) = (sum(convex(j:j+4))>3);
    end
    while P_i<k && MostlyConvex(P_i)==0
        P_i = P_i+1;
    end
end
if isempty(P_i) || P_i==k%never convex, make a linear fit
    costTerms.I = x(end);
    costTerms.Intercept = (c(end)-c(1))/(x(end)-x(1)); %beta extends to UB, fit still has a constant term
    if costTerms.Intercept<=0
        costTerms.Intercept = c(end)/x(end);
    end
    costTerms.Intercept(2:3) = 0;
else
    costTerms.I = x(P_i+1);
    alpha = x(P_i+1:end)-x(P_i+1);
    c_alpha = c(P_i+1:end)-c(P_i+1);%cost associated with alpha segment
%     c_alpha = zeros(length(alpha),1);
%     for j = 1:1:length(alpha)
%         c_alpha(j) = sum(c(P_i+1:P_i+j));
%     end
    C =[alpha .5*alpha.^2];
    maxSlopeBeta = c_alpha(end)/alpha(end);
    costTerms.Intercept = min(maxSlopeBeta,(c(P_i+1)-c(1))/(x(P_i+1)-x(1))); %slope of beta segment
    costTerms.Intercept(2:3) = lsqlin(C,c_alpha,[0 -1;-1 0;],[0;-costTerms.Intercept(1)],[C(end,1), C(end,2)],c_alpha(end),[],[],[],options);%%Quadratic fit to cost (made for alpha)
    costTerms.Intercept(3) = max(0,costTerms.Intercept(3));%ensure positive H value
end
costTerms.Intercept(4) = c(1)-x(1)*costTerms.Intercept(1); % constant term (y-intercept)

%% plot for double check
% figure(1)
% plot([0; x],[0;c],'g')
% hold on
% 
% QT = costTerms.Convex;
% beta = [linspace(0,costTerms.P) linspace(costTerms.P,costTerms.P)]';
% alpha = [linspace(0,0) linspace(0,x(end)-costTerms.P)]';
% convexFit =  QT(1)*beta + QT(2)*alpha + .5*QT(3)*alpha.^2;
% plot([0; beta+alpha],[0;convexFit],'r')
% 
% QT = costTerms.Intercept;
% beta = [linspace(0,costTerms.I) linspace(costTerms.I,costTerms.I)]';
% alpha = [linspace(0,0) linspace(0,x(end)-costTerms.I)]';
% betterFit =  QT(1)*beta + QT(2)*alpha + .5*QT(3)*alpha.^2 + QT(4);
% plot([0; beta+alpha],[0;betterFit],'b')