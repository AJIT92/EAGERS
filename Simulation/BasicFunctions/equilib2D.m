function [Out,R] = equilib2D(InFlow,T,P,H2consume,Type,PercEquilib,guess)
% performs a gradient search to find the minimum gibbs function value for
% two simultaneous reactions (methane reforming and water gas shift)
% tolerance is a relative tolence, so that it is more accurat if only a
% small % change in the CH4 and CO is occuring
%considers only 2 reactions:
% CH4+H2O <--> CO + 3H2
% CO + H2O <--> CO2 + H2
n = length(T);
specInterest = {'CH4','CO','CO2','H2','H2O'};
h = enthalpy(T,specInterest);
s = entropy(T,specInterest);

spec = fieldnames(InFlow);
spec = spec(~strcmp('T',spec));
Out = InFlow; %handles all inert gases (not in specinterest list)
for k = 1:1:n
    XnP=0;
    for i = 1:1:length(spec)
        if ~ismember(spec{i},specInterest)
            XnP = XnP + InFlow.(spec{i})(k);
        else
            g0.(specInterest{i}) = h.(specInterest{i})(k)-T(k).*s.(specInterest{i})(k);
            Inlet.(specInterest{i}) = InFlow.(specInterest{i})(k);
        end
    end
    CH4max = min(Inlet.CH4,Inlet.H2O);
    CH4min = -min(Inlet.CO,Inlet.H2);
    Span1 = CH4max-CH4min;
    if ~isempty(guess)
        x0 = guess(k);
        R.CH4(k,1) = (Span1*x0+CH4min);
        
        X.CH4 = InFlow.CH4(k) - R.CH4(k);
        X.CO = InFlow.CO(k) + R.CH4(k);
        X.CO2 = InFlow.CO2(k);
        X.H2 = InFlow.H2(k) + 3*R.CH4(k) - H2consume(k);%hydrogen consumed
        X.H2O = InFlow.H2O(k) - R.CH4(k) + H2consume(k);% water produced
        if strcmp(Type,'MCFC')
            X.CO2 = X.CO2 + H2consume(k); % CO2 brought over
        end
        for i = 1:1:length(spec)
            if ~ismember(spec{i},specInterest)
                X.(spec{i}) = InFlow.(spec{i})(k);
            end
        end
        y0 = Newton1D(.5,X,-min(X.CO2,X.H2),min(X.H2O,X.CO(k)),T(k),P,specInterest,1e-6,'GibbVal');
        Tol = max(1e-10,x0*1e-7);
    else
        Tol = 1e-6;
        x0 = 0.85;
        y0 = 0.5;
    end

    [x,y] = Newton2D(Inlet,CH4min,CH4max,H2consume(k),Type,T(k),P,g0,XnP,x0,y0,Tol);
    R.CH4(k,1) = max(CH4min,PercEquilib*(Span1*x+CH4min)); %prevents a low % equilibrium from causing problems with enough H2 for electrochemistry
    
    X.CH4 = InFlow.CH4(k) - R.CH4(k);
    X.CO = InFlow.CO(k) + R.CH4(k);
    X.CO2 = InFlow.CO2(k);
    X.H2 = InFlow.H2(k) + 3*R.CH4(k) - H2consume(k);%hydrogen consumed
    X.H2O = InFlow.H2O(k) - R.CH4(k) + H2consume(k);% water produced
    if strcmp(Type,'MCFC')
        X.CO2 = X.CO2 + H2consume(k); % CO2 brought over
    end
    for i = 1:1:length(spec)
        if ~ismember(spec{i},specInterest)
            X.(spec{i}) = InFlow.(spec{i})(k);
        end
    end
    R_COmin = -min(X.CO2,X.H2);
    R_COmax = min(X.H2O,X.CO(k));
    if PercEquilib<1
        y = Newton1D(y0,X,R_COmin,R_COmax ,T(k),P,specInterest,1e-6,'GibbVal');
    end
    R.WGS(k,1) = y*R_COmax + (1-y)*R_COmin;
    OutFlow = FlowsOut(Inlet,R.CH4(k,1),R.WGS(k,1),H2consume(k),Type);
    
    for i = 1:1:length(specInterest)
        Out.(specInterest{i})(k) = OutFlow.(specInterest{i});
    end
end
Out.T = T;

function X = FlowsOut(Inlet,CH4,WGS,H2consume,Type)
X.CH4 = Inlet.CH4 - CH4;
X.CO = Inlet.CO + CH4 - WGS;
X.CO2 = Inlet.CO2 + WGS;
X.H2 = Inlet.H2 + 3* CH4 + WGS - H2consume;%hydrogen consumed
X.H2O = Inlet.H2O - CH4 - WGS + H2consume;% water produced
if strcmp(Type,'MCFC')
    X.CO2 = X.CO2 + H2consume; % CO2 brought over
end
spec = fieldnames(Inlet);
specInterest = {'CH4','CO','CO2','H2','H2O'};
for i = 1:1:length(spec)
    if ~ismember(spec{i},specInterest)
        X.(spec{i}) = Inlet.(spec{i});
    end
end

function [x0,y0] = Newton2D(Inlet,R_CH4min,R_CH4max,H2consume,Type,T,P,g0,XnP,x0,y0,Tol)
error = 1e-4;
count = 0;
while error>Tol && count<15

    e_x = max(.01*error,1e-6);
    if x0+2*e_x>=1
        e_x = .1*(x0-1);
    end
    if x0<=1e-5 %low temperature, almost no CH4 conversion
        e_x = 1e-13;
    end
    e_y = max(.01*error,1e-6);
    if y0+2*e_y>=1
        e_y = .1*(y0-1);
    end
    G11 = GibbVal2(x0,y0,Inlet,R_CH4min,R_CH4max,H2consume,Type,T,P,g0,XnP);
    a = abs(G11);
    G21 = GibbVal2(x0+e_x,y0,Inlet,R_CH4min,R_CH4max,H2consume,Type,T,P,g0,XnP)/a;
    G31 = GibbVal2(x0+2*e_x,y0,Inlet,R_CH4min,R_CH4max,H2consume,Type,T,P,g0,XnP)/a;
    G12 = GibbVal2(x0,y0+e_y,Inlet,R_CH4min,R_CH4max,H2consume,Type,T,P,g0,XnP)/a;
    G13 = GibbVal2(x0,y0+2*e_y,Inlet,R_CH4min,R_CH4max,H2consume,Type,T,P,g0,XnP)/a;
    G22 = GibbVal2(x0+e_x,y0+e_y,Inlet,R_CH4min,R_CH4max,H2consume,Type,T,P,g0,XnP)/a;
    G11 = G11/a;
    
    dGdx1 = (G21-G11)/e_x;
    dGdx2 = (G31-G21)/e_x;
    dGdx3 = (G22-G12)/e_x;
    dGdxdx = (dGdx2-dGdx1)/e_x;
    dGdxdy = (dGdx3-dGdx1)/e_y;
    
    dGdy1 = (G12-G11)/e_y;
    dGdy2 = (G13-G12)/e_y;
    dGdy3 = (G22-G21)/e_y;
    dGdydy = (dGdy2-dGdy1)/e_y;
    dGdydx = (dGdy3-dGdy1)/e_x;
    J = [dGdxdx dGdxdy; dGdydx dGdydy;];
    f = -[dGdx1;dGdy1];
    
    if dGdydy==0 && dGdxdx~=0 || abs(dGdxdx)>1e3*abs(dGdydy) %1D problem
        deltaX(1) = -dGdx1/dGdxdx;
        deltaX(2) = 0;
    elseif dGdxdx==0 && dGdydy~=0%1D problem
        deltaX(1) = 0;
        deltaX(2) = -dGdy1/dGdydy;
    elseif any(J~=0)
        deltaX = J\f;
    else deltaX= [0;0];
    end
    a = x0 + deltaX(1);
    if a>1
        scale = .75*(1-x0)/deltaX(1);%take a smaller step if iteration takes it beyond 0 or 1
    elseif a<0
        scale = .75*(x0)/(-deltaX(1));
    else scale = 1; 
    end
    x0 = x0 + deltaX(1)*scale;
    b = y0 + deltaX(2);
    if b>1
        scale = .75*(1-y0)/deltaX(2); %take a smaller step if iteration takes it beyond 0 or 1
    elseif b<0
        scale = .75*(y0)/(-deltaX(2));
    else scale = 1; 
    end
    y0 = y0 + deltaX(2)*scale;
    error = max(abs(deltaX));
    if x0<=1e-9 %low temperature, no CH4 conversion
        error = 0;
        disp('low temperature, no CH4 conversion')
    end
    count = count+1;
end

function G = GibbVal2(x,y,Inlet,R_CH4min,R_CH4max,H2consume,Type,T,P,g0,XnP)
global Ru
if isempty(R_CH4min)
    CH4 = x;
else
    CH4 = x*R_CH4max + (1-x)*R_CH4min;
end
R_COmin = -min(Inlet.CO2,Inlet.H2 + 3*CH4 - H2consume);
R_COmax = min((Inlet.H2O-CH4),(Inlet.CO+CH4));
WGS = R_COmin + y*(R_COmax -R_COmin);
X = FlowsOut(Inlet,CH4,WGS,H2consume,Type);
G = 0;
spec = fieldnames(g0);
sumX = 0;
for i = 1:1:length(spec)
    sumX = sumX + X.(spec{i});
end
sumX = sumX + XnP;
for i = 1:1:length(spec)
    G = G+X.(spec{i}).*(g0.(spec{i})/(Ru*T)+log(X.(spec{i})./sumX));
    %     G = G+X.(spec{i})*(g0.(spec{i})/(Ru*T)+log(X.(spec{i})/sumX*P));
    %     G = G+X.(spec{i})*(g0.(spec{i}) + Ru*T*log(P)+Ru*T*log(X.(spec{i})/sumX));
end