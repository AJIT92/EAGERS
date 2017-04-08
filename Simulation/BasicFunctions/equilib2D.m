function [Out,R] = equilib2D(InFlow,T,P,H2consume,COconsume,Type,PercEquilib,guess)
% performs a gradient search to find the minimum gibbs function value for
% two simultaneous reactions (methane reforming and water gas shift)
% tolerance is a relative tolence, so that it is more accurate if only a
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
for k = 1:1:n %treat each node independently
    for i = 1:1:length(specInterest)
        g0.(specInterest{i}) = h.(specInterest{i})(k)-T(k).*s.(specInterest{i})(k);
    end
    for i = 1:1:length(spec)
        Inlet.(spec{i}) = InFlow.(spec{i})(k);
    end
    CH4max = min(InFlow.CH4(k),InFlow.H2O(k)+H2consume(k));
    CH4min = -min(InFlow.CO(k)-COconsume,(InFlow.H2(k)+InFlow.CO(k)-COconsume-H2consume(k))/4);
    Span1 = CH4max-CH4min;
    if ~isempty(guess)
        x0 = guess(k);
        R.CH4(k,1) = (Span1*x0+CH4min);
        
        X.CH4 = InFlow.CH4(k) - R.CH4(k);
        X.CO = InFlow.CO(k) + R.CH4(k) - COconsume(k);
        X.CO2 = InFlow.CO2(k) + COconsume(k);
        X.H2 = InFlow.H2(k) + 3*R.CH4(k) - H2consume(k);%hydrogen consumed
        X.H2O = InFlow.H2O(k) - R.CH4(k) + H2consume(k);% water produced
        if any(strcmp(Type,{'MCFC';'MCEC'}))
            X.CO2 = X.CO2 + H2consume(k); % CO2 brought over
        end
        for i = 1:1:length(spec)
            if ~ismember(spec{i},specInterest)
                X.(spec{i}) = InFlow.(spec{i})(k);
            end
        end
        y0 = Newton1D(.5,X,-min(X.CO2,X.H2),min(X.H2O,X.CO),T(k),P,specInterest,1e-6,'GibbVal');
        Tol = max(1e-10,x0*1e-7);
    else
        Tol = 1e-6;
        x0 = 0.85;
        y0 = 0.5;
    end

    [x,y] = Newton2D(Inlet,CH4min,CH4max,H2consume(k),COconsume(k),Type,T(k),P,g0,x0,y0,Tol);
    R.CH4(k,1) = max(CH4min,PercEquilib*(Span1*x+CH4min)); %prevents a low % equilibrium from causing problems with enough H2 for electrochemistry
    X.CH4 = InFlow.CH4(k) - R.CH4(k);
    X.CO = InFlow.CO(k) + R.CH4(k) - COconsume(k);
    X.CO2 = InFlow.CO2(k) + COconsume(k);
    X.H2 = InFlow.H2(k) + 3*R.CH4(k) - H2consume(k);%hydrogen consumed
    X.H2O = InFlow.H2O(k) - R.CH4(k) + H2consume(k);% water produced
    if any(strcmp(Type,{'MCFC';'MCEC'}))
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
    
    Out.CH4(k) = InFlow.CH4(k) - R.CH4(k);
    Out.CO(k) = InFlow.CO(k) + R.CH4(k) - R.WGS(k) - COconsume(k);
    Out.CO2(k) = InFlow.CO2(k) + R.WGS(k) + COconsume(k);
    Out.H2(k) = InFlow.H2(k) + 3*R.CH4(k) + R.WGS(k) - H2consume(k);
    Out.H2O(k) = InFlow.H2O(k) - R.CH4(k) - R.WGS(k) +H2consume(k);
end
Out.T = T;


function [x0,y0] = Newton2D(Inlet,R_CH4min,R_CH4max,H2consume,COconsume,Type,T,P,g0,x0,y0,Tol)
%standard 2-D newtonian method (gradient search) for a minimum Gibbs energy
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
    %%G11
    [X,WGS] = FlowOut(Inlet,x0,y0,R_CH4min,R_CH4max,H2consume,COconsume,Type);
    G11 = GibbVal(X,WGS,T,P,g0);
    a = abs(G11);
    [X,WGS] = FlowOut(Inlet,x0+e_x,y0,R_CH4min,R_CH4max,H2consume,COconsume,Type);
    G21 = GibbVal(X,WGS,T,P,g0)/a;
    
    [X,WGS] = FlowOut(Inlet,x0+2*e_x,y0,R_CH4min,R_CH4max,H2consume,COconsume,Type);
    G31 = GibbVal(X,WGS,T,P,g0)/a;
    
    [X,WGS] = FlowOut(Inlet,x0,y0+e_y,R_CH4min,R_CH4max,H2consume,COconsume,Type);
    G12 = GibbVal(X,WGS,T,P,g0)/a;
    
    [X,WGS] = FlowOut(Inlet,x0,y0+2*e_y,R_CH4min,R_CH4max,H2consume,COconsume,Type);
    G13 = GibbVal(X,WGS,T,P,g0)/a;
    
    [X,WGS] = FlowOut(Inlet,x0+e_x,y0+e_y,R_CH4min,R_CH4max,H2consume,COconsume,Type);
    G22 = GibbVal(X,WGS,T,P,g0)/a;
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

function [X,WGS] = FlowOut(Inlet,x,y,R_CH4min,R_CH4max,H2consume,COconsume,Type)
%Calculates outlet concentrations without accounting for WGS reaction (done in GibbVal)
if isempty(R_CH4min)
    CH4 = x;
else
    CH4 = x*R_CH4max + (1-x)*R_CH4min;
end
R_COmin = -min(Inlet.CO2-COconsume,Inlet.H2 + 3*CH4 - H2consume);
R_COmax = min((Inlet.H2O + H2consume - CH4),(Inlet.CO+CH4-COconsume));
WGS = R_COmin + y*(R_COmax -R_COmin);
%find exit species not considering WGS
spec = fieldnames(Inlet);
for i = 1:1:length(spec)
    X.(spec{i}) = Inlet.(spec{i});
end
X.CH4 = Inlet.CH4 - CH4;
X.CO = Inlet.CO + CH4 - COconsume;
X.CO2 = Inlet.CO2 + COconsume;
X.H2 = Inlet.H2 + 3* CH4 - H2consume;%hydrogen consumed
X.H2O = Inlet.H2O - CH4 + H2consume;% water produced
if any(strcmp(Type,{'MCFC';'MCEC'}))
    X.CO2 = X.CO2 + H2consume; % CO2 brought over
end