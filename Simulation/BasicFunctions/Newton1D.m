function y0 = Newton1D(y0,Flow,R_COmin,R_COmax,T,P,specInterest,Tol,func)
h = enthalpy(T,specInterest);
s = entropy(T,specInterest);
for i = 1:1:length(specInterest)
    g0.(specInterest{i}) = h.(specInterest{i})-T.*s.(specInterest{i});
end
if any((Flow.H2 + R_COmax)<0)
    y0 = ones(length(T),1);
else
    dY = zeros(length(T),1);
    scale = dY;
    error = dY + max(10*Tol,1e-6);
    count = 0;
    while max(error)>Tol && count<12
        e_y = max(.001*error,1e-6);
        if any(y0+2*e_y>=1)
            e_y(y0+2*e_y>=1) = -.1*abs(e_y(y0+2*e_y>=1));
        end
        WGS = y0.*R_COmax + (1-y0).*R_COmin;
        WGS_plus = (y0+e_y).*R_COmax + (1-(y0+e_y)).*R_COmin;
        WGS_minus = (y0-e_y).*R_COmax + (1-(y0-e_y)).*R_COmin;
        G11 = feval(func,Flow,WGS,T,P,g0);
        a = abs(G11);
        G12 = feval(func,Flow,WGS_plus,T,P,g0)./a;
        G13 = feval(func,Flow,WGS_minus,T,P,g0)./a;
        G11 = G11./a;

        dGdy1 = (G12-G11)./e_y;
        dGdy2 = (G11-G13)./e_y;
        dGdydy = (dGdy1-dGdy2)./e_y;
        
        dY = -dGdy1./dGdydy;
        dY(dGdydy==0) =0;
        b = y0 + dY;
        scale(:) = 1;
        if any(b>1) || any(b<0) %take a smaller step if iteration takes it beyond 0 or 1
            scale(b>1) = .75*(1-y0(b>1))./dY(b>1);
            scale(b<0) = .75*y0(b<0)./(-dY(b<0));
        end
        y0 = y0 +dY.*scale; 
        error = abs(dY);
        count = count+1;
    end
end
% disp(strcat('Newton1D count is:',num2str(count)))