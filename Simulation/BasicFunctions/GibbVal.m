function G = GibbVal(X,WGS,T,P,g0)
global Ru
X.CO = X.CO - WGS;
X.CO2 = X.CO2 + WGS;
X.H2 = X.H2  + WGS;
X.H2O = X.H2O - WGS;

spec = fieldnames(X);
sumX = 0;
for i = 1:1:length(spec)
    sumX = sumX + X.(spec{i});
end
spec = fieldnames(g0);
G = 0;
for i = 1:1:length(spec)
    G = G+X.(spec{i}).*(g0.(spec{i})./(Ru*T)+log(X.(spec{i})./sumX));
    %     G = G+X.(spec{i})*(g0.(spec{i})/(Ru*T)+log(X.(spec{i})/sumX*P));
    %     G = G+X.(spec{i})*(g0.(spec{i}) + Ru*T*log(P)+Ru*T*log(X.(spec{i})/sumX));
end