function SSnew = changeTimestep(SS,Tnew,Toriginal)
SSmodel = ss(SS.A, SS.B,SS.C,SS.D,Toriginal);
SSmodel = d2d(SSmodel,Tnew);
r = size(SSmodel);
SSnew.A = SSmodel(1,1).a;
q = size(SSnew.A);
SSnew.B = zeros(q(1),r(2));
for i = 1:1:r(2)
    SSnew.B(:,i) = SSmodel(1,i).b;
end
SSnew.C = zeros(r(1),q(2));
for i = 1:1:r(1)
    SSnew.C(i,:) = SSmodel(i,1).c;
end