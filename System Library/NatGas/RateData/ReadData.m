stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';'USavg';};
for j = 1:1:300
GasRate.Date(j,1) = datenum(1989, j,01);
end
for j = 1:1:156
GasRate.Date(j,2) = datenum(2001, j,01);
end
for j = 1:1:144
GasRate.Date(j,3) = datenum(2002, j,01);
end
A = CityGateRate;
for i = 1:1:51
    for j = 1:1:max(size(A))
        if isnan(A(i,j))
            if j>1
                k = i;
                while k+1<max(size(A))  && isnan(A(i,k))
                    k = k+1;
                end
                if isnan(A(i,k))
                    A(i,j) = A(i,j-1);
                else A(i,j) = (A(i,j-1)+A(i,k))/2;
                end
            else k =1;
               while k+1<max(size(A))  && isnan(A(i,k))
                    k = k+1;
                end
                if isnan(A(i,k))
                    A(i,j) = 0;
                else A(i,j) = A(i,k);
                end
            end
        end
    end
end
CityGateRate = A;

A = ResidentialRate;
for i = 1:1:51
    for j = 1:1:max(size(A))
        if isnan(A(i,j))
            if j>1
                k = i;
                while k+1<max(size(A))  && isnan(A(i,k))
                    k = k+1;
                end
                if isnan(A(i,k))
                    A(i,j) = A(i,j-1);
                else A(i,j) = (A(i,j-1)+A(i,k))/2;
                end
            else k =1;
               while k+1<max(size(A))  && isnan(A(i,k))
                    k = k+1;
                end
                if isnan(A(i,k))
                    A(i,j) = 0;
                else A(i,j) = A(i,k);
                end
            end
        end
    end
end
ResidentialRate = A;


A = CommercialRate;
for i = 1:1:51
    for j = 1:1:max(size(A))
        if isnan(A(i,j))
            if j>1
                k = i;
                while k+1<max(size(A))  && isnan(A(i,k))
                    k = k+1;
                end
                if isnan(A(i,k))
                    A(i,j) = A(i,j-1);
                else A(i,j) = (A(i,j-1)+A(i,k))/2;
                end
            else k =1;
               while k+1<max(size(A))  && isnan(A(i,k))
                    k = k+1;
                end
                if isnan(A(i,k))
                    A(i,j) = 0;
                else A(i,j) = A(i,k);
                end
            end
        end
    end
end
CommercialRate = A;

A = IndustrialRate;
for i = 1:1:51
    for j = 1:1:max(size(A))
        if isnan(A(i,j))
            if j>1
                k = i;
                while k+1<max(size(A))  && isnan(A(i,k))
                    k = k+1;
                end
                if isnan(A(i,k))
                    A(i,j) = A(i,j-1);
                else A(i,j) = (A(i,j-1)+A(i,k))/2;
                end
            else k =1;
               while k+1<max(size(A))  && isnan(A(i,k))
                    k = k+1;
                end
                if isnan(A(i,k))
                    A(i,j) = 0;
                else A(i,j) = A(i,k);
                end
            end
        end
    end
end
IndustrialRate = A;

A = ForElectricGenRate;
for i = 1:1:51
    if nnz(isnan(A(i,:)))>30
        for j = 1:1:max(size(A))
            A(i,j) = IndustrialRate(i,j+12);
        end
    end
    for j = 1:1:max(size(A))
        if isnan(A(i,j))
            if j>1
                k = i;
                while k+1<max(size(A))  && isnan(A(i,k))
                    k = k+1;
                end
                if isnan(A(i,k))
                    A(i,j) = A(i,j-1);
                else A(i,j) = (A(i,j-1)+A(i,k))/2;
                end
            else k =1;
               while k+1<max(size(A))  && isnan(A(i,k))
                    k = k+1;
                end
                if isnan(A(i,k))
                    A(i,j) = 0;
                else A(i,j) = A(i,k);
                end
            end
        end
    end
end
ForElectricGenRate = A;

    
for i = 1:1:51
    state = char(stateAbrev(i));
    GasRate.(state).CityGate = CityGateRate(i,:);
    GasRate.(state).Residential = ResidentialRate(i,:);
    GasRate.(state).Commercial = CommercialRate(i,:);
    GasRate.(state).Industrial = IndustrialRate(i,:);
    GasRate.(state).ForElectricGeneration = ForElectricGenRate(i,:);
end