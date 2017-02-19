stateAbrev = {'AL';'AK';'AZ';'AR';'CA';'CO';'CT';'DE';'FL';'GA';'HI';'ID';'IL';'IN';'IA';'KS';'KY';'LA';'ME';'MD';'MA';'MI';'MN';'MS';'MO';
              'MT';'NE';'NV';'NH';'NJ';'NM';'NY';'NC';'ND';'OH';'OK';'OR';'PA';'RI';'SC';'SD';'TN';'TX';'UT';'VT';'VA';'WA';'WV';'WI';'WY';'USavg';};

% ElecRate.Year = linspace(2012,1990,23);
% for j = 1:1:156
% ElecRate.Date(j) = datenum(2001, j,01);
% end

% A = AllUserRate;
% for i = 1:1:51
%     for j = 1:1:max(size(A))
%         if isnan(A(i,j))
%             if j>1
%                 k = i;
%                 while k+1<max(size(A))  && isnan(A(i,k))
%                     k = k+1;
%                 end
%                 if isnan(A(i,k))
%                     A(i,j) = A(i,j-1);
%                 else A(i,j) = (A(i,j-1)+A(i,k))/2;
%                 end
%             else k =1;
%                while k+1<max(size(A))  && isnan(A(i,k))
%                     k = k+1;
%                 end
%                 if isnan(A(i,k))
%                     A(i,j) = 0;
%                 else A(i,j) = A(i,k);
%                 end
%             end
%         end
%     end
% end
% AllUserRate = A;
% 
% A = ResidentialRate;
% for i = 1:1:51
%     for j = 1:1:max(size(A))
%         if isnan(A(i,j))
%             if j>1
%                 k = i;
%                 while k+1<max(size(A))  && isnan(A(i,k))
%                     k = k+1;
%                 end
%                 if isnan(A(i,k))
%                     A(i,j) = A(i,j-1);
%                 else A(i,j) = (A(i,j-1)+A(i,k))/2;
%                 end
%             else k =1;
%                while k+1<max(size(A))  && isnan(A(i,k))
%                     k = k+1;
%                 end
%                 if isnan(A(i,k))
%                     A(i,j) = 0;
%                 else A(i,j) = A(i,k);
%                 end
%             end
%         end
%     end
% end
% ResidentialRate = A;
% 
% 
% A = CommercialRate;
% for i = 1:1:51
%     for j = 1:1:max(size(A))
%         if isnan(A(i,j))
%             if j>1
%                 k = i;
%                 while k+1<max(size(A))  && isnan(A(i,k))
%                     k = k+1;
%                 end
%                 if isnan(A(i,k))
%                     A(i,j) = A(i,j-1);
%                 else A(i,j) = (A(i,j-1)+A(i,k))/2;
%                 end
%             else k =1;
%                while k+1<max(size(A))  && isnan(A(i,k))
%                     k = k+1;
%                 end
%                 if isnan(A(i,k))
%                     A(i,j) = 0;
%                 else A(i,j) = A(i,k);
%                 end
%             end
%         end
%     end
% end
% CommercialRate = A;
% 
% A = IndustrialRate;
% for i = 1:1:51
%     for j = 1:1:max(size(A))
%         if isnan(A(i,j))
%             if j>1
%                 k = i;
%                 while k+1<max(size(A))  && isnan(A(i,k))
%                     k = k+1;
%                 end
%                 if isnan(A(i,k))
%                     A(i,j) = A(i,j-1);
%                 else A(i,j) = (A(i,j-1)+A(i,k))/2;
%                 end
%             else k =1;
%                while k+1<max(size(A))  && isnan(A(i,k))
%                     k = k+1;
%                 end
%                 if isnan(A(i,k))
%                     A(i,j) = 0;
%                 else A(i,j) = A(i,k);
%                 end
%             end
%         end
%     end
% end
% IndustrialRate = A;


    
for i = 1:1:51
    state = char(stateAbrev(i));
%     ElecRate.(state).AllSector = AllUserRate(i,:);
%     ElecRate.(state).Residential = ResidentialRate(i,:);
%     ElecRate.(state).Commercial = CommercialRate(i,:);
%     ElecRate.(state).Industrial = IndustrialRate(i,:);
    ElecRate.(state).AllSectorAnnual = AllSectorAnnual(i,:);
    ElecRate.(state).ResidentialAnnual = ResidentialAnnual(i,:);
    ElecRate.(state).CommercialAnnual = CommercialAnnual(i,:);
    ElecRate.(state).IndustrialAnnual = IndustrialAnnual(i,:);
end