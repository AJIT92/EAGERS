function [ map ] = RTdownscaling( input_args )
mapLen = length(SurfArea);%number of nodes in the original map
xlen = mapLen/FC.rows;%length of each node in the final map
ylen = mapLen/FC.columns;%height of each node in the final map
nodes = FC.rows*FC.columns;
Sigma = 5.670367e-11;%kW/(m^2 K^4)


map = zeros(nodes,nodes);
for i = 1:mapLen
    for j = 1:mapLen
        nncoord = ((i-1)*mapLen+j);%stacked node position in original coordinates
        nposxr = i/xlen;%node coordinates
        nmodxr = mod(i,xlen);
        nposxl = (i-1)/xlen;
        nmodxl = mod(i-1,xlen);
        nposyu = j/ylen;
        nmodyu = mod(j,ylen);
        nposyd = (j-1)/ylen;
        nmodyd = mod(j-1,ylen);

        if (nposxr-nmodxr)>(nposxl-nmodxl) %crosses a boundary in x direction
            nx0 = ((nposxr-nmodxr)-nposxl)/(nposxr-nposxl);
            nx1 = nmodxr/(nposxr-nposxl);
        else%doesn't cross boundary
            nx0 = 1;
            nx1 = 0;
        end
        if (nposyu-nmodyu)>(nposyd-nmodyd) %crosses boundary in y direction
            ny0 = ((nposyu-nmodyu)-nposyd)/(nposyu-nposyd);
            ny1 = nmodyu/(nposyu-nposxd);
        else%doesn't cross boundary
            ny0 = 1;
            ny1 = 0;
        end
        nlowx = round(nposxl-nmodxl);%need to make sure these are ints
        nhighx = round(nposxr-nmodxr);
        nlowy = round(nposyd-nmodyd);
        nhighy = round(nposyu-nmodyu);
        
        nval(1) = nx0*ny0;
        nval(2) = nx1*ny0;
        nval(3) = nx0*ny1;
        nval(4) = nx1*ny1;
        
        ncoord(1) = (nlowx-1)*FC.columns + nlowy;%LowLow
        ncoord(2) = (nhighx-1)*FC.columns + nlowy;%HighLow
        ncoord(3) = (nlowx-1)*FC.columns + nhighy;%LowHigh
        ncoord(4) = (nhighx-1)*FC.columns + nhighy;%HighHigh
            
        for k = 1:mapLen
            for m = 1:mapLen
                ttcoord = (k-1)*mapLen+m;%stacked target position in original coordinates
                tposxr = k/xlen;%target coordinates
                tmodxr = mod(k,xlen);
                tposxl = (k-1)/xlen;
                tmodxl = mod(k-1,xlen);
                tposyu = m/ylen;
                tmodyu = mod(m,ylen);
                tposyd = (m-1)/ylen;
                tmodyd = mod(m-1,ylen);
                
                if (tposxr-tmodxr)>(tposxl-tmodxl) %crosses a boundary in x direction
                    tx0 = ((tposxr-tmodxr)-tposxl)/(tposxr-tposxl);
                    tx1 = tmodxr/(tposxr-tposxl);
                else%doesn't cross boundary
                    tx0 = 1;
                    tx1 = 0;
                end
                if (tposyu-tmodyu)>(tposyd-tmodyd) %crosses boundary in y direction
                    ty0 = ((tposyu-tmodyu)-tposyd)/(tposyu-tposyd);
                    ty1 = tmodyu/(tposyu-tposxd);
                else%doesn't cross boundary
                    ty0 = 1;
                    ty1 = 0;
                end
                tlowx = round(nposxl-nmodxl);%need to make sure these are ints
                thighx = round(nposxr-nmodxr);
                tlowy = round(nposyd-nmodyd);
                thighy = round(nposyu-nmodyu);
                
                tval(1) = tx0*ty0;
                tval(2) = tx1*ty0;
                tval(3) = tx0*ty1;
                tval(4) = tx1*ty1;
                
                tcoord(1) = (tlowx-1)*FC.columns + tlowy;%LowLow
                tcoord(2) = (thighx-1)*FC.columns + tlowy;%HighLow
                tcoord(3) = (tlowx-1)*FC.columns + thighy;%LowHigh
                tcoord(4) = (thighx-1)*FC.columns + thighy;%HighHigh
                
                for nc = 1:length(ncoord)
                    for tc = 1:length(tcoord)
                        map(ncoord(nc),tcoord(tc)) = map(ncoord(nc),tcoord(tc)) + nval(nc)*tval(tc)*Sigma*SurfArea(nncoord)*ViewFactor(nncoord,ttcoord);
                    end
                end
            end
        end

    end
end

end

