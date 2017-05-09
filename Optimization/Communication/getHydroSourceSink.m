function SourceSink = getHydroSourceSink(Date,node)
global Plant
x1 = max(1,nnz(Plant.Data.Hydro.Timestamp<Date(1))); 
equip = Plant.subNet.Hydro(node).Equipment;%all equipment at the node specified in the first part of the line name
for k = 1:1:length(equip)
    I = find(strcmp(equip(k),Plant.Data.Hydro.Equipment));%column index of this dam in the stored matrices of Data.Hydro
    if length(Date) == 1
        r = (Date - Plant.Data.Hydro.Timestamp(x1))/(Plant.Data.Hydro.Timestamp(x1+1) - Plant.Data.Hydro.Timestamp(x1));
        SourceSink = (1-r)*Plant.Data.Hydro.SourcesandSinks(x1,I) + r*Plant.Data.Hydro.SourcesandSinks(x1+1,I);
    else
        x2 = nnz(Plant.Data.Hydro.Timestamp<Date(end))+1; 
        n = nnz(Date<=Plant.Data.Hydro.Timestamp(x1));%take care of any initial conditions before data exists
        SourceSink(1:n) = Plant.Data.Hydro.Outflow(x1);
        SourceSink(n+1:end) = interp1(Plant.Data.Hydro.Timestamp(x1:x2),Plant.Data.Hydro.SourcesandSinks(x1:x2+1,I),Date);
    end
end