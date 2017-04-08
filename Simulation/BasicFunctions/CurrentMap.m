function CurrentMap(Y,block,FigNum)
% plots the fuel cell or electrolyzer current distribution
columns = block.columns;
rows = block.rows;
nodes = block.nodes;
switch block.direction
    case 'coflow'
        Xlabel = 'Fuel/Air Flow Direction -->';
        Ylabel = 'Cell Width (parallel channels)';
    case 'counterflow'
        Xlabel = 'Fuel Flow --> & Air Flow <--';
        Ylabel = 'Cell Width (parallel channels)';
    case 'crossflow'
        Xlabel = 'Cell Length (m)';
        Ylabel = 'Cell Width (m)';
end

wC = block.L_Cell/columns;
wR = block.W_Cell/rows;
FCstates = Y(block.States).*block.Scale';
CurrentDen = FCstates(end-nodes-1:end-2)/block.A_Node/1e4;%current density in A/cm^2
if rows ==1
    X = [0 linspace(.5*wC,block.L_Cell -.5*wC,columns) block.L_Cell];
    Y = [0 block.W_Cell];
    Z(1,1:columns+2) = 0;
    C(1,1:columns+2) =  [CurrentDen(1) CurrentDen(1:columns) CurrentDen(columns)];
    Z(2,1:columns+2) = 0;
    C(2,1:columns+2) =  [CurrentDen(1) CurrentDen(1:columns) CurrentDen(columns)];
    %Draw line coordinates
    for j = 1:1:rows+1
        HorizLinesX(3*j-2:3*j) = [0 block.L_Cell 0];
        HorizLinesY(3*j-2:3*j) = (j-1)*wR;
        HorizLinesZ(3*j-2:3*j) = 0;
    end
    for j = 1:1:columns+1
        VertLinesX(3*j-2:3*j) =(j-1)*wC;
        VertLinesY(3*j-2:3*j) = [0 block.W_Cell 0];
        VertLinesZ(3*j-2:3*j) = 0;
    end
elseif rows>1
    %Organize Temperature & node coordinates
    X = [0 linspace(.5*wC,block.L_Cell -.5*wC,columns) block.L_Cell];
    Y = [0 linspace(.5*wR,block.L_Cell -.5*wR,rows) block.W_Cell];
    for j = 1:1:rows
        Z(j+1,2:columns+1) = 0;
        C(j+1,2:columns+1) = CurrentDen(1+(j-1)*columns:j*columns);
    end
    Z(1,2:columns+2) = 0;
    C(1,2:columns+2) =  [CurrentDen(1:columns) CurrentDen(columns)];
    Z(2:rows+2,columns+2) = 0;
    C(2:rows+2,columns+2) = [CurrentDen(linspace(columns,rows*columns,rows)) CurrentDen(rows*columns)];
    Z(rows+2,linspace(columns+1,1,columns+1)) = 0;
    C(rows+2,linspace(columns+1,1,columns+1)) = [CurrentDen(linspace(rows*columns,(rows-1)*columns+1,columns)) CurrentDen((rows-1)*columns+1)];
    Z(linspace(rows+1,1,rows+1),1) = 0;
    C(linspace(rows+1,1,rows+1),1) = [CurrentDen(linspace((rows-1)*columns+1,1,rows)) CurrentDen(1)];
    %Draw line coordinates
    for j = 1:1:rows+1
        HorizLinesX(3*j-2:3*j) = [0 block.L_Cell 0];
        HorizLinesY(3*j-2:3*j) = (j-1)*wR;
        HorizLinesZ(3*j-2:3*j) = 0;
    end
    for j = 1:1:columns+1
        VertLinesX(3*j-2:3*j) =(j-1)*wC;
        VertLinesY(3*j-2:3*j) = [0 block.W_Cell 0];
        VertLinesZ(3*j-2:3*j) = 0;
    end
end

figure(FigNum);
surf(X, Y, Z, C,'FaceColor','interp','EdgeColor','none')
caxis([min(CurrentDen),max(CurrentDen)]);
hold on
plot3(HorizLinesX,HorizLinesY,HorizLinesZ,'k-','LineWidth',1);
plot3(VertLinesX,VertLinesY,VertLinesZ,'k-','LineWidth',1)
hold off
axis([0 block.L_Cell 0 block.W_Cell])
view([0 90])
xlabel(Xlabel,'FontSize',14)
ylabel(Ylabel,'FontSize',14)
c = colorbar;
% c.Label.String = 'Current Density(A/cm^2)';