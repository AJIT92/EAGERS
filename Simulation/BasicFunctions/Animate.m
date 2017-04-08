function Animate(Time,Y,block,FigNum)
%% Animate, produces several graphs and videos of transients
%%Static Background & Coordinates
global Model_dir TagInf
FCstates_alltime = zeros(length(Y(:,1)),length(block.States));
columns = block.columns;
rows = block.rows;
nodes = block.nodes;
L_Cell = block.L_Cell;
W_Cell = block.W_Cell;
if block.direction ==1
    Xlabel = 'Fuel/Air Flow Direction -->';
    Ylabel = 'Cell Width (parallel channels)';
elseif block.direction ==2
    Xlabel = 'Fuel Flow --> & Air Flow <--';
    Ylabel = 'Cell Width (parallel channels)';
elseif block.direction ==3
    Xlabel = 'Fuel Flow Direction -->';
    Ylabel = 'Air Flow Direction --> ';
end
wC = L_Cell/columns;
wR = W_Cell/rows;
if rows ==1
    Xcoord = [0 linspace(.5*wC,L_Cell -.5*wC,columns) L_Cell];
    Ycoord = [0 W_Cell];
    Zcoord(1,1:columns+2) = 0;
    Zcoord(2,1:columns+2) = 0;
    %Draw line coordinates
    for j = 1:1:rows+1
        HorizLinesX(3*j-2:3*j) = [0 L_Cell 0];
        HorizLinesY(3*j-2:3*j) = (j-1)*wR;
        HorizLinesZ(3*j-2:3*j) = 0;
    end
    for j = 1:1:columns+1
        VertLinesX(3*j-2:3*j) =(j-1)*wC;
        VertLinesY(3*j-2:3*j) = [0 W_Cell 0];
        VertLinesZ(3*j-2:3*j) = 0;
    end
    %draw box coordinates
    Box1X = [0 L_Cell L_Cell 0 0];
    Box1Y = [0 0 W_Cell W_Cell 0];
    Box1Z = [0 0 0 0 0];
elseif Manifold== 0 && rows>1
    %Organize Temperature & node coordinates
    Xcoord = [0 linspace(.5*wC,L_Cell -.5*wC,columns) L_Cell];
    Ycoord = [0 linspace(.5*wR,L_Cell -.5*wR,rows) W_Cell];
    for j = 1:1:rows
        Zcoord(j+1,2:columns+1) = 0;
    end
    Zcoord(1,2:columns+2) = 0;
    Zcoord(2:rows+2,columns+2) = 0;
    Zcoord(rows+2,linspace(columns+1,1,columns+1)) = 0;
    Zcoord(linspace(rows+1,1,rows+1),1) = 0;
    %Draw line coordinates
    for j = 1:1:rows+1
        HorizLinesX(3*j-2:3*j) = [0 L_Cell 0];
        HorizLinesY(3*j-2:3*j) = (j-1)*wR;
        HorizLinesZ(3*j-2:3*j) = 0;
    end
    for j = 1:1:columns+1
        VertLinesX(3*j-2:3*j) =(j-1)*wC;
        VertLinesY(3*j-2:3*j) = [0 W_Cell 0];
        VertLinesZ(3*j-2:3*j) = 0;
    end
    %draw box coordinates
    Box1X = [0 L_Cell L_Cell 0 0];
    Box1Y = [0 0 W_Cell W_Cell 0];
    Box1Z = [0 0 0 0 0];
elseif Manifold == 1
    wE1 = (L_Plate-L_Cell)/2;
    wE2 = (W_Plate-W_Cell)/2;
    %Organize Temperature & node coordinates
    Xcoord = [0 .5*wE1 linspace(wE1+.5*wC,L_Plate -(wE1+.5*wC),columns) L_Plate-.5*wE1 L_Plate];
    Ycoord = [0 .5*wE2 linspace(wE2+.5*wR,W_Plate -(wE2+.5*wR),rows) W_Plate-.5*wE2 W_Plate];
    for j = 1:1:rows
        Zcoord(j+2,3:columns+2) = 0;
    end
    %Edges
    Zcoord(2,3:columns+3) = 0;
    Zcoord(3:rows+3,columns+3) = 0;
    Zcoord(rows+3,linspace(columns+2,2,columns+1)) = 0;
    Zcoord(linspace(rows+2,2,rows+1),2) = 0;

    Zcoord(1,2:columns+4) = 0;
    Zcoord(2:rows+4,columns+4) = 0;
    Zcoord(rows+4,linspace(columns+3,1,columns+3)) = 0;
    Zcoord(linspace(rows+3,1,rows+3),1) = 0;

    %Draw line coordinates
    for j = 1:1:rows+3
        HorizLinesX(3*j-2:3*j) = [0 L_Plate 0];
        HorizLinesY(3*j-2:3*j) = wE2+(j-2)*wR;
        HorizLinesZ(3*j-2:3*j) = 0;
    end
    HorizLinesY(1:3) = 0;
    HorizLinesY(3*j-2:3*j) = W_Plate;
    for j = 1:1:columns+3
        VertLinesX(3*j-2:3*j) = wE1+(j-2)*wC;
        VertLinesY(3*j-2:3*j) = [0 W_Plate 0];
        VertLinesZ(3*j-2:3*j) = 0;
    end
    VertLinesX(1:3) = 0;
    VertLinesX(3*j-2:3*j) = L_Plate;
    %draw box coordinates
    Box1X = [0 L_Plate L_Plate 0 0];
    Box1Y = [0 0 W_Plate W_Plate 0];
    Box1Z = [0 0 0 0 0];
    Box2X = [wE1 wE1+L_Cell wE1+L_Cell wE1 wE1];
    Box2Y = [wE2 wE2 wE2+W_Cell wE2+W_Cell wE2];
    Box2Z = [0 0 0 0 0];
end

%% Dynamic
VideoName = strcat(component,datestr(now),'.avi');
filename = fullfile(Model_dir,'graphics','Videos',VideoName);
vidObj = VideoWriter(filename);

points = length(Time);
PENavgT = zeros(points,1);
Voltage = zeros(points,1);
Current = zeros(points,1);
Power = zeros(points,1);
for t = 1:1:points
    FCstates_alltime(t,:) = FCstates;
    FCstates = Y(t,block.States).*block.Scale';
    PENavgT(t) = mean(FCstates(2*nodes+1:3*nodes));
    Voltage(t) = TagInf.(block.name).Voltage(t);
    Current(t) = sum(FCstates(end-nodes-1:end-2));
    Power(t) = Voltage(t)*Current(t);
end
Frames = 100;
dt = ((time2+100)-(time1-50))/Frames;

point1 = 1;
while Time(point1)<time1-50
    point1 = point1+1;
end
point1 = point1-1;
point2 = point1+2;
while Time(point2)<time2+100 && point2<points
    point2 = point2+1;
end
TimePlot = (Time(1:point2)-(time1-50))/minute;
strTime= 'Time (Minutes)';
%% Time frame for animation
posttime=1;
for f = 1:1:Frames+1
    while Time(posttime)<(time1-50)+(f-1)*dt
        posttime = posttime+1;
    end
    pretime = posttime -1;
    stepSize = Time(posttime)-Time(pretime);
    a = (Time(posttime)-((time1-50)+(f-1)*dt))/stepSize;
    %% Interpolate data for all frams from start of transient to t + 100s
    Yt = a*FCstates_alltime(pretime,:)+(1-a)*FCstates_alltime(posttime,:);
    CellT = Yt(2*nodes+1:3*nodes);
    C = Zcoord*0;
    if Manifold== 0
        for j = 1:1:rows
            C(j+1,2:columns+1) = CellT(1+(j-1)*columns:j*columns);
        end
        C(1,2:columns+2) =  [CellT(1:columns) CellT(columns)];
        C(2:rows+2,columns+2) = [CellT(linspace(columns,rows*columns,rows)) CellT(rows*columns)];
        C(rows+2,linspace(columns+1,1,columns+1)) = [CellT(linspace(rows*columns,(rows-1)*columns+1,columns)) CellT((rows-1)*columns+1)];
        C(linspace(rows+1,1,rows+1),1) = [CellT(linspace((rows-1)*columns+1,1,rows)) CellT(1)];
    elseif Manifold == 1
        EdgeT = a*EdgeStates(pretime,2*(rows+columns+2)+1:4*(rows+columns+2)+1)+(1-a)*EdgeStates(pretime,2*(rows+columns+2)+1:4*(rows+columns+2)+1);
        for j = 1:1:rows
            C(j+2,3:columns+2) = CellT(1+(j-1)*columns:j*columns);
        end
        %Edges
        C(2,3:columns+3) =  EdgeT(1:columns+1);
        C(3:rows+3,columns+3) = EdgeT(columns+2:rows+columns+2);
        C(rows+3,linspace(columns+2,2,columns+1)) = EdgeT(rows+columns+3:rows+2*columns+3);
        C(linspace(rows+2,2,rows+1),2) = EdgeT(rows+2*columns+4:2*rows+2*columns+4);

        C(1,2:columns+4) =  [EdgeT(1) EdgeT(1:columns+1) EdgeT(columns+1)];
        C(2:rows+4,columns+4) = [ EdgeT(columns+2) EdgeT(columns+2:rows+columns+2) EdgeT(rows+columns+2)];
        C(rows+4,linspace(columns+3,1,columns+3)) = [EdgeT(rows+columns+3) EdgeT(rows+columns+3:rows+2*columns+3) EdgeT(rows+2*columns+3)];
        C(linspace(rows+3,1,rows+3),1) = [EdgeT(rows+2*columns+4) EdgeT(rows+2*columns+4:2*rows+2*columns+4) EdgeT(2*rows+2*columns+4)];
    end

       %%%%%%%Begin Graphing Portion
    if FCtype ==1
        Tmin= 800;
        Tmax=1000;
    elseif FCtype ==0
        Tmin = 900;
        Tmax = 1100;
    end
    time = [(f-1)*dt, (f-1)*dt]; 
    one = [-inf,inf];

    %%%%%%Electrolyte
    hidden on
    h=figure(FigNum);
    subplot('Position',[.1 .125 .4 .775]), surf(Xcoord, Ycoord, Zcoord, C,'FaceColor','interp','EdgeColor','none')
    caxis([Tmin Tmax])
    colorbar
    hold on
    plot3(HorizLinesX,HorizLinesY,HorizLinesZ,'k-','LineWidth',1)
    plot3(VertLinesX,VertLinesY,VertLinesZ,'k-','LineWidth',1)
    plot3(Box1X,Box1Y,Box1Z,'k-','LineWidth',4)
    if Manifold ==1
        plot3(Box2X,Box2Y,Box2Z,'k-','LineWidth',3)
    end
    if Manifold ==0
        axis([0 L_Cell 0 W_Cell])
    else
        axis([0 L_Plate 0 W_Plate])
    end
    view([180 90])
    xlabel(Xlabel)
    ylabel(Ylabel)
    hold off
    title('Electrolyte Temperature Profile','fontsize',14)
    xlabel(Xlabel,'Fontsize',14)
    ylabel(Ylabel,'Fontsize',14)
    subplot('Position',[.6 .125 .3 .775]), [Ax,H1,H2] =plotyy(TimePlot(point1:point2),PENavgT(point1:point2),TimePlot(point1:point2),Power(point1:point2));
    set(H1, 'LineWidth',2,'Color','red','LineStyle','-.')
    set(H2, 'LineWidth',2,'Color','cyan','LineStyle','-')
    set(Ax(1),'XColor','k','YColor','k','Fontsize',12)
    set(Ax(2),'XColor','k','YColor','k','Fontsize',12)
    hold(Ax(1));
    plot(Ax(1),time/minute,one,'b-','LineWidth',3)
    hold off
    xlabel(strTime,'Fontsize',16,'Color',[0 0 0])
    axis(Ax(1),[0 TimePlot(length(TimePlot)) Tmin Tmax])
    axis(Ax(2),[0 TimePlot(length(TimePlot)) 0 10*round(0.11*max(Power(point1:point2)))])
    ylabel(Ax(1),'PEN Avg Temp (K)','Fontsize',16,'Color',[0 0 0])
    ylabel(Ax(2),'Stack Power (%)','Fontsize',16,'Color',[0 0 0])
    legend(Ax(1),'Stack Power','PEN Avg Temp','Location','NorthEast')
    set(Ax(1),'YTick',Tmin:20:Tmax)
    set(Ax(2),'YTick',0:ceil(0.11*max(Power(point1:point2))):10*ceil(0.11*max(Power(point1:point2))))
    %% Create movie frame
    writeVideo(vidObj,currFrame);
%         Electrolyte(f)=getframe(h,[25 0 975 440]);
end %ends frame loop      
close(vidObj);    
    
% movie2avi(Electrolyte, 'SOFCPowerSweep')%%Make frames into a movie
% movie(h,Electrolyte,1);