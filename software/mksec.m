function mksec(v1,v2,v3,dem,s,m)
load color/color_rainbow.mat
[A,R1] = readgeoraster(dem);
A=double(A);
[lata,lona] = geographicGrid(R1);

[B1,~] = readgeoraster(v1);
[B2,~] = readgeoraster(v2);
[B3,R2] = readgeoraster(v3);
B1=double(B1);
B2=double(B2);
B3=double(B3);
[lat,lon] = geographicGrid(R2);

CR=R1;
C1 = zeros(CR.RasterSize);
C1 = C1- Inf;
C1 = imbedm(lat,lon,B1,C1,CR);


C2 = zeros(CR.RasterSize);
C2 = C2- Inf;
C2 = imbedm(lat,lon,B2,C2,CR);

C3 = zeros(CR.RasterSize);
C3 = C3- Inf;
C3 = imbedm(lat,lon,B3,C3,CR);


imagesc(C2);
colormap(color_rainbow)
hold on 
[x,y]=getline;
close all
line=GetLinePixels([x(1), y(1)], [x(2), y(2)]);
for i=1:length(line)
    C11(i)=C1(line(i,2),line(i,1));
    C22(i)=C2(line(i,2),line(i,1));
    C33(i)=C3(line(i,2),line(i,1));
    A11(i)=A(line(i,2),line(i,1));
end
CC=sqrt(C11.^2+C22.^2);

lat1=lata(round(y(1)),round(x(1)));
lon1=lona(round(y(1)),round(x(1)));
lat2=lata(round(y(2)),round(x(2)));
lon2=lona(round(y(2)),round(x(2)));
lat1 = deg2rad(lat1);
lon1 = deg2rad(lon1);
lat2 = deg2rad(lat2);
lon2 = deg2rad(lon2);

radius = 6371;  
dlat = lat2 - lat1;
dlon = lon2 - lon1;
a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
distance = radius * c;
x1=1:1:length(line);

numPoints = 100000;
xi = linspace(min(x1), max(x1), numPoints); 
yi = interp1(x1, A11, xi, 'spline');
xi1=[xi max(xi) 1];
yi1=[yi 0 0];
patch(xi1,yi1,[0.925490196 0.925490196 0.925490196],'EdgeColor','none');
hold on
plot(xi,yi,'LineWidth',1.5,'Color','black');
% plot(x1,A11);
hold on
sc=((max(A11)-min(A11))/(max(x1)-min(x1)))*5;
quiver(x1,A11,CC/sc,C33,s,'LineWidth',0.8,'Color','r','MaxHeadSize',m);
set(gcf,'unit','centimeters','position',[10 10 20 6])
set(gca,'linewidth',1.5,'Fontname','Airl','fontsize',12); 
xlabel('Distance(km)');
ylabel('Altitude(m)');
xlim([0,max(x1)+5]);
ylim([min(A11)-(max(A11)-min(A11))/10,max(A11)+(max(A11)-min(A11))/10]);
set(gca,'XTick',[0,round(max(x1)/3),round(2*max(x1)/3),max(x1)], ...
    'xticklabel',{'0',num2str(round(distance/3,2)),num2str(round(2*distance/3,2)),num2str(round(distance,2))})
box off

end
function pixels = GetLinePixels(point1, point2)
x1 = round(point1(1));
y1 = round(point1(2));
x2 = round(point2(1));
y2 = round(point2(2));

dx = abs(x2 - x1);
dy = abs(y2 - y1);
sx = sign(x2 - x1);
sy = sign(y2 - y1);
err = dx - dy;

pixels = [x1, y1];

while ~(x1 == x2 && y1 == y2)
    e2 = 2 * err;
    if e2 > -dy
        err = err - dy;
        x1 = x1 + sx;
    end
    if e2 < dx
        err = err + dx;
        y1 = y1 + sy;
    end
    pixels = [pixels; x1, y1];
end
end