function mkfigure_v(obj,v1,v2,v3,dem,scale,p,q,varargin)
load color/color_rainbow.mat
if length(varargin)==1
    shp=varargin{1,1};
    D=shaperead(shp);
    D=struct2table(D);
end

[B1,~] = readgeoraster(v1);
[B2,~] = readgeoraster(v2);
[B3,R2] = readgeoraster(v3);
B1=double(B1);
B2=double(B2);
B3=double(B3);
if strcmpi('insar',obj.tech)==1
    B1=B1*1000;
    B2=B2*1000;
    B3=B3*1000;
end
[lat,lon] = geographicGrid(R2);

if strcmpi('no',dem)==1
    load mkfigure/color_rainbow2.mat
    figW = 0.9;
    figH = (size(B3,1)/size(B3,2))*0.9;
    Z1=B3;
    Z1(:)=8000;
    Z2=B3;
    Z2(:)=0;
    indh=1:p:size(B3,1);
    indl=1:p:size(B3,2);
    A1=B3;
    A1(isnan(A1))=0;B1(isnan(B1))=0;B2(isnan(B2))=0;
    A1(indh,indl)=0;
    B1(~A1==0)=0;
    B2(~A1==0)=0;
    B3(B3==0)=nan;
    figure
    ax1=subplot(1,2,1);
    ax=worldmap(R2.LatitudeLimits,R2.LongitudeLimits);
    setm(ax,"FontSize",12,"FontName","Airal")
    meshm(B3,R2);
    colormap(color_rainbow)
    if length(varargin)==1
        for j=1:size(D,1)
            Z=cell2mat(D.Y(j));
            Z(:)=8000;
            plot3m(cell2mat(D.Y(j)),cell2mat(D.X(j)),Z,"Color","black");
            clear Z
            hold on
        end
    end
    quiver3m(lat,lon,Z1,B1,B2,Z2,"black",q);
    set(gca,"FontSize",12,"FontName","Airal");
    c=colorbar(ax1,'Position',[figW-0.057 0.12 0.02 figH/3],"FontSize",12);
    if strcmpi('insar',obj.tech)==1
        c.Label.String = 'velocity(mm/a)';
    else
        c.Label.String = 'velocity(m/a)';
    end
    ax2=subplot(1,2,2);
    ax=worldmap(R2.LatitudeLimits,R2.LongitudeLimits);
    setm(ax,"FontSize",12,"FontName","Airal",'FontColor','none','GColor','none')
    scaleruler('units','km')
    setm(handlem('scaleruler1'),'RulerStyle',"lines",'TickDir','down',"FontSize",12,"LineWidth",1,"FontName","Airal",'MajorTick',0:scale/2:scale,'MinorTick',0:scale/4:scale/2)
    northarrow('latitude',R2.LatitudeLimits(2)-diff(R2.LatitudeLimits)/10,'longitude',R2.LongitudeLimits(2)-diff(R2.LongitudeLimits)/20)
    framem('off')
    set(ax1,"Position",[0.12 0.12 figW figH]);
    set(ax2,"Position",[0.1 0.1 figW figH]);
else

    [A,R1] = readgeoraster(dem);
    A=double(A);
    [lat1,lon1] = geographicGrid(R1);

    CR=R1;
    C1 = zeros(CR.RasterSize);
    C1 = C1- Inf;
    C1 = imbedm(lat,lon,B1,C1,CR);
    C1(isnan(C1))= -Inf;
    C1(C1==0)= -Inf;

    C2 = zeros(CR.RasterSize);
    C2 = C2- Inf;
    C2 = imbedm(lat,lon,B2,C2,CR);
    C2(isnan(C2))= -Inf;
    C2(C2==0)= -Inf;

    C3 = zeros(CR.RasterSize);
    C3 = C3- Inf;
    C3 = imbedm(lat,lon,B3,C3,CR);
    C3(isnan(C3))= -Inf;
    C3(C3==0)= -Inf;

    Z1=A;
    Z1(:)=8000;
    Z2=A;
    Z2(:)=0;

    C1(C1==-inf)=0;
    C2(C2==-inf)=0;
    indh=1:p:size(A,1);
    indl=1:p:size(A,2);
    A1=A;
    A1(indh,indl)=0;
    C1(~A1==0)=0;
    C2(~A1==0)=0;


    figW = 0.9;
    figH = (size(A,1)/size(A,2))*0.9;

    figure
    ax1=subplot(1,2,1);
    ax=worldmap(R1.LatitudeLimits,R1.LongitudeLimits);
    setm(ax,"FontSize",12,"FontName","Airal")
    h = meshm(A,R1,size(A),A);
    %plot3m(cell2mat(D.Y),D.X,D.Z,"Color","black");
    if length(varargin)==1
        for j=1:size(D,1)
            Z=cell2mat(D.Y(j));
            Z(:)=8000;
            plot3m(cell2mat(D.Y(j)),cell2mat(D.X(j)),Z,"Color","black");
            clear Z
            hold on
        end
    end
    quiver3m(lat1,lon1,Z1,C1,C2,Z2,"black",q);
    h.CData = C3;
    h.FaceColor = 'texturemap';
    colormap(color_rainbow)
    material dull %dull shiny
    camlight
    lighting gouraud
    axis tight
    box off
    set(gca,"FontSize",12,"FontName","Airal");
    c=colorbar(ax1,'Position',[figW-0.057 0.12 0.02 figH/3],"FontSize",12);
    if strcmpi('insar',obj.tech)==1
        c.Label.String = 'velocity(mm/a)';
    else
        c.Label.String = 'velocity(m/a)';
    end
    ax2=subplot(1,2,2);
    ax=worldmap(R1.LatitudeLimits,R1.LongitudeLimits);
    setm(ax,"FontSize",12,"FontName","Airal",'FontColor','none','GColor','none')
    scaleruler('units','km')
    setm(handlem('scaleruler1'),'RulerStyle',"lines",'TickDir','down',"FontSize",12,"LineWidth",1,"FontName","Airal",'MajorTick',0:scale/2:scale,'MinorTick',0:scale/4:scale/2)
    northarrow('latitude',R1.LatitudeLimits(2)-diff(R1.LatitudeLimits)/10,'longitude',R1.LongitudeLimits(2)-diff(R1.LongitudeLimits)/20)
    framem('off')
    set(ax1,"Position",[0.12 0.12 figW figH]);
    set(ax2,"Position",[0.1 0.1 figW figH]);
end