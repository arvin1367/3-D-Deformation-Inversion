function mkgifs(obj,v,dem,scale,colrange,name,varargin)
load color/color_rainbow.mat
mkdir picture
imds=imageDatastore(v,"FileExtensions",".tif","IncludeSubfolders",true);
if length(varargin)==1
    shp=varargin{1,1};
    D=shaperead(shp);
    D=struct2table(D);
end

for i=1:length(imds.Files)

    [B,R2] = readgeoraster(imds.Files{i});
    B=double(B);
    [lat,lon] = geographicGrid(R2);

    if strcmpi('insar',obj.tech)==1
        B=B*1000;
    end

    ind=strfind(imds.Files{i},'\');
    str=imds.Files{i};
    nstr=str(ind(end)+1:end);
    nstr1=extractBefore(nstr,'.');
    nstr2=extractBefore(nstr,'_');
    nstr3=extractAfter(nstr,'_');
    tstr=[nstr2,'--',nstr3];
    nstr11=['picture\',nstr1];

    if strcmpi('no',dem)==1
        %load mkfigure/color_rainbow_up.mat
        figW = 0.9;
        figH = (size(B,1)/size(B,2))*0.9;

        B(B==0)=nan;
        figure
        ax1=subplot(1,2,1);
        ax=worldmap(R2.LatitudeLimits,R2.LongitudeLimits);
        setm(ax,"FontSize",12,"FontName","Airal")
        meshm(B,R2);
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
        set(gca,"FontSize",12,"FontName","Airal");
        c=colorbar(ax1,'Position',[figW-0.057 0.12 0.02 figH/3],"FontSize",12);
        caxis(colrange)
        if strcmpi('insar',obj.tech)==1
            c.Label.String = 'velocity(mm/d)';
        else
            c.Label.String = 'velocity(m/d)';
        end
        title(tstr)

        ax2=subplot(1,2,2);
        ax=worldmap(R2.LatitudeLimits,R2.LongitudeLimits);
        setm(ax,"FontSize",12,"FontName","Airal",'FontColor','none','GColor','none')
        scaleruler('units','km')
        setm(handlem('scaleruler1'),'RulerStyle',"lines",'TickDir','down',"FontSize",12,"LineWidth",1,"FontName","Airal",'MajorTick',0:scale/2:scale,'MinorTick',0:scale/4:scale/2)
        northarrow('latitude',R2.LatitudeLimits(2)-diff(R2.LatitudeLimits)/10,'longitude',R2.LongitudeLimits(2)-diff(R2.LongitudeLimits)/20)
        framem('off')
        set(ax1,"Position",[0.12 0.12 figW figH]);
        set(ax2,"Position",[0.1 0.1 figW figH]);

        print(nstr11,'-dpng','-r600');

    else
        [A,R1] = readgeoraster(dem);
        A=double(A);
        CR=R1;
        figW = 0.9;
        figH = (size(A,1)/size(A,2))*0.9;

        C = zeros(CR.RasterSize);
        C = C- Inf;
        C = imbedm(lat,lon,B,C,CR);
        C(isnan(C))= -Inf;
        C(C==0)= -Inf;

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
        h.CData = C;
        h.FaceColor = 'texturemap';
        colormap(color_rainbow)
        material dull %dull shiny
        camlight
        lighting gouraud
        axis tight
        box off
        set(gca,"FontSize",12,"FontName","Airal");
        c=colorbar(ax1,'Position',[figW-0.057 0.12 0.02 figH/3],"FontSize",12);
        caxis(colrange)
        if strcmpi('insar',obj.tech)==1
            c.Label.String = 'velocity(mm/d)';
        else
            c.Label.String = 'velocity(m/d)';
        end
        title(tstr)

        ax2=subplot(1,2,2);
        ax=worldmap(R1.LatitudeLimits,R1.LongitudeLimits);
        setm(ax,"FontSize",12,"FontName","Airal",'FontColor','none','GColor','none')
        scaleruler('units','km')
        setm(handlem('scaleruler1'),'RulerStyle',"lines",'TickDir','down',"FontSize",12,"LineWidth",1,"FontName","Airal",'MajorTick',0:scale/2:scale,'MinorTick',0:scale/4:scale/2)
        northarrow('latitude',R1.LatitudeLimits(2)-diff(R1.LatitudeLimits)/10,'longitude',R1.LongitudeLimits(2)-diff(R1.LongitudeLimits)/20)
        framem('off')
        set(ax1,"Position",[0.12 0.12 figW figH]);
        set(ax2,"Position",[0.1 0.1 figW figH]);

        print(nstr11,'-dpng','-r600');
    end
    clear B C ind str nstr nstr11 
    close all
end

imds1=imageDatastore('picture',"FileExtensions",".png","IncludeSubfolders",true);
nFrames=length(imds.Files);
mov(1:nFrames) = struct('cdata', [],'colormap', []);
%namell=dir([ '*' '.png']);
set(gcf,'color',[1 1 1]);

for k=1:length(imds.Files)
    Img=imread (imds1.Files{k}); 
    imshow(Img,[]);
    frame=getframe(gcf);
    im=frame2im(frame);
    [I,map]=rgb2ind(Img,256);
    I=uint8(I);
    mov(i)=getframe(gcf);
    if k == 1
        imwrite(I,map,name,'gif', 'Loopcount',Inf,'delaytime',0.2);  
    else
        imwrite(I,map,name,'gif','WriteMode','Append','delaytime',0.2);  
    end
end
rmdir('picture','s')
end