classdef myReprocess
    properties (SetAccess=public)
        sta % Uncertainty statistical indicators
        rect % Index of stable region
    end

    methods
        function obj=myReprocess(obj1,obj2)
            if strcmpi('POT',obj1.tech)==1
                disp('Please determine the stable area yourself in ArcGIS and draw a rectangle based on the provided image:')
                load mkfigure/color_rainbow.mat
                subplot(1,2,1)
                h1=imagesc(obj1.nan_mask);
                subplot(1,2,2)
                h2=imagesc(obj1.m.ascending_azimuth(:,:,1));
                colormap(color_rainbow)
                [~,obj.rect]=imcrop(h1);
                disp('Please close the image to continue......');
                obj.rect=round(obj.rect);
                close all

                disp('Generating statistical chart:');
                [obj.sta.std.ascending_azimuth.each,obj.sta.std.ascending_azimuth.ave]=compute_uncertainty_from_SAR(obj,obj1.m.ascending_azimuth,obj2.aspect.a,'Ascending azimuth',obj1.date.ascending,obj1);
                [obj.sta.std.ascending_range.each,obj.sta.std.ascending_range.ave]=compute_uncertainty_from_SAR(obj,obj1.m.ascending_range,obj2.aspect.a,'Ascending range',obj1.date.ascending,obj1);
                [obj.sta.std.descending_azimuth.each,obj.sta.std.descending_azimuth.ave]=compute_uncertainty_from_SAR(obj,obj1.m.descending_azimuth,obj2.aspect.d,'Descending azimuth',obj1.date.descending,obj1);
                [obj.sta.std.descending_range.each,obj.sta.std.descending_range.ave]=compute_uncertainty_from_SAR(obj,obj1.m.descending_range,obj2.aspect.d,'Descending range',obj1.date.descending,obj1);

            elseif strcmpi('insar',obj1.tech)==1
                load mkfigure/color_rainbow.mat
                disp('Please determine the stable area yourself in ArcGIS and draw a rectangle based on the provided image:')
                subplot(1,2,1)
                h1=imagesc(obj1.nan_mask);
                subplot(1,2,2)
                h2=imagesc(obj1.m.ascending_los(:,:,1));
                colormap(color_rainbow)
                [~,obj.rect]=imcrop(h1);
                disp('Please close the image to continue......');
                obj.rect=round(obj.rect);
                close all

                disp('Generating statistical chart:');
                [obj.sta.std.ascending_los.each,obj.sta.std.ascending_los.ave]=compute_uncertainty_from_SAR(obj,obj1.m.ascending_los,obj2.aspect.a,'Ascending los',obj1.date.ascending,obj1);
                [obj.sta.std.descending_los.each,obj.sta.std.descending_los.ave]=compute_uncertainty_from_SAR(obj,obj1.m.descending_los,obj2.aspect.d,'Descending los',obj1.date.descending,obj1);

            end

            obj=compute_uncertainty_from_3D(obj,obj1);

        end

        function [temp_std,ave_std]=compute_uncertainty_from_SAR(obj,data,len,str,date,obj1)
            temp=zeros(obj.rect(4),obj.rect(3));
            temp_std=ones(len,1);
            dd=days(date(:,2)-date(:,1));
            for k=1:len
                for i=1:obj.rect(4)
                    for j=1:obj.rect(3)
                        temp(i,j)=100*data(i+obj.rect(2),j+obj.rect(1),k)/dd(k);
                    end
                end
                stdk=std(temp,0,'all');
                temp_std(k)=stdk;
            end
            ave_std=mean(temp_std);
            if strcmpi('POT',obj1.tech)==1
                up=ceil(max(temp_std));
                down=floor(min(temp_std));
                x=down:0.1:up;
                for g=1:length(x)-1
                    y(g)=(length(find(temp_std>x(g) & temp_std<x(g+1)))/length(temp_std))*100;
                end
                xy=(down+0.05):0.1:(up-0.05);
            else
                up=max(temp_std);
                down=min(temp_std);
                x=down:0.001:up;
                for g=1:length(x)-1
                    y(g)=(length(find(temp_std>x(g) & temp_std<x(g+1)))/length(temp_std))*100;
                    xy(g)=x(g)+0.0005;
                end
            end
          
            figure
            bar(xy,y);
            set(gca, 'FontSize',12,'linewidth',1.5);
            xlabel('Standard deviation(cm/day)');
            ylabel('Precentage(%)');
            title([str,'   Mean  ',num2str(ave_std)]);
            clear up down i j k g xy x y dd
        end

        function obj=compute_uncertainty_from_3D(obj,obj1)

            %error propagation law y=k1*x1+k2*x2, dy=sqrt((k1*dx1)^2+(k2*dx2)^2)
            if obj1.islocal==1
                phia=mean(obj1.sup.lv_phi_a,'all');
                phid=mean(obj1.sup.lv_phi_d,'all');
                thetaa=mean(obj1.sup.lv_theta_a,'all');
                thetad=mean(obj1.sup.lv_theta_d,'all');

                if strcmpi('POT',obj1.tech)==1
                    A=zeros(4,3);

                    A(1,1)=-cos(phia);
                    A(1,2)=sin(phia);
                    A(1,3)=0;

                    A(2,1)=sin(phia)*cos(thetaa);
                    A(2,2)=cos(phia)*cos(thetaa);
                    A(2,3)=sin(thetaa);

                    A(3,1)=-cos(phid);
                    A(3,2)=sin(phid);
                    A(3,3)=0;

                    A(4,1)=sin(phid)*cos(thetad);
                    A(4,2)=cos(phid)*cos(thetad);
                    A(4,3)=sin(thetad);

                    AA=pinv(A);
                    obj.sta.uncertainty.vn=sqrt((AA(1,1)*obj.sta.std.ascending_azimuth.ave)^2+ ...
                        (AA(1,2)*obj.sta.std.ascending_range.ave)^2+(AA(1,3)*obj.sta.std.descending_azimuth.ave)^2+ ...
                        (AA(1,4)*obj.sta.std.descending_range.ave)^2);
                    obj.sta.uncertainty.ve=sqrt((AA(2,1)*obj.sta.std.ascending_azimuth.ave)^2+ ...
                        (AA(2,2)*obj.sta.std.ascending_range.ave)^2+(AA(2,3)*obj.sta.std.descending_azimuth.ave)^2+ ...
                        (AA(2,4)*obj.sta.std.descending_range.ave)^2);
                    obj.sta.uncertainty.vu=sqrt((AA(3,1)*obj.sta.std.ascending_azimuth.ave)^2+ ...
                        (AA(3,2)*obj.sta.std.ascending_range.ave)^2+(AA(3,3)*obj.sta.std.descending_azimuth.ave)^2+ ...
                        (AA(3,4)*obj.sta.std.descending_range.ave)^2);

                    clear A AA phid phia thetaa thetad
                elseif strcmpi('insar',obj1.tech)==1
                    A=zeros(2,3);

                    A(1,1)=sin(phia)*cos(thetaa);
                    A(1,2)=cos(phia)*cos(thetaa);
                    A(1,3)=sin(thetaa);

                    A(2,1)=sin(phid)*cos(thetad);
                    A(2,2)=cos(phid)*cos(thetad);
                    A(2,3)=sin(thetad);

                    AA=pinv(A);
                    obj.sta.uncertainty.vn=sqrt((AA(1,1)*obj.sta.std.ascending_los.ave)^2+ ...
                        (AA(1,2)*obj.sta.std.descending_los.ave)^2);
                    obj.sta.uncertainty.ve=sqrt((AA(2,1)*obj.sta.std.ascending_los.ave)^2+ ...
                        (AA(2,2)*obj.sta.std.descending_los.ave)^2);
                    obj.sta.uncertainty.vu=sqrt((AA(3,1)*obj.sta.std.ascending_los.ave)^2+ ...
                        (AA(3,2)*obj.sta.std.descending_los.ave)^2);

                    clear A AA phid phia thetaa thetad
                end

            else
                if strcmpi('POT',obj1.tech)==1
                    A=zeros(4,3);

                    A(1,1)=cosd(obj1.angle.a.azimuth);
                    A(1,2)=sind(obj1.angle.a.azimuth);
                    A(1,3)=0;

                    A(2,1)=sind(obj1.angle.a.range)*sind(obj1.angle.a.azimuth);
                    A(2,2)=-sind(obj1.angle.a.range)*cosd(obj1.angle.a.azimuth);
                    A(2,3)=cosd(obj1.angle.a.range);

                    A(3,1)=cosd(obj1.angle.d.azimuth);
                    A(3,2)=sind(obj1.angle.d.azimuth);
                    A(3,3)=0;

                    A(4,1)=sind(obj1.angle.d.range)*sind(obj1.angle.d.azimuth);
                    A(4,2)=-sind(obj1.angle.d.range)*cosd(obj1.angle.d.azimuth);
                    A(4,3)=cosd(obj1.angle.d.range);

                    AA=pinv(A);
                    obj.sta.uncertainty.vn=sqrt((AA(1,1)*obj.sta.std.ascending_azimuth.ave)^2+ ...
                        (AA(1,2)*obj.sta.std.ascending_range.ave)^2+(AA(1,3)*obj.sta.std.descending_azimuth.ave)^2+ ...
                        (AA(1,4)*obj.sta.std.descending_range.ave)^2);
                    obj.sta.uncertainty.ve=sqrt((AA(2,1)*obj.sta.std.ascending_azimuth.ave)^2+ ...
                        (AA(2,2)*obj.sta.std.ascending_range.ave)^2+(AA(2,3)*obj.sta.std.descending_azimuth.ave)^2+ ...
                        (AA(2,4)*obj.sta.std.descending_range.ave)^2);
                    obj.sta.uncertainty.vu=sqrt((AA(3,1)*obj.sta.std.ascending_azimuth.ave)^2+ ...
                        (AA(3,2)*obj.sta.std.ascending_range.ave)^2+(AA(3,3)*obj.sta.std.descending_azimuth.ave)^2+ ...
                        (AA(3,4)*obj.sta.std.descending_range.ave)^2);

                    clear A AA

                elseif strcmpi('insar',obj1.tech)==1
                    A=zeros(2,3);

                    A(1,1)=sind(obj1.angle.a.range)*sind(obj1.angle.a.azimuth);
                    A(1,2)=-sind(obj1.angle.a.range)*cosd(obj1.angle.a.azimuth);
                    A(1,3)=cosd(obj1.angle.a.range);

                    A(2,1)=sind(obj1.angle.d.range)*sind(obj1.angle.d.azimuth);
                    A(2,2)=-sind(obj1.angle.d.range)*cosd(obj1.angle.d.azimuth);
                    A(2,3)=cosd(obj1.angle.d.range);

                    AA=pinv(A);
                    obj.sta.uncertainty.vn=sqrt((AA(1,1)*obj.sta.std.ascending_los.ave)^2+ ...
                        (AA(1,2)*obj.sta.std.descending_los.ave)^2);
                    obj.sta.uncertainty.ve=sqrt((AA(2,1)*obj.sta.std.ascending_los.ave)^2+ ...
                        (AA(2,2)*obj.sta.std.descending_los.ave)^2);
                    obj.sta.uncertainty.vu=sqrt((AA(3,1)*obj.sta.std.ascending_los.ave)^2+ ...
                        (AA(3,2)*obj.sta.std.descending_los.ave)^2);

                    clear A AA

                end
            end
        end

    end
end