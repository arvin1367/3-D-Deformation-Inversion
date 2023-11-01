classdef myPreprocess
    properties (SetAccess=public)
        files
        geo
        m
        sy
        index
        date
        nan_mask
        sup
        islocal
        angle
        tech
        strsup
    end
    methods
        function obj=myPreprocess(dirc,method,varargin)
            obj.tech=method;
            if ~isempty(varargin)
                obj.angle.a.range=varargin{1};
                obj.angle.a.azimuth=varargin{2};
                obj.angle.d.range=varargin{3};
                obj.angle.d.azimuth=varargin{4};
                obj.islocal=0;
            else
                obj.islocal=1; %flag=0 no local incidence angle ,flag=1 have local incidence angle
            end

            if isunix==1
                obj.sy='/';
            else
                obj.sy='\';
            end
            obj=creat_folder(obj,dirc);
            if strcmpi('POT',obj.tech)==1
                %Define file path
                imds(1).path=imageDatastore(obj.files(1).folder,"FileExtensions",".tif","IncludeSubfolders",true);
                imds(2).path=imageDatastore(obj.files(2).folder,"FileExtensions",".tif","IncludeSubfolders",true);
                imds(3).path=imageDatastore(obj.files(3).folder,"FileExtensions",".tif","IncludeSubfolders",true);
                imds(4).path=imageDatastore(obj.files(4).folder,"FileExtensions",".tif","IncludeSubfolders",true);
                %Check data
                ishasdata(obj,imds(1).path);
                ishasdata(obj,imds(3).path);
                ishasdata(obj,imds(2).path);
                ishasdata(obj,imds(4).path);
            elseif strcmpi('insar',obj.tech)==1
                %Define file path
                imds(1).path=imageDatastore(obj.files(1).folder,"FileExtensions",".tif","IncludeSubfolders",true);
                imds(2).path=imageDatastore(obj.files(2).folder,"FileExtensions",".tif","IncludeSubfolders",true);
                %Check data
                ishasdata(obj,imds(1).path);
                ishasdata(obj,imds(2).path);
            end
            disp("Continue");
            disp('Start reading data:...');

            %load data
            tstart=tic;
            if strcmpi('POT',obj.tech)==1
                [~,obj.geo.R]=readgeoraster(imds(1).path.Files{1});
                obj.geo.info=geotiffinfo(imds(1).path.Files{1});
                ascending_azimuth=load_data(obj,imds(1).path,dirc);
                save('myFile.mat',"ascending_azimuth",'-v7.3');
                clear ascending_azimuth
                obj.m = matfile('myFile.mat','Writable',true);

                ascending_range=load_data(obj,imds(2).path,dirc);
                obj.m.ascending_range=ascending_range;
                clear ascending_range
                descending_azimuth=load_data(obj,imds(3).path,dirc);
                obj.m.descending_azimuth=descending_azimuth;
                clear descending_azimuth
                descending_range=load_data(obj,imds(4).path,dirc);
                obj.m.descending_range=descending_range;
                clear descending_range

                tnow=toc(tstart);
                express=['Data read completed !!!','Completion time:',num2str(tnow),'s'];
                disp(express);
                clear tstart tnow
            elseif strcmpi('insar',obj.tech)==1
                [~,obj.geo.R]=readgeoraster(imds(1).path.Files{1});
                obj.geo.info=geotiffinfo(imds(1).path.Files{1});
                ascending_los=load_data(obj,imds(1).path,dirc);
                save('myFile.mat',"ascending_los",'-v7.3');
                clear ascending_los
                obj.m = matfile('myFile.mat','Writable',true);

                descending_los=load_data(obj,imds(2).path,dirc);
                obj.m.descending_los=descending_los;
                clear descending_los

                tnow=toc(tstart);
                express=['Data read completed !!!','Completion time:',num2str(tnow),'s'];
                disp(express);
                clear tstart tnow
            end

            %mask index
            tstart=tic;
            if strcmpi('POT',obj.tech)==1
                indexaa=extract_mask_zeros(obj,obj.m.ascending_azimuth,imds(1).path,dirc);
                indexar=extract_mask_zeros(obj,obj.m.ascending_range,imds(2).path,dirc);
                indexda=extract_mask_zeros(obj,obj.m.descending_azimuth,imds(3).path,dirc);
                indexdr=extract_mask_zeros(obj,obj.m.descending_range,imds(4).path,dirc);
                obj.index=unique(cat(1,indexaa,indexar,indexda,indexdr));
            elseif strcmpi('insar',obj.tech)==1
                indexa=extract_mask_zeros(obj,obj.m.ascending_los,imds(1).path,dirc);
                indexd=extract_mask_zeros(obj,obj.m.descending_los,imds(2).path,dirc);
                obj.index=unique(cat(1,indexa,indexd));
            end

            %Storage mask tif
            if strcmpi('POT',obj.tech)==1
                temp=obj.m.descending_range(:,:,1);
            elseif strcmpi('insar',obj.tech)==1
                temp=obj.m.descending_los(:,:,1);
            end
            temp(obj.index)=nan;
            obj.nan_mask=temp;
            temp(~isnan(temp))=1;
            temp(isnan(temp))=0;
            geotiffwrite('nan_mask.tif',temp,obj.geo.R,'GeoKeyDirectoryTag',obj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);
            clear temp
            tnow=toc(tstart);
            express=['Extract mask completed !!!','Completion time:',num2str(tnow),'s'];
            disp(express);

            %load date
            disp('Calculate the time:...');
            d1=[dirc,obj.sy,'date_ascending.xls'];d2=[dirc,obj.sy,'date_descending.xls'];
            obj.date.ascending=datetime(readmatrix(d1),"ConvertFrom",'yyyyMMdd');
            obj.date.descending=datetime(readmatrix(d2),"ConvertFrom",'yyyyMMdd');
            obj.date.time=unique([obj.date.ascending(:,1);obj.date.ascending(:,2);obj.date.descending(:,1);obj.date.descending(:,2)]);
            obj.date.delta_t=days(diff(obj.date.time));
            clear d1 d2

            %Import supporting files
            obj.strsup=[dirc,obj.sy,'supporting_documents'];
            if isfolder(obj.strsup)==1
                disp('Import supporting files:...');
                sup=[dirc,obj.sy,'supporting_documents'];
                if strcmpi('POT',obj.tech)==1
                    imds(5).path=imageDatastore(sup,"FileExtensions",".tif","IncludeSubfolders",true);
                    obj=load_sup_data(obj,imds(5).path);
                elseif strcmpi('insar',obj.tech)==1
                    imds(3).path=imageDatastore(sup,"FileExtensions",".tif","IncludeSubfolders",true);
                    obj=load_sup_data(obj,imds(3).path);
                end
            end
            disp('Preprocessing completed !!!!!!!!!!!!!!!');

        end

        function obj=creat_folder(obj,dirc)
            if strcmpi('POT',obj.tech)==1
                files1(1).name='ascending_azimuth';
                files1(1).folder=[dirc,obj.sy,'ascending_azimuth'];
                files1(2).name='ascending_range';
                files1(2).folder=[dirc,obj.sy,'ascending_range'];
                files1(3).name='descending_azimuth';
                files1(3).folder=[dirc,obj.sy,'descending_azimuth'];
                files1(4).name='descending_range';
                files1(4).folder=[dirc,obj.sy,'descending_range'];
                obj.files=files1;
            elseif strcmpi('insar',obj.tech)==1
                files1(1).name='ascending_los';
                files1(1).folder=[dirc,obj.sy,'ascending_los'];
                files1(2).name='descending_los';
                files1(2).folder=[dirc,obj.sy,'descending_los'];
                obj.files=files1;
            end
        end

        function ishasdata(obj,data)
            if hasdata(data)==0
                disp("error: No data in folder");
            end
        end

        function data=load_data(obj,imds,dirc)
            dirc1=[dirc,obj.sy];
            h=waitbar(0,'Start reading data.....');
            [eg,~]=readgeoraster(imds.Files{1});
            data=single(zeros(size(eg,1),size(eg,2),length(imds.Files)));
            k=length(imds.Files);
            for i=1:k
                [tt,~]=readgeoraster(imds.Files{i});
                tt(tt>10000)=0;tt(tt<-10000)=0;
                data(:,:,i)=tt;
                str=['Exporting  ',extractAfter(imds.Folders{1,1},dirc1),',  current progress is ',num2str(i),'/',num2str(k)];
                waitbar(i/k,h,str);
                clear tt
            end
            disp(['Exporting  ',extractAfter(imds.Folders{1,1},dirc1)]);
            close(h)
            clear h i k str eg
        end

        function indexout=extract_mask_zeros(obj,a,imds,dirc)
            h=waitbar(0,'Start extract masking data.....');
            index1=find(a(:,:,1)==0);
            dirc1=[dirc,obj.sy];
            for i=1:size(a,3)
                index2=find(a(:,:,i)==0);
                index3=cat(1,index1,index2);
                index1=unique(index3);
                str=['Extract masking ',extractAfter(imds.Folders{1,1},dirc1),',  current progress is ',num2str(i),'/',num2str(size(a,3))];
                waitbar(i/size(a,3),h,str);
                clear index2 index3
            end
            indexout=index1;
            close(h)
            disp(['Exporting  ',extractAfter(imds.Folders{1,1},dirc1)]);
            clear i index1 a
        end
     
        function obj=load_sup_data(obj,imds)
            str1=['supporting_documents',obj.sy];
            str2='.tif';
            if obj.islocal==1
                for i=1:length(imds.Files)
                    str_name=extractAfter(imds.Files{i},str1);
                    str_name=extractBefore(str_name,str2);
                    if strcmp(str_name,'dem_ddew')==1
                        obj.sup.dem_ddew=readgeoraster(imds.Files{i});
                    elseif strcmp(str_name,'dem_ddns')==1
                        obj.sup.dem_ddns=readgeoraster(imds.Files{i});
                    elseif strcmp(str_name,'lv_phi_a')==1
                        obj.sup.lv_phi_a=readgeoraster(imds.Files{i});
                    elseif strcmp(str_name,'lv_phi_d')==1
                        obj.sup.lv_phi_d=readgeoraster(imds.Files{i});
                    elseif strcmp(str_name,'lv_theta_a')==1
                        obj.sup.lv_theta_a=readgeoraster(imds.Files{i});
                    else
                        obj.sup.lv_theta_d=readgeoraster(imds.Files{i});
                    end
                    clear str_name
                end
            else
                for i=1:length(imds.Files)
                    str_name=extractAfter(imds.Files{i},str1);
                    str_name=extractBefore(str_name,str2);
                    if strcmp(str_name,'dem_ddew')==1
                        obj.sup.dem_ddew=readgeoraster(imds.Files{i});
                    elseif strcmp(str_name,'dem_ddns')==1
                        obj.sup.dem_ddns=readgeoraster(imds.Files{i});
                    end
                    clear str_name
                end
            end
        end

    end
end