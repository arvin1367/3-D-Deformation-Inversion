classdef myProcess
    properties (SetAccess=public)
        aspect
        m2
        lamda
        normv
        flag
        method
    end

    methods
        function obj=myProcess(preobj,varargin)
            %Compute
            disp('Bulid storage output:...')
            if ~isempty(varargin)
                if strcmp(varargin{1},'Tikhonov')==1
                    obj.method=0;
                    obj.lamda=varargin{2};
                    obj.flag=varargin{3};
                elseif strcmp(varargin{1},'TSVD')==1
                    obj.method=1;
                    obj.lamda=varargin{2};
                else
                    error('Input parameter error!!!');
                end
            end
            
            obj.aspect.length=size(preobj.nan_mask,1);
            obj.aspect.width=size(preobj.nan_mask,2);
            obj.aspect.phase=size(preobj.date.delta_t,1);
            obj.aspect.a=size(preobj.date.ascending,1);
            obj.aspect.d=size(preobj.date.descending,1);
            if isfolder(preobj.strsup)==1
                obj.aspect.spf=length(preobj.date.delta_t);
            else
                obj.aspect.spf=0;
            end
            obj.aspect.tikhonov=(length(preobj.date.delta_t)-1)*3;
            if strcmpi('POT',preobj.tech)==1
                if obj.method==0
                    if obj.flag==0
                        obj.aspect.pairs=obj.aspect.a*2+obj.aspect.d*2+obj.aspect.tikhonov+obj.aspect.spf+3;
                    elseif obj.flag==1
                        obj.aspect.pairs=obj.aspect.a*2+obj.aspect.d*2+obj.aspect.tikhonov+obj.aspect.spf;
                    elseif obj.flag==2
                        obj.aspect.pairs=obj.aspect.a*2+obj.aspect.d*2+obj.aspect.tikhonov+obj.aspect.spf-3;
                    end
                elseif obj.method==1
                    obj.aspect.pairs=obj.aspect.a*2+obj.aspect.d*2+obj.aspect.spf;
                end
            elseif strcmpi('insar',preobj.tech)==1
                if obj.method==0
                    obj.aspect.pairs=obj.aspect.a+obj.aspect.d+obj.aspect.tikhonov+obj.aspect.spf;
                    if obj.flag==0
                        obj.aspect.pairs=obj.aspect.pairs+3;
                    elseif obj.flag==2
                        obj.aspect.pairs=obj.aspect.pairs-3;
                    end
                elseif obj.method==1
                    obj.aspect.pairs=obj.aspect.a+obj.aspect.d+obj.aspect.spf;
                end
            end

            ve=single(zeros(obj.aspect.length,obj.aspect.width,obj.aspect.phase));
            vn=single(zeros(obj.aspect.length,obj.aspect.width,obj.aspect.phase));
            vu=single(zeros(obj.aspect.length,obj.aspect.width,obj.aspect.phase));
            if strcmpi('POT',preobj.tech)==1
                a.aa=preobj.m.ascending_azimuth;a.ar=preobj.m.ascending_range;
                a.da=preobj.m.descending_azimuth;a.dr=preobj.m.descending_range;
            elseif strcmpi('insar',preobj.tech)==1
                a.a=preobj.m.ascending_los;a.d=preobj.m.descending_los;
            end

            tstart=tic;
            disp('Starting 3D calculation may take a long time, please be patient and wait:...');
            count=1;all_count=obj.aspect.length;lable=0;str='Completed  0%';
            h=waitbar(0,'Please be patient and wait.....');
            waitbar(0,h,str);
            norm_x=0;countt=0;norm_axy=0;

            for i=1:obj.aspect.length
                for j=1:obj.aspect.width
                    if ~isnan(preobj.nan_mask(i,j)) %Determine whether it is a null value
                        if obj.method==0
                            b=make_b(obj,preobj,a,i,j); %b is deformation data pairs
                            A=make_A(obj,preobj,i,j); %A is time coefficient matrix
                            temp=pinv(A)*b; %SVD
                        elseif obj.method==1
                            b=make_b(obj,preobj,a,i,j); %b is deformation data pairs
                            A=make_A(obj,preobj,i,j); %A is time coefficient matrix
                            [U,S,V]=svd(A);
                            S(S<obj.lamda)=0;
                            temp=pinv(U*S*V')*b;
                        end
                        for k=1:obj.aspect.phase
                            vn(i,j,k)=temp((k-1)*3+1);
                            ve(i,j,k)=temp((k-1)*3+2);
                            vu(i,j,k)=temp((k-1)*3+3);
                        end
                        norm_t=squeeze(cat(3,vn(i,j,:),ve(i,j,:),vu(i,j,:)));
                        if  norm(norm_t)<10000
                            norm_x=norm_x+norm(norm_t);
                            norm_axy=norm_axy+norm(A*temp);
                            countt=countt+1;
                        end
                        clear temp A b k norm_t
                    end
                end
                count_a=100*count/all_count;
                count_b=floor(count_a);
                count_c=count_a-count_b;
                if lable>count_c
                    str=['Completed  ',num2str(count_b),'%'];
                    disp(str);
                    waitbar(i/all_count,h,str);
                end
                lable=count_c;
                count=count+1;
            end
            
            obj.normv.x=norm_x/countt;
            obj.normv.axy=norm_axy/countt;
            disp(['Computed ||x|| and ||Ax-Y|| norms: ',num2str(obj.normv.x),'    ',num2str(obj.normv.axy)]);
            close(h)
            clear h lable count_c count_a all_count count_b str count i j q a

            disp('Storing result as mat........myFile_out-----velocity');
            for k=1:obj.aspect.phase
                ve(:,:,k)=del_std(obj,ve(:,:,k));
                vn(:,:,k)=del_std(obj,vn(:,:,k));
                vu(:,:,k)=del_std(obj,vu(:,:,k));
            end
            save('myFile_out.mat',"ve",'-v7.3');
            obj.m2 = matfile('myFile_out.mat','Writable',true);
            obj.m2.vn=vn;
            obj.m2.vu=vu;

            tnow=toc(tstart);
            disp(['Completion time:',num2str(tnow/3600),'h']);
            clear tstart tnow

            %Create cumulative time series
            disp('Start building cumulative time series output.....');
            tstart=tic;
            ve_cm=single(zeros(obj.aspect.length,obj.aspect.width,obj.aspect.phase));
            vn_cm=single(zeros(obj.aspect.length,obj.aspect.width,obj.aspect.phase));
            vu_cm=single(zeros(obj.aspect.length,obj.aspect.width,obj.aspect.phase));
            for k=1:obj.aspect.phase
                if k==1
                    ve_cm(:,:,k)=ve(:,:,k).*preobj.date.delta_t(k);
                    ve_cm(:,:,k)=del_std(obj,ve_cm(:,:,k));
                    vn_cm(:,:,k)=vn(:,:,k).*preobj.date.delta_t(k);
                    vn_cm(:,:,k)=del_std(obj,vn_cm(:,:,k));
                    vu_cm(:,:,k)=vu(:,:,k).*preobj.date.delta_t(k);
                    vu_cm(:,:,k)=del_std(obj,vu_cm(:,:,k));
                else
                    ve_cm(:,:,k)=ve(:,:,k).*preobj.date.delta_t(k)+ve_cm(:,:,k-1);
                    ve_cm(:,:,k)=del_std(obj,ve_cm(:,:,k));
                    vn_cm(:,:,k)=vn(:,:,k).*preobj.date.delta_t(k)+vn_cm(:,:,k-1);
                    vn_cm(:,:,k)=del_std(obj,vn_cm(:,:,k));
                    vu_cm(:,:,k)=vu(:,:,k).*preobj.date.delta_t(k)+vu_cm(:,:,k-1);
                    vu_cm(:,:,k)=del_std(obj,vu_cm(:,:,k));
                end
            end
            obj.m2.ve_cm=ve_cm;
            obj.m2.vn_cm=vn_cm;
            obj.m2.vu_cm=vu_cm;

            ve_cm_v=(ve_cm(:,:,end)./days(preobj.date.time(end)-preobj.date.time(1))).*365;
            ve_cm_v=del_std(obj,ve_cm_v);
            vn_cm_v=(vn_cm(:,:,end)./days(preobj.date.time(end)-preobj.date.time(1))).*365;
            vn_cm_v=del_std(obj,vn_cm_v);
            vu_cm_v=(vu_cm(:,:,end)./days(preobj.date.time(end)-preobj.date.time(1))).*365;
            vu_cm_v=del_std(obj,vu_cm_v);
            obj.m2.ve_cm_v=ve_cm_v;
            obj.m2.vn_cm_v=vn_cm_v;
            obj.m2.vu_cm_v=vu_cm_v;
            tnow=toc(tstart);
            disp(['completion time: ',num2str(tnow),'s']);
            clear tstart tnow k

            %Export as tif
            disp('Export as tif write to disk.....');
            tstart=tic;
            if isunix==1
                mkdir results_velocity/ve
                mkdir results_velocity/vn
                mkdir results_velocity/vu
                mkdir results_cumulative/ve
                mkdir results_cumulative/vn
                mkdir results_cumulative/vu
            else
                mkdir results_velocity\ve
                mkdir results_velocity\vn
                mkdir results_velocity\vu
                mkdir results_cumulative\ve
                mkdir results_cumulative\vn
                mkdir results_cumulative\vu
            end
            for k=1:obj.aspect.phase
                str_ve_v=['results_velocity',preobj.sy,'ve',preobj.sy,datestr(preobj.date.time(k),'yyyymmdd'),'_',datestr(preobj.date.time(k+1),'yyyymmdd'),'.tif'];
                geotiffwrite(str_ve_v,ve(:,:,k),preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);

                str_vn_v=['results_velocity',preobj.sy,'vn',preobj.sy,datestr(preobj.date.time(k),'yyyymmdd'),'_',datestr(preobj.date.time(k+1),'yyyymmdd'),'.tif'];
                geotiffwrite(str_vn_v,vn(:,:,k),preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);

                str_vu_v=['results_velocity',preobj.sy,'vu',preobj.sy,datestr(preobj.date.time(k),'yyyymmdd'),'_',datestr(preobj.date.time(k+1),'yyyymmdd'),'.tif'];
                geotiffwrite(str_vu_v,vu(:,:,k),preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);

                str_ve_cm=['results_cumulative',preobj.sy,'ve',preobj.sy,datestr(preobj.date.time(k),'yyyymmdd'),'_',datestr(preobj.date.time(k+1),'yyyymmdd'),'.tif'];
                geotiffwrite(str_ve_cm,ve_cm(:,:,k),preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);

                str_vn_cm=['results_cumulative',preobj.sy,'vn',preobj.sy,datestr(preobj.date.time(k),'yyyymmdd'),'_',datestr(preobj.date.time(k+1),'yyyymmdd'),'.tif'];
                geotiffwrite(str_vn_cm,vn_cm(:,:,k),preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);

                str_vu_cm=['results_cumulative',preobj.sy,'vu',preobj.sy,datestr(preobj.date.time(k),'yyyymmdd'),'_',datestr(preobj.date.time(k+1),'yyyymmdd'),'.tif'];
                geotiffwrite(str_vu_cm,vu_cm(:,:,k),preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);
            end
            geotiffwrite('ve_ave_year.tif',ve_cm_v,preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);
            geotiffwrite('vn_ave_year.tif',vn_cm_v,preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);
            geotiffwrite('vu_ave_year.tif',vu_cm_v,preobj.geo.R,'GeoKeyDirectoryTag',preobj.geo.info.GeoTIFFTags.GeoKeyDirectoryTag);
            tnow=toc(tstart);
            disp(['Completion time:',num2str(tnow),'s']);
            clear ve vn vu vu_cm_v ve_cm_v vn_cm_v vu_cm  vn_cm ve_cm str tstart tnow str_ve_v str_vn_v str_vu_v str_ve_cm str_vn_cm str_vu_c
        end

        function A=make_A(obj,preobj,i,j)
            %Cteate time coefficient matrix
            A=zeros(obj.aspect.pairs,length(preobj.date.delta_t)*3);
            time=preobj.date.time;
            delta_t=preobj.date.delta_t;

            if isfolder(preobj.strsup)==1
                ddew=preobj.sup.dem_ddew;
                ddns=preobj.sup.dem_ddns;
            end

            if preobj.islocal==1
                phia=preobj.sup.lv_phi_a;
                phid=preobj.sup.lv_phi_d;
                thetaa=preobj.sup.lv_theta_a;
                thetad=preobj.sup.lv_theta_d;
            end

            if strcmpi('POT',preobj.tech)==1
                %1:a_a
                for k=1:obj.aspect.a
                    index1=find(time==preobj.date.ascending(k,1));
                    index2=find(time==preobj.date.ascending(k,2))-1;
                    if preobj.islocal==1
                        A(k,3*(index1-1)+1)=-delta_t(index1)*cos(phia(i,j));
                        A(k,3*(index1-1)+2)=delta_t(index1)*sin(phia(i,j));
                        A(k,3*(index1-1)+3)=0;

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=-delta_t(index1)*cos(phia(i,j));
                            A(k,3*(index2-1)+2)=delta_t(index1)*sin(phia(i,j));
                            A(k,3*(index2-1)+3)=0;
                            index2=index2-1;
                        end
                    else
                        A(k,3*(index1-1)+1)=delta_t(index1)*cosd(preobj.angle.a.azimuth);
                        A(k,3*(index1-1)+2)=delta_t(index1)*sind(preobj.angle.a.azimuth);
                        A(k,3*(index1-1)+3)=0;

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=delta_t(index1)*cosd(preobj.angle.a.azimuth);
                            A(k,3*(index2-1)+2)=delta_t(index1)*sind(preobj.angle.a.azimuth);
                            A(k,3*(index2-1)+3)=0;
                            index2=index2-1;
                        end
                    end
                    clear index1 index2
                end

                %a_a:a_r
                for k=obj.aspect.a+1:2*obj.aspect.a
                    index1=find(time==preobj.date.ascending(k-obj.aspect.a,1));
                    index2=find(time==preobj.date.ascending(k-obj.aspect.a,2))-1;
                    if preobj.islocal==1
                        A(k,3*(index1-1)+1)=delta_t(index1)*sin(phia(i,j))*cos(thetaa(i,j));
                        A(k,3*(index1-1)+2)=delta_t(index1)*cos(phia(i,j))*cos(thetaa(i,j));
                        A(k,3*(index1-1)+3)=delta_t(index1)*sin(thetaa(i,j));

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=delta_t(index1)*sin(phia(i,j))*cos(thetaa(i,j));
                            A(k,3*(index2-1)+2)=delta_t(index1)*cos(phia(i,j))*cos(thetaa(i,j));
                            A(k,3*(index2-1)+3)=delta_t(index1)*sin(thetaa(i,j));
                            index2=index2-1;
                        end
                    else
                        A(k,3*(index1-1)+1)=delta_t(index1)*sind(preobj.angle.a.range)*sind(preobj.angle.a.azimuth);
                        A(k,3*(index1-1)+2)=-delta_t(index1)*sind(preobj.angle.a.range)*cosd(preobj.angle.a.azimuth);
                        A(k,3*(index1-1)+3)=delta_t(index1)*cosd(preobj.angle.a.range);

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=delta_t(index1)*sind(preobj.angle.a.range)*sind(preobj.angle.a.azimuth);
                            A(k,3*(index2-1)+2)=-delta_t(index1)*sind(preobj.angle.a.range)*cosd(preobj.angle.a.azimuth);
                            A(k,3*(index2-1)+3)=delta_t(index1)*cosd(preobj.angle.a.range);
                            index2=index2-1;
                        end
                    end
                    clear index1 index2
                end

                %a_r:d_a
                for k=obj.aspect.a*2+1:obj.aspect.d+obj.aspect.a*2
                    index1=find(time==preobj.date.descending(k-obj.aspect.a*2,1));
                    index2=find(time==preobj.date.descending(k-obj.aspect.a*2,2))-1;
                    if preobj.islocal==1
                        A(k,3*(index1-1)+1)=-delta_t(index1)*cos(phid(i,j));
                        A(k,3*(index1-1)+2)=delta_t(index1)*sin(phid(i,j));
                        A(k,3*(index1-1)+3)=0;

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=-delta_t(index1)*cos(phid(i,j));
                            A(k,3*(index2-1)+2)=delta_t(index1)*sin(phid(i,j));
                            A(k,3*(index2-1)+3)=0;
                            index2=index2-1;
                        end
                    else
                        A(k,3*(index1-1)+1)=delta_t(index1)*cosd(preobj.angle.d.azimuth);
                        A(k,3*(index1-1)+2)=delta_t(index1)*sind(preobj.angle.d.azimuth);
                        A(k,3*(index1-1)+3)=0;

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=delta_t(index1)*cosd(preobj.angle.d.azimuth);
                            A(k,3*(index2-1)+2)=delta_t(index1)*sind(preobj.angle.d.azimuth);
                            A(k,3*(index2-1)+3)=0;
                            index2=index2-1;
                        end
                    end
                    clear index1 index2
                end

                %d_a:d_r
                for k=obj.aspect.a*2+1+obj.aspect.d:obj.aspect.d*2+obj.aspect.a*2
                    index1=find(time==preobj.date.descending(k-obj.aspect.a*2-obj.aspect.d,1));
                    index2=find(time==preobj.date.descending(k-obj.aspect.a*2-obj.aspect.d,2))-1;
                    if preobj.islocal==1
                        A(k,3*(index1-1)+1)=delta_t(index1)*sin(phid(i,j))*cos(thetad(i,j));
                        A(k,3*(index1-1)+2)=delta_t(index1)*cos(phid(i,j))*cos(thetad(i,j));
                        A(k,3*(index1-1)+3)=delta_t(index1)*sin(thetad(i,j));

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=delta_t(index1)*sin(phid(i,j))*cos(thetad(i,j));
                            A(k,3*(index2-1)+2)=delta_t(index1)*cos(phid(i,j))*cos(thetad(i,j));
                            A(k,3*(index2-1)+3)=delta_t(index1)*sin(thetad(i,j));
                            index2=index2-1;
                        end
                    else
                        A(k,3*(index1-1)+1)=delta_t(index1)*sind(preobj.angle.d.range)*sind(preobj.angle.d.azimuth);
                        A(k,3*(index1-1)+2)=-delta_t(index1)*sind(preobj.angle.d.range)*cosd(preobj.angle.d.azimuth);
                        A(k,3*(index1-1)+3)=delta_t(index1)*cosd(preobj.angle.d.range);

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=delta_t(index1)*sind(preobj.angle.d.range)*sind(preobj.angle.d.azimuth);
                            A(k,3*(index2-1)+2)=-delta_t(index1)*sind(preobj.angle.d.range)*cosd(preobj.angle.d.azimuth);
                            A(k,3*(index2-1)+3)=delta_t(index1)*cosd(preobj.angle.d.range);
                            index2=index2-1;
                        end
                    end
                    clear index1 index2
                end

                %d_r:SPF
                if isfolder(preobj.strsup)==1
                    for k=obj.aspect.a*2+1+obj.aspect.d*2:obj.aspect.d*2+obj.aspect.a*2+obj.aspect.spf
                        q=k-obj.aspect.a*2-obj.aspect.d*2;
                        A(k,3*(q-1)+1)=ddew(i,j);
                        A(k,3*(q-1)+2)=ddns(i,j);
                        A(k,3*(q-1)+3)=-1;
                        clear q
                    end
                end

                if obj.method==0
                    %SPF:Tikhonov
                    if obj.flag==0
                        for k=obj.aspect.a*2+1+obj.aspect.d*2+obj.aspect.spf:obj.aspect.pairs
                            q=k-obj.aspect.a*2-obj.aspect.d*2-obj.aspect.spf;
                            A(k,q)=obj.lamda;
                            clear q
                        end
                    elseif obj.flag==1
                        for k=obj.aspect.a*2+1+obj.aspect.d*2+obj.aspect.spf:obj.aspect.pairs
                            q=k-obj.aspect.a*2-obj.aspect.d*2-obj.aspect.spf;
                            A(k,q)=-obj.lamda;
                            A(k,q+3)=obj.lamda;
                            clear q
                        end
                    elseif obj.flag==2
                        for k=obj.aspect.a*2+1+obj.aspect.d*2+obj.aspect.spf:obj.aspect.pairs
                            q=k-obj.aspect.a*2-obj.aspect.d*2-obj.aspect.spf;
                            A(k,q)=obj.lamda;
                            A(k,q+3)=-2*obj.lamda;
                            A(k,q+6)=obj.lamda;
                            clear q
                        end
                    end
                end

                clear k time delta_t ddew ddns phia phid thetaa thetad

            elseif strcmpi('insar',preobj.tech)==1
                %1:a
                for k=1:obj.aspect.a
                    index1=find(time==preobj.date.ascending(k,1));
                    index2=find(time==preobj.date.ascending(k,2))-1;
                    if preobj.islocal==1
                        A(k,3*(index1-1)+1)=delta_t(index1)*sin(phia(i,j))*cos(thetaa(i,j));
                        A(k,3*(index1-1)+2)=delta_t(index1)*cos(phia(i,j))*cos(thetaa(i,j));
                        A(k,3*(index1-1)+3)=delta_t(index1)*sin(thetaa(i,j));

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=delta_t(index1)*sin(phia(i,j))*cos(thetaa(i,j));
                            A(k,3*(index2-1)+2)=delta_t(index1)*cos(phia(i,j))*cos(thetaa(i,j));
                            A(k,3*(index2-1)+3)=delta_t(index1)*sin(thetaa(i,j));
                            index2=index2-1;
                        end
                    else
                        A(k,3*(index1-1)+1)=delta_t(index1)*sind(preobj.angle.a.range)*sind(preobj.angle.a.azimuth);
                        A(k,3*(index1-1)+2)=-delta_t(index1)*sind(preobj.angle.a.range)*cosd(preobj.angle.a.azimuth);
                        A(k,3*(index1-1)+3)=delta_t(index1)*cosd(preobj.angle.a.range);

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=delta_t(index1)*sind(preobj.angle.a.range)*sind(preobj.angle.a.azimuth);
                            A(k,3*(index2-1)+2)=-delta_t(index1)*sind(preobj.angle.a.range)*cosd(preobj.angle.a.azimuth);
                            A(k,3*(index2-1)+3)=delta_t(index1)*cosd(preobj.angle.a.range);
                            index2=index2-1;
                        end
                    end
                    clear index1 index2
                end

                %a:d
                for k=obj.aspect.a+1:obj.aspect.d+obj.aspect.a
                    index1=find(time==preobj.date.descending(k-obj.aspect.a,1));
                    index2=find(time==preobj.date.descending(k-obj.aspect.a,2))-1;
                    if preobj.islocal==1
                        A(k,3*(index1-1)+1)=delta_t(index1)*sin(phid(i,j))*cos(thetad(i,j));
                        A(k,3*(index1-1)+2)=delta_t(index1)*cos(phid(i,j))*cos(thetad(i,j));
                        A(k,3*(index1-1)+3)=delta_t(index1)*sin(thetad(i,j));

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=delta_t(index1)*sin(phid(i,j))*cos(thetad(i,j));
                            A(k,3*(index2-1)+2)=delta_t(index1)*cos(phid(i,j))*cos(thetad(i,j));
                            A(k,3*(index2-1)+3)=delta_t(index1)*sin(thetad(i,j));
                            index2=index2-1;
                        end
                    else
                        A(k,3*(index1-1)+1)=delta_t(index1)*sind(preobj.angle.d.range)*sind(preobj.angle.d.azimuth);
                        A(k,3*(index1-1)+2)=-delta_t(index1)*sind(preobj.angle.d.range)*cosd(preobj.angle.d.azimuth);
                        A(k,3*(index1-1)+3)=delta_t(index1)*cosd(preobj.angle.d.range);

                        while (index2-index1)~=0
                            A(k,3*(index2-1)+1)=delta_t(index1)*sind(preobj.angle.d.range)*sind(preobj.angle.d.azimuth);
                            A(k,3*(index2-1)+2)=-delta_t(index1)*sind(preobj.angle.d.range)*cosd(preobj.angle.d.azimuth);
                            A(k,3*(index2-1)+3)=delta_t(index1)*cosd(preobj.angle.d.range);
                            index2=index2-1;
                        end
                    end
                    clear index1 index2
                end

                %d:SPF
                if isfolder(preobj.strsup)==1
                    for k=obj.aspect.a+1+obj.aspect.d:obj.aspect.d+obj.aspect.a+obj.aspect.spf
                        q=k-obj.aspect.a-obj.aspect.d;
                        A(k,3*(q-1)+1)=ddew(i,j);
                        A(k,3*(q-1)+2)=ddns(i,j);
                        A(k,3*(q-1)+3)=-1;
                        clear q
                    end
                end

                %SPF:Tikhonov
                if obj.method==0
                    if obj.flag==0
                        for k=obj.aspect.a+1+obj.aspect.d+obj.aspect.spf:obj.aspect.pairs
                            q=k-obj.aspect.a-obj.aspect.d-obj.aspect.spf;
                            A(k,q)=obj.lamda;
                            clear q
                        end
                    elseif obj.flag==1
                        for k=obj.aspect.a+1+obj.aspect.d+obj.aspect.spf:obj.aspect.pairs
                            q=k-obj.aspect.a-obj.aspect.d-obj.aspect.spf;
                            A(k,q)=-obj.lamda;
                            A(k,q+3)=obj.lamda;
                            clear q
                        end
                    elseif obj.flag==2
                        for k=obj.aspect.a+1+obj.aspect.d+obj.aspect.spf:obj.aspect.pairs
                            q=k-obj.aspect.a-obj.aspect.d-obj.aspect.spf;
                            A(k,q)=obj.lamda;
                            A(k,q+3)=-2*obj.lamda;
                            A(k,q+6)=obj.lamda;
                            clear q
                        end
                    end
                end

                clear k time delta_t ddew ddns phia phid thetaa thetad
            end
        end

        function b=make_b(obj,preobj,a,i,j)
            b=zeros(obj.aspect.pairs,1);
            if strcmpi('POT',preobj.tech)==1
                for k=1:obj.aspect.a
                    b(k)=a.aa(i,j,k);
                end
                for k=1:obj.aspect.a
                    b(k+obj.aspect.a)=a.ar(i,j,k);
                end
                for k=1:obj.aspect.d
                    b(k+obj.aspect.a*2)=a.da(i,j,k);
                end
                for k=1:obj.aspect.d
                    b(k+obj.aspect.a*2+obj.aspect.d)=a.dr(i,j,k);
                end
                clear k
            elseif strcmpi('insar',preobj.tech)==1
                for k=1:obj.aspect.a
                    b(k)=a.a(i,j,k);
                end
                for k=1:obj.aspect.d
                    b(k+obj.aspect.a)=a.d(i,j,k);
                end
                clear k
            end
        end

        function e=del_std(obj,a)
            b=std(a,1,'all','omitnan' );
            c=mean(a,'all','omitnan' );
            d=[c-10*b c+10*b];
            a(a>d(2))=nan;
            a(a<d(1))=nan;
            e=a;
            clear b c d a
        end


    end
end