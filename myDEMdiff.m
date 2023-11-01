classdef myDEMdiff
    properties
        dem;
        geoInfo;
    end
    
    methods
        function obj = myDEMdiff(demtif,objout)
            
            [dem,geoInfo] = readgeoraster(demtif);
            dem = double(dem);
            obj.dem = dem;
            obj.geoInfo = geoInfo;
            [~,crop_geoInfo] = readgeoraster(objout);
            
            %filter
            kernel_size = [17,17] ;%Filter Window
            demFilter = filtering(obj,dem,kernel_size,0);
            
            % first-order diff
            demFilter_resolution = 30;% demFilter resolution
            [demGrad_ew,demGrad_ns] = demGrad_ew_ns(obj,demFilter,demFilter_resolution); % Perform first order difference in north-south and east-west directions
            
            % resample
            pixel_length = crop_geoInfo.CellExtentInLatitude; 
            pixel_wideth = crop_geoInfo.CellExtentInLongitude; 
            demGrad_geoInfo = geoInfo;
            
            % Output as sampled grid data and its geographic coordinate information
            [resampledDEM_ew,resampledDEM_geoInfo] = demResampling(obj,demGrad_ew,demGrad_geoInfo,pixel_length,pixel_wideth,'nearest');
            [demCropResult_ew,demCropResult_geoInfo_ew] = cropByMask(obj,resampledDEM_ew,crop_geoInfo);
            geotiffwrite('ew.tif', single(demCropResult_ew), demCropResult_geoInfo_ew);
            
            [resampledDEM_ns,resampledDEM_geoInfo] = demResampling(obj,demGrad_ns,demGrad_geoInfo,pixel_length,pixel_wideth,'nearest');
            [demCropResult_ns,demCropResult_geoInfo_ns] = cropByMask(obj,resampledDEM_ns,crop_geoInfo);
            geotiffwrite('ns.tif', -single(demCropResult_ns), demCropResult_geoInfo_ns);
            
        end
        
        function demFilter = filtering(obj,demtiff,kernel_size, inputArg2)
            %   raw_resolution: Expected resolution for input DEM,lowSampleing_resolution: Expected DEM resolution to be output
            %   inputArg1 for filters: 'mean filtering', 'median filtering', 'Gaussian filtering'         
                        
            if (inputArg2 == 0)
                filter = fspecial('average',kernel_size);
                demFilter = imfilter(demtiff,filter,'replicate');
            elseif (inputArg2 == 1)
                demFilter = medfilt2(demtiff,kernel_size);
            else
                filter = fspecial('gaussian',kernel_size);
                demFilter = imfilter(demtiff,filter);
            end
        end

        function [demGrad_ew,demGrad_ns] = demGrad_ew_ns(obj,demFilter,resolution)
            [demGrad_ew,demGrad_ns] = gradient(demFilter,resolution);
        end

        function [demResampled,demResampled_geoInfo] = demResampling(obj,demGrad,demGrad_geoInfo,pixel_length,pixel_wideth,inputArg1)
            %   demGrad_ew,demGrad_ns Should be DEM data after first order difference
            %   raw_resolution,upSample_resolution,sampleRate are DEM image resolution and sampling rate for input and output, respectively
            %   inputArg1 is Interpolation Method of Resampling Algorithm：'nearest' 'linear' 'cubic' 'makima' 'spline'
            %sampleRate = double(upSample_resolution/raw_resolution);
            
            %Original DEM coordinate grid
            raw_lat = demGrad_geoInfo.LatitudeLimits(1)+0.5*demGrad_geoInfo.CellExtentInLatitude : demGrad_geoInfo.CellExtentInLatitude : demGrad_geoInfo.LatitudeLimits(2)-0.5*demGrad_geoInfo.CellExtentInLatitude; % 计算出X坐标
            raw_lon = demGrad_geoInfo.LongitudeLimits(2)-0.5*demGrad_geoInfo.CellExtentInLongitude : -demGrad_geoInfo.CellExtentInLongitude : demGrad_geoInfo.LongitudeLimits(1)+0.5*demGrad_geoInfo.CellExtentInLongitude; % 计算出Y坐标
            raw_lat = fliplr(raw_lat);
            raw_lon = fliplr(raw_lon);
            [DEM_rawlon,DEM_rawlat] = meshgrid(raw_lon,raw_lat);
            %Create a new coordinate grid
            Sample_lat = demGrad_geoInfo.LatitudeLimits(1)+0.5*pixel_length: pixel_length : demGrad_geoInfo.LatitudeLimits(2)-0.5*pixel_length;
            Sample_lon = demGrad_geoInfo.LongitudeLimits(2)-0.5*pixel_wideth : -pixel_wideth : demGrad_geoInfo.LongitudeLimits(1)+0.5*pixel_wideth;
            Sample_lat = fliplr(Sample_lat);
            Sample_lon = fliplr(Sample_lon);
            [DEM_Samplelon,DEM_Samplelat] = meshgrid(Sample_lon,Sample_lat);
            
            demResampled = interp2(DEM_rawlon, DEM_rawlat, demGrad, DEM_Samplelon,DEM_Samplelat , inputArg1);
            %Preserve DEM geographic information after resampling
            demResampled_geoInfo= georefcells(demGrad_geoInfo.LatitudeLimits,demGrad_geoInfo.LongitudeLimits,size(demResampled));
            demResampled_geoInfo.ColumnsStartFrom = 'north';
            %demResampled_ew = imresize(demGrad_ew,sampleRate,'bilinear');
            %demResampled_ns = imresize(demGrad_ns,sampleRate,'bilinear');

            %geotiffwrite('demResampled_ew.tif', output_resampledDEM, resampledDEM_geoInfo);
        end

        function [demCropResult,demCropResult_geoInfo] = cropByMask(obj,DEM,standardImg_geoInfo)
            %CropDEM is the DEM to be cropped, standardImg_ GeoInfo is the geographic coordinate information of standard images (including geographic coordinates)

            %Load geographic coordinate information of DEM to be cropped
            DEM_geoInfo = georefcells(obj.geoInfo.LatitudeLimits,obj.geoInfo.LongitudeLimits,size(DEM));
            DEM_geoInfo.ColumnsStartFrom = 'north';
            %Construct a latitude and longitude network for the DEM to be cropped
            DEM_lat = DEM_geoInfo.LatitudeLimits(1)+0.5*DEM_geoInfo.CellExtentInLatitude:DEM_geoInfo.CellExtentInLatitude:DEM_geoInfo.LatitudeLimits(2)-0.5*DEM_geoInfo.CellExtentInLatitude;
            DEM_lon = DEM_geoInfo.LongitudeLimits(2)-0.5*DEM_geoInfo.CellExtentInLongitude:-DEM_geoInfo.CellExtentInLongitude:DEM_geoInfo.LongitudeLimits(1)+0.5*DEM_geoInfo.CellExtentInLongitude;
            DEM_lat = fliplr(DEM_lat);
            DEM_lon = fliplr(DEM_lon);
            [DEM_lon,DEM_lat] = meshgrid(DEM_lon,DEM_lat);

            %Construct a latitude and longitude grid within the clipping range
            cropDEM_lat = standardImg_geoInfo.LatitudeLimits(1)+0.5*standardImg_geoInfo.CellExtentInLatitude:standardImg_geoInfo.CellExtentInLatitude:standardImg_geoInfo.LatitudeLimits(2)-0.5*standardImg_geoInfo.CellExtentInLatitude;
            cropDEM_lon = standardImg_geoInfo.LongitudeLimits(2)-0.5*standardImg_geoInfo.CellExtentInLongitude:-standardImg_geoInfo.CellExtentInLongitude:standardImg_geoInfo.LongitudeLimits(1)+0.5*standardImg_geoInfo.CellExtentInLongitude;
            cropDEM_lat = fliplr(cropDEM_lat);
            cropDEM_lon = fliplr(cropDEM_lon);
            [cropDEM_lon,cropDEM_lat] = meshgrid(cropDEM_lon,cropDEM_lat);
            %Interpolate the DEM to be cropped to the cropped range
            demCropResult = interp2(DEM_lon,DEM_lat, DEM,cropDEM_lon,cropDEM_lat , "nearest");
            %Geographic coordinate information of DEM after loading and cropping
            demCropResult_geoInfo = georefcells(standardImg_geoInfo.LatitudeLimits,standardImg_geoInfo.LongitudeLimits,size(demCropResult));
            demCropResult_geoInfo.ColumnsStartFrom = 'north';
            %Save cropped DEM
            %geotiffwrite('dem_clip.tif', demCropResult, demCropResult_geoInfo);
            
            
            %转换为像空间坐标
            %[lower_left_col,lower_left_row] = geographicToIntrinsic(DEM_geoInfo,position(1,1),position(2,1));
            %[upper_right_col,upper_right_row] = geographicToIntrinsic(DEM_geoInfo,position(1,2),position(2,2));
            %像素取整--根据matlab像素坐标的定义，采用四舍五入（round）
            %lower_left_row = round(lower_left_row);
            %lower_left_col = round(lower_left_col);
            %upper_right_row = round(upper_right_row);
            %upper_right_col = round(upper_right_col);
            %取整后转为经纬度
            %[crop_start_lat,crop_start_lon] = intrinsicToGeographic(DEM_geoInfo,lower_left_col-0.5,lower_left_row+0.5);
            %[crop_end_lat,crop_end_lon] = intrinsicToGeographic(DEM_geoInfo,upper_right_col+0.5,upper_right_row-0.5);
            %像空间坐标下进行裁剪
            %demCropResult = DEM(upper_right_row:lower_left_row, lower_left_col:upper_right_col);
            %加载裁剪后DEM的地理坐标信息
            %demCropResult_geoInfo = georefcells([crop_start_lat,crop_end_lat],[crop_start_lon,crop_end_lon],size(demCropResult));
            %demCropResult_geoInfo.ColumnsStartFrom = 'north';
            %保存裁剪后的DEM
            %geotiffwrite('dem_clip.tif', demCropResult, demCropResult_geoInfo);
        end
    end
end

