# 3-D-Deformation-Inversion
<img width="416" alt="image" src="https://github.com/arvin1367/3-D-Deformation-Inversion/assets/49364261/5d3bd8bd-1613-4d0a-ad6c-36321716a209">

## What is this repository for?

A MATLAB toolbox for automatically calculating SAR-derived 3-D deformation maps of glacier, landslide, and land subsidence.


## How do I get set up?
Copy the above folder to MATLAB 2022a or higher. We provide real-time scripts (main_3d_glacier.mlx, main_3d_landslide.mlx, main_3d_land_subsidence.mlx) for different research subjects. You will have a clear understanding of how the software operates.

## Usage
We have provided an example of a glacier from ([Dem Glacier Data](https://data.mendeley.com/datasets/nj7956xtbm/1)) to download the data.
Here is the real-time script for the glacier, which is the usage steps：
## 3D inversion for glaciers (SBAS-POT results)
Software and system preparation: Supports Windows, Mac, and Linux systems, requiring Matlab 2022a or higher.
#### Data preparation:
Add this code package to the path first. Contains a parent folder and five sub folders, as well as a time table of image pairs for ascending and descending.
#### example:
* ascending_azimuth      (dir) 
 
* ascending_range        (dir)

* descending_azimuth     (dir)
 
* descending_range       (dir)

* supporting_documents   (dir)

* date_ascending.xls     (table)

* date_descending.xls    (table)

Notice: Please ensure that the folder and file names are consistent, except for the sbas-pot processing results, which are the files in the first four folders.

supporting_documents include the local evalution angle(lv_theta.tif) and orientation angle(lv_phi.tif) obtained from GAMMA software(look_vector),dem_ddew.tif and dem_ddns.tif obtained from myDEMdiff(). lv_theta_a.tif and lv_theta_d.tif are ascending and descending,respectively. lv_phi_a.tif and lv_phi_d.tif are ascending and descending,respectively.

To calculate the first order difference between the east-west and north-south directions of the DEM, it is necessary to provide larger DEM data and reference data than the research area. Final obtained dem_ddew.tif and dem_ddns.tif, place it in supporting_documents folder.
####  step 1: myDEMdiff()  
Start data preprocessing, mainly including reading data, constructing zero mask, and calculating time.
####  step 2: myPreprocess()
Start data processing, mainly including constructing coefficient matrices and constant matrices, conducting three-dimensional inversion, saving data(TIF and MAT formats, speed, cumulative, and annual average), and calculating norms.
####  step 3: myProcess()
Result analysis, mainly calculating statistical indicators of uncertainty, and drawing and displaying them.
####  step 4: myReprocess()
```m
clear
clc
```
#### example:
##### (1) Calculate the first-order difference of DEM using functions (myDEMdiff) in east-west (dem_ddew.tif) and north-south (dem_ddns.tif) directions for SPF constraint. Input parameter one: slightly larger than the DEM of the research area; Output parameter 2: The size and geographic information tif that needs to be output.
```m
myDEMdiff('other/njbw_dem.tif','njbw/ascending_range/20221004_20221016r.tif') 
```
##### (2) dirc is the parent floder, this example is 'njbw'.
```m
dirc='njbw';
```
##### (3) Start data preprocessing by entering the name of the parent folder. If you do not have a local incidence angle file, you can enter the incidence angle of ascending in the second parameter, the azimuth of ascending in the third parameter, the incidence angle of descending in the fourth parameter, and the azimuth of descending in the fifth parameter.     eg. obj_arvin=myPreprocess_3dglacier(dirc,36.9156,347.2974,37.8561,192.6945);POT InSAR
```m
obj_arvin=myPreprocess(dirc,'POT'); 
save('obj_arvin.mat');
```
##### (4) To start 3D inversion, input data preprocessing is required for the structure. Including two methods for solving ill conditioned equations: Tikhonov and TSVD.For Tikhonov: obj_data=myProcess_3dglacier(obj_arvin,'Tikhonov',0.1,1);The second parameter is the regularization parameter, recommended 0.1 (large deformation, e.g. glacier), and the third parameter is the regularization order, usually 0, 1, 2 (recommended 1). For TSVD: obj_data=myProcess_3dglacier(obj_arvin,'TSVD',0.5);The second parameter is the singular value truncation value, recommended 0.5 (large deformation, e.g. glacier).
```m
obj_data=myProcess(obj_arvin,'TSVD',0.5); 
save('obj_data.mat');
```
##### (5) Start the post-processing process, interactively select a stable region (draw a box in the stable region) as the random error, and then calculate the uncertainty.If an error is reported, it indicates that there are a large number of empty values in the stable region and a new selection is needed.
```m
obj_un=myReprocess(obj_arvin,obj_data); 
save('obj_un.mat');
```
##### (6) Draw deformation maps in various directions, parameter 1: calculated tif; parameter 2: DEM of the area, which needs to be consistent with the deformation range. If not, you can use 'no' instead; parameter 3: scale factor, set according to the user's drawing range; parameter 4: select input, if any, it is the line shp file of the deformation contour, such as the glacier range.
```m
mkfigure(obj_arvin,'dem_t1_01/vu.tif','other/demkk.tif',2,'other/njbw.shp'); 
```
##### (7) Draw dynamic deformation maps in various directions, parameter 1: calculated tif, rate between each time period; parameter 2: DEM of the area, which needs to be consistent with the deformation range, If not, you can use 'no' instead; parameter 3: scale factor, set according to the user's drawing range; parameter 4: color bar range, set according to the deformation size, parameter 5: select input, and some are the deformation contour line shp files, such as glacier range.
```m
mkgifs('njbw_results/results_velocity/ve','other/demkk.tif',2,[-1.5 1.5],'ve.gif','other/njbw.shp'); 
```
##### (8) Draw a deformation rate map with vector arrows, with the base color representing the upper and lower shape variables, and the arrow length representing the horizontal shape variables. The first three parameters represent the solving tif in the north, east, and vertical directions, respectively. Parameter 4: DEM of the area, which needs to be consistent with the deformation range, If not, you can use 'no' instead. Parameter 5: Scale factor, set according to the user's drawing range. Parameter 6: Vector arrows use spacing, reflecting arrow density. Parameter 7: Arrow scaling factor. Parameter 8: Select input. If there are any, it will be the line shp file of the deformation contour, such as the glacier range.
```m
mkfigure_v('njbw_results/vn_mask.tif','njbw_results/ve_mask.tif','njbw_results/vu_mask.tif','other/demkk.tif',2,10,10,'other/njbw.shp'); 
```
##### (9) Draw a vector arrow profile. Interactive selection of the profile to be drawn, supporting a straight line. The first three parameters represent the solving tif in the north, east, and vertical directions, respectively. Parameter 4: The DEM of the region needs to be consistent with the deformation range.
```m
mksec('njbw_results/vn_mask.tif','njbw_results/ve_mask.tif','njbw_results/vu_mask.tif','other/demkk.tif',0.05); 
```
##### (10) Draw a timing chart. Interactive selection of points and drawing of their long-term cumulative deformation map.Parameter 1: results; Parameter 2: time.
```m
mkpt('njbw_results/myFile_out',obj_arvin.date.time)
```
