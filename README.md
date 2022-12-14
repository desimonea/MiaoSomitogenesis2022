# MiaoSomitogenesis2022

Somitoids and segmentoids image processing. Sample code from:
Miao et al., Reconstruction and deconstruction of human somitogenesis in vitro, Nature, 2022

Copyright (C) 2022  The authors

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

The code requires MATLAB (MathWorks).
The code runs in MATLAB_R2022a. The code runs on a 
MacBook Pro (14-inch, 2021) macOS Monterey
Apple M1 Pro
Memory 32 GB 

INSTALLATION AND USAGE - Kimographs, nematic order and spatial correlation

The user can run the code by changing 
MATLAB's present working directory to 'MiaoSomitogenesis2022-main'.
MS Windows users have to adapt paths to MS Windows sintax.
 
Data files are too large to include within this Github repo and must be
downloaded separately from:
https://drive.google.com/file/d/15o-6nOMArxW6rLGYv5WqQHi1q7tN9ZGz/view?usp=share_link

The downloaded folder must be unzipped and tne 'data' folder must be
placed in 'MiaoSomitogenesis2022-main'.

MATLAB Toolboxes required are:
- Image Processing Toolbox 
- Statistics and Machine Learning Toolbox 
- Signal Processing Toolbox 
- Curve Fitting Toolbox 

The code requires external MATLAB functions:
loadtiff.m, saveastiff.m, nanconv.m, brewermap.m
Those functions be downloaded from MathWorks. We attach them to the code together with
their licenses.

The main image processing routine are:
1) mainSegmentoidKymos.m - It generates Mesp2 and Hes7 kimographs in segmentoids and calculates their oscillation properties.
2) mainNematic_segmentoid.m - It calculates the nematic order of Mesp2 signal in segmentoids. To run afte the previous.
3) mainSomitoid_Mesp2Uncx.m - It calculates the spatial correlation of Mesp2, Uncx and Mesp2-Uncx signals in somitoids.

The script can run altogether or each section sequentially.

A sample dataset is provided is provided (see Installation). Expected output for each step is provided. 8Gb RAM is required. Data sample runtime: minutes.

INSTALLATION AND USAGE - SORTING FLOWS

You current Folder needs to be "MiaoSomitogenesis2022-main/sortingFlows/" while running this pipeline.

Data files are too large to include within this Github repo and must be
downloaded separately from:
https://drive.google.com/file/d/190IquISjvHIV3okzDnBo805sy3sXpJPW/view?usp=share_link
The downloaded folder must be unzipped and the "data_tif" and "cell_struct" folders must be
placed in 'MiaoSomitogenesis2022-main/sortingFlows/'.

1. run experiment_*.m
2. run analysis_*.m except analysis_*_merge.m

These two steps take input from ilastik results and original movies, which are large files and not included. 
The outputs of these two steps are in the "cell_struct" folder from Google Drive.

The following scripts analyze different replicates of somitoid/segmentoid and are similar to each other.

3. run revision2_*_tracking.m and revision2_somitoid_disp.m to produce output figures into "MiaoSomitogenesis2022-main/sortingFlows_figures/figure_revision2_segmentoid" or "MiaoSomitogenesis2022-main/sortingFlows_figures/figure_revision2_somitoid"
This step takes input from "sortingFlows\cell_struct". The user will be able to run these scripts (for somitoids and segmentoids respectively)

Sample output can be downloaded from: 
https://drive.google.com/file/d/1wAUL9sT1bplaU9MkZxJR7p-szKSrFHop/view?usp=share_link

4. For PIV, revision2_segmentoid_piv.m and revision2_somitoid_piv.m are used for generating preliminary PIV outputs and plots. They take input from ilastik results and original movies, which are large files that are not included. 
5. revision2_segmentoid_piv_find_offset.m find offsets for removing MESP2- cells, as described in Methods. Output of this step is already stored in the "cell_struct" folder.
6. revision2_somitoid_piv_drawROI.m is used to manually draw ROIs for statistical tests. The generated ROIs are already stored in "cell_struct". 
7. revision2_segmentoid_piv_AD16oct22_1031.m generate final figures for PIV into "MiaoSomitogenesis2022-main/sortingFlows_figures/figure_revision2_segmentoid_piv", taking inputs from "cell_struct" and "data_tif".
8. revision2_somitoid_piv_AD16oct22.m generate final figures for PIV into "MiaoSomitogenesis2022-main/sortingFlows_figures/figure_revision2_somitoid_piv", taking inputs from "cell_struct" and "data_tif".

Finally, source data of all figures from "sortingFlows" are generated and stored in "MiaoSomitogenesis2022-main/excel_data/"
