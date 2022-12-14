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
