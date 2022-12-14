0. You Curent Folder needs to be "\code_revision2_alvin\scripts_final" while running this pipeline.
1. run experiment_*.m
2. run analysis_*.m except analysis_*_merge.m

These two steps take input from ilastik results and original movies, which are extremely large files that are not included. 
The outputs of these two steps are in "code_revision2_alvin\scripts_final\cell_struct"

These scripts for different replicates of somitoid/segmentoid are more or less duplicates of eachothers.

3. run revision2_*_tracking.m and revision2_somitoid_disp.m to produce output figures into "code_revision2_alvin\figure_revision2_segmentoid" or "somitoid"
This step takes input from "code_revision2_alvin\scripts_final\cell_struct". You should be able to run these scripts (for somitoids and segmentoids respectively)
Feel free to delete figures in the figure-folder and reproduce them.

4. For PIV, revision2_segmentoid_piv.m and revision2_somitoid_piv.m are used for generating preliminary PIV outputs and plots. They take input from ilastik results and original movies, which are large files that are not included. 
5. revision2_segmentoid_piv_find_offset.m find offsets for removing MESP2- cells, as described in Methods.
6. revision2_somitoid_piv_drawROI.m is used to manually raw ROIs for statistical tests. The generated ROIs are already stored in "code_revision2_alvin\scripts_final\cell_struct".
7. revision2_segmentoid_piv_AD16oct22_1031.m generate final figures for PIV into "code_revision2_alvin\figure_revision2_segmentoid_piv", taking inputs from "cell_struct" and "data_tif".
8. revision2_somitoid_piv_AD16oct22.m generate final figures for PIV into "code_revision2_alvin\figure_revision2_somitoid_piv", taking inputs from "cell_struct" and "data_tif".