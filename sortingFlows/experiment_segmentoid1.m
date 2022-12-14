%% reinitialize
close all;
clear;
clc;
addpath(genpath('./functions/'))
%% define paths *
paths = [];
paths.directory = 'D:\BIO\PhD\ditalia\somitoid\data\image\Segmentoid_H2B spiking\';
paths.dataname = 'ROI_H2B spiking_Mesp2_Segmentoid_Day4 start_6min interval_2022_04_03__16_42_48_Subset-data';
xy_pxsize = 0.692;
z_pxsize = 8;
t_size = 6; % min
%% load track matrix % object center is decimal, center_of_the_object is integer
track_csv=readtable([paths.directory,paths.dataname, '_CSV-Table.csv']);
% trim table for long enough lineages
column_names = ["frame" "labelimageId" "trackId" "lineageId" "parentTrackId" "mergerLabelId"...
    "Center_of_the_object_0" "Center_of_the_object_1" "Center_of_the_object_2"...
    "Size_in_pixels_0" "Mean_Intensity_0" "Mean_Intensity_1" "Object_Center_2"];
track_use = track_csv(:,column_names);
track_use.x = track_use.Center_of_the_object_0.*xy_pxsize;
track_use.y = track_use.Center_of_the_object_1.*xy_pxsize;
track_use.z = track_use.Center_of_the_object_2.*z_pxsize;
track_use = track_use(:,[14:16,2:13,1]);

[track_use.Properties.VariableNames,['lineagdID_haem']];

%% track with the haematopoetic package

%%
param.mem = 2;
param.dim = 3;
param.good = 20;
param.quiet = 0;
track_output = track(table2array(track_use),20*xy_pxsize,param);
track_use_haem = array2table(track_output,'VariableNames',[track_use.Properties.VariableNames,['lineageId_haem']]);

%% create struct for individual cells from 2channel track csv ** segmentoidID = somitoidID+100
% (cells after divisions are considered as different cells, no division detection for segmentoids)
% ilastik exporting results can change each time you switch output type! because
% it's non-deterministic
% intensity_0: sparse label (green)
% intensity_1: signal (red)
cell_struct = [];
cell_count = 0;
somitoidID = 101;
lineage_id = unique(track_use_haem.lineageId_haem);
lineage_id = lineage_id(lineage_id>=0);
for i = 1:numel(lineage_id)
    lineage = track_use_haem(track_use_haem.lineageId_haem==lineage_id(i),:);
    cell_struct_i = [];
    % loop over dividion tracks and separate them into individual cells
    division = lineage(lineage.parentTrackId>0,:);
    if ~isempty(division)
        merge_tracks = tracks_to_merge(division);
        % construct structure for individual cells
        for k = 1:numel(merge_tracks)
            lineage_singlecell = lineage(ismember(lineage.trackId,merge_tracks{k}),:);
            cell_struct_i = create_one_struct(lineage_singlecell,cell_count,somitoidID,division);
            
            % update cell count and add to cell struct
            cell_count = cell_count+1;
            cell_struct = [cell_struct,cell_struct_i];
    
        end
    else
        % for lineage with no division
        lineage_singlecell = lineage;
        cell_struct_i = create_one_struct(lineage_singlecell,cell_count,somitoidID,division);
        % update cell count and add to cell struct
        cell_count = cell_count+1;
        cell_struct = [cell_struct,cell_struct_i];
    end
end
%% save data
save(['./cell_struct/',paths.dataname,'_cell_struct.mat'], "cell_struct")
disp(['./cell_struct/',paths.dataname,'_cell_struct.mat' ' Saved']);



