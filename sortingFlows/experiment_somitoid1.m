%% reinitialize
close all;
clear;
clc;
addpath(genpath('./functions/'))
%% define paths *
paths = [];
paths.directory = 'D:\BIO\PhD\ditalia\somitoid\data\image\Somitoid_single cell tracking\';
paths.dataname = 'Somitoid1_ 72-84hr_6min interval_2021_09_10__10_15_25-data';

%% load track matrix,CVS-table output from ilastik % object_center is float, center_of_the_object is integer
track_csv=readtable([paths.directory,paths.dataname, '_CSV-Table.csv']);
column_names = ["frame" "labelimageId" "trackId" "lineageId" "parentTrackId" "mergerLabelId"...
    "Center_of_the_object_0" "Center_of_the_object_1" "Center_of_the_object_2"...
    "Size_in_pixels_0" "Mean_Intensity_0" "Mean_Intensity_1" "Object_Center_2"];
track_use = track_csv(:,column_names);



%% create struct for individual cells from 2channel track csv **
% (cells after divisions are considered as different cells)
% ilastik exporting results can change each time you switch output type! perhaps because
% it's non-deterministic
% intensity_0: sparse label (green) H2B-GFP
% intensity_1: signal (red) MESP2-mCherry
cell_struct = [];
cell_count = 0;
somitoidID = 1;
lineage_id = unique(track_use.lineageId);
for i = 1:numel(lineage_id)
    lineage = track_use(track_use.lineageId==lineage_id(i),:);
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
