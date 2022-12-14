%% reinitialize
close all;
clear;
clc;
addpath(genpath('./functions/'))
%% define paths 
paths = [];
paths.directory = 'D:\BIO\PhD\ditalia\somitoid\data\image\Segmentoid_H2B spiking\';
paths.expname = 'ROID_2022_05_12__23_33_17_Subset-D';
paths.dataname = [paths.expname,'-data'];
paths.plotFolder = 'D:\BIO\PhD\ditalia\somitoid\github\ditalia-somitoid\figure\segmentoid';
% define pixelsize
xy_pxsize = 0.692;
z_pxsize = 8;
t_size = 6; % min
%% load data
load(['./cell_struct/',paths.dataname,'_cell_struct.mat'])








%% calculate per frame per z gfp value
%% load raw c,(y,x),z,t; transposed compared with tracking-results *
h5disp([paths.directory,paths.expname,'.h5'])
raw_image = h5read([paths.directory,paths.expname,'.h5'], '/data');
raw_image_h2b = squeeze(raw_image(1,:,:,:,:));
%% load tracking-results *
tracking_results = h5read([paths.directory,paths.dataname,'_Tracking-Result.h5'], '/exported_data');
cell_mask = squeeze(tracking_results)>0;

%% t_by_z matrix
t_by_z_h2b = nan(size(cell_mask,[4,3]));
for t = 1:size(cell_mask,[4])
    for z = 1:size(cell_mask,[3])
        cell_mask_zt = cell_mask(:,:,z,t);
        raw_h2b_zt = raw_image_h2b(:,:,z,t);
        t_by_z_h2b(t,z) = sum(raw_h2b_zt(cell_mask_zt),'all')./sum(cell_mask_zt,'all');


    end
end







%% identify segment number and correct time
%%
segment1_x = uint16([46,426]);
segment1_y = uint16([36,200]);
segment1_t0 = 70;
segment2_x = uint16([26,431]);
segment2_y = uint16([205,365]);
segment2_t0 = 100;
segment3_x = uint16([39,407]);
segment3_y = uint16([396,515]);
segment3_t0 = 170;

%% note the coordinates are not transposed, because the x,y coor for ilastik output are relative to the original image (as viewed in FIji)
cell_struct_segments = [];
for i = [cell_struct.cellId]
    cell_struct_i = cell_struct([cell_struct.cellId]==i);
    if ~isempty(find(cell_struct_i.frame==segment1_t0,1)) && ...
            all(cell_struct_i.coordinates(find(cell_struct_i.frame==segment1_t0),[1,2]) < [segment1_x(2),segment1_y(2)],"all") && ...
            all(cell_struct_i.coordinates(find(cell_struct_i.frame==segment1_t0),[1,2]) > [segment1_x(1),segment1_y(1)],"all")
        cell_struct_i.segment_num = 1;
        cell_struct_i.segment_t0 = segment1_t0;
        align_logical = cell_struct_i.frame >=segment1_t0;

    elseif ~isempty(find(cell_struct_i.frame==segment2_t0,1)) && ...
            all(cell_struct_i.coordinates(find(cell_struct_i.frame==segment2_t0),[1,2]) < [segment2_x(2),segment2_y(2)],"all") && ...
            all(cell_struct_i.coordinates(find(cell_struct_i.frame==segment2_t0),[1,2]) > [segment2_x(1),segment2_y(1)],"all")
        cell_struct_i.segment_num = 2;
        cell_struct_i.segment_t0 = segment2_t0;
        align_logical = cell_struct_i.frame >=segment2_t0;

    elseif ~isempty(find(cell_struct_i.frame==segment3_t0,1)) && ...
            all(cell_struct_i.coordinates(find(cell_struct_i.frame==segment3_t0),[1,2]) < [segment3_x(2),segment3_y(2)],"all") && ...
            all(cell_struct_i.coordinates(find(cell_struct_i.frame==segment3_t0),[1,2]) > [segment3_x(1),segment3_y(1)],'all')
        cell_struct_i.segment_num = 3;
        cell_struct_i.segment_t0 = segment3_t0;
        align_logical = cell_struct_i.frame >=segment3_t0;
        
    else
        continue

    end
    cell_struct_i.frame = cell_struct_i.frame(align_logical);
    cell_struct_i.Mean_Intensity_0 = cell_struct_i.Mean_Intensity_0(align_logical);
    cell_struct_i.Mean_Intensity_1 = cell_struct_i.Mean_Intensity_1(align_logical);
    cell_struct_i.Size_in_pixels_0 = cell_struct_i.Size_in_pixels_0(align_logical);
    cell_struct_i.coordinates = cell_struct_i.coordinates(align_logical,:);
    cell_struct_segments = [cell_struct_segments cell_struct_i];
end









%% pre-processing
%% trim cell_struct to remove dim cells and short time courses
frame_coverage = 0.9;
max_frame = max(cat(1,cell_struct_segments.frame))+1; % start with frame 0 so +1
good_cell = [];
% image_size = size(raw_image_h2b,1);
for i = [cell_struct_segments.cellId]
    cell_struct_i = cell_struct_segments([cell_struct_segments.cellId]==i);
    % coverage
    if numel(cell_struct_i.frame) >= 70%(frame_coverage * (max_frame-cell_struct_i.segment_t0))
        if numel(cell_struct_i.division) <= 2
            good_cell = [good_cell i];
        end
    end
end

cell_struct_trim = cell_struct_segments(ismember([cell_struct_segments.cellId],good_cell));

% add intensity ratio red(dense)/green(sparse)
ratios_old = arrayfun(@(x)x.Mean_Intensity_1./x.Mean_Intensity_0,cell_struct_trim,'UniformOutput',false);
[cell_struct_trim.intensity_ratio_old] = ratios_old{:};

% add average pixel intensity per t/z
meanh2b_t_z = arrayfun(@(x)t_by_z_h2b(sub2ind(size(t_by_z_h2b),x.frame+1,x.coordinates(:,3)+1)),...
    cell_struct_trim,'UniformOutput',false);
[cell_struct_trim.meanh2b_t_z] = meanh2b_t_z{:};

% add intensity ratio red(dense)/green(sparse)
ratios = arrayfun(@(x)x.Mean_Intensity_1./x.meanh2b_t_z,cell_struct_trim,'UniformOutput',false);
[cell_struct_trim.intensity_ratio] = ratios{:};
%%



%% calculate intensity normalized by average for each segment
smallest_final = min(arrayfun(@(x)x.frame(end)-x.segment_t0,cell_struct_trim));
largest_final = max(arrayfun(@(x)x.frame(end)-x.segment_t0,cell_struct_trim));

for i = 1:max([cell_struct_trim.segment_num])
mean_segment{i} = arrayfun(@(y)mean(cell2mat(arrayfun(@(x)x.intensity_ratio(x.frame-x.segment_t0==y)',...
    cell_struct_trim([cell_struct_trim.segment_num]==i),'UniformOutput',false))),...
    0:largest_final);
end


% intensity_ratio_norm = arrayfun(@(x)x.intensity_ratio./mean_segment{x.segment_num}(1:numel(x.intensity_ratio))', cell_struct_trim, 'UniformOutput', false);
intensity_ratio_norm = arrayfun(@(x)x.intensity_ratio./mean_segment{x.segment_num}(x.frame-x.segment_t0+1)', cell_struct_trim, 'UniformOutput', false);

[cell_struct_trim.intensity_ratio_norm] = intensity_ratio_norm{:};

%%
% IDD = num2cell([cell_struct_trim.somitoidID]+1);
% [cell_struct_trim.somitoidID] = IDD{:};


%% save data
save(['./cell_struct/',paths.dataname,'_cell_struct_trim.mat'], "cell_struct_trim")
disp(['./cell_struct/',paths.dataname,'_cell_struct_trim.mat' ' Saved']);
