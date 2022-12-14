%% reinitialize
close all;
clear;
clc;
addpath(genpath('./functions/'))
%% define paths 
paths = [];
paths.directory = 'D:\BIO\PhD\ditalia\somitoid\data\image\Somitoid_single cell tracking\';
paths.dataname = 'Somitoid1_ 72-84hr_6min interval_2021_09_10__10_15_25-data';
paths.plotFolder = 'D:\BIO\PhD\ditalia\somitoid\github\ditalia-somitoid\figure\';
% define pixelsize
xy_pxsize = 0.692;
z_pxsize = 2.24;
t_size = 6; % min
%% load data (output from experiment_somitoid1)
load(['./cell_struct/',paths.dataname,'_cell_struct.mat'])






%% calculate per frame per z gfp value for calibration
%% load raw c,(y,x),z,t; transposed compared with tracking-results *
h5disp([paths.directory,'Somitoid1_ 72-84hr_6min interval_2021_09_10__10_15_25.h5'])
raw_image = h5read([paths.directory,'Somitoid1_ 72-84hr_6min interval_2021_09_10__10_15_25.h5'], '/data');
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









%%
%% pre-processing
%% trim cell_struct and calculate normalized intensities
frame_coverage = 0.9;
max_frame = max(cat(1,cell_struct.frame))+1; % start with frame 0 so +1
good_cell = [];
image_size = size(raw_image_h2b,1);
for i = [cell_struct.cellId]
    cell_struct_i = cell_struct([cell_struct.cellId]==i);
    % coverage
    if numel(cell_struct_i.frame) >= (frame_coverage * max_frame)
        % remove cells from the edges
        nonedge_filter = all(cell_struct_i.coordinates(:,1:2)<=0.9*image_size,2) & ...
            all(cell_struct_i.coordinates(:,1:2)>=0.1*image_size,2);
        if sum(nonedge_filter) >= (0.8 * numel(cell_struct_i.frame))
            % forbid more than two divisions each cell track
            if numel(cell_struct_i.division) <= 2
                good_cell = [good_cell i];
            end
        end
    end
end

cell_struct_trim = cell_struct(ismember([cell_struct.cellId],good_cell));

% add intensity ratio red(dense)/green(sparse)
ratios_old = arrayfun(@(x)x.Mean_Intensity_1./x.Mean_Intensity_0,cell_struct_trim,'UniformOutput',false);
[cell_struct_trim.intensity_ratio_old] = ratios_old{:};

% add average pixel intensity per t/z
meanh2b_t_z = arrayfun(@(x)t_by_z_h2b(sub2ind(size(t_by_z_h2b),x.frame+1,x.coordinates(:,3)+1)),...
    cell_struct_trim,'UniformOutput',false);
[cell_struct_trim.meanh2b_t_z] = meanh2b_t_z{:};

% add intensity ratio normalized by t_by_z_h2b
ratios = arrayfun(@(x)x.Mean_Intensity_1./x.meanh2b_t_z,cell_struct_trim,'UniformOutput',false);
[cell_struct_trim.intensity_ratio] = ratios{:};


%%
%%
%% surrounding intensity
%% load raw c,(y,x),z,t; transposed compared with tracking-results *
h5disp([paths.directory,'Somitoid1_ 72-84hr_6min interval_2021_09_10__10_15_25.h5'])
raw_image = h5read([paths.directory,'Somitoid1_ 72-84hr_6min interval_2021_09_10__10_15_25.h5'], '/data');

% max projection for raw mesp2 (c2)
% raw_proj = squeeze(max(raw_image(2,:,:,:,:),[],4));

%% surrounding intensity calculate ----- old (I call it old because I had something 'new' but didn't work
%% max proj
raw_image_mesp2 = squeeze(raw_image(2,:,:,:,:));
raw_proj = squeeze(max(raw_image_mesp2(:,:,:,:),[],3));

%% guassfilt to create surrounding background
% resolution = 62 pixels
raw_proj_filt = imgaussfilt(raw_proj(:,:,:),51);
%% calculate surrounding signal
% mean_wholeimage_raw = squeeze(mean(raw_proj, [1,2]));
surrounding_intensity_cell = {};
for i = 1:numel(cell_struct_trim)
    cell_struct_i = cell_struct_trim(i);
    surrounding_intensity = nan([numel(cell_struct_i.frame),1]);
    for j = 1:numel(cell_struct_i.frame)
        frame = cell_struct_i.frame(j)+1;
        xy_coor = cell_struct_i.coordinates(j,1:2)+1;
        surrounding_intensity(j) = double(raw_proj_filt(xy_coor(1),xy_coor(2),frame));%./mean_wholeimage_raw(frame);
    end
    surrounding_intensity_cell = [surrounding_intensity_cell, surrounding_intensity];
end

[cell_struct_trim.surrounding_intensity_old] = surrounding_intensity_cell{:};



%% save data
save(['./cell_struct/',paths.dataname,'_cell_struct_trim.mat'], "cell_struct_trim")
disp(['./cell_struct/',paths.dataname,'_cell_struct_trim.mat' ' Saved']);




%% save surrounding last
somitoid1_surrounding_last = raw_proj_filt(:,:,end);
save(['./cell_struct/','somitoid1_surrounding_last.mat'], "somitoid1_surrounding_last")
disp(['./cell_struct/','somitoid1_surrounding_last.mat' ' Saved']);
