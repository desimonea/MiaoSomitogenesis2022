%%
% segmentoid1-101, segmentoidD-102, segmentoidA-103
%% reinitialize
close all;
clear;
clc;
addpath(genpath('./functions/'))
%% define paths 
paths = [];
paths.directory = 'D:\BIO\PhD\ditalia\somitoid\data\image\Somitoid_single cell tracking\piv\';
paths.expname1 = 'Somitoid1_mesp2_ 72-84hr_6min interval_2021_09_10__10_15_25_MAX stitch';
paths.expname2 = 'Somitoid2_mesp2_ 72-84hr_6min interval_2021_09_18__15_24_26_MAX stitch';

% paths.plotFolder = 'D:\BIO\PhD\ditalia\somitoid\github\ditalia-somitoid\code_alvin\figure_scripts_final_segmentoid';
% define pixelsize
xy_pxsize = 0.692;
z_pxsize = 2.24;
t_size = 6; % min



%% somitoid1
segmentoid_proj = loadtiff([paths.directory paths.(['expname' num2str(1)]) '.tif']);
imgEnd = segmentoid_proj(:,:,120);

for i = 1:5
    ROI_somitoid1_mesp2_high_i = roipoly(imadjust(imgEnd));
    save(['./cell_struct/','ROI_somitoid1_mesp2_high_' num2str(i) '.mat'], "ROI_somitoid1_mesp2_high_i")
end

for i = 1:5
    ROI_somitoid1_mesp2_low_i = roipoly(imadjust(imgEnd));
    save(['./cell_struct/','ROI_somitoid1_mesp2_low_' num2str(i) '.mat'], "ROI_somitoid1_mesp2_low_i")
end

%% somitoid2
segmentoid_proj = loadtiff([paths.directory paths.(['expname' num2str(2)]) '.tif']);
imgEnd = segmentoid_proj(:,:,120);

for i = 1:5
    ROI_somitoid2_mesp2_high_i = roipoly(imadjust(imgEnd));
    save(['./cell_struct/','ROI_somitoid2_mesp2_high_' num2str(i) '.mat'], "ROI_somitoid2_mesp2_high_i")
end

for i = 1:5
    ROI_somitoid2_mesp2_low_i = roipoly(imadjust(imgEnd));
    save(['./cell_struct/','ROI_somitoid2_mesp2_low_' num2str(i) '.mat'], "ROI_somitoid2_mesp2_low_i")
end

% 
% %% correct name somitoid1
% for i = 1:5
%     load(['./cell_struct/','ROI_somitoid1_mesp2_high_' num2str(i) '.mat'], "ROI_somitoid1_mesp2_high_i")
%     figure;
%     imshow(ROI_somitoid1_mesp2_high_i);
% %     ROI_somitoid2_mesp2_high_i = ROI_somitoid1_mesp2_high_i;
% %     save(['./cell_struct/','ROI_somitoid2_mesp2_high_' num2str(i) '.mat'], "ROI_somitoid2_mesp2_high_i")
% end
% 
% for i = 1:5
%     load(['./cell_struct/','ROI_somitoid1_mesp2_low_' num2str(i) '.mat'], "ROI_somitoid1_mesp2_low_i")
%     figure;
%     imshow(ROI_somitoid1_mesp2_low_i);
% %     ROI_somitoid2_mesp2_low_i = ROI_somitoid1_mesp2_low_i;
% %     save(['./cell_struct/','ROI_somitoid2_mesp2_low_' num2str(i) '.mat'], "ROI_somitoid2_mesp2_low_i")
% end
% %% correct name somitoid1
% for i = 1:5
%     load(['./cell_struct/','ROI_somitoid2_mesp2_high_' num2str(i) '.mat'], "ROI_somitoid2_mesp2_high_i")
%     figure;
%     imshow(ROI_somitoid2_mesp2_high_i);
% %     ROI_somitoid2_mesp2_high_i = ROI_somitoid1_mesp2_high_i;
% %     save(['./cell_struct/','ROI_somitoid2_mesp2_high_' num2str(i) '.mat'], "ROI_somitoid2_mesp2_high_i")
% end
% 
% for i = 1:5
%     load(['./cell_struct/','ROI_somitoid2_mesp2_low_' num2str(i) '.mat'], "ROI_somitoid2_mesp2_low_i")
%     figure;
%     imshow(ROI_somitoid2_mesp2_low_i);
% %     ROI_somitoid2_mesp2_low_i = ROI_somitoid1_mesp2_low_i;
% %     save(['./cell_struct/','ROI_somitoid2_mesp2_low_' num2str(i) '.mat'], "ROI_somitoid2_mesp2_low_i")
% end
% 

