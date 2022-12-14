%%
% segmentoid1-101, segmentoidD-102, segmentoidA-103
%% reinitialize
close all;
clear;
clc;
addpath(genpath('./functions/'))
%% define paths 
paths = [];
paths.directory = 'D:\BIO\PhD\ditalia\somitoid\data\image\Segmentoid_H2B spiking\';
paths.expname1 = 'ROI_H2B spiking_Mesp2_Segmentoid_Day4 start_6min interval_2022_04_03__16_42_48_Subset';
paths.expname2 = 'ROID_2022_05_12__23_33_17_Subset-D';
paths.expname3 = 'ROIA_2022_05_12__23_33_17_Subset-A';

paths.dataname1 = [paths.expname1,'-data'];
paths.dataname2 = [paths.expname2,'-data'];
paths.dataname3 = [paths.expname3,'-data'];
paths.plotFolder = 'D:\BIO\PhD\ditalia\somitoid\github\ditalia-somitoid\code_revision2_alvin\figure_scripts_final_segmentoid';
% define pixelsize
xy_pxsize = 0.692;
z_pxsize = 8;
t_size = 6; % min

%% load and merge data
%% load data
%% load data
% cell_struct_trim1 = load(['./cell_struct_mid/',paths.dataname1,'_cell_struct_trim.mat']);
% cell_struct_trim2 = load(['./cell_struct_mid/',paths.dataname2,'_cell_struct_trim.mat']);
% cell_struct_trim3 = load(['./cell_struct_mid/',paths.dataname3,'_cell_struct_trim.mat']);
% 
cell_struct_trim1 = load(['./cell_struct/',paths.dataname1,'_cell_struct_trim.mat']);
cell_struct_trim2 = load(['./cell_struct/',paths.dataname2,'_cell_struct_trim.mat']);
cell_struct_trim3 = load(['./cell_struct/',paths.dataname3,'_cell_struct_trim.mat']);


cell_struct_trim = cat(2,cell_struct_trim1.cell_struct_trim, cell_struct_trim2.cell_struct_trim, cell_struct_trim3.cell_struct_trim);


%% 101
cell_struct_trim_101 = cell_struct_trim([cell_struct_trim.somitoidID]==101);

first_10_frames = cell2mat(arrayfun(@(s)s.Mean_Intensity_1(1:10),cell_struct_trim_101,'UniformOutput',0));

offset_101 = prctile(first_10_frames(:),25)

histogram(first_10_frames,20)

%% 102
cell_struct_trim_102 = cell_struct_trim([cell_struct_trim.somitoidID]==102);

first_10_frames = cell2mat(arrayfun(@(s)s.Mean_Intensity_1(1:10),cell_struct_trim_102,'UniformOutput',0));

offset_102 = prctile(first_10_frames(:),25)

histogram(first_10_frames,20)

%% 103
cell_struct_trim_103 = cell_struct_trim([cell_struct_trim.somitoidID]==103);

first_10_frames = cell2mat(arrayfun(@(s)s.Mean_Intensity_1(1:10),cell_struct_trim_103,'UniformOutput',0));

offset_103 = prctile(first_10_frames(:),25)

histogram(first_10_frames,20)
