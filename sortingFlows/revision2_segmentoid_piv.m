%%
% segmentoid1-101, segmentoidD-102, segmentoidA-103
%% reinitialize
close all;
clear;
clc;
addpath(genpath('./functions/'))
%% define paths 
paths = [];
paths.directory = 'D:\BIO\PhD\ditalia\somitoid\data\image\Segmentoid_H2B spiking\to_alessandro_and_piv\';
paths.expname1 = 'segmentoid1_ROT_t35_MAX_ROI_H2B spiking';
paths.expname2 = 'segmentoidD_REG_t1_MAX_ROID_2022_05_12__23_33_17_Subset-D';
paths.expname3 = 'segmentoidA_REG_t100_MAX_ROIA_2022_05_12__23_33_17_Subset-A';

paths.dataname1 = [paths.expname1,'-data'];
paths.dataname2 = [paths.expname2,'-data'];
paths.dataname3 = [paths.expname3,'-data'];
paths.plotFolder = 'D:\BIO\PhD\ditalia\somitoid\github\ditalia-somitoid\code_revision2_alvin\figure_scripts_final_segmentoid';
% define pixelsize
xy_pxsize = 0.692; % um
z_pxsize = 8; % um
t_size = 6; % min

optsSave=[];
optsSave.compress='lzw';
optsSave.overwrite=true;

% roiStack = loadtiff([roipath]);
% saveastiff(s2,[paths.outFolder filesep st_dir(i).name(1:end)],optsSave);














%% segmentoid1
%% load image 
segmentoid1_proj = loadtiff([paths.directory paths.expname1 '.tif']);
%% save as separate
mkdir([paths.directory paths.expname1])

%% mask for tissue
figure
segmentoid1_mask = false(size(segmentoid1_proj));
for i = 1:size(segmentoid1_proj,3)
im = segmentoid1_proj(:,:,i);
% figure;
imshow(im-offset_101)
filt_size = 30;
im_ave = imfilter(im,ones(filt_size)./filt_size);
% imshow(im_ave)

thres_otsu = graythresh(im_ave);
im_otsu = imbinarize(im_ave,thres_otsu);
% imshow(im_otsu)
segmentoid1_mask(:,:,i) = (im_otsu);


hold on;
B = bwboundaries(segmentoid1_mask(:,:,i));
boundary = B{1};
plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1);
title(num2str(i));
hold off;
drawnow;
saveastiff(im-offset_101,[[paths.directory paths.expname1] filesep 'frame' sprintf('%03d',i) '.tif']);
end


%% save piv
mask = segmentoid1_mask;
offset = offset_101;
save(['./cell_struct/piv_offset_' paths.expname1 '.mat'],"u_original","v_original","x","y","mask","offset");
clear divergence
%% load
load(['./cell_struct/piv_offset_' paths.expname1 '.mat'])
%% post-process piv
u=[];
v=[];
filt_size = 10;
% for i=1:80
for i=40:113

maski = imresize(mask(:,:,i),0.2);
maski = maski(2:end-1,2:end-1);

ui = nanconv(u_original{i},ones(filt_size)./filt_size,'nanout');
vi = nanconv(v_original{i},ones(filt_size)./filt_size,'nanout');
ui(~maski)=nan;
vi(~maski)=nan;

u=cat(3,u,ui);
v=cat(3,v,vi);
end
uAll = nansum(u,3);
vAll = nansum(v,3);
div = divergence(uAll,vAll)./t_size;
figure;
imagesc(div,[-5 5]);
axis equal


%% mean divergence in y
figure();
plot(mean(div,2),'-')
ylabel('average divergence in y (1/min)')
xlabel('distance from top')

%% smoothening with csaps
x = {1:size(div,1),1:size(div,2)};
[xx,yy] = ndgrid(x{1},x{2});
% y = peaks(xx, yy);
% figure
% surf(xx,yy,div)
% axis off
[sval,p] = csaps(x,div,0.01,x);
figure;
imagesc(sval,[-5 5]);
axis equal

%% velocity field
spacing = 5;
hold on;
quiver(yy(1:spacing:end,1:spacing:end),xx(1:spacing:end,1:spacing:end),uAll(1:spacing:end,1:spacing:end),vAll(1:spacing:end,1:spacing:end),'k')
hold off;

















%% segmentoidD
%% load image 
segmentoidD_proj = loadtiff([paths.directory paths.expname2 '.tif']);
%% save as separate
mkdir([paths.directory paths.expname2])

%% mask for tissue
figure
segmentoidD_mask = false(size(segmentoidD_proj));
for i = 1:size(segmentoidD_proj,3)
im = segmentoidD_proj(:,:,i);
% figure;
imshow(im-offset_102)
filt_size = 30;
im_ave = imfilter(im,ones(filt_size)./filt_size);
% imshow(im_ave)

thres_otsu = graythresh(im_ave);
im_otsu = imbinarize(im_ave,thres_otsu);
% imshow(im_otsu)
segmentoidD_mask(:,:,i) = (im_otsu);


hold on;
B = bwboundaries(segmentoidD_mask(:,:,i));
boundary = B{1};
plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1);
title(num2str(i));
hold off;
drawnow;
saveastiff(im-offset_102,[[paths.directory paths.expname2] filesep 'frame' sprintf('%03d',i) '.tif']);
end

%% save piv
mask = segmentoidD_mask;
offset = offset_102;
save(['./cell_struct/piv_offset_' paths.expname2 '.mat'],"u_original","v_original","x","y","mask","offset");
clear divergence
%% load
load(['./cell_struct/piv_offset_' paths.expname2 '.mat']);
%% post-process piv
u=[];
v=[];
filt_size = 10;
for i=70:150 %1:size(v_original)
% for i=100:180 %1:size(v_original)
% for i=170:250 %1:size(v_original)

maski = imresize(mask(:,:,i),0.2);
maski = maski(2:end-1,2:end-1);

ui = nanconv(u_original{i},ones(filt_size)./filt_size,'nanout');
vi = nanconv(v_original{i},ones(filt_size)./filt_size,'nanout');
ui(~maski)=nan;
vi(~maski)=nan;

u=cat(3,u,ui);
v=cat(3,v,vi);
end
uAll = nansum(u,3);
vAll = nansum(v,3);
div = divergence(uAll,vAll)./t_size;
% div(~maski)=nan;
figure;
imagesc(div,[-5 5]);
axis equal


%% mean divergence in y
figure();
plot(mean(div,2),'-')
ylabel('average divergence in y (1/min)')
xlabel('distance from top')


% %% pre-process before csaps
% div = div(~all(isnan(div),2),~all(isnan(div),1))';
% % div = div(:,1);
% % div(isnan(div))=0;
% % div(1,1)=nan;
% div = div(10:end-10,:);
% %%
% [sval,p] = csaps(1:size(div,1),div,0.01,1:size(div,1));
% figure;
% imagesc(sval,[-5 5]);
% axis equal
%% smoothening with csaps
x = {1:size(div,1),1:size(div,2)};
[xx,yy] = ndgrid(x{1},x{2});
% y = peaks(xx, yy);
% figure
% surf(xx,yy,div)
% axis off
[sval,p] = csaps(x,div,0.01,x);
figure;
imagesc(sval,[-5 5]);
axis equal


%% velocity field
spacing = 5;
hold on;
quiver(yy(1:spacing:end,1:spacing:end),xx(1:spacing:end,1:spacing:end),uAll(1:spacing:end,1:spacing:end),vAll(1:spacing:end,1:spacing:end),'k')
hold off;















%% segmentoidA
%% load image 
segmentoidA_proj = loadtiff([paths.directory paths.expname3 '.tif']);
%% save as separate
mkdir([paths.directory paths.expname3])

%% mask for tissue
figure
segmentoidA_mask = false(size(segmentoidA_proj));
for i = 1:size(segmentoidA_proj,3)
im = segmentoidA_proj(:,:,i);
% figure;
imshow(im-offset_103)
filt_size = 30;
im_ave = imfilter(im,ones(filt_size)./filt_size);
% imshow(im_ave)

thres_otsu = graythresh(im_ave);
im_otsu = imbinarize(im_ave,thres_otsu);
% imshow(im_otsu)
segmentoidA_mask(:,:,i) = (im_otsu);


hold on;
B = bwboundaries(segmentoidA_mask(:,:,i));
boundary = B{1};
plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1);
title(num2str(i));
hold off;
drawnow;
saveastiff(im-offset_103,[[paths.directory paths.expname3] filesep 'frame' sprintf('%03d',i) '.tif']);
end


%% save piv
mask = segmentoidA_mask;
offset = offset_103;
save(['./cell_struct/piv_offset_' paths.expname3 '.mat'],"u_original","v_original","x","y","mask","offset");
clear divergence
%% load
load(['./cell_struct/piv_' paths.expname3 '.mat']);
%% post-process piv
u=[];
v=[];
filt_size = 10;
% for i=1:80
for i=50:130 %1:size(v_original)
maski = imresize(mask(:,:,i),0.2);
maski = maski(2:end-1,2:end-1);

ui = nanconv(u_original{i},ones(filt_size)./filt_size,'nanout');
vi = nanconv(v_original{i},ones(filt_size)./filt_size,'nanout');
ui(~maski)=nan;
vi(~maski)=nan;

u=cat(3,u,ui);
v=cat(3,v,vi);
end
uAll = nansum(u,3);
vAll = nansum(v,3);
div = divergence(uAll,vAll)./t_size;
figure;
imagesc(div,[-5 5]);
axis equal


%% mean divergence in y
figure();
plot(mean(div,2),'-')
ylabel('average divergence in y (1/min)')
xlabel('distance from top')

%% smoothening with csaps
x = {1:size(div,1),1:size(div,2)};
[xx,yy] = ndgrid(x{1},x{2});
% y = peaks(xx, yy);
% figure
% surf(xx,yy,div)
% axis off
[sval,p] = csaps(x,div,0.01,x);
figure;
imagesc(sval,[-5 5]);
axis equal

%% velocity field
spacing = 5;
hold on;
quiver(yy(1:spacing:end,1:spacing:end),xx(1:spacing:end,1:spacing:end),uAll(1:spacing:end,1:spacing:end),vAll(1:spacing:end,1:spacing:end),'k')
hold off;
