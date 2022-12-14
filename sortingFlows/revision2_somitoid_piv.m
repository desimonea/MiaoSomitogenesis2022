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


paths.plotFolder = 'D:\BIO\PhD\ditalia\somitoid\github\ditalia-somitoid\code_revision2_alvin\figure_revision2_somitoid_piv';
mkdir(paths.plotFolder);
% define pixelsize
xy_pxsize = 0.692; % um
z_pxsize = 2.24; % um
t_size = 6; % min

optsSave=[];
optsSave.compress='lzw';
optsSave.overwrite=true;

% roiStack = loadtiff([roipath]);
% saveastiff(s2,[paths.outFolder filesep st_dir(i).name(1:end)],optsSave);














%% segmentoid1
%% load image 
somitoid1_proj = loadtiff([paths.directory paths.expname1 '.tif']);
%% save as separate
mkdir([paths.directory filesep paths.expname1])

%% mask for tissue
figure
somitoid1_mask = false(size(somitoid1_proj));
for i = 1:size(somitoid1_proj,3)
im = somitoid1_proj(:,:,i);
% figure;
imshow(im,[0,max(im,[],'all')])
filt_size = 400;
im_ave = imfilter(im,ones(filt_size)./filt_size);
% imshow(im_ave)

thres_otsu = graythresh(im_ave);
im_otsu = imbinarize(im_ave,thres_otsu);
% imshow(im_otsu)
somitoid1_mask(:,:,i) = (im_otsu);


hold on;
B = bwboundaries(somitoid1_mask(:,:,i));
boundary = B{1};
plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1);
title(num2str(i));
hold off;
drawnow;
% saveastiff(im,[[paths.directory paths.expname1] filesep 'frame' sprintf('%03d',i) '.tif']);
end


%% save piv
mask = somitoid1_mask;
save(['./cell_struct/piv_' paths.expname1 '.mat'],"u_original","v_original","x","y","mask");
clear divergence
%% load
load(['./cell_struct/piv_' paths.expname1 '.mat'])
%% post-process piv
u=[];
v=[];
filt_size = 10;
% for i=1:80
for i=1:120

maski = imresize(mask(:,:,i),size(v_original{i}));

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

















%% somitoid2
%% load image 
somitoid2_proj = loadtiff([paths.directory paths.expname2 '.tif']);
%% save as separate
mkdir([paths.directory paths.expname2])

%% mask for tissue
figure
somitoid2_mask = false(size(somitoid2_proj));
for i = 1:size(somitoid2_proj,3)
im = somitoid2_proj(:,:,i);
% figure;
imshow(im,[0,max(im,[],'all')])
filt_size = 400;
im_ave = imfilter(im,ones(filt_size)./filt_size);
% imshow(im_ave)

thres_otsu = graythresh(im_ave);
im_otsu = imbinarize(im_ave,thres_otsu);
% imshow(im_otsu)
somitoid2_mask(:,:,i) = (im_otsu);


hold on;
B = bwboundaries(somitoid2_mask(:,:,i));
boundary = B{1};
plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1);
title(num2str(i));
hold off;
drawnow;
saveastiff(im,[[paths.directory paths.expname2] filesep 'frame' sprintf('%03d',i) '.tif']);
end

%% save piv
mask = somitoid2_mask;
save(['./cell_struct/piv_' paths.expname2 '.mat'],"u_original","v_original","x","y","mask");
clear divergence
%% load
load(['./cell_struct/piv_' paths.expname2 '.mat']);
%% post-process piv
u=[];
v=[];
filt_size = 10;
for i=1:70 %1:size(v_original)
% for i=100:180 %1:size(v_original)
% for i=170:250 %1:size(v_original)

maski = imresize(mask(:,:,i),size(v_original{i}));

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













