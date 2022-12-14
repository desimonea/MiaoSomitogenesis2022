%%
% segmentoid1-101, segmentoidD-102, segmentoidA-103
%% reinitialize
close all;
clear;
clc;
addpath(genpath('./functions/'))
addpath(genpath('./data/'))

%% define paths - treat segments separetely

paths = [];
paths.directory = './data_tif/';
paths.expname1 = 'segmentoid1_ROT_t35_MAX_ROI_H2B spiking';
paths.expname2 = 'segmentoid1_ROT_t35_MAX_ROI_H2B spiking';
paths.expname3 = 'segmentoidD_REG_t1_MAX_ROID_2022_05_12__23_33_17_Subset-D';
paths.expname4 = 'segmentoidD_REG_t1_MAX_ROID_2022_05_12__23_33_17_Subset-D';
paths.expname5 = 'segmentoidA_REG_t100_MAX_ROIA_2022_05_12__23_33_17_Subset-A';
paths.expname6 = 'segmentoidA_REG_t100_MAX_ROIA_2022_05_12__23_33_17_Subset-A';

paths.dataname1 = [paths.expname1,'-data'];
paths.dataname2 = [paths.expname2,'-data'];
paths.dataname3 = [paths.expname3,'-data'];
paths.dataname4 = [paths.expname4,'-data'];
paths.dataname5 = [paths.expname5,'-data'];
paths.dataname6 = [paths.expname6,'-data'];

runOffset = false;
if(runOffset)
    labelOffset = 'offset_';
else
    labelOffset = '';
end

paths.plotFolder = ['D:\BIO\PhD\ditalia\somitoid\github\ditalia-somitoid\code_revision2_alvin\figure_revision2_segmentoid_piv' labelOffset    '/'];
mkdir(paths.plotFolder);
mkdir([paths.plotFolder '/tests']);

% define pixelsize
dx = 0.692; % um
dz = 8; % um
dt = 6/60; % h-1

optsSave=[];
optsSave.compress='lzw';
optsSave.overwrite=true;

% roiStack = loadtiff([roipath]);
% saveastiff(s2,[paths.outFolder filesep st_dir(i).name(1:end)],optsSave);
colormap([0 0 0; brewermap(1000,'RdYlGn')]);

%tRange = [40 113; 70 120; 50 90]; %Ale all movie
%tRange = [40 113; 70 130; 50 130]; %Alvin

tRange = [20 70; 50 113; 35 100; 70 120; 1 75; 50 90];
DT = (diff(tRange').*dt);

ROIlabel = [1 2 1 2 1 2];

% process all segmentoids

plus = []; minus=[];
Mesp2Plus =[]; Mesp2Minus = [];

for n=1:6
    
    % load image and data
    segmentoid_proj = loadtiff([paths.directory paths.(['expname' num2str(n)]) '.tif']);
    load(['./cell_struct/piv_' labelOffset paths.(['expname' num2str(n)])  '.mat']);

    % post-process piv
    u=[]; v=[]; maskResize=[];
    filt_sizeV = 10;
    filtDiv = fspecial('gaussian',10,5);

    for i=tRange(n,1):tRange(n,2)
       
        uHere = u_original{i}; vHere = v_original{i};
        maski = imresize(mask(:,:,i),size(uHere));
        %maski = maski(2:end-1,2:end-1); %Alvin used this one, unsure why
         
        ui = nanconv(uHere,ones(filt_sizeV)./filt_sizeV,'same','noedge','nanout');
        vi = nanconv(vHere,ones(filt_sizeV)./filt_sizeV,'same','noedge','nanout');
        ui(~maski)=nan;
        vi(~maski)=nan;
        
        u=cat(3,u,ui);
        v=cat(3,v,vi);
        maskResize = cat(3,maskResize,maski);

    end
    uAll = nanmean(u,3); %undecided if we should do nansum here
    vAll = nanmean(v,3); %undecided if we should do nansum here
   
    maskAll = any(maskResize(:,:,end),3); %before "all"
    maskAll = imclose(maskAll,strel('disk',5));
    maskAll = bwmorph(maskAll,'majority');
    maskAll2 = bwareaopen(maskAll,1000);
    maskAllShrink = imerode(maskAll2,strel('disk',5));

    uAll(~maskAllShrink)=NaN;
    vAll(~maskAllShrink)=NaN;

    divNotNaN = divergence(x{1},y{1},uAll,vAll)./dt;
    div = divNotNaN;
    div(~maskAllShrink)=NaN;

    %smoothening with averaging filter
    divAv = nanconv(div,filtDiv,'same','nanout'); %AD changed from nanout

    %smoothening with csaps
    divAvForCsaps = nanconv(div,filtDiv,'same','noedge'); %AD changed from nanout
    divForCsaps = div; 
    divForCsaps(isnan(divForCsaps))= divAvForCsaps(isnan(divForCsaps));
    divForCsaps(isnan(divForCsaps))= divNotNaN(isnan(divForCsaps));

%     x = {1:size(div,1),1:size(div,2)};
%     [xx,yy] = ndgrid(x{1},x{2});
%     % y = peaks(xx, yy);
%     % figure
%     % surf(xx,yy,div)
%     % axis off
%     [divCsaps,p] = csaps(x,divForCsaps,0.01,x);
%     divCsaps(isnan(div))=NaN'

    %plotting for comparison Smoothing
    figure;
    colormap([0 0 0; brewermap(1000,'RdYlGn')]);

    tiledlayout(1,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
    
    nexttile
    imagesc([div],[-1 1]);
    title('No smooth');
    axis equal
    axis tight
    xticks('');   yticks('');
    %waitforbuttonpress

    nexttile
    imagesc([divAv],[-1 1]);
    title('Gauss Filter');
    axis equal
    axis tight
    xticks('');   yticks('');
    %waitforbuttonpress

    if(~isempty(paths.plotFolder))
        print([paths.plotFolder 'tests/' paths.(['expname' num2str(n)]) '_divSmooths'],'-r300','-dpng');   
        print([paths.plotFolder 'tests/' paths.(['expname' num2str(n)]) '_divSmooths'],'-depsc');   
    end

%     nexttile
%     imagesc([divCsaps],[-0.07 0.07]);
%     title('Csaps');
%     axis equal
%     axis tight
%     xticks('');   yticks('');
%     %waitforbuttonpress

     % plot overlay
    imgStart = segmentoid_proj(:,:,tRange(n,1));
    imgEnd = segmentoid_proj(:,:,tRange(n,2));
    imgEnd(~maskAllShrink)=NaN;

    % divergence mask
%     thres_otsu = graythresh(divAv);
%     mask_div = imbinarize(divAv,thres_otsu);
    thres = 0.15;
    mask_div = divAv>thres;
    mask_div_large = imresize(mask_div,size(imgEnd));
%     imshow(mask_div_large,[0,max(mask_div_large,[],'all')])

    % plot onto imgEnd
    f = figure();
    cm = brewermap(11,'RdYlGn');
    colorHere = cm(7,:);
%   colorHere = 'c';

    %imshow(imgEnd,[0,max(imgEnd,[],'all')/2]) %ALVIN
    imshow(imadjust(imgEnd,[0 0.3])); %ALVIN
    hold on;
    B = bwboundaries(mask_div_large(:,:));
    for b = 1:numel(B)
        boundary = B{b};
        plot(boundary(:,2), boundary(:,1), 'color',colorHere, 'LineWidth', 5);
    end
%     contour(imresize(divAv,size(imgEnd))','r')
    hold off;

    exportgraphics(f,[paths.plotFolder  paths.(['expname' num2str(n)]) '_overlay_' num2str(n) '.eps'],'BackgroundColor','none')
    exportgraphics(f,[paths.plotFolder  paths.(['expname' num2str(n)]) '_overlay_' num2str(n) '.png'],'BackgroundColor','none')
    
    % stats
    ROI1plus = loadtiff([paths.directory 'ROIsSegmSeparated/' paths.(['expname' num2str(n)])  '_ROI' num2str(ROIlabel(n)) '_Mesp2plus.tif']);
    ROI1minus = loadtiff([paths.directory 'ROIsSegmSeparated/' paths.(['expname' num2str(n)]) '_ROI' num2str(ROIlabel(n)) '_Mesp2minus.tif']);

    ROI1plusResized = round(imresize(ROI1plus,size(div),'nearest'))>0;
    ROI1minusResized  = round(imresize(ROI1minus,size(div),'nearest'))>0;

    ROI1plusResized = imerode(ROI1plusResized,strel('disk',1));
    ROI1minusResized = imerode(ROI1minusResized,strel('disk',1));

    imgStart = segmentoid_proj(:,:,tRange(n,1));
    imgEnd = segmentoid_proj(:,:,tRange(n,2));
    imgStartResized = imresize(imgStart,size(div));
    imgEndResized  =  imresize(imgEnd,size(div));

%     figure;
%     ROIs = zeros(size(ROI1plusResized));
%     ROIs(ROI1plusResized>0)=imgEndResized(ROI1plusResized>0);
%     ROIs(ROI1minusResized>0)=imgEndResized(ROI1minusResized>0);
%     imshow(ROIs,[],'InitialMagnification',300);

    Mesp2all = nanmean(imgEnd(ROI1plus|ROI1minus));
    Mesp2Plus = vertcat(Mesp2Plus,nanmean(imgEnd(ROI1plus>0))./Mesp2all-1);
    Mesp2Minus = vertcat(Mesp2Minus,nanmean(imgEnd(ROI1minus>0))./Mesp2all-1);

    plus = vertcat(plus,nanmean(div(ROI1plusResized )));
    minus = vertcat(minus,nanmean(div(ROI1minusResized )));

%     figure;
%     divs = nan(size(div));
%     divs(ROI1plusResized>0)=divAv(ROI1plusResized>0);
%     divs(ROI1minusResized>0)=divAv(ROI1minusResized>0);

    figure;
    colormap([0 0 0; brewermap(1000,'RdYlGn')]);
    divAvPlot = divAv; 
    divAvPlot(divAvPlot>1)=  1;
    divAvPlot(divAvPlot<-1)=-1;
  
    imagesc(divAv,[-1 1]);
    axis equal
    xlim1 = xlim;
    ylim1 = ylim;

    % velocity field
    spacing = 5;
    hold on;
    uAllPlot = uAll; 
    vAllPlot = vAll;
    uAllPlot(isnan(divAvPlot)) = NaN; 
    vAllPlot(isnan(divAvPlot)) = NaN;

    [xx,yy] = ndgrid(1:size(div,1),1:size(div,2));
    quiver(yy(1:spacing:end,1:spacing:end),xx(1:spacing:end,1:spacing:end),1.5*uAllPlot(1:spacing:end,1:spacing:end),1.5*vAllPlot(1:spacing:end,1:spacing:end),0,'b','LineWidth',1);
    hold off;
    axis tight
    axis off
    ax=gca;
    ax.Position=[0 0 1 1];
    set(gcf, 'InvertHardCopy', 'off');
    xlim(xlim1);
    ylim(ylim1);
    mygcf = gcf;

    set(gcf, 'InvertHardCopy', 'off');
    if(~isempty(paths.plotFolder))

        print([paths.plotFolder  paths.(['expname' num2str(n)]) '_divAV_' num2str(n)],'-r300','-dpng');   
        print([paths.plotFolder  paths.(['expname' num2str(n)]) '_divAv_' num2str(n)],'-depsc'); 

    end
    %% Generate csv dataframe for plot
    df_out = cell(size(divAv)+1);
    df_out{1,1} = 'Divergence Map';
    df_out{2,1} = 'y';
    df_out{1,2} = 'x';
    % [df_out{3:end,1}] = time_column_cell{:};
    df_out(2:end,2:end) = num2cell(divAv);
    if n==2
        writecell(df_out, '../excel_data/Fig4g_left.xls')
    elseif n==5
        writecell(df_out, '../excel_data/ExtFig9g_top.xls')
    elseif n==4
        writecell(df_out, '../excel_data/ExtFig9g_bottom.xls')
    end

end

% 
% figure;
% plot(Mesp2Plus,plus,'o')
% hold on;
% plot( Mesp2Minus,minus,'o')
% [h,p]=corr([Mesp2Minus(:) ;Mesp2Plus(:)],[minus(:); plus(:)],'Type','Pearson')

figure;
opts =[];
opts.color = [0 0 0];
% plotBarsPaper(1,minus([1 3 5]),opts);
% plotBarsPaper(1,minus([2 4 6]),opts);
% plotBarsPaper(2,plus([1 3 5]),opts);
% plotBarsPaper(2,plus([2 4 6]),opts);
opts.symbol ='o';
plotBarsPaper(1,minus(:),opts);
plotBarsPaper(2,plus(:),opts);

xlim([0.5 2.5]);
ylim([-0.3 0.3]);
ylabel('Divergence velocity field (h-1)');
xticks([1 2])
xticklabels({'Mesp2-' 'Mesp2+'})
box

standardizePlotAle(gcf,gca,[paths.plotFolder 'divStat']);
[~,p]=ttest2(minus,plus)
DT
%% Generate csv dataframe for plot
df_out = cell(size(minus,1)+1,2);
df_out{1,1} = 'MESP-';
df_out{1,2} = 'MESP+';
df_out(2:end,1) = num2cell(minus);
df_out(2:end,2) = num2cell(plus);


writecell(df_out, '../excel_data/Fig4g_right.xls')
%%

function myquiver(Xi,Yi,vx,vy,quiverStep)

   lmean1=3;
   K1 = ones(lmean1)./lmean1^2;
   myvx =  nanconv(vx,K1,'same','nanout');
   myvy =  nanconv(vy,K1,'same','nanout');
   %myvx =  imfilter(vx,K1);
   %myvy =  imfilter(vy,K1);
   lmean1=3;
   K1 = ones(lmean1)./lmean1^2;
   myvx =  nanconv(vx,K1,'same','nanout');
   myvy =  nanconv(vy,K1,'same','nanout');

   
  % quiver(imresize(Xi,1/quiverStep,'bilinear'),imresize(Yi,1/quiverStep,'bilinear'),imresize(vxMPlot,1/quiverStep,'bilinear').*scaleFactor.*quiverStep./dx,imresize(vyMPlot,1/quiverStep,'bilinear').*scaleFactor.*quiverStep./dy,0,'color','b','LineWidth',2);
   quiver(Xi(1:quiverStep:end),Yi(1:quiverStep:end),myvx(1:quiverStep:end).*quiverStep,myvy(1:quiverStep:end).*quiverStep,0,'color','b','LineWidth',2);

end

