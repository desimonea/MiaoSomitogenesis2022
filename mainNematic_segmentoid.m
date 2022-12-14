addpath(genpath([pwd '/functions']));

%% calculate ROIs on shiftet for nematic - set 1

% this reads objects that must be calculated with mainElongationKymos

paths =[];
paths.masterFolder = [pwd '/data/segmentoid/wt/'];% run also for Hes7null
%paths.masterFolder = [pwd '/data/segmentoid/Hes7null/'];% run also for Hes7null
paths.inFolder = [paths.masterFolder 'shifted/'];
paths.objFolder = [paths.masterFolder 'objects/'];
dirs  = dir([paths.inFolder '*Mesp2*interval.tif']);

for n = 1:numel(dirs)
    
    disp(n);

    paths.basename = dirs(n).name;
    paths.savename = replace(dirs(n).name,'Mesp2','ROI');

    img = loadtiff([paths.inFolder paths.basename]);
    load([paths.objFolder strrep(paths.basename,'.tif','.mat')]);
    
    dx = myPSM.dx;
    
    img1 = imgaussfilt(img,round(7/dx)); 
    thresh = graythresh(img1(img1>0))*255;
    bw= img1>thresh;
    bw2= imclose(bw,strel('disk',15));
    bw4 = bw2;
    
    for i=1:size(img,3)
        bw3Here = bwmorph(imfill(bw2(:,:,i),'holes'),'majority');
        bw4(:,:,i) = bwareaopen(bw3Here,round(1900/(dx*dx)),4);
    end    
       
    projROI = nanmean(bw4,3);
    th=graythresh(projROI);
    xMax = find(squeeze(sum(projROI,1))>th,1,'last')-round(70/dx);
       
    if(contains(paths.masterFolder,'wt'))
        imgHU = loadtiff([paths.inFolder replace(paths.basename,'Mesp2','Hes7_Uncx')]);
    elseif(contains(paths.masterFolder,'Hes7null'))
        imgHU = loadtiff([paths.inFolder replace(paths.basename,'Mesp2','pseudoHes7_Uncx')]);
    else
         error('Wrong filename');
    end
     
    s1 = img; 
    s2= imgHU;
    
    tot1 = sum(sum(s1,1),2)./(size(s1,1)*size(s1,2));
    tot2 = sum(sum(s2,1),2)./(size(s1,1)*size(s1,2));
    m1 = repmat(tot1,[size(s1,1) size(s1,2)]);
    m2 = repmat(tot2,[size(s2,1) size(s2,2)]);
    
    img = ((double(s1)./m1+double(s2)./m2).*(m1+m2)./2);
    img(:,xMax:end,:) = 0;
    
    img1 = uint8(imgaussfilt(img,5));
    thresh = graythresh((img1(img1>0)))*255;
    bw= img1>thresh;
    bw2= imclose(bw,strel('disk',15));
    
    for i=1:size(img,3)
        bw3Here = bwmorph(imfill(bw2(:,:,i),'holes'),'majority');
        bw4(:,:,i) = bwareaopen(bw3Here,1000,4);
    end   
    bw4(:,xMax:end,:) = 0;

    saveOptions.compress = 'lzw';
    saveOptions.overwrite = true;
    saveastiff(uint8(bw4),[paths.inFolder paths.savename],saveOptions);
end

%% calculate nematic order wt 

colors = cool(100);

paths =[];
paths.masterFolder = [pwd '/data/segmentoid/wt/'];

paths.inFolder = [paths.masterFolder 'shifted/'];
paths.resFolder = [paths.masterFolder 'Qresults/'];
paths.metaFolder = [paths.masterFolder 'raw/'];
paths.objFolder = [paths.masterFolder 'objects/'];
paths.objexpFolder = [paths.masterFolder 'objects/'];

dirs  = dir([paths.inFolder '*Mesp2*interval.tif']);
mkdir(paths.resFolder);
mkdir([paths.resFolder '/png/']);
mkdir(paths.objFolder);

opts = [];
opts.verbose = true;
opts.dxStepNematicUm = 14;
opts.LNematicWindowUm = 124; 
opts.rhoRangeUm = [0.01445 Inf]; 
opts.minFractionInRoi = 0.9;
opts.printStepSizeH = 5; %in h
opts.xRange = [-543 -181]; %in px
opts.tBins = 0:10:50; %in h
opts.toPrint = true;

QWT=[];
tExactWT=[];
QmeanTbinnedWT =[];

for i = 1:numel(dirs)

    load([paths.objexpFolder strrep(dirs(i).name,'.tif','.mat')]);
    paths.name1 = dirs(i).name;
    paths.name2 = replace(paths.name1,'Mesp2','Hes7_Uncx');
    paths.nameROI = replace(paths.name1,'Mesp2','ROI');

    % comment if do not want to to recalculate
    myPSM = nematicElongation(paths,opts,myPSM);
    save([paths.objFolder strrep(dirs(i).name,'.tif','.mat')],'myPSM');
    
    QWT{i}= myPSM.nematic.Q;
    tExactWT{i} = myPSM.nematic.tExact;
    QmeanTbinnedWT = cat(2,QmeanTbinnedWT,myPSM.nematic.QmeanTbinned(:));
    ylim([0 0.5]);
    myPSM.dx
end

%% calculate nematic order Hes7null

paths =[];
paths.masterFolder = [pwd '/data/segmentoid/Hes7null/'];

paths.objexpFolder = [paths.masterFolder 'objects/'];
                 
paths.inFolder = [paths.masterFolder 'shifted/'];
paths.resFolder = [paths.masterFolder 'Qresults/'];
paths.metaFolder = [paths.masterFolder 'raw/'];
paths.objFolder = [paths.masterFolder 'objects/'];
dirs  = dir([paths.inFolder '*Mesp2*interval.tif']);
mkdir(paths.resFolder);
mkdir([paths.resFolder '/png/']);
mkdir(paths.objFolder);

%for opts read above to not risk them being different

QHN =[];
tExactHN =[];
QmeanTbinnedHN = [];

for i = 1:numel(dirs)

    load([paths.objexpFolder strrep(dirs(i).name,'.tif','.mat')]);
    paths.name1 = dirs(i).name;
    paths.name2 = replace(paths.name1,'Mesp2','pseudoHes7_Uncx'); %pseudoHes7
    paths.nameROI = replace(paths.name1,'Mesp2','ROI');
      
    % comment if do not want to to recalculate
    myPSM = nematicElongation(paths,opts,myPSM);
    save([paths.objFolder strrep(dirs(i).name,'.tif','.mat')],'myPSM');
        
    QHN{i}= myPSM.nematic.Q;
    tExactHN{i} = myPSM.nematic.tExact;
    QmeanTbinnedHN = cat(2,QmeanTbinnedHN,myPSM.nematic.QmeanTbinned(:));
end


%% plot average time points

% Here we print nematic order for the sample data
% We read the entire dataset, not only the test data

paths.resFolder = [pwd '/data/segmentoid/Qresults/'];
mkdir(paths.resFolder);
mkdir([paths.resFolder '/png/']);

load([pwd '/data/segmentoid/wt/objects/QAll.mat']);
load([pwd '/data/segmentoid/Hes7null/objects/QAll.mat'])

figure;
opts.symbol = 'o';
colors = lines(2);

opts = [];
opts.meanWidth = 1;
opts.CapSize = 10;
opts.MarkerSize=10;

%tBins need to be the same above!
tBins =  0:10:50; %in h
dtBins = tBins(2)-tBins(1);

RStest = nan(numel(tBins,1));
Ttest  = nan(numel(tBins,1));


for idxT=1:numel(tBins)
    
    opts.color = colors(1,:);
    plotBarsPaper(tBins(idxT)-dtBins/8+dtBins/2,QmeanTbinnedWT(idxT,:),opts);
    
    opts.color = colors(2,:);
    plotBarsPaper(tBins(idxT)+dtBins/8+dtBins/2,QmeanTbinnedHN(idxT,:),opts);

    if(all(isnan(QmeanTbinnedWT(idxT,:)))| all(isnan(QmeanTbinnedHN(idxT,:))))
        RStest(idxT) =  NaN;
        Ttest(idxT) =   NaN;
        logNtest(idxT) = NaN;
    else
        RStest(idxT)    =  ranksum(QmeanTbinnedWT(idxT,:),QmeanTbinnedHN(idxT,:));
        [~,Ttest(idxT)] =  ttest2(QmeanTbinnedWT(idxT,:),QmeanTbinnedHN(idxT,:));
        [~,logNtest(idxT)] =  ttest2(log(QmeanTbinnedWT(idxT,:)),log(QmeanTbinnedHN(idxT,:)));
    end
    
    text(tBins(idxT)-2+dtBins/2,0.08,num2str(RStest(idxT),['%.' num2str(floor(abs(log10(RStest(idxT)))))+1 'f']),'FontSize',24);
end

set(gca,'box','on','FontSize',24,'LineWidth',3)
ylabel('Condition')
xtickangle(0)
axis([0 50 0 0.09]);
xticks(5:10:50);
xticklabels({'5','15','25','35','45'});
xlabel('Time (h)');
ylabel('Nematic order Q');

%legend({'','control','','','Hes7 -/-'},'Location','southwest','box','off');

pbaspect([1 1 1]);

print([paths.resFolder 'QmeanTbinned'],'-depsc','-loose','-painters');
print([paths.resFolder  '/png/' 'QmeanTbinned'],'-dpng','-loose','-painters');