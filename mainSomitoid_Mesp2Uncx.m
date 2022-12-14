%% calculate ROI1 - set 1

setnames = {'sub14'};


for n=1:numel(setnames)
    paths.masterFolder = [pwd '/data/somitoid/mesp2uncx/'];
    paths.inFolder = [paths.masterFolder 'raw/'];
    paths.basename = ['Mesp2_24pt5hr start_20min interval_2020_11_20__15_34_56_' setnames{n}];
    paths.savename = strrep(paths.basename,'Mesp2','ROI');
    paths.objFolder = [paths.masterFolder 'objects/'];

    img = imread([paths.inFolder paths.basename '.tif'],1);

    img1 = imgaussfilt(img,5);
    bw=imbinarize(img1);
    bw2= imclose(bw,strel('disk',10));
    bw3 = bwmorph(imfill(bw2,'holes'),'majority');
    bw4 = bwareaopen(bw3,1000);

    erosionROI = 50; %in px
    warning('Erosion is in px');
    
    roi = imerode(bw4,strel('disk',erosionROI));

    saveOptions.compress = 'lzw';
    saveOptions.overwrite = true;
    saveastiff(uint8(roi),[paths.inFolder strrep(paths.savename,'T1ForROI','14Test') '.tif'],saveOptions); %for compatibility
    imshowpair(img,roi)
end

%% merge mesp2 and uncx

setnames = {'sub14'};

for n=1:numel(setnames)
    paths=[];
    paths.masterFolder = [pwd '/data/somitoid/mesp2uncx/'];
    paths.inFolder = [paths.masterFolder 'raw/'];
    
    stdir=dir([paths.inFolder  'Uncx*' setnames{n} '*']);
    paths.basename1 = strrep(stdir(1).name,'.tif','');
    paths.basename2 = strrep(paths.basename1,'Uncx','Mesp2');
    paths.basenameMerged = strrep(paths.basename1,'Uncx','Merged');
    paths.basenameROI= strrep(paths.basename1,'Uncx','ROI');

    paths.objFolder = [paths.masterFolder 'objects/'];

    s1 = loadtiff([paths.inFolder paths.basename1 '.tif']);
    s2 = loadtiff([paths.inFolder paths.basename2 '.tif']);

    roi = loadtiff([paths.inFolder paths.basenameROI '.tif'])>0;
    s1(~roi)=NaN;
    s2(~roi)=NaN;
    
    tot1 = nansum(nansum(s1,1),2)./nansum(nansum(roi,2),1);
    tot2 = nansum(nansum(s2,1),2)./nansum(nansum(roi,2),1);
    %plot(squeeze(tot1)); hold on; plot(squeeze(tot2));

    m1 = repmat(tot1,[size(s1,1) size(s1,2)]);
    m2 = repmat(tot2,[size(s2,1) size(s2,2)]);

    s3 = ((double(s1)./m1-double(s2)./m2).*(m1+m2)./2)+(2^16-1)./2;

    saveOptions.compress = 'lzw';
    saveOptions.overwrite = true;
    saveastiff(uint16(s3),[paths.inFolder paths.basenameMerged '.tif'],saveOptions);
end

%% process all sets

setnames = {'sub14'}; 

for n=1:numel(setnames)
    % create somitoid
    paths=[];
    paths.masterFolder = 'data/somitoid/mesp2uncx/';
    paths.inFolder = [paths.masterFolder 'raw/'];
    paths.basename = ['Mesp2_24pt5hr start_20min interval_2020_11_20__15_34_56_' setnames{n}];
    paths.ROIname =  strrep(paths.basename,'Mesp2','ROI');
    paths.resFolder = [paths.masterFolder 'plotSpCorr/'];
    paths.objFolder = [paths.masterFolder 'objects/'];
    mkdir(paths.resFolder);
    mkdir([paths.resFolder '/png/']);
    mkdir(paths.objFolder);

    [myPSM0]=createSomitoid(paths);

    % process somitoid
    opts=[];
    opts.tStep= 1;
    opts.tStart = round((24.5-myPSM0.t0)./myPSM0.dt)+1; 
    opts.tEnd = round(90./myPSM0.dt-myPSM0.t0./myPSM0.dt)+1; 
    opts.tEnd = myPSM0.tMax;
    opts.tStepPrint = round(6/myPSM0.dt); %hours/dt
    opts.tStepPrintFine = round(3/myPSM0.dt); %hours/dt
    opts.tEndFine = round((50-myPSM0.t0)./myPSM0.dt)+1;
    opts.tPrint =  [opts.tStart:opts.tStepPrintFine:opts.tEndFine opts.tEndFine:opts.tStepPrint:opts.tEnd]; 
    opts.barPlotLim = [20 100];

    opts.resFactor = 5.5/myPSM0.dx;
    opts.dStep = 15; % 15 default in um
    opts.xLimUm = [0 400];
    opts.xMaxAutocorrUm = 200;
    opts.gaussFiltParams = [25 50 280 280]; %in um, first two are smoothening, the other two for bk removal
    %opts.gaussFiltParams = [2 4 280 280]; %in um, first two are smoothening, the other two for bk removal

    opts.noGca = 1;
    opts.sizeMultipliers = [1.2 1.1];

    opts.verbose = true;
    opts.toPrint = true;

    opts.fractionMinMin = 0.95;
    opts.valueMinMin = -0.02;
    opts.valueMinStd = 2;
% 
    paths.basename = strrep(paths.ROIname,'ROI','Merged');
    opts.label = 'Merged';
    opts.labelPrint = 'Mesp2 and Uncx';
    [myPSM,figs]=spCorrSomitoid(myPSM0,paths,opts);
    save([paths.objFolder paths.basename],'myPSM');
    myPSM_merged{n} = myPSM;
% 
    paths.basename = strrep(paths.ROIname,'ROI','Mesp2');
    opts.label = 'Mesp2';
    opts.labelPrint = 'Mesp2';
    [myPSM,figs]=spCorrSomitoid(myPSM0,paths,opts);
    save([paths.objFolder paths.basename],'myPSM');
    myPSM_mesp2{n} = myPSM;

    paths.basename = strrep(paths.ROIname,'ROI','Uncx');
    opts.label = 'Uncx';
    opts.labelPrint = 'Uncx';
    [myPSM,figs]=spCorrSomitoid(myPSM0,paths,opts);
    save([paths.objFolder paths.basename],'myPSM');
    myPSM_uncx_{n} = myPSM;
end

%% recalculate and plot correlations with tBin

setnames = {'sub14'}; 
labels = {'Merged','Mesp2','Uncx'};
labelsPrint = {'Mesp2 and Uncx','Mesp2','Uncx'};
colors = lines(10);
colorsPrint = [0.5 0.5 0.5  ;colors([4 5],:)];

paths=[];
paths.masterFolder = 'data/somitoid/mesp2uncx/';
paths.objFolder = [paths.masterFolder 'objects/'];
paths.resFolder = [paths.masterFolder 'sprCorrPlotTbin/'];

mkdir(paths.resFolder);
mkdir([paths.resFolder '/png/']);
mkdir(paths.objFolder);

opts=[];
opts.barPlotLim = [20 100];
opts.posyPlotLim = [80 160];
opts.tBinH = 3;

opts.fractionMinMin = 0.95;
opts.valueMinMin = -0.02;
opts.valueMinStd = 2.5;

opts.verbose = true;
opts.toPrint = true;
opts.xLimUm = [0 400];
opts.xMaxAutocorrUm = 200;

for n=1:numel(setnames)    
    
    % process somitoids
    for i = 1:numel(labels)
        stdir=dir([paths.objFolder labels{i} '*' setnames{n} '.mat']);
        load([paths.objFolder stdir(1).name]);
        opts.tStep = 1;
        opts.tStart = round((24.5-myPSM.t0)./myPSM.dt)+1; %70 - in frames
        opts.tEnd = round(90./myPSM.dt-myPSM.t0./myPSM.dt)+1; %to be checked
        opts.tStepPrint = round(6/myPSM.dt); %hours/dt
        opts.tStepPrintFine = round(3/myPSM.dt); %hours/dt
        opts.tEndFine = round((50-myPSM.t0)./myPSM.dt)+1;
        opts.tPrint =  [opts.tStart:opts.tStepPrintFine:opts.tEndFine];
        opts.tPrint =  [opts.tPrint opts.tPrint(end)+opts.tStepPrint:opts.tStepPrint:opts.tEnd]; 
        myPSM.label = labels{i};
        myPSM.labelPrint = labelsPrint{i};
        myPSM.basename = stdir(1).name;
        myPSM=spCorrTBinRecal(myPSM,paths,opts);
        myPSM.colorPrint = colorsPrint(i,:);
        save([paths.objFolder myPSM.basename],'myPSM');
        myPSM_merged{n} = myPSM;
    end
end