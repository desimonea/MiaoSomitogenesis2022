addpath(genpath([pwd '/functions']));

%% process wt
paths=[];

paths.masterFolder = [pwd '/data/segmentoid/wt/'];
paths.inFolder = [paths.masterFolder 'raw/'];
paths.resultsFolder = [paths.masterFolder 'resultsKymo/'];
paths.objFolder = [paths.masterFolder 'objects/'];
paths.shiftFolder = [paths.masterFolder 'shifted/'];

paths.label1='Mesp2';
paths.label2='Hes7_Uncx';

dirs = dir([paths.inFolder '*' paths.label1 '*.tif']);
mkdir(paths.resultsFolder);
mkdir([paths.resultsFolder '/png/']);

% parameters
opts.x0 = 10;
opts.toFlip = [0];
opts.thresholdsP = [NaN];
opts.bwthresh = [0.7 0.7]/10;
opts.toSkh = [0 ];
opts.colorsKymo = [[1 0 1];[0 1 0]]; 
opts.rangeKymos = [NaN NaN NaN NaN; ]; %range1 range2
opts.rangeKymos = repmat([0 10 0 2],[1 1]);

optsOsc.label1 ='Mesp2';
optsOsc.label2 = 'Hes7';
optsOsc.verbose = 1;
optsOsc.print = 1;
optsOsc.distHes7 = 150; %where to calculate Hes7
optsOsc.tspanHours = 4.5; % for averaging oscillations
optsOsc.f0 = 1/4.5;
optsOsc.df = (0.1/4.5)*optsOsc.f0;
colors = lines(5);
optsOsc.colors = colors([5 4],:);
optsOsc.randomizeM = false;
optsOsc.simulated = false;
optsOsc.spacelag = 0;
opts.optsOsc = optsOsc;

myPSMAll=[];

for i = 1:size(dirs,1)
    
    disp(i);

    optsHere = opts;
    optsHere.x0 = opts.x0;
    optsHere.toFlip = opts.toFlip(i);
    optsHere.thresholdsP = opts.thresholdsP(i);
    optsHere.bwthresh = opts.bwthresh(i,:);
    optsHere.toSkh = opts.toSkh(i);
    optsHere.number = i;
    optsHere.rangeKymos = opts.rangeKymos(i,:);
    
    paths.name = dirs(i).name;
    myPSM = kymoElongation(paths,optsHere);
    myPSMAll{i} = myPSM;
    save([paths.objFolder strrep(dirs(i).name,'.tif','.mat')],'myPSM');

end

%we do not save myPSMAll here, because the provided file contains results
%from the entire dataset

%% process Hes7null
close all;
paths=[];

paths.masterFolder = [pwd '/data/segmentoid/Hes7null/']; 
paths.inFolder = [paths.masterFolder 'raw/'];
paths.resultsFolder = [paths.masterFolder 'resultsKymo/'];
paths.objFolder = [paths.masterFolder 'objects/'];
paths.shiftFolder = [paths.masterFolder 'shifted/'];

paths.label1='Mesp2';
paths.label2='pseudoHes7_Uncx';

dirs = dir([paths.inFolder '*' paths.label1 '*.tif']);
mkdir(paths.resultsFolder);
mkdir([paths.resultsFolder '/png/']);

% parameters
opts.x0 = 10;
opts.toFlip = [0 ];
opts.thresholdsP = [NaN ];
opts.bwthresh = [0.7 0.7]/10;
opts.toSkh = [0 ];
opts.colorsKymo = [0.7*[1 0 1];[0 1 0]];
opts.rangeKymos = repmat([0 10 0 2],[1 1]);

optsOsc.label1 ='Mesp2';
optsOsc.label2 = 'Hes7';
optsOsc.verbose = 1;
optsOsc.print = 1;
optsOsc.distHes7 = 150; %where to calculate Hes7
optsOsc.tspanHours = 4.5; % for averaging oscillations
optsOsc.f0 = 1/4.5;
optsOsc.df = (0.1/4.5)*optsOsc.f0;
colors = lines(5);
optsOsc.colors = colors([5 4],:);
optsOsc.randomizeM = false;
optsOsc.simulated = false;
optsOsc.spacelag = 0;
opts.optsOsc = optsOsc;

myPSMAllHes7null=[]
for i = 1:size(dirs,1)
    
    disp(i);

    optsHere = opts;
    optsHere.x0 = opts.x0;
    optsHere.toFlip = opts.toFlip(i);
    optsHere.thresholdsP = opts.thresholdsP(i);
    optsHere.bwthresh = opts.bwthresh(i,:);
    optsHere.toSkh = opts.toSkh(i);
    optsHere.number = i;
    optsHere.rangeKymos = opts.rangeKymos(i,:);

    paths.name = dirs(i).name;
    myPSM = kymoElongation(paths,optsHere);
    myPSMAllHes7null{i} = myPSM; 
    save([paths.objFolder strrep(dirs(i).name,'.tif','.mat')],'myPSM');

end

%we do not save myPSMAll here, because the provided file contains results
%from the entire dataset

%% plots autocorr in wild-type and Hes7 

% we read the entire dataset, not only the test version
load([pwd '/data/segmentoid/wt/objects/PSMkymoall.mat']);
load([pwd '/data/segmentoid/Hes7null/objects/PSMkymoall.mat']);

paths.resultFolder = [pwd '/data/segmentoid/freq/'];
mkdir(paths.resultFolder);

colors = lines(5);

figure;
ccm = []; cch = [];
for i=1:numel(myPSMAll)
    ccm = horzcat(ccm,[myPSMAll{i}.xcorrm(1:100,2)]);
    cch = horzcat(cch,[myPSMAll{i}.xcorrh(1:100,2)]);
    lags = horzcat(myPSMAll{i}.xcorrh(1:100,1));
end

T4c = table;

ccmMean = mean(ccm,2);

plot(lags,mean(cch,2),'-','DisplayName','Hes7','LineWidth',3,'Color',colors(5,:));hold on;
plot(lags,mean(ccm,2),'-','DisplayName','Mesp2','LineWidth',3,'Color',colors(4,:));hold on;

T4c = [T4c table(lags,'VariableNames',{'Lags'})];
T4c = [T4c table(mean(cch,2),'VariableNames',{'Auto-corr Hes7  wt'})];
T4c = [T4c table(mean(ccm,2),'VariableNames',{'Auto-corr Mesp2 wt'})];

ccmHN = []; cchHN = [];
for i=1:numel(myPSMAllHes7null)
    ccmHN = horzcat(ccmHN,[myPSMAllHes7null{i}.xcorrm(1:100,2)]);
    cchHN = horzcat(cchHN,[myPSMAllHes7null{i}.xcorrh(1:100,2)]);
end

plot(lags,mean(cchHN,2),'--','DisplayName','pseudoHes7','LineWidth',3,'Color',colors(5,:));hold on;
plot(lags,mean(ccmHN,2),'--','DisplayName','Mesp2 (Hes7 -/-)','LineWidth',3,'Color',colors(4,:));hold on;

T4c = [T4c table(mean(cchHN,2),'VariableNames',{'Auto-corr pseudoHes7  hes7 -/-'})];
T4c = [T4c table(mean(ccmHN,2),'VariableNames',{'Auto-corr Mesp2 hes7 -/-'})];

set(gca,'fontname','arial','FontSize',24,'LineWidth',3);
xlabel('Time period (h)','fontname','arial','FontSize',24);
ylabel('Time auto-correlation','fontname','arial','FontSize',24);
legend('Show','Location','SouthEast');
legend boxoff
ylim([-1.1 1.1])

[~,locs]=findpeaks(ccmMean);
xp = lags(locs(1));
yp = ccmMean(locs(1));
idxFit = lags>4 & lags<6;
fPeakCCM = fit(lags(idxFit),ccmMean(idxFit),'a*(x-x0).^2+b');
%plot(xp,1, 'v', 'Color',colors(5,:),'LineWidth',2, 'MarkerSize', 10); % Plot triangles atop the peaks.
%hold on; plot(fPeakCCM);
phaseCC = yp;
stdx0M = diff(confint(fPeakCCM,0.68))/2;
%title('Average auto-correlation');
%disp(['Freq Mesp2 in wt is ' num2str(fPeakCCM.x0) ' pm ' num2str(stdx0M(3))]);
pbaspect([1 1 1]);

mkdir([paths.resultFolder '/png/']);
print([paths.resultFolder 'autoCAverages'],'-depsc');
print([paths.resultFolder '/png/' 'autoCAverages'],'-dpng');

% print source data
fout = [paths.resultFolder 'sourceKymoAutoCorr_Wt_Hes7.txt'];
FID = fopen(fout, 'w');
writetable(T4c, fout,'WriteVariableNames',true,'WriteRowNames',false,'Delimiter','tab');
fclose(FID);


%% calculate base frequency in wt Hes7 - run after previous - USED FOR PAPER

T_8a = table;
T_8a = [T_8a table(lags,'VariableNames',{'Lags (h)'})];

figure;
Th = [];
for i=1:size(cch,2)
    cchHere = cch(:,i);
    plot(lags,cchHere,'-','LineWidth',2,'Color',colors(5,:));hold on;
    T_8a = [T_8a table(cchHere,'VariableNames',{['Auto-corr Hes7 ' num2str(i)]})];
    [~,locs]=findpeaks(cchHere);
    xp = lags(locs(1));
    yp = cchHere(locs(1));
    plot(xp,1, 'v', 'Color',colors(5,:),'LineWidth',2, 'MarkerSize', 10); % Plot triangles atop the peaks.
    Th(i) = xp;
end
set(gca,'FontSize',24,'LineWidth',3,'fontname','arial');
xlabel('Time period (h)','fontname','arial');
ylabel('Time auto-correlation','fontname','arial');
ylim([-1.1 1.1])
legend({'Hes7 (wt)'},'Location','SouthEast'); legend boxoff;
%title('Freq Hes7 individual movies');
disp(['Freq Hes7 in wt is ' num2str(nanmean(Th)) ' pm ' num2str(nanstd(Th))]);
pbaspect([1 1 1]);

print([paths.resultFolder 'autoCindividualHes7Wt'],'-depsc');
print([paths.resultFolder '/png/' 'autoCindividualHes7Wt'],'-dpng');

% print source data 8a
fout = [paths.resultFolder 'sourceKymoAutoCorr_Individual_WT_Hes7_Fig8a.txt'];
FID = fopen(fout, 'w');
if FID == -1, error('Cannot open file %s',fout); end
writetable(T_8a, fout,'WriteVariableNames',true,'WriteRowNames',false,'Delimiter','tab');
fclose(FID);

%% calculate base frequency in wt Mesp2 - run after previous - USED FOR PAPER

T_8b = table;
T_8b = [T_8b table(lags,'VariableNames',{'Lags (h)'})];

figure;
Tm = [];
for i=1:size(ccm,2)
    ccmHere = smooth(ccm(:,i),10);
    plot(lags,ccmHere,'-','LineWidth',2,'Color',colors(4,:)); hold on;
    T_8b = [T_8b table(ccmHere,'VariableNames',{['Auto-corr Mesp2 ' num2str(i)]})];

    [~,locs]=findpeaks(ccmHere);
    if(isempty(locs))
        Tm(i) = NaN;
    else
        xp = lags(locs(1));
        yp = ccmHere(locs(1));
        plot(xp,1, 'v', 'Color',colors(4,:),'LineWidth',2, 'MarkerSize', 10); % Plot triangles atop the peaks.
        Tm(i) = xp;
    end
end
set(gca,'FontSize',24,'LineWidth',3);
xlabel('Time period (h)','fontname','arial');
ylabel('Time auto-correlation','fontname','arial');
ylim([-1.1 1.1])
legend({'Mesp2 (wt)'},'Location','SouthEast'); legend boxoff;
title('Freq Mesp2 individual movies');
disp(['Freq Mesp2 in wt is ' num2str(nanmean(Tm)) ' pm ' num2str(nanstd(Tm))]);
pbaspect([1 1 1]);

print([paths.resultFolder 'autoCindividualMesp2Wt'],'-depsc');
print([paths.resultFolder '/png/' 'autoCindividualMesp2Wt'],'-dpng');

% print source data 8b
fout = [paths.resultFolder 'sourceKymoAutoCorr_Individual_WT_Mesp2_Fig8b.txt'];
FID = fopen(fout, 'w');
if FID == -1, error('Cannot open file %s',fout); end
writetable(T_8b, fout,'WriteVariableNames',true,'WriteRowNames',false,'Delimiter','tab');
fclose(FID);
