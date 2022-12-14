function [myPSM,opts] = kymoElongation(paths,opts,myPSM)

    if(nargin<3)
        myPSM = [];
    end
   
    x0=opts.x0; %space between posterior and end image after shifting
    toFlip=opts.toFlip; % does the organoid need to be flipped
    toSkh=opts.toSkh;
    thresholdsP=opts.thresholdsP; %threshold for posterior (if NaN it is calculated)
    bwthresh=opts.bwthresh;
    
    % colormap range for kymos
    if(isfield(opts,'rangeKymos'))
        rangeKymos = opts.rangeKymos;
    else
        rangeKymos = nan(1,4);
    end
    
    %labels
    label1 = replace(paths.label1,'_','_');
    label2=  replace(paths.label2,'_','_');
    
    % read stacks
    sHU = loadtiff([paths.inFolder replace(paths.name, label1,label2)]);
    sM = loadtiff([paths.inFolder paths.name]);

    % projections
    pHU = max(sHU,[],3);
    pM =  max(sM ,[],3);

    % calculate dx
    tMax = size(sHU,3);
    info=imfinfo([paths.inFolder paths.name]);
    dx=1./info(1).XResolution;
    dy=1./info(1).YResolution;

    % parse dt
    fname = strrep(info(1).Filename,'pt','.');
    delimiters = regexp(fname,'_');
    posMin =     regexp(fname,'min');
    delimitersMin = delimiters(find(delimiters<posMin,1,'last'));
    dt = str2num(fname(delimitersMin+1:posMin-1))/60; %in h
 
    % movie internal label
    i = opts.number;
    
    %% ROI calculation
    
    [XXimg,YYimg] = meshgrid(1:size(pM,2),1:size(pM,1));
    
    % smoothening and thresholding (just for ROI calculation)
    pHU2 = imgaussfilt(pHU,5,'FilterSize',11);
    pM2 =  imgaussfilt(pM,5,'FilterSize',11);

    roiHU=imbinarize(pHU2,bwthresh(1));
    roiM = imbinarize(pM2,bwthresh(2));

    figure
    imshowpair(pHU,roiHU);

    figure
    imshowpair(pM,roiM*0.2);

    % sum the two stacks and morph operations
    tot = logical(roiHU+roiM);
    warning('Morphological operations are not adapted for different px size - 1.38um standard');
    tot = imclose(tot, strel('disk',5)); %to do: adapt these 
    tot = bwmorph(tot,'fill');
    tot = imopen(tot,strel('disk',20));

    % select biggest object
    cc = bwconncomp(tot);
    pp=regionprops(cc,'Orientation','Area');
    [~,idxMaxArea]=max([pp.Area]);
    th1= pp(idxMaxArea).Orientation;
    tot = zeros(size(tot),'logical');
    tot(cc.PixelIdxList{idxMaxArea})=1;
    
    %plot ROI
    figure;
    imshow(pM2+pHU2,[0 50]); hold on;
    perimTot = bwperim(tot);
    plot(XXimg(perimTot(:)),YYimg(perimTot(:)),'o');
   
    % rotate first using main axis, then using fit. Calculate total rotation
    % angle and apply
    tot2 = imrotate(tot,-th1);
    
    [XX,YY] = meshgrid(1:size(tot2,2),1:size(tot2,1));
    fLine = fit(XX(tot2(:)),YY(tot2(:)),'poly1');
    th2 = atan(fLine.p1);
    th = -th1+th2;
   
    % flip if needed
     if(toFlip==1)
        disp('This stack is rotated');
        totRot = imrotate(tot,th+180);
        pHUrot = imrotate(pHU2,th+180);
        pMrot = imrotate(pM2,th+180);
        sHUrot = imrotate(sHU,th+180);
        sMrot =  imrotate(sM,th+180);
     else
        totRot = imrotate(tot,th);
        pHUrot = imrotate(pHU2,th);
        pMrot = imrotate(pM2,th);
        sHUrot = imrotate(sHU,th);
        sMrot =  imrotate(sM,th);
        perimTotRot = imrotate(perimTot,th);
     end

    figure;
    imshow(totRot); 
        
    % decide if to use median line or skeleton
    if(toSkh)
        % create skeleton
        skh=bwskel(totRot);
        [row,col]=ind2sub(size(skh),find(bwmorph(skh, 'endpoints')));
        [~,idxMax]=max(pdist([row col]));
        combos = nchoosek(1:numel(row),2);
        maxx = combos(idxMax,:);
        D1 = bwdistgeodesic(skh, col(maxx(1)), row(maxx(1)), 'quasi-euclidean');
        D2 = bwdistgeodesic(skh, col(maxx(2)), row(maxx(2)), 'quasi-euclidean');
        D = D1 + D2;
        D = round(D * 8) / 8;
        D(isnan(D)) = inf;
        skeleton_path = imregionalmin(D);
        xsk = 1:size(pMrot,2);
        ysk=   sum(double(skeleton_path).*YY,1)./sum(double(skeleton_path),1);
    else
        % take midline
        [XX,YY] = meshgrid(1:size(pMrot,2),1:size(pMrot,1));
        xsk = 1:size(pMrot,2);
        ysk=   sum(double(totRot).*YY,1)./sum(double(totRot),1);
    end
    
    % smoothen skeleton
    ysk2 = movmean(ysk,50); % not adapted to 
    ysk2(isnan(ysk2(:)))=ysk(isnan(ysk2(:)));
    idxs=~isnan(ysk2);
    xsk = xsk(idxs);
    ysk2 =ysk2(idxs);

    % plot skeleton
    figure;
    imshow(pMrot+pHUrot,[]); hold on;
    imshow(totRot);
    plot(xsk,ysk2);
    
    % rename for simplicity
    roiC = totRot;
    MC = pMrot;
    HUC = pHUrot;

    % mask skeleton
    skForAv = zeros(size(roiC));
    idxSk = sub2ind([size(roiC)],round(ysk2),round(xsk));
    skForAv(idxSk) = 1;

    % dilate mask skeleton for Hes7 
    skForAv = imdilate(skForAv,strel('line',75,90)); %75 %thickness stripe
    skForAv(~roiC(:))= 0; 
    
    % show mask skeleton
    figure;
    imshowpair(roiC,skForAv);
    
    % dilate mask skeleton for Hes7 
    skForAvM = imdilate(skForAv,strel('line',75,90)); %75 %thickness stripe
    skForAvM(~roiC(:))= 0;
    
    % make stack out of img
    skForAv = repmat(skForAv,[1 1 size(sHU,3)]);
    skForAvM = repmat(skForAvM,[1 1 size(sHU,3)]);

    % select skeleton regions
    sHU3 = sHUrot; sHU3(~skForAv)=NaN;
    sM3 =   sMrot; sM3(~skForAvM)=NaN;

    % create kimos and smoothen in space
    kimoHU = nanmean(sHU3,1); kimoHU= permute(kimoHU,[3 2 1]);
    kimoM  = nanmean(sM3,1);   kimoM= permute(kimoM,[3 2 1]);
    kimoHUf = imfilter(kimoHU,ones(1,50)./50); % this is not adapted for different px sizes
    kimoMf = imfilter(kimoM,ones(1,50)./50);
    
    %% calculate positions posterior and shift kimos
    
    if(isnan(thresholdsP))
        [~,idxM]=max(kimoHUf(1,:));
        threshold = nanmean(kimoHUf(1,idxM+100:idxM+200))*1.2;
    else
        threshold = thresholdsP;
    end
    
    % calculation xEdge
    xEdge= nan(tMax,1);
    for r=1:tMax
        idx = find(kimoHUf(r,:)>threshold,1,'last');
        if(~isempty(idx))
            xEdge(r)= idx;
        end
    end
    xEdge = round(smooth(xEdge,9/dt));
    
    %colormaps 
    if(isnan(rangeKymos(1)))
        rangeKymos(1) = quantile(kimoMf(:),0);
    end
    
    if(isnan(rangeKymos(2)))
        rangeKymos(2) = quantile(kimoMf(:),1);
    end
    
    if(isnan(rangeKymos(3)))
        rangeKymos(3) = quantile(kimoHUf(:),0);
    end
    
    if(isnan(rangeKymos(4)))
        rangeKymos(4) = quantile(kimoHUf(:),1);
    end
    opts.rangeKymos = rangeKymos;
    
    %% plot channel2 without alignment

    N = size(kimoHUf,2);
    
    figure
    imagesc(kimoHUf,rangeKymos(3:4)); hold on; %,[0.07 0.7]
    plot(xEdge,(1:tMax)','-ow','LineWidth',3);
    c=colorbar;
    c.Label.FontSize = 24;
    c.Label.FontName = 'Arial';    
    c.Label.String = strrep([label2 ' (a.u.)'] ,'_',' ');    
    c.LineWidth= 3;

    % adapt tick size
    mygca = gca;
    xt = mygca.XTick;    % Original 'XTick' Values
    realTicks = (0:500:(N.*dx)); realTicks = sort(realTicks);
    dLX = mygca.XLim(2)-mygca.XLim(1);
    posTicks = dLX*(realTicks/dx)./N+mygca.XLim(1);
    set(gca, 'XTick',posTicks, 'XTickLabel',realTicks, 'XTickLabelRotation',0)   % Label Ticks
  
    yt = mygca.YTick;    % Original 'XTick' Values
    realTicks = 0:10:(tMax.*dt);
    dLY = mygca.YLim(2)-mygca.YLim(1);
    posTicks = dLY*(realTicks/dt)./tMax+mygca.YLim(1);
    set(gca,'fontname','arial', 'FontSize',24,'LineWidth',3,'YTick',posTicks, 'YTickLabel',realTicks, 'YTickLabelRotation',0)   % Label Ticks
    xlabel('AP Position (um)');
    ylabel('Time (h)');
    %title(label2);
    pbaspect([1 1 1]);

    print([paths.resultsFolder '/png/' 'wt' num2str(i) '_kymo_' label2 '_noalign'],'-dpng','-loose','-painters');
    print([paths.resultsFolder 'wt' num2str(i) '_kymo_' label2 '_noalign'],'-depsc','-loose','-painters');

    %% plot channel1 without alignment

    figure
    p=imagesc(kimoMf,rangeKymos(1:2)); hold on; %,[0.07 5]
    plot(xEdge,(1:tMax)','-ow','LineWidth',3);
    c=colorbar;
    c.Label.FontSize = 24;
    c.Label.FontName = 'Arial';    
    c.Label.String = strrep([label1 ' (a.u.)'] ,'_',' ');    
    c.LineWidth= 3;
    
    mygca = gca;
    xt = mygca.XTick;    % Original 'XTick' Values
    realTicks = (0:500:(N.*dx));realTicks = sort(realTicks);
    dLX = mygca.XLim(2)-mygca.XLim(1);
    posTicks = dLX*(realTicks/dx)./N+mygca.XLim(1);
    set(gca, 'XTick',posTicks, 'XTickLabel',realTicks, 'XTickLabelRotation',0)   % Label Ticks
  
    yt = mygca.YTick;    % Original 'XTick' Values
    realTicks = 0:10:(tMax.*dt);
    dLY = mygca.YLim(2)-mygca.YLim(1);
    posTicks = dLY*(realTicks/dt)./tMax+mygca.YLim(1);
    set(gca,'fontname','arial', 'FontSize',24,'LineWidth',3,'YTick',posTicks, 'YTickLabel',realTicks, 'YTickLabelRotation',0)   % Label Ticks
    xlabel('AP Position (um)');
    ylabel('Time (h)');
    %title(label1);
    pbaspect([1 1 1]);

    print([paths.resultsFolder '/png/' 'wt' num2str(i) '_kymo_' label1 '_noalign'],'-dpng','-loose','-painters');
    print([paths.resultsFolder 'wt' num2str(i) '_kymo_' label1 '_noalign'],'-depsc','-loose','-painters');

    %% shift kimo usig positions posterior

    % kimos shifted
    kimoHUf_t= zeros(size(kimoHUf));
    kimoMf_t=  zeros(size(kimoMf));
    
    % stacks shifted
    sHUrot_t = zeros(size(sHUrot));
    sMrot_t =  zeros(size(sMrot));

    xShift = -xEdge+N-x0;

    for r=1:tMax %rows
        if(xShift>0)
            
            % shift kymo
            kimoHUf_t(r,1+xShift(r):N) = kimoHUf(r,1:N-xShift(r)); 
             kimoMf_t(r,1+xShift(r):N) = kimoMf(r,1:N-xShift(r)); 
            
            % shift stack
            sHUrot_t(:,1+xShift(r):N,r) = sHUrot(:,1:N-xShift(r),r); 
             sMrot_t(:,1+xShift(r):N,r) =   sMrot(:,1:N-xShift(r),r); 

        else
            
            % shift kymo
            kimoHUf_t(r,1:N-xShift(r)) = kimoHUf(r,1-xShift(r):N); 
            kimoMf_t(r,1:N-xShift(r)) =  kimoMf(r,1-xShift(r):N); 
            
            % shift stack
            sHUrot_t(:,1:N-xShift(r),r) =  sHUrot(:,1-xShift(r):N,r); 
             sMrot_t(:,1:N-xShift(r),r) =   sMrot(:,1-xShift(r):N,r); 
            
        end
    end
    
    saveOptions =[];
    saveOptions.overwrite = true;
    saveOptions.compress = 'lzw';
    
    saveastiff(uint8(sHUrot_t),[paths.shiftFolder replace(paths.name, label1,label2)],saveOptions);
    saveastiff(uint8(sMrot_t),[paths.shiftFolder paths.name],saveOptions);

%% plot kymo shifted channel1

    figure
    imagesc(kimoMf_t,rangeKymos(1:2)); hold on;  %,[0.07 0.7]
    c=colorbar;
    c.Label.FontSize = 24;
    c.Label.FontName = 'Arial';    
    c.Label.String = strrep([label1 ' (a.u.)'] ,'_',' ');    
    c.LineWidth= 3;

    mygca = gca;
    xt = mygca.XTick;    % Original 'XTick' Values
    realTicks = -(0:500:(N.*dx));realTicks = sort(realTicks);
    dLX = mygca.XLim(2)-mygca.XLim(1);
    posTicks = dLX*(realTicks/dx+N-x0)./N+mygca.XLim(1);
    set(gca, 'XTick',posTicks, 'XTickLabel',abs(realTicks), 'XTickLabelRotation',0)   % Label Ticks
  
    yt = mygca.YTick;    % Original 'XTick' Values
    realTicks = 0:10:(tMax.*dt);
    dLY = mygca.YLim(2)-mygca.YLim(1);
    posTicks = dLY*(realTicks/dt)./tMax+mygca.YLim(1);
    set(gca,'fontname','arial', 'FontSize',24,'LineWidth',3,'YTick',posTicks, 'YTickLabel',realTicks, 'YTickLabelRotation',0)   % Label Ticks
    xlabel('Distance from posterior (um)');
    ylabel('Time (h)');
    %title(label1);
    pbaspect([1 1 1]);

    print([paths.resultsFolder '/png/' 'wt' num2str(i) '_kymoShifted_' label1],'-dpng','-loose','-painters');
    print([paths.resultsFolder 'wt' num2str(i) '_kymoShifted_' label1],'-depsc','-loose','-painters');

    
    %% plot kymo shifted channel2

    figure
    imagesc(kimoHUf_t,rangeKymos(3:4)); hold on;  %,[0.07 0.7]
    c=colorbar;
    c.Label.FontSize = 24;
    c.Label.FontName = 'Arial';    
    c.Label.String = strrep([label2 ' (a.u.)'] ,'_',' ');    
    c.LineWidth= 3;

    mygca = gca;
    xt = mygca.XTick;    % Original 'XTick' Values
    realTicks = -(0:500:(N.*dx));realTicks = sort(realTicks);
    dLX = mygca.XLim(2)-mygca.XLim(1);
    posTicks = dLX*(realTicks/dx+N-x0)./N+mygca.XLim(1);
    set(gca, 'XTick',posTicks, 'XTickLabel',abs(realTicks), 'XTickLabelRotation',0)   % Label Ticks
  
    yt = mygca.YTick;    % Original 'XTick' Values
    realTicks = 0:10:(tMax.*dt);
    dLY = mygca.YLim(2)-mygca.YLim(1);
    posTicks = dLY*(realTicks/dt)./tMax+mygca.YLim(1);
    set(gca,'fontname','arial', 'FontSize',24,'LineWidth',3,'YTick',posTicks, 'YTickLabel',realTicks, 'YTickLabelRotation',0)   % Label Ticks
    xlabel('Distance from posterior (um)');
    ylabel('Time (h)');
    %title(label2);
    pbaspect([1 1 1]);
    mygcf = gcf;
    
    print([paths.resultsFolder '/png/' 'wt' num2str(i) '_kimoShifted_' label2],'-dpng','-loose','-painters');
    print([paths.resultsFolder 'wt' num2str(i) '_kimoShifted_' label2],'-depsc','-loose','-painters');

    
%% plot kymo shifted merged
 
   figure;
   color1 = 255*opts.colorsKymo(1,:);
   color2 = 255*opts.colorsKymo(2,:);

   km1 = (kimoMf_t-rangeKymos(1))./(rangeKymos(2)-rangeKymos(1));
   km2 = (kimoHUf_t-rangeKymos(3))./(rangeKymos(4)-rangeKymos(3));

   kimoRGB = cat(3,color1(1)*km1+color2(1)*km2,...
                   color1(2)*km1+color2(2)*km2,...
                   color1(3)*km1+color2(3)*km2);
                           
   kimoRGB=imresize(kimoRGB,[size(kimoRGB,1) size(kimoRGB,1)]);
   imshow(uint8(kimoRGB),'Border','Tight','InitialMagnification',500);
   
   set(gcf,'Position',mygcf.Position);
   set(gca,'Position',mygca.Position);
   set(gcf,'Position',mygcf.Position);
   gcaShow = gca;
   gcfShow = gcf;
   
    mygca = gca;
    xt = mygca.XTick;    % Original 'XTick' Values
    realTicks = -(0:500:(N.*dx));realTicks = sort(realTicks);
    dLX = mygca.XLim(2)-mygca.XLim(1);
    posTicks = dLX*(realTicks/dx+N-x0)./N+mygca.XLim(1);
    set(gca, 'XTick',posTicks, 'XTickLabel',abs(realTicks), 'XTickLabelRotation',0)   % Label Ticks
  
    yt = mygca.YTick;    % Original 'XTick' Values
    realTicks = 0:10:(tMax.*dt);
    dLY = mygca.YLim(2)-mygca.YLim(1);
    posTicks = dLY*(realTicks/dt)./tMax+mygca.YLim(1);
    set(gca,'fontname','arial', 'FontSize',24,'LineWidth',3,'YTick',posTicks, 'YTickLabel',realTicks, 'YTickLabelRotation',0)   % Label Ticks
    xlabel('Distance from posterior (um)');
    ylabel('Time (h)');
   axis on

   print([paths.resultsFolder '/png/' 'wt' num2str(i) '_kymoMerged' ],'-dpng','-loose','-painters');
   print([paths.resultsFolder 'wt' num2str(i) '_kymoMerged'],'-depsc','-loose','-painters');
   saveOptions = [];
   saveOptions.color = true;
   saveOptions.overwrite = true;
   saveOptions.compress = 'lzw';
   saveastiff(uint8(kimoRGB),[paths.resultsFolder 'wt' num2str(i) '_kymoMerged.tif'],saveOptions);
   
    %% Calculate and plot oscillations Hes7 multiple positions
    
    fOsc=figure;
    for x= ((N-x0-(0:50:250)/dx))
        kPlot = kimoHUf_t(:,round(x));
        kPlot2 = smooth(kPlot,round(0.5/dt)); 
        kPlot3 = sgolayfilt(kPlot2,3,2*round(1.8/dt)+1);
        plot((1:tMax)*dt,kPlot3,'LineWidth',3,'DisplayName',[num2str(abs(round((x-N+x0)*dx))) ' um']); hold on;
    end
    xlabel('Time (h)','FontSize',24);
    set(gca,'fontname','Arial','FontSize',24,'LineWidth',3);
    ylabel('Intensity (a.u.)');
    legend('show','Location','NorthEast');
    legend boxoff
    xlim([0 50]);
    pbaspect([1 1 1]);

    print([paths.resultsFolder '/png/' 'wt' num2str(i) '_' label2  '_intensity'],'-dpng');
    print([paths.resultsFolder 'wt' num2str(i) '_' label2  '_intensity'],'-depsc');

    
    %% Calculate and plot oscillations Mesp2
    
    fOscM=figure;
    for x= ((N-x0-(200:50:400)/dx))
        
        kMPlot = kimoMf_t(:,round(x));
        kMPlot2 = smooth(kMPlot,round(0.5/dt)); 
        kMPlot3 = sgolayfilt(kMPlot2,3,2*round(1.8/dt)+1);
        plot((1:tMax)*dt,kMPlot3,'LineWidth',3,'DisplayName',[num2str(abs(round((x-N+x0)*dx))) ' um']); hold on;
    end
    xlabel('Time (h)','FontSize',24);
    set(gca,'fontname','arial','FontSize',24,'LineWidth',3);
    ylabel('Intensity (a.u.)','Fontname','Arial');
    legend('show','Location','NorthEast');
    legend boxoff
    xlim([0 50]);
    pbaspect([1 1 1]);
    
    print([paths.resultsFolder '/png/' 'wt' num2str(i) '_' label1  '_intensity'],'-dpng','-loose','-painters');
    print([paths.resultsFolder 'wt' num2str(i) '_' label1  '_intensity'],'-depsc','-loose','-painters');
   
   %%

   myPSM.kimoMf_t = kimoMf_t;
   myPSM.kimoHUf_t = kimoHUf_t;
   myPSM.dt = dt;
   myPSM.dx = dx; %pxSize
   
   opts.optsOsc.number = opts.number;
   
   %% Calculate oscillations Hes7 and Mesp2
   
   myPSM.opts = opts;
   myPSM = kymoOscElong(myPSM,opts.optsOsc,paths);

   %%
   
end

