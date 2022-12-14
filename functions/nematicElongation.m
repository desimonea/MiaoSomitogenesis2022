function [myPSM,opts] = nematicElongation(paths,opts,myPSM)
% calculate nematic index for segmentoids 

    colorLevels = 100;
    colorTMaxH = 50; %in h
    colors = cool(colorLevels);

    dx = myPSM.dx;
    dt = myPSM.dt;
    dxStepNematic = round(opts.dxStepNematicUm/dx); %(conversion from um to dx)
    Lw = 2*round(opts.LNematicWindowUm/(dx*2)); %to make it even
    rhorange = opts.rhoRangeUm*dx;
    printStepSize = round(opts.printStepSizeH/dt);

    minFractionInRoi = opts.minFractionInRoi;
    
    % load stacks and ROIs
    number = myPSM.opts.number;
    s1 = loadtiff([paths.inFolder paths.name1]);
    rois = loadtiff([paths.inFolder paths.nameROI]);
    s2 = loadtiff([paths.inFolder paths.name2]);
    
    % normalize stacks and merge them
    s1forAv = s1; s2forAv = s2;
    s1forAv(rois==0)=NaN; s2forAv(rois==0)=NaN; %hfaving the ROI adds just a factor
    tot1 = nansum(nansum(s1forAv,1),2)./nansum(nansum(rois,1),2);
    tot2 = nansum(nansum(s2forAv,1),2)./nansum(nansum(rois,1),2);
    m1 = repmat(tot1,[size(s1,1) size(s1,2)]);
    m2 = repmat(tot2,[size(s2,1) size(s2,2)]);
    s = ((double(s1)./m1-double(s2)./m2).*(m1+m2)./2)+255/2; 
    
    % temporal idxs
    tIdx = 1:size(s,3); %idx
    tExact = (tIdx-1).*dt; %hours exact from idx

    Q= [];
    Qmean = [];
    
    disp('Elongatoid dt tMaxExact tMax (h)');
    disp([number dt tExact(end) size(s,3).*dt]);
    
    if(opts.verbose)
        figure;
    end
    
    for t = tIdx
        
        if(t>size(s,3))
            break
        end
        
        % read time frame
        img = double(s(:,:,t));
        roi = rois(:,:,t);
        
        % calculate yPosition sliding window
        [~,YY] = meshgrid(1:size(img,2),1:size(img,1));
        ww=YY.*logical(roi); ww(~roi)=NaN;
        yPos = floor(nanmean(ww,1));

        % create image nematic
        QHere = nan(size(img,2),3);
        for x=1:dxStepNematic:size(img,2)-Lw
            ymin = yPos(x)-Lw/2;
            ymax = yPos(x)+Lw/2;
            if(ymin<1 || ymax > size(img,1) || isnan(ymin*ymax))
               continue;
            end
                imgHere = img(ymin:ymax,x:x+Lw);
                roiHere = roi(ymin:ymax,x:x+Lw);
                if(~(size(imgHere,1)==size(imgHere,2)))
                    error('Sliding window is not squared');
                end
                if(sum(roiHere(:))<minFractionInRoi*numel(roiHere))
                    QHere(x,:) = [NaN NaN NaN];
                    continue
                end
               [Qimg,Qxyimg]= nematicImgColF(imgHere,rhorange); %to be checked
               QHere(x,:) = [-(size(img,2)-x-Lw/2) yPos(x) sqrt(Qimg.^2+Qxyimg.^2)];
        end
        
        Q = cat(3,Q,QHere);
            
        if(opts.verbose & mod(t,printStepSize)==1)
            colorIdx = floor(colorLevels*(t-1)*dt./(colorTMaxH-1)+1);
            idxPlot = ~isnan(QHere(:,1));
            if(~any(idxPlot))
                continue
            end
            plot(QHere(idxPlot,1).*dx,smooth(QHere(idxPlot,3),0.2,'moving'), ...
                 '-','Color',colors(colorIdx,:),'LineWidth',3,'DisplayName',[num2str(floor((t-1)*dt)) ' h']);
            hold on;
            drawnow
        end
        
        %calculate averages
        idxX = QHere(:,1)>opts.xRange(1) & QHere(:,1)<opts.xRange(2);
        Qmean = vertcat(Qmean,[number tExact(t) nanmean(QHere(idxX,3),1)]);
    end
    
    if(opts.verbose)
        xlabel('Position (um)');
        ylabel('Nematic order Q');
        xlim([-1000 0]);
        ylim([0 0.15]);
        xticks([-1000 -500 0]);
        yticks(0:0.05:0.15);
        set(gca,'FontSize',24,'FontName','Arial','LineWidth',4);
        box on
        legend('show','Location','northeast');
        legend boxoff;
        pbaspect([1 1 1]);
        
        if(opts.toPrint)
            print([paths.resFolder 'Q_' paths.name1(1:end-4)],'-depsc','-loose','-painters');
            print([paths.resFolder  '/png/' 'Q_' paths.name1(1:end-4)],'-dpng','-loose','-painters');
        end
    end
    
    % calculate t-binned averages
    idxT=discretize(Qmean(:,2),opts.tBins);
    QmeanTbinned= accumarray(idxT,Qmean(:,3),size(opts.tBins'),@nanmean,NaN);
    
    myPSM.nematic.Q = Q;
    myPSM.nematic.Qmean = Qmean;
    myPSM.nematic.QmeanTbinned = QmeanTbinned;
    myPSM.nematic.tExact = tExact;
    myPSM.nematic.nematicOpts = opts;

end