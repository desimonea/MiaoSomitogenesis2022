function [myPSM,figs]=spCorrSomitoid(myPSM,paths,opts)

    % this function calculates the spatial correlation of Mesp2, Uncx and
    % Mesp2-Uncx signal

    s=loadtiff([paths.inFolder paths.basename '.tif']);
    colorPlots = lines(2);
    
    %% read options
    
    t0 = myPSM.t0;
    tStart = opts.tStart;
    tEnd = opts.tEnd;
    tStep= opts.tStep;

    resFactor = opts.resFactor;
    dx = myPSM.dx;
    dt = myPSM.dt;
    
    if(isfield(myPSM,'ROI'))
        ROI = opts.ROI;
    end
    
    verbose = opts.verbose;
    toPrint = opts.toPrint;
    
    tStepPrint = opts.tStepPrint;
    
    pxSizeRes = dx*resFactor;
    dStepPxRes = opts.dStep/pxSizeRes; % in resized pixels
    
    label = opts.label;

    %% crop and smoothen
    
    if(isfield(paths,'ROIname'))
        roi = loadtiff([paths.inFolder paths.ROIname '.tif'])>0;
        roiRes = imresize(roi,1/resFactor,'method','nearest');
        sRes = double(imresize(s,1/resFactor,'method','bilinear'));
        if(size(roiRes,3)==1)
            roiRes = repmat(roiRes,[1 1 size(sRes,3)]);
            roi = repmat(roi,[1 1 size(sRes,3)]);
        end
        sCrop = double(s); 
        sCrop(~roi)=NaN;
    else
        sCrop = s(ROI(3):ROI(4),ROI(1):ROI(2),:);
        roiRes = ones(size(sRes));
        sRes = double(imresize(sCrop,1/resFactor,'method','bilinear'));
    end

    % size nuclei ~ 20 px, 5 px after resizing with px = 1.38
   sRes = imfilter(sRes,fspecial('gaussian',ceil(opts.gaussFiltParams(1)/pxSizeRes),2*ceil(opts.gaussFiltParams(2)/pxSizeRes/2)+1));
   bk =   imfilter(sRes,fspecial('gaussian',ceil(opts.gaussFiltParams(3)/pxSizeRes),2*ceil(opts.gaussFiltParams(4)/pxSizeRes/2)+1));
   sRes2= sRes./bk;
  
    if(isfield(opts,'xRange'))
            xStart = size(sRes2,2)+round(opts.xRange(1)/pxSizeRes);
            xEnd   = size(sRes2,2)+round(opts.xRange(2)/pxSizeRes);
            sRes2(:,[1:xStart xEnd:end],:) = NaN;
    end
      
    %% run time-points
    
    nColors = 100;
    colors = cool(nColors);

    corrL =[];
    maxPS =[];
    minsCorr =[];

    xPSall = [];
    PSall = [];

    xCorrAll = [];
    corrAll = [];

    if(verbose | toPrint)
        f1 = figure;  
    end

    tPrint = opts.tPrint; % decide if to remove tEnd
    
    tAllIdx = (tStart:tStep:tEnd);
    tAllH =   tAllIdx*dt+t0;
    
    f5 = figure;
    
    for i = tAllIdx
        
        if(i>size(sRes,3))
            break
        end

        rp = regionprops(roiRes(:,:,i));
        [~,idxRP] = max([rp.Area]);
        bb = rp(idxRP).BoundingBox;
        bb = [ceil(bb(1:2)) 2.*floor(bb(3:4)/2)];
        img2 = sRes2(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3),i);
        img2((roiRes(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3) ,i)==0))=NaN;
        img3 = img2 - nanmean(img2(:));
        nSidey = size(img3,1)-1;
        nSidex = size(img3,2)-1;
        
        idxNaN = isnan(img3);
        img3WZeros = img3; img3WZeros(idxNaN)=0;
        cc = xcorr2(img3WZeros)./(xcorr2(double(~idxNaN)));
        cc = cc./((nanstd(img3(:))).^2);
        cc(cc==0|isinf(cc))=NaN;
        
        [XXd,YYd]=meshgrid(-nSidex:nSidex,-nSidey:nSidey);
        d = sqrt(XXd.^2+YYd.^2);
        idxCCNotNan = ~isnan(cc(:));
        corrList = ([d(idxCCNotNan(:)) cc(idxCCNotNan(:))]);
        xBin = [0 0.01:dStepPxRes:300*dStepPxRes];
        [ii,xBin]=discretize(corrList(:,1),xBin);
        iiNotNaN = ~isnan(ii);
        
        disp(i);
        xMean = accumarray(ii(iiNotNaN),corrList(iiNotNaN,1),size(xBin'),@nanmean,NaN)';
        yMean = accumarray(ii(iiNotNaN),corrList(iiNotNaN,2),size(xBin'),@nanmean,NaN)';
        yStd =  accumarray(ii(iiNotNaN),corrList(iiNotNaN,2),size(xBin'),@nanstd,NaN)';
        yN =    accumarray(ii(iiNotNaN),corrList(iiNotNaN,2),size(xBin'),@(x) sum(~isnan(x)),NaN)';
        ySem  =  yStd./sqrt(yN);
        Y1 = [xMean' yMean' yStd' yN' ySem'];

        tIndex = round(nColors*(i*dt+myPSM.t0-opts.barPlotLim(1)+2)./(opts.barPlotLim(2)-opts.barPlotLim(1)+2));

        if(verbose | toPrint)
            if(any(tPrint==i))
                figure(f1);
                errorbar(xMean.*pxSizeRes,yMean,ySem,'Color',colors(tIndex,:),'LineWidth',3,'DisplayName',[num2str(round(i*dt+myPSM.t0)) ' h']);
                set(gca,'LineWidth',3,'FontSize',24);
                xlabel('Distance (um)'); ylabel('Spatial auto-correlation');
                hold on;
                xlim(opts.xLimUm);
            end
        end

        % fit minimum autocorrelation
        [xMeanF,yMeanF,ySemF]=prepareCurveData(xMean,yMean,ySem);
        idxNotInf = ~isinf(xMeanF.*yMeanF.*ySemF) & ySemF>0;
        xMeanF= xMeanF(idxNotInf);
        yMeanF= yMeanF(idxNotInf);
        ySemF=  ySemF(idxNotInf);

        f=fit(xMeanF,yMeanF,'smoothingspline','Weights',(ySemF.^(-2)));
        dxFit = 0.1;
        xMeanFit = 0:dxFit:max(xMean);
        yMeanFit=feval(f,xMeanFit);
        [valueMin,idxMin] = findpeaks(-yMeanFit,'MinPeakProminence',0.01,'MinPeakDistance',dStepPxRes*10);
        valueMin = -valueMin; 
        
        if(~isempty(idxMin))
            posMin = xMeanFit(idxMin(1));
        else
            posMin = NaN;
            valueMin = NaN;
        end
        
        if(valueMin(1)< opts.valueMinMin & ...
                         valueMin(1) < -opts.valueMinStd*nanstd(yMeanFit((xMeanFit> posMin+ 50/pxSizeRes) & (xMeanFit<posMin+300/pxSizeRes))));
            minsCorr = vertcat(minsCorr,[i*dt+myPSM.t0 posMin.*pxSizeRes feval(f,posMin)]);
        else
            minsCorr = vertcat(minsCorr,[NaN NaN NaN]);
        end
        disp([minsCorr(end,:)]);


        corrAll = cat(3,corrAll,[xMean(:).*pxSizeRes yMean(:) ySem(:)]);
    end

    if(toPrint)
        figure(f1);
        colormap(cool(100));
        c=colorbar;
        c.FontSize = 24;
        c.TickLabels = (opts.barPlotLim(2)-opts.barPlotLim(1)).*c.Ticks+opts.barPlotLim(1);
        c.Label.String = 'Time (h)';
        c.LineWidth = 3;
        ylim([-0.2 1.05]);
        print([paths.resFolder paths.basename '_autocorr'],'-depsc','-loose','-painters');
        print([paths.resFolder '/png/' paths.basename '_autocorr'],'-dpng','-loose','-painters');
    end
    
    if(verbose | toPrint)
       
        % print average and std
        f3=figure;
        yyaxis left
        sInt = []; 
        for t=tAllIdx(tAllIdx<=size(sRes,3))
            imgToMean = double(sCrop(:,:,t));
            sInt = vertcat(sInt,[nanmean(imgToMean(:)) nanstd(imgToMean(:))]);
        end

        plot(tAllH(tAllIdx<=size(sRes,3)),sInt(:,1),'o-','LineWidth',3);
        xlabel('Time (h)');
        xlim(opts.barPlotLim);
        ylabel(['Average ' opts.labelPrint ' signal']);
        set(gca,'LineWidth',3,'FontSize',24);
        yyaxis right
        ax = gca;
        ax.YAxis(2).Color = colorPlots(2,:);
        plot(tAllH(tAllIdx<=size(sRes,3)),sInt(:,2),'o-','LineWidth',3);
        ylabel(['Std ' label ' signal']);

        if(toPrint)
            print([paths.resFolder paths.basename '_signal'],'-depsc','-loose','-painters');
            print([paths.resFolder '/png/' paths.basename '_signal'],'-dpng','-loose','-painters');
        end

        f4=figure;
        yyaxis left

        toPlot = minsCorr(:,1) < opts.xMaxAutocorrUm; %Inf
        plot(minsCorr(toPlot,1),minsCorr(toPlot,2),'o-','LineWidth',3);
        xlabel('Time (h)');
        ylabel('Pos min auto-correlation (um)');
        set(gca,'LineWidth',3,'FontSize',24);
        xlim(opts.barPlotLim);
        ylim([80 160]);
        yyaxis right
        ax = gca;
        ax.YAxis(2).Color = colorPlots(2,:);
        plot(minsCorr(:,1),minsCorr(:,3),'o-','LineWidth',3);
        ylabel('Value min auto-correlation (um)');
        ylim([-0.2 0]);
        print([paths.resFolder paths.basename '_minAutoCorr.eps'],'-depsc','-loose','-painters');
        print([paths.resFolder '/png/' paths.basename '_minAutoCorr.png'],'-dpng' ,'-loose','-painters');
    end
    
    %% saving
    
    myPSM.spCorr.tAllIdx = tAllIdx;
    myPSM.spCorr.tAllH = tAllH;
    myPSM.spCorr.corr = corrAll;
    myPSM.spCorr.sInt = sInt;
    myPSM.spCorr.minCorr = [ minsCorr];
    
    figs = {f1,f3,f4};

end