function [myPSM,figs]=spCorrTBinRecal(myPSM,paths,opts)

    % this function t-bins spatial correlations and calculates minima
    colorPlots = lines(2);
    
    %% read options
    
    t0 = myPSM.t0;
    tStart = opts.tStart;
    tEnd = opts.tEnd;

    dx = myPSM.dx;
    dt = myPSM.dt;
    verbose = opts.verbose;
    toPrint = opts.toPrint;
    
    tPrint = opts.tPrint;

    if(isfield(opts,'yLimAutoCorr'))
       yLimAutoCorr = opts.yLimAutoCorr; 
    else
       yLimAutoCorr = [-0.2 1.05]; 
    end
    
    
    if(isfield(opts,'tOffset'))
       tOffset = opts.tOffset; 
    else
       tOffset = 0; 
    end
    
%% read correlations and average
       
    tBinFrame = round(opts.tBinH/myPSM.dt);
    
    cc = myPSM.spCorr.corr;
    cc(:,3,:) = cc(:,3,:).^2;
    ccTbin = movmean(cc,tBinFrame,3,'EndPoints','Discard');
    tIdxs = (opts.tStart-myPSM.spCorr.tAllIdx(1)+1):tBinFrame:size(ccTbin,3); 

    ccTbin(:,3,:) =  sqrt(ccTbin(:,3,:));
    ccTbin = ccTbin(:,:,tIdxs);
    tAllHTbin = movmean(myPSM.spCorr.tAllH,tBinFrame,'EndPoints','Discard');
    tAllHTbin=tAllHTbin(tIdxs);
    tAllIdxTbin = movmean(myPSM.spCorr.tAllIdx,tBinFrame,'EndPoints','Discard');
    tAllIdxTbin=tAllIdxTbin(tIdxs);
    
%% plot time-points
    
    nColors = 100;
    colors = cool(nColors);
  
    if(verbose | toPrint)
        f1 = figure;  
        T1 = table;
    end
    
    minsCorrTbin = [];
    tPrinted = [];
    
    for i = 1:size(ccTbin,3)
                
        tHere = tAllHTbin(i)-tOffset;
        tIndex = round(nColors*(tHere-opts.barPlotLim(1)+2)./(opts.barPlotLim(2)-opts.barPlotLim(1)+2));
        tIndex(tIndex>100) = 100;
        tIndex(tIndex<1) = 1;
        
        xMean = ccTbin(:,1,i);
        yMean = ccTbin(:,2,i);
        ySem =  ccTbin(:,3,i);
        
        if(verbose | toPrint)
            idxT=find(tAllIdxTbin(i)>=tPrint);
            if(~isempty(idxT))
                figure(f1);
                errorbar(xMean,yMean,ySem,'Color',colors(tIndex,:),'LineWidth',3,'DisplayName',[num2str(round(tHere)) ' h']);
                T1 = [T1 table(xMean, yMean, ySem,'VariableNames',{[num2str(round(tHere)) ' h d (um)'],[num2str(round(tHere)) ' h Auto-corr Mean'],[num2str(round(tHere)) ' h Auto-corr SEM']})];               
                hold on;
                set(gca,'LineWidth',3,'FontSize',24);
                xlabel('Distance (um)'); ylabel('Spatial auto-correlation');
                xlim(opts.xLimUm);
                ylim(yLimAutoCorr);
                tPrint(idxT)=NaN; %remove time-point from list to still print
                tPrinted = vertcat(tPrinted,tHere);
            end
        end

        % fit minimum autocorrelation
        [xMeanF,yMeanF,ySemF]=prepareCurveData(xMean,yMean,ySem);
        idxNotInf = ~isinf(xMeanF.*yMeanF.*ySemF) & ySemF>0;
        xMeanF= xMeanF(idxNotInf);
        yMeanF= yMeanF(idxNotInf);
        ySemF=  ySemF(idxNotInf);

        f=fit(xMeanF,yMeanF,'smoothingspline','Weights',(ySemF.^(-2)));
        dxFit = 0.1; %in um
        xMeanFit = 0:dxFit:max(xMean);
        yMeanFit=feval(f,xMeanFit);
        [valueMin,idxMin] = findpeaks(-yMeanFit,'MinPeakProminence',0.01,'MinPeakDistance',3/dxFit);
        %[valueMax,idxMax] = findpeaks( yMeanFit,'MinPeakProminence',0.01,'MinPeakDistance',dStepPxRes*10);
        valueMin = -valueMin; 
        
        if(~isempty(idxMin))
            posMin = xMeanFit(idxMin(1));
        else
            posMin = NaN;
            valueMin = NaN;
        end
    
        if(valueMin(1)< opts.valueMinMin & ...
                         valueMin(1) < -opts.valueMinStd*nanstd(yMeanFit(xMeanFit> posMin+50 &xMeanFit<posMin+300)));
                        %valueMin(1) <  opts.fractionMinMin*min(valueMin) &...
            minsCorrTbin = vertcat(minsCorrTbin,[tHere posMin feval(f,posMin)]);
        else
            minsCorrTbin = vertcat(minsCorrTbin,[NaN NaN NaN]);
        end
        disp([minsCorrTbin(end,:)]);
    end

    if(toPrint)
        figure(f1);
        %legend('Show'); legend boxoff; 
        colormap(cool(100));
        c=colorbar;
        c.FontSize = 24;
        c.TickLabels = (opts.barPlotLim(2)-opts.barPlotLim(1)).*c.Ticks+opts.barPlotLim(1);
        c.Label.String = 'Time (h)';
        c.LineWidth = 3;
        ylim(yLimAutoCorr);
        pbaspect([1 1 1]);

        % print source data
        fout = [paths.resFolder myPSM.basename '_autocorr.txt'];
        FID = fopen(fout, 'w');
        writetable(T1, fout,'WriteVariableNames',true,'WriteRowNames',false,'Delimiter','tab');
        fclose(FID);

        print([paths.resFolder myPSM.basename '_autocorr.eps'],'-depsc','-loose','-painters');
        print([paths.resFolder '/png/' myPSM.basename '_autocorr.png'],'-dpng','-loose','-painters');
        %standardizePlotAle(gcf,gca,,opts);
    end
    
    if(verbose | toPrint)
       
        f4=figure;

        T4 = table;

        yyaxis left

        toPlot = minsCorrTbin(:,2) < opts.xMaxAutocorrUm; 
        plot(minsCorrTbin(toPlot,1),minsCorrTbin(toPlot,2),'o-','LineWidth',3);


        xlabel('Time (h)');
        ylabel('Pos min auto-correlation (um)');
        set(gca,'LineWidth',3,'FontSize',24);
        xlim(opts.barPlotLim);
        ylim(opts.posyPlotLim);
        yyaxis right
        ax = gca;
        ax.YAxis(2).Color = colorPlots(2,:);
        plot(minsCorrTbin(toPlot,1),minsCorrTbin(toPlot,3),'o-','LineWidth',3);
        ylabel('Value min auto-correlation (um)');
        ylim([yLimAutoCorr(1) 0]);
        pbaspect([1 1 1]);

        T4 = table(minsCorrTbin(toPlot,1),minsCorrTbin(toPlot,2),minsCorrTbin(toPlot,3),'VariableNames',{'Time (h)','Pos min (um)','Value min'});
        % print source data
        fout = [paths.resFolder myPSM.basename '_minAutoCorr.txt'];
        FID = fopen(fout, 'w');
        writetable(T4, fout,'WriteVariableNames',true,'WriteRowNames',false,'Delimiter','tab');
        fclose(FID);

        print([paths.resFolder myPSM.basename '_minAutoCorr.eps'],'-depsc','-loose','-painters');
        print([paths.resFolder '/png/' myPSM.basename '_minAutoCorr.png'],'-dpng' ,'-loose','-painters');
        %standardizePlotAle(gcf,gca,[paths.resFolder myPSM.basename '_minAutoCorr.eps'],opts);
    end
    
    %% saving
    
    myPSM.spCorr.tAllIdxTbin = tAllIdxTbin;
    myPSM.spCorr.tAllHTbin = tAllHTbin;
    myPSM.spCorr.corrTbin = ccTbin;
    myPSM.spCorr.minCorrTbin =  minsCorrTbin;
    myPSM.spCorr.tPrintedTbin = tPrinted;
    
    figs = {f1,f4};

end