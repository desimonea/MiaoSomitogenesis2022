function [err,p1] = plotBarsPaper(x,yList,opts,fig)

% it plots data with points and errorbars in xpoint
% yList is a list of points

if(nargin == 4)
   figure(fig);
end

if(isfield(opts,'symbol'))
   symbol = opts.symbol; 
else
   symbol ='o'; 
end

if(isfield(opts,'color'))
   color = opts.color; 
else
   color =[0.5 0.5 0.5]; 
end

if(isfield(opts,'CapSize'))
   CapSize = opts.CapSize; 
else
   CapSize =30; 
end

if(isfield(opts,'MarkerSize'))
   MarkerSize = opts.MarkerSize; 
else
   MarkerSize = 16; 
end

if(isfield(opts,'meanWidth'))
   meanWidth = opts.meanWidth; 
else
   meanWidth = 0.2; 
end



if(isfield(opts,'xPlot') & strcmp(opts.xPlot,'data'))
   xPlot = x;
else
   xPlot=x+randn(size(yList)).*0.05; hold on;
end


p1 = plot(xPlot,yList,symbol,'Color',[0.5 0.5 0.5],'LineWidth',3.6,'MarkerSize',MarkerSize); hold on;

if(isfield(opts,'err') & strcmp(opts.err,'SEM'))
    err = errorbar(nanmean(x),nanmean(yList),nanstd(yList)./sqrt(sum(~isnan(yList))),'LineWidth',3.6,'CapSize',30,'MarkerSize',40,'Color',color);   
else
    err = errorbar(nanmean(x),nanmean(yList),nanstd(yList),'LineWidth',3.6,'CapSize',CapSize,'MarkerSize',40,'Color',color);
end

plot([nanmean(x)-meanWidth nanmean(x)+meanWidth],[nanmean(yList) nanmean(yList)],'LineWidth',3.6,'Color',color);


end

