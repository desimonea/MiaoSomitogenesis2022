function [] = standardizePlotAle(mygcf,mygca,printName,opts)

    if(nargin<4)
        opts=struct;
        opts.noGca=false;
    end
    
    if(isfield(opts,'sizeMultipliers'))
        sizeMultipliers = opts.sizeMultipliers;
    else
        sizeMultipliers=[1 1];
    end

    scaleFactor = 5;

    ax = mygca;
    ax.XLabel.FontSize = 7.2.*scaleFactor;
    ax.YLabel.FontSize = 7.2.*scaleFactor;
    ax.FontSize = 7.2.*scaleFactor;
    
    if(~(isfield(opts,'noGca') | opts.noGca))
        set(mygca,'LineWidth',0.6.*scaleFactor,'FontName','Helvetica');
    end
    
    set(mygcf,'units','centimeters');
    posgcf = mygcf.Position;
    
    
    if(~(isfield(opts,'noGca') | opts.noGca))
        posgca = mygca.Position;
    end

    % rectangle
    %set(mygcf,'units','centimeters','Position',[posgcf(1) posgcf(2) 4.4 3.8].*scaleFactor);
    %set(mygca,'Position',[0.13 0.17 0.85 0.75]);
    % 
    % % square
    set(mygcf,'units','centimeters','Position',[posgcf(1) posgcf(2) sizeMultipliers(1)*4.4.*scaleFactor sizeMultipliers(2)*4.4.*scaleFactor]);
    
    if(~(isfield(opts,'noGca') & opts.noGca))
        set(gca,'Position',[0.17 0.17 0.81 0.81]);
    end
    
    set(mygcf,'renderer','Painters')
    ax.Color = 'None';
    
    if(size(ax.YAxis,1)==1)
        ax.XColor = 'black';
        ax.YColor = 'black'; 
        ax.ZColor = 'black';
    end
    
    legend.FontSize=7.2.*scaleFactor;
    legend.Box='off';
    legend.Visible='on';
    
    if(strcmp(printName(end-3:end),'.eps') | strcmp(printName(end-3:end),'.png'))
        printName = printName(1:end-4);
    end
    
    if(nargin>2)
        print('-depsc','-r300',[printName '.eps'],'-loose','-painters');
        print('-dpng','-r300', [printName '.png'],'-loose','-painters');
    end

end

