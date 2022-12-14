function [myPSM] = kymoOscElong(myPSM,opts,paths)

        verbose = opts.verbose;
        spacelag = opts.spacelag;
        distHes7 = opts.distHes7; % in um, 150
        tspanHours = opts.tspanHours; %in h, 138
        x0 = myPSM.opts.x0;
        dx = myPSM.dx;
        dt=myPSM.dt;
        f0 = opts.f0;
        df = opts.df;

        kimoMf_t = myPSM.kimoMf_t;
        kimoHUf_t = myPSM.kimoHUf_t;

        [tMax,N] = size(kimoMf_t);
        
        if(verbose)
             fSomite= figure;
             yyaxis left;
        end
        
        x=round(N-x0-distHes7/dx);
        tspan = ((tspanHours)/dt);

        t = (1:tMax)*dt;
        kh = kimoHUf_t(:,x);
        
        if(isfield(opts,'simulated') & opts.simulated)
            warning('Here Hes7 is simulated. T = 7h phase = 4h');
            kh= sin((((1:tMax)*dt)/7)*2*pi)';
            myPSM.simulated = true;
        else
            myPSM.simulated = false; 
        end
        
        khOsc = (kh-smooth(kh,2*tspan,'moving'))./nanmean(kh);
        khOscPhase=(kh-smooth(kh,4*tspan,'lowess'))./nanmean(kh);

        yy = smooth(khOsc,round(1.5/dt)); 
        y = sgolayfilt(yy,3,2*round(1.8/dt)+1);
        
        Tprint = table;

        if(verbose)

            plot(t,y,'LineWidth',3,'Color',opts.colors(1,:)); hold on;

          
            fout = [paths.resultsFolder 'wt' num2str(opts.number) '_oscillHes7_Mesp2SourceData.txt'];
            FID = fopen(fout, 'w');
            
            Tprint = [Tprint table(t','VariableNames',{'Time (h)'})];
            Tprint = [Tprint table(y,'VariableNames', {'Intensity Hes7'})];

            ylabel(['Oscillation ' opts.label2   ' (a.u.)']);
            set(gca,'YColor',opts.colors(1,:));
            hold on; 
            xlabel('Time (h)','FontSize',24);
            ylim([-0.25 0.25]);
            yyaxis right
        end

        [cc,lags]=xcorr(khOsc,'normalized');
        ccp = cc(lags>=0); lagsp=lags(lags>=0);
        fccp = fit(lagsp'*dt,ccp,'smoothingspline');
        lagspHD = 0:0.1:25;
        ccpHD = feval(fccp,lagspHD);
        [~,loc]=findpeaks(ccpHD);
        myPSM.Th = lagspHD(loc(1));
        myPSM.xcorrh = [lagspHD' ccpHD];
        
        dfreq = 0.0667;
        L = round(1/(dfreq*dt));
        P  = abs(fft(kh-nanmean(kh),L)).^2;
        if(~(numel(P)==L))
           error('Problem with Fourier'); 
        end
        freqs = (0:(L/2))/(L*dt); 
        P1 = P(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        idxStart = find(freqs<f0-2*df,1,'last');
        idxEnd =   find(freqs>f0+2*df,1,'first'); 
        myPSM.Ph = sum(P1(idxStart:idxEnd))./sum(P1(freqs>0.15)); 
        
        x= floor(find(kimoMf_t(1,:)>max(kimoMf_t(1,:),[],2)/2,1,'last')+spacelag/dx);
        km = kimoMf_t(:,x);
        
        if(isfield(opts,'randomizeM') & opts.randomizeM)
            warning('Here Mesp2 is randomized');
            km= cumsum(randn(tMax,1)*0.01);
            myPSM.randomizeM = true;
        else
            myPSM.randomizeM = false; 
        end

        if(isfield(opts,'simulated') & opts.simulated)
            warning('Here Mesp2 is simulated. T = 7h phase = 4h');
            km= sin((((1:tMax)*dt)/7-1/7)*2*pi)';
            myPSM.simulated = true;
        else
            myPSM.simulated = false; 
        end
        
        
        t = (1:tMax)*dt;
        kmOsc=(km-smooth(km,2*tspan,'moving'));
        kmOscPhase=(km-smooth(km,4*tspan,'lowess'));

        P  = abs(fft(kmOsc-nanmean(kmOsc),L)).^2;
        freqs = (0:(L/2))/(L*dt); 
        P1 = P(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        idxStart = find(freqs<f0-2*df,1,'last');
        idxEnd =   find(freqs>f0+2*df,1,'first'); 
        myPSM.Pm = sum(P1(idxStart:idxEnd))./sum(P1(freqs>0.15));         
        
        yy =  smooth(kmOsc,round(1.5/dt)); %round(1.5/dt)
        y = sgolayfilt(yy,3,2*round(1.8/dt)+1);

        if(verbose)

            plot(t,y,'LineWidth',3,'Color',opts.colors(2,:)); hold on;

            Tprint = [Tprint table(t','VariableNames', {'Time2 (h)'})];
            Tprint = [Tprint table(y,'VariableNames', {'Intensity Mesp2'})];
            writetable(Tprint, fout,'WriteVariableNames',true,'WriteRowNames',false,'Delimiter','tab');

            ylabel(['Oscillation '  opts.label1 ' (a.u.)']);
            set(gca,'fontname','arial','FontSize',24,'LineWidth',3);
            xlim([0 70]);
            ylim([-1.2 1.2]);
            set(gca,'YColor',opts.colors(2,:));
            grid off
            pbaspect([1 1 1]);

            if(isfield(opts,'print') & opts.print)
                print([paths.resultsFolder 'wt' num2str(opts.number) '_oscillHes7_Mesp2'],'-depsc','-loose','-painters');
                print([paths.resultsFolder '/png/' 'wt' num2str(opts.number) '_oscillHes7_Mesp2'],'-dpng','-loose','-painters');
                fclose(FID);
            end
        end

        [cc,lags]=xcorr(kmOsc,'normalized');
        ccp = cc(lags>=0); lagsp=lags(lags>=0);
        fccp = fit(lagsp'*dt,ccp,'smoothingspline');
        lagspHD = 0:0.1:25;
        ccpHD = feval(fccp,lagspHD);
        [~,loc]=findpeaks(ccpHD);
        myPSM.Th = lagspHD(loc(1));
        myPSM.xcorrm = [lagspHD' ccpHD];
        
        % cross-correlation Hes7 and Mesp2
        
        [cc,lags]=xcorr(khOscPhase,kmOscPhase,'normalized');
        ccp = cc(lags>=0); lagsp=lags(lags>=0);
        fccp = fit(lagsp'*dt,ccp,'smoothingspline');
        lagspHD = 0:0.1:25;
        ccpHD = feval(fccp,lagspHD);
        myPSM.xcrosscorr = [lagspHD' ccpHD];
     
end

