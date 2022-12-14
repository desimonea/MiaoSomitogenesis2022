function [Q,Qxy] = nematicImgColF(img,rhorange,sizes)

    if(any(isnan(img(:))))
       Q=NaN; 
    else
        ff=fftshift(fft2(img));
        G = abs(ff).^2;               

        [XX,YY] = (meshgrid(1:size(img,2),...
                            1:size(img,1)));
        fXX = (XX-floor(size(XX,2)/2)-1)./size(XX,2);
        fYY = (YY-floor(size(YY,1)/2)-1)./size(XX,1);
        [~,idxMax] = max(G(:));

        G((fXX==0)&(fYY==0))=NaN; %remove the mean
        if(~isnan(G(idxMax)))
             warning('Max is not at fXX==0 and fYY==0, you may had a problem with Fourier freqs');
        end

        [th,rho]=cart2pol(fXX(:),fYY(:));

        idxq = rho(:)>rhorange(1) & rho(:)<rhorange(2);
        Q =     nansum(G(idxq).*(((cos(th(idxq))).^2)-0.5))./nansum(G(idxq));
        Qxy = - nansum(G(idxq).*cos(th(idxq)).*sin(th(idxq)))./nansum(G(idxq));

   end
end

