function [myPSM]=createSomitoid(paths,opts)

    % this function creates myPSM for somitoids
    %% read options
    
    % calculate dx
    info=imfinfo([paths.inFolder paths.basename '.tif']);
    dx=1./info(1).XResolution;
    dy=1./info(1).YResolution;
    tMax = numel(info);

    % parse dt
    fname = strrep(info(1).Filename,'pt','.');
    delimiters = regexp(fname,'_');
    posMin =     regexp(fname,'min');
    delimitersMin = delimiters(find(delimiters<posMin,1,'last'));
    dt = str2num(fname(delimitersMin+1:posMin-1))/60; %in h
    
    % parse t0
    delimiters = regexp(fname,'_');
    posT0 =     regexp(fname,'hr start');
    delimitersT0 = delimiters(find(delimiters<posT0,1,'last'));
    t0 = str2num(fname(delimitersT0+1:posT0-1)); %in h
    
    %% saving
   myPSM = [];
   myPSM.dt   = dt;
   myPSM.dx   = dx; 
   myPSM.tMax = tMax; 
   myPSM.t0   = t0;    

end