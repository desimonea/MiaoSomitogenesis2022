function [f] = config_plot(f)
%CONFIG_PLOT Summary of this function goes here
%   Detailed explanation goes here
f.OuterPosition = [100,100,700,700];
set(f.CurrentAxes, 'fontsize', 24,'linewidth',3','Position' ,[0.1921    0.1618    0.6553    0.689]);
box on;
end

