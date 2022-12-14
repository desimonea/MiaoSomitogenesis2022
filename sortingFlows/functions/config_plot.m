function [f] = config_plot(f,c)
arguments
    f
    c = [];
end
%CONFIG_PLOT Summary of this function goes here
%   Detailed explanation goes here

f.OuterPosition = [100,100,700,700];
set(f.CurrentAxes, 'fontsize', 24,'linewidth',3','Position' ,[0.1921    0.1618    0.6553    0.689]);


if ~isempty(c)
    f.OuterPosition = [100,100,800,700];
    set(f.CurrentAxes, 'fontsize', 24,'linewidth',3','Position' ,[0.1921    0.1618    0.5553    0.689]);
    c.LineWidth = 3;
    c.FontSize = 24;
end
box on;
end

