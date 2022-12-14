%% reinitialize
close all;
clear;
clc;
addpath(genpath('./functions/'))
%% define paths 
paths = [];
paths.expname1 = 'Somitoid1_ 72-84hr_6min interval_2021_09_10__10_15_25';
paths.expname2 = 'Somitoid2_ 72-84hr_6min interval_2021_09_18__15_24_26';
paths.dataname1 = [paths.expname1,'-data'];
paths.dataname2 = [paths.expname2,'-data'];
paths.plotFolder = '../sortingFlows_figures\figure_revision2_somitoid';
mkdir(paths.plotFolder);
mkdir('../excel_data')
% define pixelsize
xy_pxsize = 0.692;
z_pxsize = 2.24;
t_size = 6; % min

%% load and merge data
%% load data
%% load data
cell_struct_trim1 = load(['./cell_struct/',paths.dataname1,'_cell_struct_trim.mat']);
cell_struct_trim2 = load(['./cell_struct/',paths.dataname2,'_cell_struct_trim.mat']);
cell_struct_trim = cat(2,cell_struct_trim1.cell_struct_trim, cell_struct_trim2.cell_struct_trim);



%% mannual curation
bad_lineage1 = [7,8,12,14,16,54,57,70,71,128,144,147,150,195,232,963,966];
bad_lineage2 = [6,7,8,15,23,25,27,40,58,60,64,65,93,110,112,118,119,147,166,169,74,177,179,187,209,225,226,259,274,334];

curation_logical = (~ismember([cell_struct_trim.lineageId],bad_lineage1)&[cell_struct_trim.somitoidID]==1) | ...
    (~ismember([cell_struct_trim.lineageId],bad_lineage2)&[cell_struct_trim.somitoidID]==2);

cell_struct_trim = cell_struct_trim(curation_logical);

%% calculate intensity_ratio_norm
largest_final1 = max(arrayfun(@(x)x.frame(end),cell_struct_trim1.cell_struct_trim));
mean_by_time1 = arrayfun(@(y)mean(cell2mat(arrayfun(@(x)x.intensity_ratio(x.frame==y)',...
    cell_struct_trim1.cell_struct_trim,'UniformOutput',false))),...
    0:largest_final1);
largest_final2 = max(arrayfun(@(x)x.frame(end),cell_struct_trim2.cell_struct_trim));
mean_by_time2 = arrayfun(@(y)mean(cell2mat(arrayfun(@(x)x.intensity_ratio(x.frame==y)',...
    cell_struct_trim2.cell_struct_trim,'UniformOutput',false))),...
    0:largest_final2);
intensity_ratio_norm_cell = {};
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    if s.somitoidID ==1
        intensity_ratio_norm_cell{i} = s.intensity_ratio./mean_by_time1(s.frame+1)';    
    elseif s.somitoidID ==2
        intensity_ratio_norm_cell{i} = s.intensity_ratio./mean_by_time2(s.frame+1)'; 
    end
end
[cell_struct_trim.intensity_ratio_norm] = intensity_ratio_norm_cell{:};
%% backup cell_struct_trim
cell_struct_trim_backup = cell_struct_trim;










%% plot for both somitoids
cell_struct_trim = cell_struct_trim_backup;
%% plot intensityratio vs time 
% small - cyan, large - magenta
ratio_first = arrayfun(@(x)x.intensity_ratio_norm(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);
ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(end-4:end)),cell_struct_trim);

percent = 50;
logical_lower = ratio_first<prctile(ratio_first,percent);
logical_higher = ratio_first>prctile(ratio_first,100-percent);

% calculate average per cell intensity for each frame


f = figure('visible','on');
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    if logical_lower(i)
        pb = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-c','LineWidth',1);
        pcolor = pb;
        pcolor.Color = [pcolor.Color,0.5];
    elseif logical_higher(i)
        pr = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-m','LineWidth',1);
        pcolor = pr;
        pcolor.Color = [pcolor.Color,0.5];
    end

    hold on;

end

% plot median for each group
median_lower = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_lower),'UniformOutput',false))),...
    0:120);
median_higher = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_higher),'UniformOutput',false))),...
    0:120);

pb=plot(0:0.1:12,sgolayfilt(median_lower,5,21),'c-','LineWidth',10);
pr=plot(0:0.1:12,sgolayfilt(median_higher,5,21),'m-','LineWidth',10);

xlabel('Time (hr)')
ylabel('Normalized Mesp2 intensity')
xlim([0,12.2])
% ylim([0,2])
legend([pb,pr], {['lower ' num2str(percent) '% at 1st frame'],['higher ' num2str(percent) '% at 1st frame']},'Location','best','color','none','box','off');
% title('somitoid1')
config_plot(f);

hold off;
%%
savepath = [paths.plotFolder filesep 'merge_somitoid12_intensityratio_vs_time_normalized_50prctile'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);
%% plot intensityratio vs time 
% small - cyan, large - magenta
ratio_first = arrayfun(@(x)x.intensity_ratio_norm(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);
ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(end-4:end)),cell_struct_trim);

percent = 25;
logical_lower = ratio_first<prctile(ratio_first,percent);
logical_higher = ratio_first>prctile(ratio_first,100-percent);

% calculate average per cell intensity for each frame


f = figure('visible','on');
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    if logical_lower(i)
        pb = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-c','LineWidth',1);
        pcolor = pb;
        pcolor.Color = [pcolor.Color,0.5];
    elseif logical_higher(i)
        pr = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-m','LineWidth',1);
        pcolor = pr;
        pcolor.Color = [pcolor.Color,0.5];
    end

    hold on;

end

% plot median for each group
median_lower = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_lower),'UniformOutput',false))),...
    0:120);
median_higher = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_higher),'UniformOutput',false))),...
    0:120);

pb=plot(0:0.1:12,sgolayfilt(median_lower,5,21),'c-','LineWidth',10);
pr=plot(0:0.1:12,sgolayfilt(median_higher,5,21),'m-','LineWidth',10);

xlabel('Time (hr)')
ylabel('Normalized Mesp2 intensity')
xlim([0,12.2])
% ylim([0,2])
legend([pb,pr], {['lower ' num2str(percent) '% at 1st frame'],['higher ' num2str(percent) '% at 1st frame']},'Location','best','color','none','box','off');
% title('somitoid1')
config_plot(f);

hold off;
%%
savepath = [paths.plotFolder filesep 'merge_somitoid12_intensityratio_vs_time_normalized_25prctile'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);

%% mahal
d2 = mahal([ratio_first',ratio_final'],[ratio_first',ratio_final']);
f = figure();
cm = colormap('parula');
scatter(ratio_first',ratio_final',100,d2,'o','LineWidth',2)
c = colorbar;
c.Label.String = 'Mahalanobis Distance';
ylim([0,4]);
xlim([0,4]);
xlabel('Mesp2 at start')
ylabel('Mesp2 at end')
% axis equal;
config_plot(f,c);
hold off;
%%
savepath = [paths.plotFolder filesep 'merge_somitoid12_last_by_first_normalized_mahal'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);
%% last frame vs first frame

[~,rank_first] = ismember(ratio_first(:),unique(ratio_first(:)));
[~,rank_final] = ismember(ratio_final(:),unique(ratio_final(:)));
r_sprman = corr([ratio_first(:),ratio_final(:)],'type','Spearman')

f = figure('visible','on');
cm = colormap(lines(2));
hold on;


xdata = ratio_first;
ydata = ratio_final;
% xdata = ratio_first(ratio_first<0.78);
% ydata = ratio_final(ratio_first<0.78);

trim_logical = d2<1000;



xdata = xdata(trim_logical);
ydata = ydata(trim_logical);

data = plot(xdata,ydata,'o','color',cm(1,:),'MarkerSize',10,'LineWidth',2);
plot(ratio_first(~trim_logical),ratio_final(~trim_logical),'x','color','m','MarkerSize',10,'LineWidth',2)


[fitobj,gof] = fit(xdata(:),ydata(:),fittype({'x','1'}));
x_fit = linspace(0,4,100);
y_fit = feval(fitobj,x_fit);
fitted = plot(x_fit,y_fit,'-','color',cm(2,:), 'LineWidth',2);

lm = fitlm(xdata(:),ydata(:));


xlabel('Mesp2 at start')
ylabel('Mesp2 at end')
ylim([0,4])
xlim([0,4])
legend([fitted], {['R^2 = ', num2str(round(gof.rsquare,2))]},'Location','best','color','none','box','off')


config_plot(f);
hold off;
%%
savepath = [paths.plotFolder filesep 'merge_somitoid12_last_by_first_normalized_all'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);

%% last frame vs first frame

[~,rank_first] = ismember(ratio_first(:),unique(ratio_first(:)));
[~,rank_final] = ismember(ratio_final(:),unique(ratio_final(:)));
r_sprman = corr([ratio_first(:),ratio_final(:)],'type','Spearman')

f = figure('visible','on');
cm = colormap(lines(2));
hold on;


xdata = ratio_first;
ydata = ratio_final;
% xdata = ratio_first(ratio_first<0.78);
% ydata = ratio_final(ratio_first<0.78);

trim_logical = d2<8;



xdata = xdata(trim_logical);
ydata = ydata(trim_logical);

data = plot(xdata,ydata,'o','color',cm(1,:),'MarkerSize',10,'LineWidth',2);
plot(ratio_first(~trim_logical),ratio_final(~trim_logical),'x','color','m','MarkerSize',10,'LineWidth',2)


[fitobj,gof] = fit(xdata(:),ydata(:),fittype({'x','1'}));
x_fit = linspace(0,4,100);
y_fit = feval(fitobj,x_fit);
fitted = plot(x_fit,y_fit,'-','color',cm(2,:), 'LineWidth',2);

lm = fitlm(xdata(:),ydata(:));


xlabel('Mesp2 at start')
ylabel('Mesp2 at end')
ylim([0,4])
xlim([0,4])
legend([fitted], {['R^2 = ', num2str(round(gof.rsquare,2))]},'Location','best','color','none','box','off')


config_plot(f);
hold off;
%%
savepath = [paths.plotFolder filesep 'merge_somitoid12_last_by_first_normalized_outlier'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);

%% histogram for percentage/absolute change between start and end
% f = figure();
% prct_change = (ratio_final-ratio_first)./ratio_first;
% prct_change = (ratio_final-ratio_first);
% 
% histogram(prct_change)
% xlabel('absolute change between start and last frame')
% 
% 
% sum(abs(prct_change)>0.5)
% numel(prct_change)









%%
%%
%% plot for somitoid1
cell_struct_trim = cell_struct_trim_backup;
cell_struct_trim = cell_struct_trim([cell_struct_trim.somitoidID]==1);
%% plot intensityratio vs time 
% small - cyan, large - magenta
ratio_first = arrayfun(@(x)x.intensity_ratio_norm(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);
ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(end-4:end)),cell_struct_trim);

percent = 50;
logical_lower = ratio_first<prctile(ratio_first,percent);
logical_higher = ratio_first>prctile(ratio_first,100-percent);

% calculate average per cell intensity for each frame


f = figure('visible','on');
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    if logical_lower(i)
        pb = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-c','LineWidth',1);
        pcolor = pb;
        pcolor.Color = [pcolor.Color,0.5];
    elseif logical_higher(i)
        pr = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-m','LineWidth',1);
        pcolor = pr;
        pcolor.Color = [pcolor.Color,0.5];
    end

    hold on;

end

% plot median for each group
median_lower = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_lower),'UniformOutput',false))),...
    0:120);
median_higher = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_higher),'UniformOutput',false))),...
    0:120);

pb=plot(0:0.1:12,sgolayfilt(median_lower,5,21),'c-','LineWidth',10);
pr=plot(0:0.1:12,sgolayfilt(median_higher,5,21),'m-','LineWidth',10);

xlabel('Time (hr)')
ylabel('Normalized Mesp2 intensity')
xlim([0,12.2])
% ylim([0,2])
legend([pb,pr], {['lower ' num2str(percent) '% at 1st frame'],['higher ' num2str(percent) '% at 1st frame']},'Location','best','color','none','box','off');
% title('somitoid1')
config_plot(f);

hold off;
%% Generate csv dataframe for plot
time_column = unique(cat(1,cell_struct_trim.frame)./10);
grouping = logical_lower;
df_out = cell(numel(time_column)+2,numel(cell_struct_trim)+1);
df_out{1,1} = 'Group';
df_out{2,1} = 'Time(hr)';
% [df_out{3:end,1}] = time_column_cell{:};
df_out(3:end,1) = num2cell(time_column);
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    x = s.frame./10;
    y_cell = num2cell(s.intensity_ratio_norm);
    [~,idx] = ismember(x,time_column);
    df_out(idx+2,i+1) = y_cell;
    df_out{2,i+1} = ['Cell' num2str(i)];
    if grouping(i)
        df_out{1,i+1} = 'low';
    else
        df_out{1,i+1} = 'high';
    end
end
writecell(df_out, '../excel_data/Fig2f.xls')
%%
savepath = [paths.plotFolder filesep 'somitoid1_intensityratio_vs_time_normalized_50prctile'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);
%% plot intensityratio vs time 
% small - cyan, large - magenta
ratio_first = arrayfun(@(x)x.intensity_ratio_norm(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);
ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(end-4:end)),cell_struct_trim);

percent = 25;
logical_lower = ratio_first<prctile(ratio_first,percent);
logical_higher = ratio_first>prctile(ratio_first,100-percent);

% calculate average per cell intensity for each frame


f = figure('visible','on');
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    if logical_lower(i)
        pb = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-c','LineWidth',1);
        pcolor = pb;
        pcolor.Color = [pcolor.Color,0.5];
    elseif logical_higher(i)
        pr = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-m','LineWidth',1);
        pcolor = pr;
        pcolor.Color = [pcolor.Color,0.5];
    end

    hold on;

end

% plot median for each group
median_lower = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_lower),'UniformOutput',false))),...
    0:120);
median_higher = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_higher),'UniformOutput',false))),...
    0:120);

pb=plot(0:0.1:12,sgolayfilt(median_lower,5,21),'c-','LineWidth',10);
pr=plot(0:0.1:12,sgolayfilt(median_higher,5,21),'m-','LineWidth',10);

xlabel('Time (hr)')
ylabel('Normalized Mesp2 intensity')
xlim([0,12.2])
% ylim([0,2])
legend([pb,pr], {['lower ' num2str(percent) '% at 1st frame'],['higher ' num2str(percent) '% at 1st frame']},'Location','best','color','none','box','off');
% title('somitoid1')
config_plot(f);

hold off;
%%
savepath = [paths.plotFolder filesep 'somitoid1_intensityratio_vs_time_normalized_25prctile'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);

%% mahal
d2 = mahal([ratio_first',ratio_final'],[ratio_first',ratio_final']);
f = figure();
cm = colormap('parula');
scatter(ratio_first',ratio_final',100,d2,'o','LineWidth',2)
c = colorbar;
c.Label.String = 'Mahalanobis Distance';
ylim([0,4]);
xlim([0,4]);
xlabel('Mesp2 at start')
ylabel('Mesp2 at end')
% axis equal;
config_plot(f,c);
hold off;
%%
savepath = [paths.plotFolder filesep 'somitoid1_last_by_first_normalized_mahal'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);
%% last frame vs first frame

[~,rank_first] = ismember(ratio_first(:),unique(ratio_first(:)));
[~,rank_final] = ismember(ratio_final(:),unique(ratio_final(:)));
r_sprman = corr([ratio_first(:),ratio_final(:)],'type','Spearman')

f = figure('visible','on');
cm = colormap(lines(2));
hold on;

xdata = ratio_first;
ydata = ratio_final;
% xdata = ratio_first(ratio_first<0.78);
% ydata = ratio_final(ratio_first<0.78);

trim_logical = d2<1000;


xdata = xdata(trim_logical);
ydata = ydata(trim_logical);

data = plot(xdata,ydata,'o','color',cm(1,:),'MarkerSize',10,'LineWidth',2);
plot(ratio_first(~trim_logical),ratio_final(~trim_logical),'x','color','m','MarkerSize',10,'LineWidth',2)


[fitobj,gof] = fit(xdata(:),ydata(:),fittype({'x','1'}));
x_fit = linspace(0,4,100);
y_fit = feval(fitobj,x_fit);
fitted = plot(x_fit,y_fit,'-','color',cm(2,:), 'LineWidth',2);

lm = fitlm(xdata(:),ydata(:));


xlabel('Mesp2 at start')
ylabel('Mesp2 at end')
ylim([0,4])
xlim([0,4])
legend([fitted], {['R^2 = ', num2str(round(gof.rsquare,2))]},'Location','best','color','none','box','off')


config_plot(f);
hold off;
%%
savepath = [paths.plotFolder filesep 'somitoid1_last_by_first_normalized_all'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);
%% last frame vs first frame

[~,rank_first] = ismember(ratio_first(:),unique(ratio_first(:)));
[~,rank_final] = ismember(ratio_final(:),unique(ratio_final(:)));
r_sprman = corr([ratio_first(:),ratio_final(:)],'type','Spearman')

f = figure('visible','on');
cm = colormap(lines(2));
hold on;

xdata = ratio_first;
ydata = ratio_final;
% xdata = ratio_first(ratio_first<0.78);
% ydata = ratio_final(ratio_first<0.78);

trim_logical = d2<8;


xdata = xdata(trim_logical);
ydata = ydata(trim_logical);

data = plot(xdata,ydata,'o','color',cm(1,:),'MarkerSize',10,'LineWidth',2);
plot(ratio_first(~trim_logical),ratio_final(~trim_logical),'x','color','m','MarkerSize',10,'LineWidth',2)


[fitobj,gof] = fit(xdata(:),ydata(:),fittype({'x','1'}));
x_fit = linspace(0,4,100);
y_fit = feval(fitobj,x_fit);
fitted = plot(x_fit,y_fit,'-','color',cm(2,:), 'LineWidth',2);

lm = fitlm(xdata(:),ydata(:));


xlabel('Mesp2 at start')
ylabel('Mesp2 at end')
ylim([0,4])
xlim([0,4])
legend([fitted], {['R^2 = ', num2str(round(gof.rsquare,2))]},'Location','best','color','none','box','off')


config_plot(f);
hold off;
%%
savepath = [paths.plotFolder filesep 'somitoid1_last_by_first_normalized_outlier'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);















%%
%% plot for somitoid2
cell_struct_trim = cell_struct_trim_backup;
cell_struct_trim = cell_struct_trim([cell_struct_trim.somitoidID]==2);
%% plot intensityratio vs time 
% small - cyan, large - magenta
ratio_first = arrayfun(@(x)x.intensity_ratio_norm(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);
ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(end-4:end)),cell_struct_trim);

percent = 50;
logical_lower = ratio_first<prctile(ratio_first,percent);
logical_higher = ratio_first>prctile(ratio_first,100-percent);

% calculate average per cell intensity for each frame


f = figure('visible','on');
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    if logical_lower(i)
        pb = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-c','LineWidth',1);
        pcolor = pb;
        pcolor.Color = [pcolor.Color,0.5];
    elseif logical_higher(i)
        pr = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-m','LineWidth',1);
        pcolor = pr;
        pcolor.Color = [pcolor.Color,0.5];
    end

    hold on;

end

% plot median for each group
median_lower = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_lower),'UniformOutput',false))),...
    0:120);
median_higher = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_higher),'UniformOutput',false))),...
    0:120);

pb=plot(0:0.1:12,sgolayfilt(median_lower,5,21),'c-','LineWidth',10);
pr=plot(0:0.1:12,sgolayfilt(median_higher,5,21),'m-','LineWidth',10);

xlabel('Time (hr)')
ylabel('Normalized Mesp2 intensity')
xlim([0,12.2])
% ylim([0,2])
legend([pb,pr], {['lower ' num2str(percent) '% at 1st frame'],['higher ' num2str(percent) '% at 1st frame']},'Location','best','color','none','box','off');
% title('somitoid1')
config_plot(f);

hold off;
%% Generate csv dataframe for plot
time_column = unique(cat(1,cell_struct_trim.frame)./10);
grouping = logical_lower;
df_out = cell(numel(time_column)+2,numel(cell_struct_trim)+1);
df_out{1,1} = 'Group';
df_out{2,1} = 'Time(hr)';
% [df_out{3:end,1}] = time_column_cell{:};
df_out(3:end,1) = num2cell(time_column);
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    x = s.frame./10;
    y_cell = num2cell(s.intensity_ratio_norm);
    [~,idx] = ismember(x,time_column);
    df_out(idx+2,i+1) = y_cell;
    df_out{2,i+1} = ['Cell' num2str(i)];
    if grouping(i)
        df_out{1,i+1} = 'low';
    else
        df_out{1,i+1} = 'high';
    end
end
writecell(df_out, '../excel_data/ExtFig2p.xls')
%%
savepath = [paths.plotFolder filesep 'somitoid2_intensityratio_vs_time_normalized_50prctile'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);
%% plot intensityratio vs time 
% small - cyan, large - magenta
ratio_first = arrayfun(@(x)x.intensity_ratio_norm(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);
ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(end-4:end)),cell_struct_trim);

percent = 25;
logical_lower = ratio_first<prctile(ratio_first,percent);
logical_higher = ratio_first>prctile(ratio_first,100-percent);

% calculate average per cell intensity for each frame


f = figure('visible','on');
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    if logical_lower(i)
        pb = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-c','LineWidth',1);
        pcolor = pb;
        pcolor.Color = [pcolor.Color,0.5];
    elseif logical_higher(i)
        pr = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm,5,21),'-m','LineWidth',1);
        pcolor = pr;
        pcolor.Color = [pcolor.Color,0.5];
    end

    hold on;

end

% plot median for each group
median_lower = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_lower),'UniformOutput',false))),...
    0:120);
median_higher = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm(x.frame==y)',...
    cell_struct_trim(logical_higher),'UniformOutput',false))),...
    0:120);

pb=plot(0:0.1:12,sgolayfilt(median_lower,5,21),'c-','LineWidth',10);
pr=plot(0:0.1:12,sgolayfilt(median_higher,5,21),'m-','LineWidth',10);

xlabel('Time (hr)')
ylabel('Normalized Mesp2 intensity')
xlim([0,12.2])
% ylim([0,2])
legend([pb,pr], {['lower ' num2str(percent) '% at 1st frame'],['higher ' num2str(percent) '% at 1st frame']},'Location','best','color','none','box','off');
% title('somitoid1')
config_plot(f);

hold off;
%%
savepath = [paths.plotFolder filesep 'somitoid2_intensityratio_vs_time_normalized_25prctile'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);

%% mahal
d2 = mahal([ratio_first',ratio_final'],[ratio_first',ratio_final']);
f = figure();
cm = colormap('parula');
scatter(ratio_first',ratio_final',100,d2,'o','LineWidth',2)
c = colorbar;
c.Label.String = 'Mahalanobis Distance';
ylim([0,4]);
xlim([0,4]);
xlabel('Mesp2 at start')
ylabel('Mesp2 at end')
% axis equal;
config_plot(f,c);
hold off;
%%
savepath = [paths.plotFolder filesep 'somitoid2_last_by_first_normalized_mahal'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);
%% last frame vs first frame

[~,rank_first] = ismember(ratio_first(:),unique(ratio_first(:)));
[~,rank_final] = ismember(ratio_final(:),unique(ratio_final(:)));
r_sprman = corr([ratio_first(:),ratio_final(:)],'type','Spearman')

f = figure('visible','on');
cm = colormap(lines(2));
hold on;

xdata = ratio_first;
ydata = ratio_final;
% xdata = ratio_first(ratio_first<0.78);
% ydata = ratio_final(ratio_first<0.78);

trim_logical = d2<1000;

xdata = xdata(trim_logical);
ydata = ydata(trim_logical);

data = plot(xdata,ydata,'o','color',cm(1,:),'MarkerSize',10,'LineWidth',2);
plot(ratio_first(~trim_logical),ratio_final(~trim_logical),'x','color','m','MarkerSize',10,'LineWidth',2)


[fitobj,gof] = fit(xdata(:),ydata(:),fittype({'x','1'}));
x_fit = linspace(0,4,100);
y_fit = feval(fitobj,x_fit);
fitted = plot(x_fit,y_fit,'-','color',cm(2,:), 'LineWidth',2);

lm = fitlm(xdata(:),ydata(:));


xlabel('Mesp2 at start')
ylabel('Mesp2 at end')
ylim([0,4])
xlim([0,4])
legend([fitted], {['R^2 = ', num2str(round(gof.rsquare,2))]},'Location','best','color','none','box','off')


config_plot(f);
hold off;
%%
savepath = [paths.plotFolder filesep 'somitoid2_last_by_first_normalized_all'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);
%% last frame vs first frame

[~,rank_first] = ismember(ratio_first(:),unique(ratio_first(:)));
[~,rank_final] = ismember(ratio_final(:),unique(ratio_final(:)));
r_sprman = corr([ratio_first(:),ratio_final(:)],'type','Spearman')

f = figure('visible','on');
cm = colormap(lines(2));
hold on;

xdata = ratio_first;
ydata = ratio_final;
% xdata = ratio_first(ratio_first<0.78);
% ydata = ratio_final(ratio_first<0.78);

trim_logical = d2<8;

xdata = xdata(trim_logical);
ydata = ydata(trim_logical);

data = plot(xdata,ydata,'o','color',cm(1,:),'MarkerSize',10,'LineWidth',2);
plot(ratio_first(~trim_logical),ratio_final(~trim_logical),'x','color','m','MarkerSize',10,'LineWidth',2)


[fitobj,gof] = fit(xdata(:),ydata(:),fittype({'x','1'}));
x_fit = linspace(0,4,100);
y_fit = feval(fitobj,x_fit);
fitted = plot(x_fit,y_fit,'-','color',cm(2,:), 'LineWidth',2);

lm = fitlm(xdata(:),ydata(:));


xlabel('Mesp2 at start')
ylabel('Mesp2 at end')
ylim([0,4])
xlim([0,4])
legend([fitted], {['R^2 = ', num2str(round(gof.rsquare,2))]},'Location','best','color','none','box','off')


config_plot(f);
hold off;
%%
savepath = [paths.plotFolder filesep 'somitoid2_last_by_first_normalized_outlier'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);
























%% plot all data merged
cell_struct_trim = cell_struct_trim_backup;


%% plot surround vs time, but group separately for each somitoid
% small - cyan, large - magenta
percent = 50;
ratio_first = arrayfun(@(x)x.intensity_ratio(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);

somitoidID = arrayfun(@(x)x.somitoidID,cell_struct_trim);


logical_lower = (ratio_first<prctile(ratio_first(somitoidID==1), [100-percent])&somitoidID==1) | ...
    (ratio_first<prctile(ratio_first(somitoidID==2), [100-percent])&somitoidID==2);
logical_higher = (ratio_first>=prctile(ratio_first(somitoidID==1), [percent])&somitoidID==1) | ...
    (ratio_first>=prctile(ratio_first(somitoidID==2), [percent])&somitoidID==2);


median(ratio_first)
f = figure('visible','on');
for i = 1:numel(cell_struct_trim)
    if logical_lower(i)
        pb = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).surrounding_intensity_old,5,21),'-c','LineWidth',1);
        pcolor=pb;
    elseif logical_higher(i)
        pr = plot(cell_struct_trim(i).frame./10, sgolayfilt(cell_struct_trim(i).surrounding_intensity_old,5,21),'-m','LineWidth',1);
        pcolor=pr;
    else
%         pk = plot(cell_struct_trim(i).frame, sgolayfilt(cell_struct_trim(i).surrounding_intensity_norm,20,21),'-k','LineWidth',1);
    end
    pcolor.Color = [pcolor.Color,0.5];
    hold on;
end

% plot median for each group
median_lower = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.surrounding_intensity_old(x.frame==y)',...
    cell_struct_trim(logical_lower),'UniformOutput',false))),...
    0:120);
median_higher = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.surrounding_intensity_old(x.frame==y)',...
    cell_struct_trim(logical_higher),'UniformOutput',false))),...
    0:120);

pb=plot(0:0.1:12,sgolayfilt(median_lower,5,21),'c-','LineWidth',10);
pr=plot(0:0.1:12,sgolayfilt(median_higher,5,21),'m-','LineWidth',10);

xlabel('Time (hr)')
ylabel('Surrounding Intensity')
xlim([0,12.2])
ylim([180,800])
legend([pb,pr], {['lower ',num2str(percent),'% at 1st frame'],['higher ',num2str(percent),'% at 1st frame']},'Location','best','color','none','box','off');

config_plot(f)
hold off;
%% Generate csv dataframe for plot
time_column = unique(cat(1,cell_struct_trim.frame)./10);
grouping = logical_lower;
df_out = cell(numel(time_column)+2,numel(cell_struct_trim)+1);
df_out{1,1} = 'Group';
df_out{2,1} = 'Time(hr)';
% [df_out{3:end,1}] = time_column_cell{:};
df_out(3:end,1) = num2cell(time_column);
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    x = s.frame./10;
    y_cell = num2cell(s.surrounding_intensity_old);
    [~,idx] = ismember(x,time_column);
    df_out(idx+2,i+1) = y_cell;
    df_out{2,i+1} = ['Cell' num2str(i)];
    if grouping(i)
        df_out{1,i+1} = 'low';
    else
        df_out{1,i+1} = 'high';
    end
end
writecell(df_out, '../excel_data/ExtFig2q.xls')

%%
savepath = [paths.plotFolder filesep 'merge_surrounding_by_time_group_separately_plotmedian.eps'];
print(savepath,'-depsc','-loose')
saveas(f,[paths.plotFolder,filesep,'merge_surrounding_by_time_group_separately_plotmedian.png']);
saveas(f,[paths.plotFolder,filesep,'merge_surrounding_by_time_group_separately_plotmedian.fig']);

%% plot difference between 1-end for each, but separately for two somitoids
percent = 50;

change_in_surrounding = arrayfun(@(x)(mean(x.surrounding_intensity_old(end-4:end))-mean(x.surrounding_intensity_old(1:5))),cell_struct_trim);

ratio_first = arrayfun(@(x)x.intensity_ratio(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);

somitoidID = arrayfun(@(x)x.somitoidID,cell_struct_trim);


logical_lower = (ratio_first<prctile(ratio_first(somitoidID==1), [percent])&somitoidID==1) | ...
    (ratio_first<prctile(ratio_first(somitoidID==2), [percent])&somitoidID==2);
logical_higher = (ratio_first>=prctile(ratio_first(somitoidID==1), [100-percent])&somitoidID==1) | ...
    (ratio_first>=prctile(ratio_first(somitoidID==2), [100-percent])&somitoidID==2);
% modify ratio frist to remove the smallest first_frame surrounding
% ratio_first = ratio_first(start_surrounding >= start_surrounding_sorted(10));
% change_in_surrounding = change_in_surrounding(start_surrounding >= start_surrounding_sorted(10));

group_lower = change_in_surrounding(logical_lower);%%%%%*****
group_higher = change_in_surrounding(logical_higher);
surrounding_perc = cat(1,group_lower',group_higher');
grouping_perc = cat(1, zeros(size(group_lower))', ones(size(group_higher))');

f = figure('visible','on');
h = boxplot(surrounding_perc, grouping_perc,'Labels',{['lower ',num2str(percent),'%'],...
    ['higher ',num2str(percent),'%']},'Colors',colormap(lines(1)));
set(h,{'linew'},{3})
[~,pttest] = ttest2(group_higher,...
    group_lower);
[putest,~] = ranksum(group_higher,...
    group_lower);
sigstar({[1,2]}, pttest);
ylabel('Change in surrounding intensity')
title(['t test p value =' num2str(pttest)])
ylim([-350,120])
pttest
config_plot(f);
hold off;
%%
savepath = [paths.plotFolder filesep 'merge_surrounding_difference_group_separately.eps'];
print(savepath,'-depsc','-loose')
saveas(f,[paths.plotFolder,filesep,'merge_surrounding_difference_group_separately.png']);
saveas(f,[paths.plotFolder,filesep,'merge_surrounding_difference_group_separately.fig']);
%% plot surrounding intensity at the start for each group separately
percent = 50;
start_surrounding = arrayfun(@(x)mean(x.surrounding_intensity_old(1:5)),cell_struct_trim);

ratio_first = arrayfun(@(x)x.intensity_ratio(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);

somitoidID = arrayfun(@(x)x.somitoidID,cell_struct_trim);

logical_lower = (ratio_first<prctile(ratio_first(somitoidID==1), [percent])&somitoidID==1) | ...
    (ratio_first<prctile(ratio_first(somitoidID==2), [percent])&somitoidID==2);
logical_higher = (ratio_first>=prctile(ratio_first(somitoidID==1), [100-percent])&somitoidID==1) | ...
    (ratio_first>=prctile(ratio_first(somitoidID==2), [100-percent])&somitoidID==2);


group_lower = start_surrounding(logical_lower);
group_higher = start_surrounding(logical_higher);
surrounding_25 = cat(1,group_lower',group_higher');
grouping_25 = cat(1, zeros(size(group_lower))', ones(size(group_higher))');

f = figure('visible','on');
h = boxplot(surrounding_25, grouping_25,'Labels',{'lower 50%','higher 50%'},'Colors',colormap(lines(1)));
set(h,{'linew'},{3})
[~,pttest] = ttest2(group_higher,...
    group_lower);
[putest,~] = ranksum(group_higher,...
    group_lower);
ylim([280,740])
% f = figure('visible','on');
% boxplot(start_surrounding, ratio_first>=median(ratio_first),'Labels',{'lower 50%','higher 50%'},'Colors',colormap(lines(1)))
% 
% [~,pttest] = ttest2(start_surrounding(ratio_first>=median(ratio_first)),...
%     start_surrounding(ratio_first<median(ratio_first)));
title(['t test p value =' num2str(pttest)])

sigstar({[1,2]}, pttest);
ylabel('Surrounding intensity at 1st frame')


config_plot(f)
hold off;
%%
savepath = [paths.plotFolder filesep 'merge_surrounding_start_group_separately.eps'];
print(savepath,'-depsc','-loose')
saveas(f,[paths.plotFolder,filesep,'merge_surrounding_start_group_separately.png']);
saveas(f,[paths.plotFolder,filesep,'merge_surrounding_start_group_separately.fig']);

%% plot surrounding intensity at the end for each group separately
percent = 50;

end_surrounding = arrayfun(@(x)mean(x.surrounding_intensity_old(end-4:end)),cell_struct_trim);

ratio_first = arrayfun(@(x)x.intensity_ratio(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);


somitoidID = arrayfun(@(x)x.somitoidID,cell_struct_trim);

logical_lower = (ratio_first<prctile(ratio_first(somitoidID==1), [100-percent])&somitoidID==1) | ...
    (ratio_first<prctile(ratio_first(somitoidID==2), [100-percent])&somitoidID==2);
logical_higher = (ratio_first>=prctile(ratio_first(somitoidID==1), [percent])&somitoidID==1) | ...
    (ratio_first>=prctile(ratio_first(somitoidID==2), [percent])&somitoidID==2);

% 25 doesnt mean it has to be 25%
group_lower = end_surrounding(logical_lower);
group_higher = end_surrounding(logical_higher);
surrounding_25 = cat(1,group_lower',group_higher');
grouping_25 = cat(1, zeros(size(group_lower))', ones(size(group_higher))');


f = figure('visible','on');
h = boxplot(surrounding_25, grouping_25,'Labels',{'lower 50%','higher 50%'},'Colors',colormap(lines(1)));
set(h,{'linew'},{3})
[~,pttest] = ttest2(group_higher,...
    group_lower);
% [~,pttest] = ttest2(end_surrounding(ratio_first>=median(ratio_first)),...
%     end_surrounding(ratio_first<median(ratio_first)));
[putest,~] = ranksum(group_higher,...
    group_lower);
sigstar({[1,2]}, pttest);
ylabel('Surrounding intensity at last frame')
title(['t test p value =' num2str(pttest)])
ylim([180,580])

config_plot(f)
hold off;
%%
savepath = [paths.plotFolder filesep 'merge_surrounding_end_group_separately.eps'];
print(savepath,'-depsc','-loose')
saveas(f,[paths.plotFolder,filesep,'merge_surrounding_end_group_separately.png']);
saveas(f,[paths.plotFolder,filesep,'merge_surrounding_end_group_separately.fig']);



%% 
close all;
