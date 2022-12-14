%% reinitialize
close all;
clear;
clc;
addpath(genpath('./functions/'))
%% define paths 
paths = [];
paths.directory = 'D:\BIO\PhD\ditalia\somitoid\data\image\Segmentoid_H2B spiking\';
paths.expname1 = 'ROI_H2B spiking_Mesp2_Segmentoid_Day4 start_6min interval_2022_04_03__16_42_48_Subset';
paths.expname2 = 'ROID_2022_05_12__23_33_17_Subset-D';
paths.expname3 = 'ROIA_2022_05_12__23_33_17_Subset-A';

paths.dataname1 = [paths.expname1,'-data'];
paths.dataname2 = [paths.expname2,'-data'];
paths.dataname3 = [paths.expname3,'-data'];
paths.plotFolder = 'D:\BIO\PhD\ditalia\somitoid\github\ditalia-somitoid\code_revision2_alvin\figure_revision2_segmentoid';
mkdir(paths.plotFolder);
% define pixelsize
xy_pxsize = 0.692;
z_pxsize = 8;
t_size = 6; % min

%% load and merge data
%% load data
%% load data
% cell_struct_trim1 = load(['./cell_struct_mid/',paths.dataname1,'_cell_struct_trim.mat']);
% cell_struct_trim2 = load(['./cell_struct_mid/',paths.dataname2,'_cell_struct_trim.mat']);
% cell_struct_trim3 = load(['./cell_struct_mid/',paths.dataname3,'_cell_struct_trim.mat']);
% 
cell_struct_trim1 = load(['./cell_struct/',paths.dataname1,'_cell_struct_trim.mat']);
cell_struct_trim2 = load(['./cell_struct/',paths.dataname2,'_cell_struct_trim.mat']);
cell_struct_trim3 = load(['./cell_struct/',paths.dataname3,'_cell_struct_trim.mat']);


cell_struct_trim = cat(2,cell_struct_trim1.cell_struct_trim, cell_struct_trim2.cell_struct_trim, cell_struct_trim3.cell_struct_trim);


%% just look at one somitoid (may not be useful)
% cell_struct_trim = cell_struct_trim1.cell_struct_trim;
% cell_struct_trim = cell_struct_trim2.cell_struct_trim;

%% mannual curation
bad_lineage2 = [273,373,380,396,398,565,761,1063,760,770,771,774];
bad_lineage3 = [505,660,867];
% remove segment3 from segmentoid1 because time window too short
curation_logical = (~ismember([cell_struct_trim.segment_num],[3])&[cell_struct_trim.somitoidID]==101) | ...
    (~ismember([cell_struct_trim.lineageId],bad_lineage2)&[cell_struct_trim.somitoidID]==102) | ...
    (~ismember([cell_struct_trim.lineageId],bad_lineage3)&[cell_struct_trim.somitoidID]==103);

cell_struct_trim = cell_struct_trim(curation_logical);



19,28,17









%% plot intensity normalized by average for each segment
% first
smallest_final = min(arrayfun(@(x)x.frame(end)-x.segment_t0,cell_struct_trim));
ratio_first = arrayfun(@(x)x.intensity_ratio_norm(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);

% final
smallest_final = min(arrayfun(@(x)x.frame(end)-x.segment_t0,cell_struct_trim));

ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(end-4:end)),cell_struct_trim);
ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(smallest_final-4:smallest_final)),cell_struct_trim);
ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(...
    (x.frame-x.segment_t0<=smallest_final)&(x.frame-x.segment_t0>=smallest_final-4)...
    )),cell_struct_trim);

% start plotting
percent = 50;
logical_lower = ratio_first<prctile(ratio_first,percent);
logical_higher = ratio_first>prctile(ratio_first,100-percent);

f = figure('visible','on');
for i = 1:numel(cell_struct_trim)
%     if ratio_first(i) < median(ratio_first)
%         pb = plot((cell_struct_trim(i).frame(1:end)-cell_struct_trim(i).frame(1))./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm(1:end),5,9),'-c','LineWidth',1);
%     else
%         pr = plot((cell_struct_trim(i).frame(1:end)-cell_struct_trim(i).frame(1))./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm(1:end),5,9),'-m','LineWidth',1);
%     end

% restrict end point
    s = cell_struct_trim(i);
    [~,idx_final] = min(abs(cell_struct_trim(i).frame-cell_struct_trim(i).segment_t0-smallest_final));
    if logical_lower(i)
        pb = plot((cell_struct_trim(i).frame(s.frame-s.segment_t0<=smallest_final)-cell_struct_trim(i).segment_t0)./10,...
            sgolayfilt(cell_struct_trim(i).intensity_ratio_norm(s.frame-s.segment_t0<=smallest_final),5,21),'-c','LineWidth',1);
        pb.Color = [pb.Color,0.5];
    elseif logical_higher(i)
        pr = plot((cell_struct_trim(i).frame(s.frame-s.segment_t0<=smallest_final)-cell_struct_trim(i).segment_t0)./10,...
            sgolayfilt(cell_struct_trim(i).intensity_ratio_norm(s.frame-s.segment_t0<=smallest_final),5,21),'-m','LineWidth',1);
        pr.Color = [pr.Color,0.5];
    end


    hold on;

end

% plot median for each group
median_lower = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm((x.frame-x.frame(1))==y)',...
    cell_struct_trim(logical_lower),'UniformOutput',false))),...
    0:67);
median_higher = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm((x.frame-x.frame(1))==y)',...
    cell_struct_trim(logical_higher),'UniformOutput',false))),...
    0:67);
pb=plot(0:0.1:6.7,sgolayfilt(median_lower,5,21),'c-','LineWidth',10);
pr=plot(0:0.1:6.7,sgolayfilt(median_higher,5,21),'m-','LineWidth',10);

xlabel('Time (hr)')
ylabel('Normalized Mesp2 intensity')
xlim([0,72./10])
% ylim([0,1.5])
legend([pb,pr], {['lower ' num2str(percent) '% at 1st frame'],['higher ' num2str(percent) '% at 1st frame']},'Location','best','color','none','box','off');
% title('segment1')
config_plot(f)

hold off;
%% Generate csv dataframe for plot
time_column = (0:smallest_final)'/10;
grouping = logical_lower;
df_out = cell(numel(time_column)+2,numel(cell_struct_trim)+1);
df_out{1,1} = 'Group';
df_out{2,1} = 'Time(hr)';
% [df_out{3:end,1}] = time_column_cell{:};
df_out(3:end,1) = num2cell(time_column);
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    x = (s.frame(s.frame-s.segment_t0<=smallest_final)-s.segment_t0)./10;
    ydata = s.intensity_ratio_norm(s.frame-s.segment_t0<=smallest_final);
    y_cell = num2cell(ydata);
    [~,idx] = ismember(x,time_column);
    df_out(idx+2,i+1) = y_cell;
    df_out{2,i+1} = ['Cell' num2str(i)];
    if grouping(i)
        df_out{1,i+1} = 'low';
    else
        df_out{1,i+1} = 'high';
    end
end
writecell(df_out, '../excel_data/Fig4f.xls')
%%
savepath = [paths.plotFolder filesep 'merge_segmentoid123_intensityratio_vs_time_normalized_50prctile'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);
%% plot intensity normalized by average for each segment
% first
smallest_final = min(arrayfun(@(x)x.frame(end)-x.segment_t0,cell_struct_trim));
ratio_first = arrayfun(@(x)x.intensity_ratio_norm(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);

% final
smallest_final = min(arrayfun(@(x)x.frame(end)-x.segment_t0,cell_struct_trim));

ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(end-4:end)),cell_struct_trim);
ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(smallest_final-4:smallest_final)),cell_struct_trim);
ratio_final = arrayfun(@(x)median(x.intensity_ratio_norm(...
    (x.frame-x.segment_t0<=smallest_final)&(x.frame-x.segment_t0>=smallest_final-4)...
    )),cell_struct_trim);

% start plotting
percent = 25;
logical_lower = ratio_first<prctile(ratio_first,percent);
logical_higher = ratio_first>prctile(ratio_first,100-percent);

f = figure('visible','on');
for i = 1:numel(cell_struct_trim)
%     if ratio_first(i) < median(ratio_first)
%         pb = plot((cell_struct_trim(i).frame(1:end)-cell_struct_trim(i).frame(1))./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm(1:end),5,9),'-c','LineWidth',1);
%     else
%         pr = plot((cell_struct_trim(i).frame(1:end)-cell_struct_trim(i).frame(1))./10, sgolayfilt(cell_struct_trim(i).intensity_ratio_norm(1:end),5,9),'-m','LineWidth',1);
%     end

% restrict end point
    s = cell_struct_trim(i);
    [~,idx_final] = min(abs(cell_struct_trim(i).frame-cell_struct_trim(i).segment_t0-smallest_final));
    if logical_lower(i)
        pb = plot((cell_struct_trim(i).frame(s.frame-s.segment_t0<=smallest_final)-cell_struct_trim(i).segment_t0)./10,...
            sgolayfilt(cell_struct_trim(i).intensity_ratio_norm(s.frame-s.segment_t0<=smallest_final),5,21),'-c','LineWidth',1);
        pb.Color = [pb.Color,0.5];
    elseif logical_higher(i)
        pr = plot((cell_struct_trim(i).frame(s.frame-s.segment_t0<=smallest_final)-cell_struct_trim(i).segment_t0)./10,...
            sgolayfilt(cell_struct_trim(i).intensity_ratio_norm(s.frame-s.segment_t0<=smallest_final),5,21),'-m','LineWidth',1);
        pr.Color = [pr.Color,0.5];
    end


    hold on;

end

% plot median for each group
median_lower = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm((x.frame-x.frame(1))==y)',...
    cell_struct_trim(logical_lower),'UniformOutput',false))),...
    0:67);
median_higher = arrayfun(@(y)median(cell2mat(arrayfun(@(x)x.intensity_ratio_norm((x.frame-x.frame(1))==y)',...
    cell_struct_trim(logical_higher),'UniformOutput',false))),...
    0:67);
pb=plot(0:0.1:6.7,sgolayfilt(median_lower,5,21),'c-','LineWidth',10);
pr=plot(0:0.1:6.7,sgolayfilt(median_higher,5,21),'m-','LineWidth',10);

xlabel('Time (hr)')
ylabel('Normalized Mesp2 intensity')
xlim([0,72./10])
% ylim([0,1.5])
legend([pb,pr], {['lower ' num2str(percent) '% at 1st frame'],['higher ' num2str(percent) '% at 1st frame']},'Location','best','color','none','box','off');
% title('segment1')
config_plot(f)

hold off;
%%
savepath = [paths.plotFolder filesep 'merge_segmentoid123_intensityratio_vs_time_normalized_25prctile'];
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
savepath = [paths.plotFolder filesep 'merge_segmentoid123_last_by_first_normalized_mahal'];
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
savepath = [paths.plotFolder filesep 'merge_segmentoid123_last_by_first_normalized_all'];
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
% xdata = rank_first;
% ydata = rank_final;

trim_logical = d2<8;

xdata = xdata(trim_logical);
ydata = ydata(trim_logical);

data = plot(xdata,ydata,'o','color',cm(1,:),'MarkerSize',10,'LineWidth',2);
plot(ratio_first(~trim_logical),ratio_final(~trim_logical),'x','color','m','MarkerSize',10,'LineWidth',2)


[fitobj,gof] = fit(xdata(:),ydata(:),fittype({'x','1'}));
x_fit = linspace(0,max(4),100);
y_fit = feval(fitobj,x_fit);
fitted = plot(x_fit,y_fit,'-','color',cm(2,:), 'LineWidth',2);
% confidence
% pred_bound = predint(fitobj,x_fit,0.68,'observation','off');
% p_bound = plot(x_fit,pred_bound,'m--','LineWidth',2);

lm = fitlm(xdata(:),ydata(:));


xlabel('Mesp2 at start')
ylabel('Mesp2 at end')
ylim([0,4])
xlim([0,4])
legend([fitted], {['R^2 = ', num2str(round(gof.rsquare,2))]},'Location','best','color','none','box','off')


config_plot(f);
hold off;

%%
savepath = [paths.plotFolder filesep 'merge_segmentoid123_last_by_first_normalized_outlier'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);
%% histogram for percentage change between start and end
f = figure();
mesp2_change = (ratio_final-ratio_first)./ratio_first;
mesp2_change = (ratio_final-ratio_first);

histogram(mesp2_change)
xlabel('absolute change between start and last frame')


sum(abs(mesp2_change)>0.5)
numel(mesp2_change)


%%
close all;