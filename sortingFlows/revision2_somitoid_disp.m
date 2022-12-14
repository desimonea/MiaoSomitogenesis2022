%% reinitialize
close all;
clear;
clc;
addpath(genpath('./functions/'))
%% define paths 
paths = [];
paths.directory = 'D:\BIO\PhD\ditalia\somitoid\data\image\Somitoid_single cell tracking\';
paths.expname1 = 'Somitoid1_ 72-84hr_6min interval_2021_09_10__10_15_25';
paths.expname2 = 'Somitoid2_ 72-84hr_6min interval_2021_09_18__15_24_26';
paths.dataname1 = [paths.expname1,'-data'];
paths.dataname2 = [paths.expname2,'-data'];
paths.plotFolder = 'D:\BIO\PhD\ditalia\somitoid\github\ditalia-somitoid\code_revision2_alvin\figure_revision2_somitoid';
mkdir(paths.plotFolder)
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























%% plot all data merged
cell_struct_trim = cell_struct_trim_backup;


%% plot surrounddiff_final, measured by filt image of last frame vs time, but group separately for each somitoid !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% calculate
load(['./cell_struct/','somitoid1_surrounding_last.mat'], "somitoid1_surrounding_last")
load(['./cell_struct/','somitoid2_surrounding_last.mat'], "somitoid2_surrounding_last")
surrounding_final_at_first = {};
for i = 1:numel(cell_struct_trim)
    s = cell_struct_trim(i);
    if s.somitoidID ==1
        image_surrounding_last = somitoid1_surrounding_last;
    elseif s.somitoidID ==2
        image_surrounding_last = somitoid2_surrounding_last;
    end
    xy_coor_first = s.coordinates(1,1:2)+1;
    surrounding_final_at_first{i} = double(image_surrounding_last(xy_coor_first(1),xy_coor_first(2)));
end
[cell_struct_trim.surrounding_final_at_first] = surrounding_final_at_first{:};


%%
cm = turbo(10);
cm = parula(10);
cm = copper(10);
cm = brewermap(10,'RdYlGn');
c_low = cm(3,:);
c_high = cm(9,:);
% small - cyan, large - magenta
percent = 50;
ratio_first = arrayfun(@(x)x.intensity_ratio(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio_norm(1:5)),cell_struct_trim);

somitoidID = arrayfun(@(x)x.somitoidID,cell_struct_trim);


logical_lower1 = (ratio_first<prctile(ratio_first(somitoidID==1), [100-percent])&somitoidID==1) | ...
    (ratio_first<prctile(ratio_first(somitoidID==2), [100-percent])&somitoidID==2);
logical_higher1 = (ratio_first>=prctile(ratio_first(somitoidID==1), [percent])&somitoidID==1) | ...
    (ratio_first>=prctile(ratio_first(somitoidID==2), [percent])&somitoidID==2);
% logical_lower1 = ratio_first<prctile(ratio_first, [100-percent]);
% logical_higher1 = ratio_first>=prctile(ratio_first, [100-percent]);


ratio_first = [cell_struct_trim.surrounding_final_at_first];
logical_lower2 = (ratio_first<prctile(ratio_first(somitoidID==1), [100-percent])&somitoidID==1) | ...
    (ratio_first<prctile(ratio_first(somitoidID==2), [100-percent])&somitoidID==2);
logical_higher2 = (ratio_first>=prctile(ratio_first(somitoidID==1), [percent])&somitoidID==1) | ...
    (ratio_first>=prctile(ratio_first(somitoidID==2), [percent])&somitoidID==2);
% logical_lower2 = ratio_first<prctile(ratio_first, [100-percent]);
% logical_higher2 = ratio_first>=prctile(ratio_first, [100-percent]);

logical_lower = logical_lower1==logical_lower2;

logical_higher = logical_lower1~=logical_lower2;

median(ratio_first)
f = figure('visible','on');
for i = 1:numel(cell_struct_trim)
    ydata_func = @(s)sqrt((s.coordinates(:,1:2)-s.coordinates(1,1:2)).^2*[1;1]);
    ydata = ydata_func(cell_struct_trim(i)).*xy_pxsize;
    if logical_lower(i)
        pb = plot(cell_struct_trim(i).frame./10, sgolayfilt(ydata,5,21),'-','color',c_low,'LineWidth',1);
        pcolor=pb;
    elseif logical_higher(i)
        pr = plot(cell_struct_trim(i).frame./10, sgolayfilt(ydata,5,21),'-','color',c_high,'LineWidth',1);
        pcolor=pr;
    else
%         pk = plot(cell_struct_trim(i).frame, sgolayfilt(cell_struct_trim(i).surrounding_intensity_norm,20,21),'-k','LineWidth',1);
    end
    pcolor.Color = [pcolor.Color,0.5];
    hold on;
end

% plot median for each group
median_lower = arrayfun(@(y)mean(cell2mat(arrayfun(@(s)sqrt((s.coordinates(s.frame==y,1:2)-s.coordinates(1,1:2)).^2*[1;1])',...
    cell_struct_trim(logical_lower),'UniformOutput',false))),...
    0:120);
median_higher = arrayfun(@(y)mean(cell2mat(arrayfun(@(s)sqrt((s.coordinates(s.frame==y,1:2)-s.coordinates(1,1:2)).^2*[1;1])',...
    cell_struct_trim(logical_higher),'UniformOutput',false))),...
    0:120);

pb=plot(0:0.1:12,sgolayfilt(median_lower.*xy_pxsize,5,21),'-','color',c_low,'LineWidth',10);
pr=plot(0:0.1:12,sgolayfilt(median_higher.*xy_pxsize,5,21),'-','color',c_high,'LineWidth',10);

xlabel('Time (hr)')
ylabel('Displacement (\mum)')
xlim([0,12.2])
% ylim([180,800])
legend([pb,pr], {['Correct at 1st frame'],['Wrong at 1st frame']},'Location','best','color','none','box','off');

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
    ydata_func = @(s)sqrt((s.coordinates(:,1:2)-s.coordinates(1,1:2)).^2*[1;1]);
    ydata = ydata_func(s).*xy_pxsize;
    y_cell = num2cell(ydata);
    [~,idx] = ismember(x,time_column);
    df_out(idx+2,i+1) = y_cell;
    df_out{2,i+1} = ['Cell' num2str(i)];
    if grouping(i)
        df_out{1,i+1} = 'correct';
    else
        df_out{1,i+1} = 'wrong';
    end
end
writecell(df_out, '../excel_data/ExtFig2s.xls')
%%
savepath = [paths.plotFolder filesep 'displacement_by_time'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);



%% plot displacement at the end for each group separately
percent = 50;

end_displacement = arrayfun(@(s)sqrt((s.coordinates(end,1:2)-s.coordinates(1,1:2)).^2*[1;1]),cell_struct_trim).*xy_pxsize;



somitoidID = arrayfun(@(x)x.somitoidID,cell_struct_trim);



% 25 doesnt mean it has to be 25%
group_lower = end_displacement(logical_lower);
group_higher = end_displacement(logical_higher);
surrounding_25 = cat(1,group_lower',group_higher');
grouping_25 = cat(1, zeros(size(group_lower))', ones(size(group_higher))');


f = figure('visible','on');
h = boxplot(surrounding_25, grouping_25,'Labels',{'Correct','Wrong'},'Colors',colormap(lines(1)));
set(h,{'linew'},{3})
[~,pttest] = ttest2(group_higher,...
    group_lower);
% [~,pttest] = ttest2(end_surrounding(ratio_first>=median(ratio_first)),...
%     end_surrounding(ratio_first<median(ratio_first)));
[putest,~] = ranksum(group_higher,...
    group_lower);
sigstar({[1,2]}, pttest);
ylabel('Displacement at last frame (\mum)')
title(['t test p value =' num2str(pttest)])
% ylim([180,580])

config_plot(f);
hold off;
%%
savepath = [paths.plotFolder filesep 'displacement_last_frame'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);





















%%
%%
%% Mesp2 intensity

%% tracking by time
cm = turbo(10);
cm = parula(10);
cm = copper(10);
cm = brewermap(10,'RdYlGn');
c_low = cm(3,:);
c_high = cm(9,:);
% small - cyan, large - magenta
percent = 50;
ratio_first = arrayfun(@(x)x.intensity_ratio_norm(1),cell_struct_trim);
ratio_first = arrayfun(@(x)median(x.intensity_ratio(1:5)),cell_struct_trim);

somitoidID = arrayfun(@(x)x.somitoidID,cell_struct_trim);


logical_lower1 = (ratio_first<prctile(ratio_first(somitoidID==1), [100-percent])&somitoidID==1) | ...
    (ratio_first<prctile(ratio_first(somitoidID==2), [100-percent])&somitoidID==2);
logical_higher1 = (ratio_first>=prctile(ratio_first(somitoidID==1), [percent])&somitoidID==1) | ...
    (ratio_first>=prctile(ratio_first(somitoidID==2), [percent])&somitoidID==2);

ratio_first = [cell_struct_trim.surrounding_final_at_first];
logical_lower2 = (ratio_first<prctile(ratio_first(somitoidID==1), [100-percent])&somitoidID==1) | ...
    (ratio_first<prctile(ratio_first(somitoidID==2), [100-percent])&somitoidID==2);
logical_higher2 = (ratio_first>=prctile(ratio_first(somitoidID==1), [percent])&somitoidID==1) | ...
    (ratio_first>=prctile(ratio_first(somitoidID==2), [percent])&somitoidID==2);

logical_lower = logical_lower1==logical_lower2;

logical_higher = logical_lower1~=logical_lower2;

median(ratio_first)
f = figure('visible','on');
for i = 1:numel(cell_struct_trim)
    ydata_func = @(s)sqrt((s.coordinates(:,1:2)-s.coordinates(1,1:2)).^2*[1;1]);
    ydata = cell_struct_trim(i).intensity_ratio_norm;
    if logical_lower(i)
        pb = plot(cell_struct_trim(i).frame./10, sgolayfilt(ydata,5,21),'-','color',c_low,'LineWidth',1);
        pcolor=pb;
    elseif logical_higher(i)
        pr = plot(cell_struct_trim(i).frame./10, sgolayfilt(ydata,5,21),'-','color',c_high,'LineWidth',1);
        pcolor=pr;
    else
%         pk = plot(cell_struct_trim(i).frame, sgolayfilt(cell_struct_trim(i).surrounding_intensity_norm,20,21),'-k','LineWidth',1);
    end
    pcolor.Color = [pcolor.Color,0.5];
    hold on;
end

% plot median for each group
median_lower = arrayfun(@(y)mean(cell2mat(arrayfun(@(s)s.intensity_ratio_norm(s.frame==y)',...
    cell_struct_trim(logical_lower),'UniformOutput',false))),...
    0:120);
median_higher = arrayfun(@(y)mean(cell2mat(arrayfun(@(s)s.intensity_ratio_norm(s.frame==y)',...
    cell_struct_trim(logical_higher),'UniformOutput',false))),...
    0:120);

pb=plot(0:0.1:12,sgolayfilt(median_lower.*xy_pxsize,5,21),'-','color',c_low,'LineWidth',10);
pr=plot(0:0.1:12,sgolayfilt(median_higher.*xy_pxsize,5,21),'-','color',c_high,'LineWidth',10);

xlabel('Time (hr)')
ylabel('Normalized Meps2 Intensity (A.U.)')
xlim([0,12.2])
% ylim([180,800])
legend([pb,pr], {['Correct at 1st frame'],['Wrong at 1st frame']},'Location','best','color','none','box','off');

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
        df_out{1,i+1} = 'correct';
    else
        df_out{1,i+1} = 'wrong';
    end
end
writecell(df_out, '../excel_data/ExtFig2r.xls')
%%
savepath = [paths.plotFolder filesep 'correct_and_wrong_mesp2_by_time.eps'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);



%% plot mesp2 at the end for each group separately
percent = 50;

end_displacement = arrayfun(@(s)median(s.intensity_ratio_norm(end-4:end)),cell_struct_trim);



somitoidID = arrayfun(@(x)x.somitoidID,cell_struct_trim);



% 25 doesnt mean it has to be 25%
group_lower = end_displacement(logical_lower);
group_higher = end_displacement(logical_higher);
surrounding_25 = cat(1,group_lower',group_higher');
grouping_25 = cat(1, zeros(size(group_lower))', ones(size(group_higher))');


f = figure('visible','on');
h = boxplot(surrounding_25, grouping_25,'Labels',{'Correct','Wrong'},'Colors',colormap(lines(1)));
set(h,{'linew'},{3})
[~,pttest] = ttest2(group_higher,...
    group_lower);
% [~,pttest] = ttest2(end_surrounding(ratio_first>=median(ratio_first)),...
%     end_surrounding(ratio_first<median(ratio_first)));
[putest,~] = ranksum(group_higher,...
    group_lower);
sigstar({[1,2]}, pttest);
ylabel('Normalized Meps2 Last Frame')
title(['t test p value =' num2str(pttest)])
% ylim([180,580])

config_plot(f);
hold off;

%%
savepath = [paths.plotFolder filesep 'correct_and_wrong_mesp2_last_frame.eps'];
print(savepath,'-depsc','-loose')
saveas(f,[savepath '.png']);
saveas(f,[savepath '.fig']);


%% 
close all;