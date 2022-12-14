%%
%% Export excel for divergence map
df_out = cell(size(divAv)+1);
df_out{1,1} = 'Divergence Map';
df_out{2,1} = 'y';
df_out{1,2} = 'x';
% [df_out{3:end,1}] = time_column_cell{:};
df_out(2:end,2:end) = num2cell(divAv);
if n==1
    writecell(df_out, '../excel_data/ExtFig2r.xls')
elseif n==2
    writecell(df_out, '../excel_data/Fig2j.xls')
end

%% Export excel for ROI statistics
df_out = cell(size(minus,1)+1,2);
df_out{1,1} = 'MESP-';
df_out{1,2} = 'MESP+';
df_out(2:end,1) = num2cell(minus);
df_out(2:end,2) = num2cell(plus);


writecell(df_out, '../excel_data/Fig2k.xls')