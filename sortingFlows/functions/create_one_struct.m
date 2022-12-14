function [cell_struct_i] = create_one_struct(lineage_singlecell,cell_count,somitoidID,division)
%CREATE_ONE_STRUCT Summary of this function goes here
%   Detailed explanation goes here
arguments
    lineage_singlecell
    cell_count
    somitoidID
    division = [];
    
end
cell_struct_i.somitoidID = somitoidID;
cell_struct_i.cellId = cell_count;
if any(ismember(lineage_singlecell.Properties.VariableNames,'lineageId_haem'))
    cell_struct_i.lineageId = lineage_singlecell.lineageId_haem(1);
else
    cell_struct_i.lineageId = lineage_singlecell.lineageId(1);
end
cell_struct_i.frame = lineage_singlecell.frame;
cell_struct_i.Mean_Intensity_0 = lineage_singlecell.Mean_Intensity_0;
cell_struct_i.Mean_Intensity_1 = lineage_singlecell.Mean_Intensity_1;
cell_struct_i.Size_in_pixels_0 = lineage_singlecell.Size_in_pixels_0;
cell_struct_i.coordinates = lineage_singlecell{:,["Center_of_the_object_0",...
    "Center_of_the_object_1","Center_of_the_object_2"]};
% frame of division - marked as first frame with two daughter cells

cell_struct_i.division = unique(division.frame);

    

end

