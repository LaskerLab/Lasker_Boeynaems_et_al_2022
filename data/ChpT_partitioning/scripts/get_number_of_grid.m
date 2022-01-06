function [grids_sizes] = get_number_of_grid(data)

aa=extractfield(data,'PROFILE_MED');
num_cells = length(aa);
sizes=zeros(num_cells,1);
for i=1:num_cells
    sizes(i)=size(data(i).PROFILE_MED.ch(1).pixel,1);
end
grid_sizes=unique(sort(sizes));
end