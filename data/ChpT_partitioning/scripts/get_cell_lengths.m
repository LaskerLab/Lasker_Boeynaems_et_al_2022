function [lengths] = get_cell_lengths(data)
   %Save the experiment from MicrobeJ
    a=extractfield(data,'SHAPE');
    b=cell2mat(a);
    lengths = [b(:).length]';
end

