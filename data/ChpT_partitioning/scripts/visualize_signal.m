function visualize_signal(name,data,background, threshold)
%visualize_signal

aa=extractfield(data,'PROFILE_MED');
num_cells = length(aa);
sizes=zeros(num_cells,1);
for i=1:num_cells
    sizes(i)=size(data(i).PROFILE_MED.ch(2).pixel,1);
end
max_cell_length = max(sizes);
min_cell_length = min(sizes);

data_resampled1=zeros(min_cell_length-1,num_cells);

figure;
for i=1:num_cells
    subplot(ceil(num_cells/5),5,i);
    x=mean(aa{1,i}.ch(2).pixel-background,2);
    
%     if mean(x)<threshold
%         continue
%     end
    %data_(:,i) = padarray(x,[max_cell_length-size(x,1)],'post');
    rdata = resample(x,min_cell_length-1,size(x,1));
    if sum(rdata(1:5)) < sum(rdata(end-5:end))
        rdata = flip(rdata);
    end
    data_resampled1(:,i) = rdata;
    plot(rdata);
    title(mean(x));
end

end

