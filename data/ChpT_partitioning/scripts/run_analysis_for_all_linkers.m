%{
Keren Lasker
November, 2021
%}

%{
This script is used to run analysis for all linkers.
Parameters the user can set:
resampling_grid_size: how many bins per cell
num_profiles: Number of profiles to consider

Note: signal processing package should be installed.
%}
clear;

%https://www.youtube.com/watch?v=p-MGLfRzBwg
%List all the input experiments

names={'Half','Double_linker','Double_OD','Double_double'};
file_names={'6401_half_Px_img04_08.mat',...
    '6402_double_Px_segmentation_v2.mat',...
    '6452_double_OD_Px_segmentation_v2.mat',...
    '6456_double_double_Px_segmentation.mat'};

popz_background_set={2000,2100,1700,3700};
chpt_background_set={2000,3000,3500,3500};
output_filename={'6401_half','6402_double_linker',...
    '6452_double_OD',...
    '6456_double_double'};

shuffle = @(a)a(randperm(numel(a)));
%Setup the output file. It is going to a matrix with 30 pixels per cell and
%at most 250 profiles per condition.
%The output file is pc_arr_arr. pc is gor partition coefficient
resampling_grid_size=30;
num_profiles=250;
pc_arr_all=zeros(size(names,2),num_profiles);
%create a header for the table
header=[];
for k = 1:size(names,2)
    for j1 = 1:num_profiles
        header=[header , names{k}];
    end
end

cell_lengths_all_conditions = zeros(num_profiles,length(names));

% threshold: Ignore the cell if its max bin intensity (after background
% substraction) is lower than the threshold
% max_len:ignore cells of length longer than max_len
threshold=2000;
max_len=8.5;

%Iterate over conditions, (file per condition)
for k=1:size(file_names,2)
    name = names{k};
    exp_dataset=load(file_names{k}).Experiment;

    curr_len=1;
    curr_prof=1;

    signal_distribution_per_condition_popz = zeros( ...
        resampling_grid_size,num_profiles);

    signal_distribution_per_condition_chpt = zeros( ...
        resampling_grid_size,num_profiles);

     background_popz=popz_background_set{k};
     background_chpt=chpt_background_set{k};

    % Iterate over cells per condition
    for j1=1:size(exp_dataset,1)
    
        exp=exp_dataset(j1).Bacteria;
        bacteria_values = exp;
        lengths = get_cell_lengths(bacteria_values);

        lengths = shuffle(sort(lengths));
        num_profiles_to_use=min([num_profiles size(lengths,1)]);
        cell_lengths_all_conditions(curr_len:curr_len+num_profiles_to_use-1,k)=...
            lengths(1:num_profiles_to_use);

        curr_len=curr_len+num_profiles_to_use;
        get_number_of_grid(bacteria_values);


        [signal_distribution_popz,signal_distribution_chpt] = ...
            run_analysis_for_linker_two_channels(name,...
            bacteria_values,...
            background_popz,...
            background_chpt,threshold,...
            max_len,...
            resampling_grid_size, 2,3);


        n=size(signal_distribution_popz,2);
        signal_distribution_per_condition_popz(:,...
            curr_prof:curr_prof+n-1) = signal_distribution_popz;
        signal_distribution_per_condition_chpt(:,...
            curr_prof:curr_prof+n-1) = signal_distribution_chpt;
        curr_prof=curr_prof+n;
    end %j


    num_cells_per_condition(k)=size(cell_lengths_all_conditions,1);


    % figure;
    % plot(mat2gray(mean(signal_distribution_all_conditions_popz(:,1:size(Experiment,1)),2)));
    % hold;
    % plot(mat2gray(mean(signal_distribution_all_conditions_chpt(:,1:size(Experiment,1)),2)));
    % hold;


    %Now, calculate partition coefficient for this condition
    %First, align the cells to the old pole, by find the regions of large PopZ.

    % Cells are added sequentially, so we are looking for columns that do
    % not have values, i.e. everything is zero.
    % Then, we expect that after the first one of those 
    % cols_with_all_zeros(1), all columns are fully zero and before they
    % are with values. That is the reason we are running the loop from 1 to
    % cols_with_all_zeros(1).

    cols_with_all_zeros = find(all(signal_distribution_per_condition_popz==0));

    pc_arr=zeros(1,cols_with_all_zeros(1));
    for i=1:cols_with_all_zeros(1)
        if i==15
            i
        end
        popz_profile_per_cell=[0;signal_distribution_per_condition_popz(:,i);0];
        chpt_profile_per_cell=[0;signal_distribution_per_condition_chpt(:,i);0];
        %figure;
        [pks,locs, widths, prominences]=findpeaks(popz_profile_per_cell,'SortStr','descend');
        if size(pks,1)==0
            continue
        end
        m1=mean(popz_profile_per_cell);
        r_in=0.0;
        num_inds_used=0;
        %figure;
        %plot(popz_profile_per_cell);
        %hold on;
        for j1=1:size(locs,1)
            if popz_profile_per_cell(locs(j1))/m1 > 2.5 % if the value of PopZ at the peak is at least 2 times higher than the average value of PopZ in the cell
                %find full width at hald max

                maxValue = pks(j1);
                start_ind=max(1,locs(j1)-ceil(widths(j1)/2));
                end_ind=min(size(popz_profile_per_cell,1),locs(j1)+ceil(widths(j1)/2));
%                 if j1==1 
%                     start_ind=1;
%                 else
%                     start_ind=max(1,locs(j1-1));
%                 end
%                 if j1==size(locs,1)
%                     end_ind=size(popz_profile_per_cell,1);
%                 else
%                     end_ind=min(locs(j1+1),size(popz_profile_per_cell,1));
%                 end

                % Find where it's more than half the max.
                aboveHalfMax = popz_profile_per_cell(start_ind:end_ind) > maxValue/2;
                % Get the first and last index where it's more than the half max.
                first_ind = find(aboveHalfMax, 1, 'first');
                last_ind = find(aboveHalfMax, 1, 'last');
                % Draw lines there

                
%                 line([first_ind first_ind], [0 maxValue/2], 'Color', 'r', 'LineWidth', 2);
%                 line([last_ind last_ind], [0 maxValue/2], 'Color', 'r', 'LineWidth', 2);
%                 line([first_ind last_ind], [maxValue/2 maxValue/2], 'Color', 'r', 'LineWidth', 2);

                r_in=r_in+sum(chpt_profile_per_cell(first_ind:last_ind));
                
                num_inds_used = num_inds_used + last_ind-first_ind+1;
            end
        end
        cell_size=size(popz_profile_per_cell,1);
        vol_condensate=(1.0*num_inds_used-1)/(1.0*cell_size);
        vol_rest=(1.0*(cell_size-num_inds_used+1))/(1.0*cell_size);
%         up_v=r_in/num_inds_used;
%         low_v=(sum(chpt_profile_per_cell)-r_in)/(size(popz_profile_per_cell,1)-num_inds_used);

%          up_v=r_in;
%          low_v=(sum(chpt_profile_per_cell)-r_in);

         up_v=r_in/vol_condensate;
         low_v=(sum(chpt_profile_per_cell)-r_in)/vol_rest;

%           up_v=sum(chpt_profile_per_cell)/size(popz_profile_per_cell,1)+r_in/num_inds_used;
%         low_v=sum(chpt_profile_per_cell)/size(popz_profile_per_cell,1)+(sum(chpt_profile_per_cell)-r_in)/(size(popz_profile_per_cell,1)-num_inds_used);
       
         %up_v=r_in/num_inds_used;
         %low_v=(sum(chpt_profile_per_cell)-r_in)/(size(popz_profile_per_cell,1)-num_inds_used);

         
        pc_arr(i)=(up_v)/(low_v);
        
    end
    pc_arr_all(k,1:size(pc_arr,2))=pc_arr;
    
    pc_arr(isnan(pc_arr))=0;
    
    fprintf('%s: %3.3f %d\n',names{k},mean(nonzeros(pc_arr)),size(nonzeros(pc_arr)));
    
end %k

pc_arr_all(pc_arr_all==0)= NaN;
36

