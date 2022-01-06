%{
Keren Lasker
May, 2020
%}

%{
This script is used to run analysis for all linkers.
Parameters the user can set:
resampling_grid_size: how many bins per cell
num_profiles: Number of profiles to consider

Note: signal processing package should be installed.
%}
clear;
names={'Half','Double_linker','Double_OD','Double_double'};
file_names={'6401_half_Px_img04_08.mat',...
    '6402_double_Px_segmentation.mat',...
    '6452_double_OD_Px_segmentation.mat',...
    '6456_double_double_Px_segmentation.mat'};

popz_background_set={2000,2100,1700,3700};
chpt_background_set={2000,3000,4100,3500};
output_filename={'6401_half','6402_double_linker',...
    '6452_double_OD',...
    '6456_double_double'};


shuffle = @(a)a(randperm(numel(a)));

resampling_grid_size=30;
num_profiles=250;
pc_arr_all=zeros(size(names,2),num_profiles);
%create a header for the table
header=[];
for k = 1:size(names,2)
    for j = 1:num_profiles
        header=[header , names{k}];
    end
end




for k=3:size(file_names,2)
    name = names{k};
    exp_dataset=load(file_names{k}).Experiment;

    curr_len=1;
    curr_prof=1;


    num_cells=zeros(size(names,2));

    signal_distribution_all_conditions_popz = zeros( ...
        resampling_grid_size,num_profiles);

    signal_distribution_all_conditions_chpt = zeros( ...
        resampling_grid_size,num_profiles);

    cell_lengths_all_conditions = zeros(num_profiles,length(names));


    background_popz=popz_background_set{k};
    background_chpt=chpt_background_set{k};

    threshold=2000;
    max_len=8.5;


    for j=1:size(exp_dataset,1)
    
        exp=exp_dataset(j).Bacteria;
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
        signal_distribution_all_conditions_popz(:,...
            curr_prof:curr_prof+n-1) = signal_distribution_popz;
        signal_distribution_all_conditions_chpt(:,...
            curr_prof:curr_prof+n-1) = signal_distribution_chpt;
        curr_prof=curr_prof+n;
    end %j
    num_cells(k)=size(cell_lengths_all_conditions,1);


    % figure;
    % plot(mat2gray(mean(signal_distribution_all_conditions_popz(:,1:size(Experiment,1)),2)));
    % hold;
    % plot(mat2gray(mean(signal_distribution_all_conditions_chpt(:,1:size(Experiment,1)),2)));
    % hold;


    %calculate partition coefficient:
    %First, find the regions of large PopZ.


    cols_with_all_zeros = find(all(signal_distribution_all_conditions_popz==0));

    pc_arr=zeros(1,cols_with_all_zeros(1));
    for i=1:cols_with_all_zeros(1)
        a1=[0;signal_distribution_all_conditions_popz(:,i);0];
        a2=[0;signal_distribution_all_conditions_chpt(:,i);0];
        %figure;
        [pks,locs]=findpeaks(a1,'SortStr','descend');
        if size(pks,1)==0
            continue
        end
        m1=mean(a1);
        r_in=0.0;
        [i,a1(locs(1)),m1];
        for j=1:size(locs,1)
            if a1(locs(j))/m1 > 2.0
                r_in=r_in+sum(a2((max(1,locs(j)-2)):(min(locs(j)+2,size(a2,1)))));
                [i,m1,a1(locs(j)),a1(locs(j))/m1,r_in,sum(a2)];
            end
        end
        
        pc_arr(i)=r_in/(sum(a2)-r_in);
        
    end
    pc_arr_all(k,1:size(pc_arr,2))=pc_arr;
    names{k}
    mean(nonzeros(pc_arr))
    size(pc_arr)
end
36

