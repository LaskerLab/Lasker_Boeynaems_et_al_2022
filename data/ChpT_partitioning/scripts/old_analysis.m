%{
Keren Lasker
May, 2020
%}

%{
This script is used to run analysis for all linkers.
Parameters the user can set:
resampling_grid_size: how many bins per cell
num_profiles: Number of profiles to consider
%}

clear;
load('data_per_linker.mat');

shuffle = @(a)a(randperm(numel(a)));

names = {"Half_linker";"Double_linker"; %"Doulbe_oligomerization";
    "Double_double_222";"Double_oligomerization_121"};
names={"Half_linker";"Double_linker";
    "Double_double_222";"Double_oligomerization_121";
    "PG_half";"PG_full";"DE_half";"DE_full";
    "scr_5R";"scr_5"}%;};

    
names={"wildtype";"Half_linker";"Double_linker";
    "Double_double_222";"Double_oligomerization_121"};

resampling_grid_size=30;
num_profiles = 250;

%create a header for the table
header=[];
for k = 1:length(names)
    for j = 1:num_profiles
        header=[header , names{k}];
    end
end

num_cells=zeros(size(names));

signal_distribution_all_conditions = zeros( ... 
    resampling_grid_size,size(names,1)*num_profiles);

cell_lengths_all_conditions = zeros(num_profiles,length(names));

for k=1:length(names) 
    name = names{k};

    dd = dataperlinker(strcmp(dataperlinker.name, name), :);

    exp=load(dd{1,2});
    if dd{1,3}==1
            exp=exp.Experiment.Bacteria;
    else
        exp=exp.test;
        
    end
    bacteria_values = exp;

    % Combine
    for i=2:size(dd,1)
        dd{i,2}
        exp=load(dd{i,2});
        if dd{i,3}==1
            exp=exp.Experiment.Bacteria;
        else
            exp = exp.test;
        end
        bacteria_values=[bacteria_values;exp];
    end
    if size(dd,1)==1
        i=1;
    end
    
    lengths = get_cell_lengths(bacteria_values);
    lengths = lengths(lengths<dd{i,14});
    lengths = shuffle(sort(lengths));
    num_profiles_to_use=min([num_profiles size(lengths,1)]);
    cell_lengths_all_conditions(1:num_profiles_to_use,k)=...
        lengths(1:num_profiles_to_use);
    csvwrite(append('lengths_',name,'.csv'),lengths);
    
    get_number_of_grid(bacteria_values);
    
    num_cells(k)=size(lengths,1);
    
    % Visualize signal distribution
    %visualize_signal(dd{1,1},bacteria_values,dd{1,3},dd{1,4});
    
    name = dd{1,1};
    background = dd{1,4};
    threshold = dd{1,5};
    max_len=dd{1,14};
    %background=2000;
    signal_distribution = run_analysis_for_linker(name,...
                                bacteria_values,...
                                background,threshold,max_len,...
                                resampling_grid_size);
                            
    if size(signal_distribution,2) <  num_profiles
        size(signal_distribution,2)
    end
    
    a=zeros([resampling_grid_size,num_profiles]);
    num_profiles_to_use=min([num_profiles size(signal_distribution,2)]);
    a(:,1:num_profiles_to_use)= ...
        signal_distribution(:,1:num_profiles_to_use);
    a(a<0)=0;
     signal_distribution_all_conditions(:,...
        (k-1)*num_profiles+1:(k)*num_profiles) = a;
                            
                            
    writematrix(signal_distribution,append('signal_',name,'.csv'),...
                                            'Delimiter','comma');
    signal_mean = mean(signal_distribution,2);
    figure;
    plot(signal_mean);

end
writematrix(signal_distribution_all_conditions,'all_conditions.csv',...
    'Delimiter','comma');
writematrix(header,'all_conditions_header.csv','Delimiter','comma');
csvwrite('lengths_all_conditions.csv',cell_lengths_all_conditions);
num_cells

