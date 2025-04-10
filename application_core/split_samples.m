
function [] = split_samples(input_name, total_samples, split_size, split_dir)
    load(input_name, 'total_samples_scaled');
    for i = 1:split_size:total_samples
       samples_split_scaled = total_samples_scaled(:,i:(i+split_size-1));
       save(fullfile(string(split_dir),"samples_scaled_"+string(i)+".mat"), 'samples_split_scaled', '-v7.3');
    end
end
