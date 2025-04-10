function [] = combineSampleRatios(split_dir, slot_num, split_size, output_ratios_file, output_d_ratios_file)
    load(string(split_dir)+"/mid_ratios_100_1.mat","data");
    sample_ratios = zeros(size(data,1),size(data,2),slot_num*split_size,'single');
    %sample_d_ratios = zeros(size(data,1),size(data,2),slot_num*split_size,'single');
    for slot_i = 1:slot_num
        slot_start_i = 1+(slot_i-1)*split_size;
        load(string(split_dir)+"/mid_ratios_100_"+string(slot_start_i)+".mat","data");
        if size(data,1) ~= size(sample_ratios,1)
            continue
        end
        sample_ratios(:,:,slot_start_i:(slot_start_i+split_size-1)) = data;
        %load(string(split_dir)+"/mid_d_ratios_100_"+string(slot_start_i)+".mat","data");
        %sample_d_ratios(:,:,slot_start_i:(slot_start_i+split_size-1)) = data;
    end
    save(output_ratios_file, "sample_ratios","-v7.3");
    %save(output_d_ratios_file, "sample_d_ratios","-v7.3");
end