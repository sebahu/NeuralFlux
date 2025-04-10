function [error] = workflowLearnNNs(config_script, start_mid_index, end_mid_index)
    addpath("application_core");
    addpath("configs");
    % the learning can be split up per mid_name_input
    % if only config_script is given, all NNs are learned
    % if only start_mid_index is given, only one NN is learned
    % the indices reference to mid_name_input array from the config script
    error = 0;
    if ~exist('config_script', 'var')
        disp("call with appropriate config, start and end MID index eg. workflowLearnNNs('configMinTest',1,10)");
        error = -1;
        return;
    end

    cfg = eval(config_script);
    
    if ~exist('start_mid_index', 'var')
        start_mid_index = 1;
        end_mid_index = length(cfg.mid_name_input);
    end
    if ~exist('end_mid_index', 'var')
        end_mid_index = start_mid_index;
    end
    if start_mid_index < 1 || start_mid_index > end_mid_index || end_mid_index > length(cfg.mid_name_input)
        disp("Invalid arguments for start/end_mid_index: "+string(start_mid_index)+"/"+string(end_mid_index));
        disp(" MID indices range from 1 to "+string(length(cfg.mid_name_input))+ ", start_mid_index must be <= end_mid_index.");
        return;
    end
    
    learnNNsForOneMIDFromNspLong(strjoin(cfg.mid_name_input(start_mid_index:end_mid_index),";"), cfg.model_file, ...
        cfg.scaled_fd_sample_file, cfg.sample_decomp_ratios_file, cfg.met_pools_sample_file, ...
        cfg.selected_timepoints, cfg.NN_output_dir, [], []);
end
