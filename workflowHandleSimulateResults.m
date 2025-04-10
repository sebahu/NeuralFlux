function [error] = workflowHandleSimulateResults(config_script)
    addpath("application_core");
    addpath("configs");
    error = 0;
    if ~exist('config_script', 'var')
        disp("call with appropriate config eg. workflowHandleSimulateResults('configMinTest')");
        error = -1;
        return;
    end

    cfg = eval(config_script);

    combineSampleRatios(cfg.split_dir, cfg.num_slots, cfg.slot_size, cfg.sample_ratios_file, cfg.sample_d_ratios_file)

    decompartmentalizeSampleRatios(cfg.model_file, cfg.total_samples, cfg.slot_size, cfg.split_dir, ...
        cfg.sample_ratios_file, cfg.sample_d_ratios_file, cfg.sample_decomp_ratios_file, cfg.sample_decomp_d_ratios_file);

end
