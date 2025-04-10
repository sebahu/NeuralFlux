function [error] = workflowSample(config_script, gurobi_path, cobra_path, max_input)
    addpath("application_core");
    addpath("configs");
    error = 0;
    if ~exist('cobra_path', 'var')
        disp("call with appropriate config and paths, eg. workflowSample('configMinTest','~/apps/cobratoolbox/linux64/matlab', '~/apps/cobratoolbox')");
        error = -1;
        return;
    end

    cfg = eval(config_script);
    eval(cfg.model_preparation_script+"(cfg.model_file)");
    
    % from here, all is parameterized
    for fraction_i = 1:length(cfg.opt_bio_fractions)
        sampleModel(cfg.model_file, cfg.biomass_rxn_id, cfg.opt_bio_fractions(fraction_i), ...
            cfg.num_samples_per_opt_bio_fraction, cfg.samples_dir, gurobi_path, cobra_path);
    end
    if exist('max_input', 'var')
        unite_and_scale_fds(cfg.model_file, cfg.samples_dir, cfg.prefixes, cfg.num_samples_per_opt_bio_fraction, cfg.scaled_fd_sample_file, max_input);
    else
        unite_and_scale_fds(cfg.model_file, cfg.samples_dir, cfg.prefixes, cfg.num_samples_per_opt_bio_fraction, cfg.scaled_fd_sample_file);
    end
    split_samples(cfg.scaled_fd_sample_file, cfg.total_samples, cfg.slot_size, cfg.split_dir);
    createMetConcSampling(cfg.model_file, cfg.total_samples, cfg.slot_size, cfg.split_dir, cfg.met_pools_sample_file);
end
