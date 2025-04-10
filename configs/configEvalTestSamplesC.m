function [config] = configEvalTestSamplesC()
% no function declaration, the variables are to be imported by the caller

% part 1: the model and experiment-specific configuration
config.num_samples_per_opt_bio_fraction = 1000;
config.opt_bio_fractions = [0.31, 0.37, 0.45, 0.53, 0.60, 0.67, 0.73, 0.79, 0.84, 0.88, 0.92, 0.96];
config.prefixes = string(config.opt_bio_fractions*100);

config.total_samples = 12000; %num_samples_per_opt_bio_fraction * len(opt_bio_fractions)
config.num_slots = 300;
config.slot_size = 40;

config.biomass_rxn_id = "Bio_opt";
config.logsPerHour = 20;
config.simDurationHours = 4;
config.selected_timepoints = [4,8,12,20,28,36,44,52,60,70,81];
config.atomLabelFraction = 0; % or 0.011 natural abundance or 0.008 - C3 plants have a δ13C of −33 to −24‰.

config.mid_name_input = ["Glc[c]:C#1,C#2,C#3,C#4,C#5,C#6.0", "Glc[c]:C#1,C#2,C#3,C#4,C#5,C#6.1", "Glc[c]:C#1,C#2,C#3,C#4,C#5,C#6.2", ...
    "Glc[c]:C#1,C#2,C#3,C#4,C#5,C#6.3", "Glc[c]:C#1,C#2,C#3,C#4,C#5,C#6.4", "Glc[c]:C#1,C#2,C#3,C#4,C#5,C#6.5", "Glc[c]:C#1,C#2,C#3,C#4,C#5,C#6.6", ...
    "Glc[h]:C#1,C#2,C#3,C#4,C#5,C#6.0", "Glc[h]:C#1,C#2,C#3,C#4,C#5,C#6.1", "Glc[h]:C#1,C#2,C#3,C#4,C#5,C#6.2" ...
    "Glc[h]:C#1,C#2,C#3,C#4,C#5,C#6.3", "Glc[h]:C#1,C#2,C#3,C#4,C#5,C#6.4", "Glc[h]:C#1,C#2,C#3,C#4,C#5,C#6.5", "Glc[h]:C#1,C#2,C#3,C#4,C#5,C#6.6"];

% part 2: path and file names for intermediate storage
config.base_dir = fullfile("runtime_data", "run_EvalTestSamplesC");
config.model_file = fullfile(config.base_dir,"preparedAraCore");
config.samples_dir = fullfile(config.base_dir,"samples");
config.scaled_fd_sample_file = fullfile(config.base_dir,"samples_totalfluxlimited_30_50_65_80_90_95_scaled.mat");
config.met_pools_sample_file = fullfile(config.base_dir,"sample_met_pools.mat");
config.split_dir = fullfile(config.base_dir,"split_data");
config.sample_ratios_file = fullfile(config.base_dir,"normalized_met_sample_ratios_mid_TestSamplesAraCoreC.mat");
config.sample_d_ratios_file = fullfile(config.base_dir,"normalized_met_sample_d_ratios_mid_TestSamplesAraCoreC.mat");
config.sample_decomp_ratios_file = fullfile(config.base_dir,"normalized_met_sample_decomp_ratios_mid_TestSamplesAraCoreC.mat");
config.sample_decomp_d_ratios_file = fullfile(config.base_dir,"normalized_met_sample_decomp_d_ratios_mid_TestSamplesAraCoreC.mat");
config.NN_output_dir = fullfile(config.base_dir,"NN_araCore");
config.model_preparation_script = "prepareAraCoreC";
config.all_result_data_file = "all_result_dataC.mat";
config.ci_results_dir = "results_ciC";


% make sure, all folders exist
mkdir(config.base_dir);
mkdir(config.samples_dir);
mkdir(config.split_dir);
mkdir(config.NN_output_dir);
end