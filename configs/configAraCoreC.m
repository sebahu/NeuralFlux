function [config] = configAraCoreC()
% no function declaration, the variables are to be imported by the caller

% part 1: the model and experiment-specific configuration
config.num_samples_per_opt_bio_fraction = 50000;
config.opt_bio_fractions = [0.3,0.4,0.5,0.575,0.65,0.725,0.8,0.85,0.9,0.925,0.95,0.975];
config.prefixes = string(config.opt_bio_fractions*100);

config.total_samples = 600000; %num_samples_per_opt_bio_fraction * len(opt_bio_fractions)
config.num_slots = 6000;
config.slot_size = 100;

config.biomass_rxn_id = "Bio_opt";
config.logsPerHour = 20;
config.simDurationHours = 4;
%config.selected_timepoints = [5,12,23,52,81];
config.selected_timepoints = [3,5,8,12,17,23,37,52,66,81];
config.atomLabelFraction = 0; % or 0.011 natural abundance or 0.008 - C3 plants have a δ13C of −33 to −24‰.

config.met_mid_sds = ["2PGA", 0.01
    "6PG", 0.01
    "Ala", 0.01
    "Arg", 0.01
    "Asn", 0.01
    "Asp", 0.01
    "Cit", 0.01
    "Cys", 0.01
    "F6P", 0.01
    "FBP", 0.01
    "Fum", 0.01
    "Glc", 0.01
    "Gln", 0.01
    "Glu", 0.01
    "Gly", 0.01
    "His", 0.01
    "Ile", 0.01
    "Leu", 0.01
    "Lys", 0.01
    "Met", 0.01
    "Phe", 0.01
    "Pro", 0.01
    "Ru5P", 0.01
    "RuBP", 0.01
    "Ser", 0.01
    "Thr", 0.01
    "Trp", 0.01
    "Tyr", 0.01
    "Val", 0.01];

config.mid_name_input = [    "2PGA:C#1,C#2,C#3.0"
    "2PGA:C#1,C#2,C#3.1"
    "2PGA:C#1,C#2,C#3.2"
    "2PGA:C#1,C#2,C#3.3"
    "6PG:C#1,C#2,C#3,C#4,C#5,C#6.1"
    "6PG:C#1,C#2,C#3,C#4,C#5,C#6.2"
    "6PG:C#1,C#2,C#3,C#4,C#5,C#6.3"
    "6PG:C#1,C#2,C#3,C#4,C#5,C#6.4"
    "6PG:C#1,C#2,C#3,C#4,C#5,C#6.5"
    "6PG:C#1,C#2,C#3,C#4,C#5,C#6.6"
    "Ala:C#1,C#2,C#3.0"
    "Ala:C#1,C#2,C#3.1"
    "Ala:C#1,C#2,C#3.2"
    "Ala:C#1,C#2,C#3.3"
    "Arg:C#1,C#2,C#3,C#4,C#5,C#6.0"
    "Arg:C#1,C#2,C#3,C#4,C#5,C#6.1"
    "Asn:C#1,C#2,C#3,C#4.0"
    "Asn:C#1,C#2,C#3,C#4.1"
    "Asn:C#1,C#2,C#3,C#4.2"
    "Asn:C#1,C#2,C#3,C#4.3"
    "Asp:C#1,C#2,C#3,C#4.0"
    "Asp:C#1,C#2,C#3,C#4.1"
    "Asp:C#1,C#2,C#3,C#4.2"
    "Asp:C#1,C#2,C#3,C#4.3"
    "Asp:C#1,C#2,C#3,C#4.4"
    "Cit:C#1,C#2,C#3,C#4,C#5,C#6.0"
    "Cit:C#1,C#2,C#3,C#4,C#5,C#6.1"
    "Cit:C#1,C#2,C#3,C#4,C#5,C#6.2"
    "Cit:C#1,C#2,C#3,C#4,C#5,C#6.3"
    "Cit:C#1,C#2,C#3,C#4,C#5,C#6.4"
    "Cit:C#1,C#2,C#3,C#4,C#5,C#6.5"
    "Cit:C#1,C#2,C#3,C#4,C#5,C#6.6"
    "Cys:C#1,C#2,C#3.0"
    "Cys:C#1,C#2,C#3.1"
    "Cys:C#1,C#2,C#3.2"
    "Cys:C#1,C#2,C#3.3"
    "F6P:C#1,C#2,C#3,C#4,C#5,C#6.0"
    "F6P:C#1,C#2,C#3,C#4,C#5,C#6.1"
    "F6P:C#1,C#2,C#3,C#4,C#5,C#6.2"
    "F6P:C#1,C#2,C#3,C#4,C#5,C#6.3"
    "F6P:C#1,C#2,C#3,C#4,C#5,C#6.4"
    "F6P:C#1,C#2,C#3,C#4,C#5,C#6.5"
    "F6P:C#1,C#2,C#3,C#4,C#5,C#6.6"
    "FBP:C#1,C#2,C#3,C#4,C#5,C#6.0"
    "FBP:C#1,C#2,C#3,C#4,C#5,C#6.1"
    "FBP:C#1,C#2,C#3,C#4,C#5,C#6.2"
    "FBP:C#1,C#2,C#3,C#4,C#5,C#6.3"
    "FBP:C#1,C#2,C#3,C#4,C#5,C#6.4"
    "FBP:C#1,C#2,C#3,C#4,C#5,C#6.5"
    "FBP:C#1,C#2,C#3,C#4,C#5,C#6.6"
    "Fum:C#1,C#2,C#3,C#4.0"
    "Fum:C#1,C#2,C#3,C#4.1"
    "Fum:C#1,C#2,C#3,C#4.2"
    "Fum:C#1,C#2,C#3,C#4.3"
    "Fum:C#1,C#2,C#3,C#4.4"
    "Glc:C#1,C#2,C#3,C#4,C#5,C#6.0"
    "Glc:C#1,C#2,C#3,C#4,C#5,C#6.1"
    "Glc:C#1,C#2,C#3,C#4,C#5,C#6.2"
    "Glc:C#1,C#2,C#3,C#4,C#5,C#6.3"
    "Glc:C#1,C#2,C#3,C#4,C#5,C#6.4"
    "Glc:C#1,C#2,C#3,C#4,C#5,C#6.5"
    "Glc:C#1,C#2,C#3,C#4,C#5,C#6.6"
    "Gln:C#1,C#2,C#3,C#4,C#5.0"
    "Glu:C#1,C#2,C#3,C#4,C#5.0"
    "Glu:C#1,C#2,C#3,C#4,C#5.1"
    "Gly:C#1,C#2.0"
    "Gly:C#1,C#2.1"
    "Gly:C#1,C#2.2"
    "His:C#1,C#2,C#3,C#4,C#5,C#6.0"
    "His:C#1,C#2,C#3,C#4,C#5,C#6.1"
    "His:C#1,C#2,C#3,C#4,C#5,C#6.2"
    "His:C#1,C#2,C#3,C#4,C#5,C#6.3"
    "His:C#1,C#2,C#3,C#4,C#5,C#6.4"
    "His:C#1,C#2,C#3,C#4,C#5,C#6.5"
    "His:C#1,C#2,C#3,C#4,C#5,C#6.6"
    "Ile:C#1,C#2,C#3,C#4,C#5,C#6.0"
    "Ile:C#1,C#2,C#3,C#4,C#5,C#6.1"
    "Ile:C#1,C#2,C#3,C#4,C#5,C#6.2"
    "Ile:C#1,C#2,C#3,C#4,C#5,C#6.3"
    "Ile:C#1,C#2,C#3,C#4,C#5,C#6.4"
    "Ile:C#1,C#2,C#3,C#4,C#5,C#6.5"
    "Ile:C#1,C#2,C#3,C#4,C#5,C#6.6"
    "Leu:C#1,C#2,C#3,C#4,C#5,C#6.0"
    "Leu:C#1,C#2,C#3,C#4,C#5,C#6.1"
    "Leu:C#1,C#2,C#3,C#4,C#5,C#6.2"
    "Leu:C#1,C#2,C#3,C#4,C#5,C#6.3"
    "Leu:C#1,C#2,C#3,C#4,C#5,C#6.4"
    "Lys:C#1,C#2,C#3,C#4,C#5,C#6.0"
    "Lys:C#1,C#2,C#3,C#4,C#5,C#6.1"
    "Lys:C#1,C#2,C#3,C#4,C#5,C#6.2"
    "Lys:C#1,C#2,C#3,C#4,C#5,C#6.3"
    "Lys:C#1,C#2,C#3,C#4,C#5,C#6.4"
    "Lys:C#1,C#2,C#3,C#4,C#5,C#6.5"
    "Lys:C#1,C#2,C#3,C#4,C#5,C#6.6"
    "Met:C#1,C#2,C#3,C#4,C#5.0"
    "Met:C#1,C#2,C#3,C#4,C#5.1"
    "Met:C#1,C#2,C#3,C#4,C#5.2"
    "Met:C#1,C#2,C#3,C#4,C#5.3"
    "Met:C#1,C#2,C#3,C#4,C#5.4"
    "Phe:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.0"
    "Phe:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.1"
    "Phe:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.2"
    "Phe:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.3"
    "Phe:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.4"
    "Phe:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.5"
    "Phe:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.6"
    "Phe:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.7"
    "Pro:C#1,C#2,C#3,C#4,C#5.0"
    "Ru5P:C#1,C#2,C#3,C#4,C#5.0"
    "Ru5P:C#1,C#2,C#3,C#4,C#5.1"
    "Ru5P:C#1,C#2,C#3,C#4,C#5.2"
    "Ru5P:C#1,C#2,C#3,C#4,C#5.3"
    "Ru5P:C#1,C#2,C#3,C#4,C#5.4"
    "Ru5P:C#1,C#2,C#3,C#4,C#5.5"
    "RuBP:C#1,C#2,C#3,C#4,C#5.0"
    "RuBP:C#1,C#2,C#3,C#4,C#5.1"
    "RuBP:C#1,C#2,C#3,C#4,C#5.2"
    "RuBP:C#1,C#2,C#3,C#4,C#5.3"
    "RuBP:C#1,C#2,C#3,C#4,C#5.4"
    "RuBP:C#1,C#2,C#3,C#4,C#5.5"
    "Ser:C#1,C#2,C#3.0"
    "Ser:C#1,C#2,C#3.1"
    "Ser:C#1,C#2,C#3.2"
    "Ser:C#1,C#2,C#3.3"
    "Thr:C#1,C#2,C#3,C#4.0"
    "Thr:C#1,C#2,C#3,C#4.1"
    "Thr:C#1,C#2,C#3,C#4.2"
    "Thr:C#1,C#2,C#3,C#4.3"
    "Thr:C#1,C#2,C#3,C#4.4"
    "Trp:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9,C#10,C#11.0"
    "Trp:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9,C#10,C#11.1"
    "Trp:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9,C#10,C#11.2"
    "Trp:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9,C#10,C#11.3"
    "Trp:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9,C#10,C#11.4"
    "Trp:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9,C#10,C#11.5"
    "Trp:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9,C#10,C#11.6"
    "Trp:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9,C#10,C#11.7"
    "Tyr:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.0"
    "Tyr:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.1"
    "Tyr:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.2"
    "Tyr:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.3"
    "Tyr:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.4"
    "Tyr:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.5"
    "Tyr:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.6"
    "Tyr:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.7"
    "Tyr:C#1,C#2,C#3,C#4,C#5,C#6,C#7,C#8,C#9.8"
    "Val:C#1,C#2,C#3,C#4,C#5.0"
    "Val:C#1,C#2,C#3,C#4,C#5.1"
    "Val:C#1,C#2,C#3,C#4,C#5.2"
    "Val:C#1,C#2,C#3,C#4,C#5.3"
    "Val:C#1,C#2,C#3,C#4,C#5.4"
    "Val:C#1,C#2,C#3,C#4,C#5.5"];

% part 2: path and file names for intermediate storage
config.base_dir = fullfile("runtime_data", "run_AraCoreC");
config.model_file = fullfile(config.base_dir,"preparedAraCore");
config.samples_dir = fullfile(config.base_dir,"samples");
config.scaled_fd_sample_file = fullfile(config.base_dir,"samples_totalfluxlimited_30_50_65_80_90_95_scaled.mat");
config.met_pools_sample_file = fullfile(config.base_dir,"sample_met_pools.mat");
config.split_dir = fullfile(config.base_dir,"split_data");
config.sample_ratios_file = fullfile(config.base_dir,"normalized_met_sample_ratios_mid_AraCoreC.mat");
config.sample_d_ratios_file = fullfile(config.base_dir,"normalized_met_sample_d_ratios_mid_AraCoreC.mat");
config.sample_decomp_ratios_file = fullfile(config.base_dir,"normalized_met_sample_decomp_ratios_mid_AraCoreC.mat");
config.sample_decomp_d_ratios_file = fullfile(config.base_dir,"normalized_met_sample_decomp_d_ratios_mid_AraCoreC.mat");
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
