function [config] = configAraCoreN()
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
config.selected_timepoints = [3,5,8,12,17,23,37,52,66,81];
%config.selected_timepoints = [4,8,12,20,28,36,44,52,60,70,81];
% config.selected_timepoints = [2,3,5,6,16,17,23,24,32,37,40,48,66];
% config.selected_timepoints = [2,3,4,5,6,8,12,16,17,20,23,24,28,32,36,37,40,44,48,52,60,66,70,81];
config.atomLabelFraction = 0; % or 0.011 natural abundance or 0.008 - C3 plants have a δ13C of −33 to −24‰.

config.mid_name_input = ["Ala:N#1.0"; 
    "Arg:N#1,N#2,N#3,N#4.0";
    "Asn:N#1,N#2.0";
    "Asp:N#1.0";
    "CTP:N#1,N#2,N#3.0";
    "Cys:N#1.0";
    "GABA:N#1.0";
    "GTP:N#1,N#2,N#3,N#4,N#5.0"; 
    "Gln:N#1,N#2.0";
    "Glu:N#1.0";
    "Gly:N#1.0";
    "His:N#1,N#2,N#3.0";
    "Ile:N#1.0";
    "Leu:N#1.0";
    "Lys:N#1,N#2.0";
    "Met:N#1.0";
    "Phe:N#1.0";
    "Pro:N#1.0";
    "Ser:N#1.0";
    "Thr:N#1.0";
    "Trp:N#1,N#2.0";
    "Tyr:N#1.0";
    "UTP:N#1,N#2.0";
    "Val:N#1.0";
    "urea:N#1,N#2.0";
    "Ala:N#1.1"; 
    "Arg:N#1,N#2,N#3,N#4.1";"Arg:N#1,N#2,N#3,N#4.2";"Arg:N#1,N#2,N#3,N#4.3";"Arg:N#1,N#2,N#3,N#4.4";
    "Asn:N#1,N#2.1";"Asn:N#1,N#2.2";
    "Asp:N#1.1";
    "CTP:N#1,N#2,N#3.1";"CTP:N#1,N#2,N#3.2";"CTP:N#1,N#2,N#3.3";
    "Cys:N#1.1";
    "GABA:N#1.1";
    "GTP:N#1,N#2,N#3,N#4,N#5.1"; "GTP:N#1,N#2,N#3,N#4,N#5.2"; "GTP:N#1,N#2,N#3,N#4,N#5.3"; "GTP:N#1,N#2,N#3,N#4,N#5.4"; "GTP:N#1,N#2,N#3,N#4,N#5.5"; 
    "Gln:N#1,N#2.1";"Gln:N#1,N#2.2";
    "Glu:N#1.1";
    "Gly:N#1.1";
    "His:N#1,N#2,N#3.1";"His:N#1,N#2,N#3.2";"His:N#1,N#2,N#3.3";
    "Ile:N#1.1";
    "Leu:N#1.1";
    "Lys:N#1,N#2.1";"Lys:N#1,N#2.2";
    "Met:N#1.1";
    "Phe:N#1.1";
    "Pro:N#1.1";
    "Ser:N#1.1";
    "Thr:N#1.1";
    "Trp:N#1,N#2.1";"Trp:N#1,N#2.2";
    "Tyr:N#1.1";
    "UTP:N#1,N#2.1"; "UTP:N#1,N#2.2";
    "Val:N#1.1";
    "urea:N#1,N#2.1";"urea:N#1,N#2.2"
];

config.met_mid_sds = [
    "Ala", 0.01
    "Arg", 0.01
    "Asn", 0.01
    "Asp", 0.01
    "CTP", 0.01
    "Cys", 0.01
    "GABA", 0.01
    "GTP", 0.01
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
    "Ser", 0.01
    "Thr", 0.01
    "Trp", 0.01
    "Tyr", 0.01
    "UTP", 0.01
    "Val", 0.01
    "urea", 0.01];


% part 2: path and file names for intermediate storage
config.base_dir = fullfile("runtime_data", "run_AraCoreN");
config.model_file = fullfile(config.base_dir,"preparedAraCore");
config.samples_dir = fullfile(config.base_dir,"samples");
config.scaled_fd_sample_file = fullfile(config.base_dir,"samples_totalfluxlimited_30_50_65_80_90_95_scaled.mat");
config.met_pools_sample_file = fullfile(config.base_dir,"sample_met_pools.mat");
config.split_dir = fullfile(config.base_dir,"split_data");
config.sample_ratios_file = fullfile(config.base_dir,"normalized_met_sample_ratios_mid_AraCoreN.mat");
config.sample_d_ratios_file = fullfile(config.base_dir,"normalized_met_sample_d_ratios_mid_AraCoreN.mat");
config.sample_decomp_ratios_file = fullfile(config.base_dir,"normalized_met_sample_decomp_ratios_mid_AraCoreN.mat");
config.sample_decomp_d_ratios_file = fullfile(config.base_dir,"normalized_met_sample_decomp_d_ratios_mid_AraCoreN.mat");
config.NN_output_dir = fullfile(config.base_dir,"NN_araCore");
config.model_preparation_script = "prepareAraCoreN";
config.all_result_data_file = "all_result_dataN.mat";
config.ci_results_dir = "results_ciN";


% make sure, all folders exist
mkdir(config.base_dir);
mkdir(config.samples_dir);
mkdir(config.split_dir);
mkdir(config.NN_output_dir);
end
