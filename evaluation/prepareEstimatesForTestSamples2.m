function [error] = prepareEstimatesForTestSamples2(config_script_NN, config_script_testsamples, config_script_NN2, config_script_testsamples2)
% This part: find closest samples of the original samples, used to learn
% the NNs from a test set of different samples according to their
% "measurements", i.e. decompartmentalized mids of certain metabolites at
% certain time points (defined in the configs)
% the flux/metabolite concentraions of those closest samples are then used
% as starting values for the optimization problem in the next step.
% separating the steps is due to memory efficiency, as the total samples
% mid data is very big
error = 0;
if ~exist('config_script_testsamples', 'var')
    disp("call with appropriate configs eg. prepareEstimatesForTestSamples('configAraCoreC','configEvalTestSamplesC')");
    error = -1;
    return;
end

cfgNN = eval(config_script_NN);
cfgTest = eval(config_script_testsamples);
cfgNN2 = eval(config_script_NN2);
cfgTest2 = eval(config_script_testsamples2);

load(cfgNN.model_file, "measured_emu_mids");
measured_emu_mids1 = measured_emu_mids;
measured_decomp_emu_mids1 = unique(extractBefore(measured_emu_mids1,"[")+extractAfter(measured_emu_mids1,"]"));
load(cfgNN2.model_file, "measured_emu_mids");
measured_emu_mids2 = measured_emu_mids;
measured_decomp_emu_mids2 = unique(extractBefore(measured_emu_mids2,"[")+extractAfter(measured_emu_mids2,"]"));

load(cfgTest.sample_decomp_ratios_file,"sample_decomp_ratios");
sample_decomp_ratios_test1 = sample_decomp_ratios;
load(cfgTest.met_pools_sample_file,"sample_met_pools");
sample_met_pools_test1 = sample_met_pools;
load(cfgNN.sample_decomp_ratios_file,"sample_decomp_ratios");
sample_decomp_ratios1 = sample_decomp_ratios;
load(cfgNN.met_pools_sample_file,"sample_met_pools");
sample_met_pools1 = sample_met_pools;

load(cfgTest2.sample_decomp_ratios_file,"sample_decomp_ratios");
sample_decomp_ratios_test2 = sample_decomp_ratios;
load(cfgTest2.met_pools_sample_file,"sample_met_pools");
sample_met_pools_test2 = sample_met_pools;
load(cfgNN2.sample_decomp_ratios_file,"sample_decomp_ratios");
sample_decomp_ratios2 = sample_decomp_ratios;
load(cfgNN2.met_pools_sample_file,"sample_met_pools");
sample_met_pools2 = sample_met_pools;

decomp_relevant_measured_mids_b1 = ismember(extractAfter(cfgNN.mid_name_input,"."), ...
                                            ["0", "1", "2", "3", "4", "5"]);
decomp_measured_mids1 = cfgNN.mid_name_input(decomp_relevant_measured_mids_b1);
needed_mids_ix1 = zeros(size(decomp_measured_mids1));
for  mid_i=1:length(decomp_measured_mids1)
    needed_mids_ix1(mid_i) = find(measured_decomp_emu_mids1 == decomp_measured_mids1(mid_i));
end

decomp_relevant_measured_mids_b2 = ismember(extractAfter(cfgNN2.mid_name_input,"."), ...
                                            ["0", "1", "2", "3", "4", "5"]);
decomp_measured_mids2 = cfgNN2.mid_name_input(decomp_relevant_measured_mids_b2);
needed_mids_ix2 = zeros(size(decomp_measured_mids2));
for  mid_i=1:length(decomp_measured_mids2)
    needed_mids_ix2(mid_i) = find(measured_decomp_emu_mids2 == decomp_measured_mids2(mid_i));
end

num_test_samples = size(sample_decomp_ratios_test1,3);
closest1000_samples = zeros(num_test_samples,1000,2);
for test_i = 1:num_test_samples

    test_Y = [squeeze(sample_decomp_ratios_test1(needed_mids_ix1,cfgNN.selected_timepoints,test_i));
                squeeze(sample_decomp_ratios_test2(needed_mids_ix2,cfgNN.selected_timepoints,test_i))];

    sample_diff = [sample_decomp_ratios1(needed_mids_ix1, cfgNN.selected_timepoints,:);
                sample_decomp_ratios2(needed_mids_ix2, cfgNN.selected_timepoints,:)]-test_Y;
    sample_diff_squaresum_per_mid = squeeze(sum(sample_diff.^2,2));
    sample_diff_squaresum = sum(sample_diff_squaresum_per_mid,1)';
    %sample_met_diff_squaresum = sum((100*(rel_decomp_met_samples-rel_measured_decomp_met_conc)).^2,1)';
   
    [sorted_diff, best_fit_samples] = sort(sample_diff_squaresum);
    %[~, best_fit_samples] = sort(sample_diff_squaresum+sample_met_diff_squaresum);
    closest1000_samples(test_i,:,:) = [best_fit_samples(1:1000), sorted_diff(1:1000)];
end

save("evalClosestSamples_"+config_script_NN+"_"+config_script_testsamples + "_" + ...
    config_script_NN2+"_"+config_script_testsamples2+".mat","closest1000_samples");

