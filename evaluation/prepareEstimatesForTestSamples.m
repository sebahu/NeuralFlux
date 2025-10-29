function [error] = prepareEstimatesForTestSamples(config_script_NN, config_script_testsamples)
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


load(cfgNN.model_file, "measured_emu_mids");
measured_decomp_emu_mids = unique(extractBefore(measured_emu_mids,"[")+extractAfter(measured_emu_mids,"]"));

load(cfgTest.sample_decomp_ratios_file,"sample_decomp_ratios");
sample_decomp_ratios_test = sample_decomp_ratios;
load(cfgTest.met_pools_sample_file,"sample_met_pools");
sample_met_pools_test = sample_met_pools;
load(cfgNN.sample_decomp_ratios_file,"sample_decomp_ratios");
load(cfgNN.met_pools_sample_file,"sample_met_pools");

decomp_relevant_measured_mids_b = ismember(extractAfter(cfgNN.mid_name_input,"."), ...
                                            ["0", "1", "2", "3", "4", "5"]);
decomp_measured_mids = cfgNN.mid_name_input(decomp_relevant_measured_mids_b);
needed_mids_ix = zeros(size(decomp_measured_mids));
for  mid_i=1:length(decomp_measured_mids)
    needed_mids_ix(mid_i) = find(measured_decomp_emu_mids == decomp_measured_mids(mid_i));
end

num_test_samples = size(sample_decomp_ratios_test,3);
closest1000_samples = zeros(num_test_samples,1000,2);
for test_i = 1:num_test_samples

    test_Y = squeeze(sample_decomp_ratios_test(needed_mids_ix,cfgNN.selected_timepoints,test_i));

    sample_diff = sample_decomp_ratios(needed_mids_ix, cfgNN.selected_timepoints,:) -test_Y;
    sample_diff_squaresum_per_mid = squeeze(sum(sample_diff.^2,2));
    sample_diff_squaresum = sum(sample_diff_squaresum_per_mid,1)';
    %sample_met_diff_squaresum = sum((100*(rel_decomp_met_samples-rel_measured_decomp_met_conc)).^2,1)';
   
    [sorted_diff, best_fit_samples] = sort(sample_diff_squaresum);
    %[~, best_fit_samples] = sort(sample_diff_squaresum+sample_met_diff_squaresum);
    closest1000_samples(test_i,:,:) = [best_fit_samples(1:1000), sorted_diff(1:1000)];
end

save("evalClosestSamples_"+config_script_NN+"_"+config_script_testsamples+".mat","closest1000_samples");

