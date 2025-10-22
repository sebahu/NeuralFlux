%function [error] = createFiguresForTestSamples(config_script_NN, config_script_testsamples, ...
%                                   config_script_NN2, config_script_testsamples2, rng_start, start_is, end_is)

error = 0;
if ~exist('config_script_testsamples', 'var')
    disp("call with appropriate configs eg. createFiguresForTestSamples('configAraCoreC','configEvalTestSamplesC')");
    error = -1;
    return;
end

cfgTest = eval(config_script_testsamples);
cfgTest2 = eval(config_script_testsamples2);
load(cfgTest.sample_decomp_ratios_file,"sample_decomp_ratios");
sample_decomp_ratios_test1 = sample_decomp_ratios;
load(cfgTest2.sample_decomp_ratios_file,"sample_decomp_ratios");
sample_decomp_ratios_test2 = sample_decomp_ratios;
num_test_samples = size(sample_decomp_ratios_test1,3);
load(cfgTest.met_pools_sample_file,"sample_met_pools");
sample_met_pools_test = sample_met_pools;
load(cfgTest.scaled_fd_sample_file,"total_samples_scaled");
total_samples_scaled_test = total_samples_scaled;

rng(rng_start);
rand_samples = randperm(num_test_samples);

num_eval_samples = 100;
num_best_eval_samples = 100;

sample_decomp_ratios_rnd1 = sample_decomp_ratios_test1(:,:,rand_samples(1:num_eval_samples));
sample_decomp_ratios_rnd2 = sample_decomp_ratios_test2(:,:,rand_samples(1:num_eval_samples));
sample_met_pools_rnd = sample_met_pools_test(:,rand_samples(1:num_eval_samples));
total_samples_scaled_rnd = total_samples_scaled_test(:,rand_samples(1:num_eval_samples));

cfgNN = eval(config_script_NN);
load(cfgNN.scaled_fd_sample_file,"total_samples_scaled");
load(cfgNN.met_pools_sample_file,"sample_met_pools");
load(cfgNN.model_file, "measured_emu_mids", "model", "needed_emus", "rnullspace");
measured_decomp_emu_mids1 = unique(extractBefore(measured_emu_mids,"[")+extractAfter(measured_emu_mids,"]"));
measured_decomp_mets1 = unique(extractBefore(measured_emu_mids,"["));
decomp_relevant_measured_mids1_b = ismember(extractAfter(cfgNN.mid_name_input,"."), ...
                                            ["0", "1", "2", "3", "4", "5"]);
decomp_measured_mids1 = cfgNN.mid_name_input(decomp_relevant_measured_mids1_b);
needed_mids_ix1 = zeros(size(decomp_measured_mids1));
for  mid_i=1:length(decomp_measured_mids1)
    needed_mids_ix1(mid_i) = find(measured_decomp_emu_mids1 == decomp_measured_mids1(mid_i));
end
[emu_mid_met_inx1, ~, ~, ~] = create_emu_data(model, needed_emus);

cfgNN2 = eval(config_script_NN2);
load(cfgNN2.model_file, "measured_emu_mids", "needed_emus");
measured_decomp_emu_mids2 = unique(extractBefore(measured_emu_mids,"[")+extractAfter(measured_emu_mids,"]"));
measured_decomp_mets2 = unique(extractBefore(measured_emu_mids,"["));
decomp_relevant_measured_mids2_b = ismember(extractAfter(cfgNN2.mid_name_input,"."), ...
                                            ["0", "1", "2", "3", "4", "5"]);
decomp_measured_mids2 = cfgNN2.mid_name_input(decomp_relevant_measured_mids2_b);
needed_mids_ix2 = zeros(size(decomp_measured_mids2));
for  mid_i=1:length(decomp_measured_mids2)
    needed_mids_ix2(mid_i) = find(measured_decomp_emu_mids2 == decomp_measured_mids2(mid_i));
end
[emu_mid_met_inx2, ~, ~, ~] = create_emu_data(model, needed_emus);

relevant_met_is1 = unique(emu_mid_met_inx1);
relevant_met_is2 = unique(emu_mid_met_inx2);
relevant_met_is = unique([emu_mid_met_inx1; emu_mid_met_inx2]);
    
rnsp_factor_samples = (total_samples_scaled(1:end-1,:)'/rnullspace')';
num_nspfactors = size(rnullspace,2);
num_relevant_mets = length(relevant_met_is);
minX = min([rnsp_factor_samples',sample_met_pools(relevant_met_is,:)']);
rangeX = max([rnsp_factor_samples',sample_met_pools(relevant_met_is,:)']) - minX;
rangeX(rangeX<=0)=1;


NNs_per_MID_and_tp = {};
mu_per_MID_and_tp = {};
sigma_per_MID_and_tp = {};
test_per_MID_and_tp = {};
relevant_decomp_met_names = extractBefore(string(model.mets(relevant_met_is)),"[");
met_decomp_matrix = zeros(length(relevant_met_is), length(measured_decomp_mets1));
met_decomp_min = zeros(length(measured_decomp_mets1),1);
for  comp_met_i=1:length(relevant_met_is)
    measured_decomp_met_i = find(measured_decomp_mets1 == extractBefore(string(model.mets(relevant_met_is(comp_met_i))),"["));
    met_decomp_matrix(comp_met_i, measured_decomp_met_i) = rangeX(num_nspfactors+comp_met_i);
    met_decomp_min(measured_decomp_met_i) = met_decomp_min(measured_decomp_met_i) + minX(num_nspfactors+comp_met_i);
end

relative_constraint_factor = 1.5; % could also be an array with values for
% each measured decompartmentalized mjetabolite, baed on sd of the corresponding measurements

for  mid_i=1:length(decomp_measured_mids1)
    for tp_i = 1:length(cfgNN.selected_timepoints)
        current_tp = cfgNN.selected_timepoints(tp_i);
        load(cfgNN.NN_output_dir + "/NN6_for_"+decomp_measured_mids1(mid_i)+"_tp_+"+string(current_tp)+"_withMets_3L.mat",'compactMdl',...
            'filter_i','rand_indices','mu','sigma','test_result');

        NNs_per_MID_and_tp{mid_i,tp_i} = compactMdl;
        mu_per_MID_and_tp{mid_i,tp_i} = mu;
        sigma_per_MID_and_tp{mid_i,tp_i} = sigma;
        test_per_MID_and_tp{mid_i,tp_i} = test_result;
    end
end

offset2 = length(decomp_measured_mids1);
for  mid_i=1:length(decomp_measured_mids2)
    for tp_i = 1:length(cfgNN.selected_timepoints)
        current_tp = cfgNN.selected_timepoints(tp_i);
        load(cfgNN2.NN_output_dir + "/NN6_for_"+decomp_measured_mids2(mid_i)+"_tp_+"+string(current_tp)+"_withMets_3L.mat",'compactMdl',...
            'filter_i','rand_indices','mu','sigma','test_result');
        NNs_per_MID_and_tp{mid_i+offset2,tp_i} = compactMdl;
        mu_per_MID_and_tp{mid_i+offset2,tp_i} = mu;
        sigma_per_MID_and_tp{mid_i+offset2,tp_i} = sigma;
        test_per_MID_and_tp{mid_i+offset2,tp_i} = test_result;
    end
end
x_selectors_per_mid = true(size(NNs_per_MID_and_tp,1),num_nspfactors + length(relevant_met_is));
x_selectors_per_mid(1:length(decomp_measured_mids1),num_nspfactors+1:end) = ...
    repmat(ismember(relevant_met_is, relevant_met_is1)',length(decomp_measured_mids1), 1);
x_selectors_per_mid(length(decomp_measured_mids1)+1 : end,num_nspfactors+1:end) = ...
    repmat(ismember(relevant_met_is, relevant_met_is2)', length(decomp_measured_mids2),1);


all_flux_diff_log_no_met_data = [];
all_flux_diff_log_met_constraints = [];
all_flux_diff_log_perfect_mets = [];
all_met_diff_log_no_met_data = [];
all_met_diff_log_met_constraints = [];
all_met_diff_log_perfect_mets = [];
all_pred_fluxes_no_met_data = [];
all_pred_mets_no_met_data = [];
all_pred_fluxes_met_constraints = [];
all_pred_mets_met_constraints = [];
all_pred_fluxes_perfect_mets = [];
all_pred_mets_perfect_mets = [];
all_true_fluxes = [];
all_true_mets = [];
all_true_fluxes_b = [];
all_true_mets_b = [];

for batch_i = 1:length(start_is)
    load(resultsDir+"evalResultsConstrainedMetsChi2_"+config_script_NN+"_"+config_script_testsamples+"_" ...
        +config_script_NN2+"_"+config_script_testsamples2+"_"+string(rng_start)+"_"+string(start_is(batch_i))+ ...
        "_"+string(end_is(batch_i))+".mat","flux_diff_log2","met_diff_log2", ...
            "pred_fluxes2", "pred_mets2", "true_fluxes", "true_mets");
    all_flux_diff_log_no_met_data = [all_flux_diff_log_no_met_data, flux_diff_log2];
    all_met_diff_log_no_met_data = [all_met_diff_log_no_met_data, met_diff_log2];
    all_pred_fluxes_no_met_data = [all_pred_fluxes_no_met_data, pred_fluxes2];
    all_pred_mets_no_met_data = [all_pred_mets_no_met_data, pred_mets2];
    true_fluxes = total_samples_scaled_test(1:end-1,rand_samples(start_is(batch_i):end_is(batch_i)));
    all_true_fluxes = [all_true_fluxes, true_fluxes];
    true_fluxes1 = true_fluxes;
    true_mets = sample_met_pools_test(relevant_met_is,rand_samples(start_is(batch_i):end_is(batch_i)));
    all_true_mets = [all_true_mets, true_mets];
    true_mets1 = true_mets;
    load(resultsDir+"evalResultsConstrainedMetsChi2_"+config_script_NN+"_"+config_script_testsamples+"_" ...
        +config_script_NN2+"_"+config_script_testsamples2+"_"+string(rng_start)+"_"+string(start_is(batch_i))+ ...
        "_"+string(end_is(batch_i))+".mat","flux_diff_log2","met_diff_log2", ...
            "pred_fluxes2", "pred_mets2", "true_fluxes", "true_mets");
    all_flux_diff_log_met_constraints = [all_flux_diff_log_met_constraints, flux_diff_log2];
    all_met_diff_log_met_constraints = [all_met_diff_log_met_constraints, met_diff_log2];
    all_pred_fluxes_met_constraints = [all_pred_fluxes_met_constraints, pred_fluxes2];
    all_pred_mets_met_constraints = [all_pred_mets_met_constraints, pred_mets2];
    if any(any(true_fluxes1 ~= true_fluxes)) || any(any(true_mets1 ~= true_mets))
        disp("evalResultsConstrainedMets2_"+config_script_NN+"_"+config_script_testsamples+"_" ...
            +config_script_NN2+"_"+config_script_testsamples2+"_"+string(rng_start)+"_"+string(start_is(batch_i))+ ...
        "_"+string(end_is(batch_i))+".mat differs!")
    end
    load(resultsDir+"evalResultsKnownMetsChi2_"+config_script_NN+"_"+config_script_testsamples+"_" ...
            +config_script_NN2+"_"+config_script_testsamples2+"_"+string(rng_start)+"_"+string(start_is(batch_i))+ ...
        "_"+string(end_is(batch_i))+".mat","flux_diff_log2","met_diff_log2", ...
            "pred_fluxes2", "pred_mets2", "true_fluxes", "true_mets");
    all_flux_diff_log_perfect_mets = [all_flux_diff_log_perfect_mets, flux_diff_log2];
    all_met_diff_log_perfect_mets = [all_met_diff_log_perfect_mets, met_diff_log2];
    all_pred_fluxes_perfect_mets = [all_pred_fluxes_perfect_mets, pred_fluxes2];
    all_pred_mets_perfect_mets = [all_pred_mets_perfect_mets, pred_mets2];
    if any(any(true_fluxes1 ~= true_fluxes)) || any(any(true_mets1 ~= true_mets))
        disp("evalResultsConstrainedMets_"+config_script_NN+"_"+config_script_testsamples+"_" ...
            +config_script_NN2+"_"+config_script_testsamples2+"_"+string(rng_start)+"_"+string(start_is(batch_i))+ ...
        "_"+string(end_is(batch_i))+".mat differs!")
    end
end

all_residuals_no_met_data = zeros(num_eval_samples,1);
all_residuals_met_constraints = zeros(num_eval_samples,1);
all_residuals_perfect_mets = zeros(num_eval_samples,1);
Y_len = (length(needed_mids_ix1) + length(needed_mids_ix2)) * length(cfgNN.selected_timepoints);
all_Y_diffs_no_met_data = zeros(Y_len,num_eval_samples);
all_Y_diffs_met_constraints = zeros(Y_len,num_eval_samples);
all_Y_diffs_perfect_mets = zeros(Y_len,num_eval_samples);
all_Y_values_no_met_data = zeros(Y_len,num_eval_samples);
all_Y_values_met_constraints = zeros(Y_len,num_eval_samples);
all_Y_values_perfect_mets = zeros(Y_len,num_eval_samples);
all_Y_values_true = zeros(Y_len,num_eval_samples);

test_pred_mids = zeros(length(needed_mids_ix1) + length(needed_mids_ix2), length(cfgNN.selected_timepoints), num_eval_samples);
for r_test_i = 1:num_eval_samples
    test_i = rand_samples(r_test_i);
    x_raw = [ (total_samples_scaled_test(1:end-1,test_i)'/rnullspace')'; sample_met_pools_test(relevant_met_is,test_i) ];
    x = (x_raw-minX')./rangeX';
    Y_values = forwardSigmoidnet2v2(NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp, x', x_selectors_per_mid);
    test_pred_mids(:,:,r_test_i) = reshape(Y_values,[],length(cfgNN.selected_timepoints));
end
pred_mid_diff = test_pred_mids - [sample_decomp_ratios_test1(needed_mids_ix1,cfgNN.selected_timepoints,rand_samples(1:num_eval_samples)) ; ...
                                 sample_decomp_ratios_test2(needed_mids_ix2,cfgNN.selected_timepoints,rand_samples(1:num_eval_samples)) ];

for pred_i = 1:size(all_true_fluxes,2)    
    test_Y = [squeeze(sample_decomp_ratios_rnd1(needed_mids_ix1,cfgNN.selected_timepoints, pred_i)) ; ...
                squeeze(sample_decomp_ratios_rnd2(needed_mids_ix2,cfgNN.selected_timepoints, pred_i)) ];
    test_Y_flat = reshape(test_Y',[],1);
    all_Y_values_true(:,pred_i) = test_Y_flat;
    x_raw = [ (all_pred_fluxes_no_met_data(:,pred_i)'/rnullspace')'; all_pred_mets_no_met_data(:,pred_i) ];
    x = (x_raw-minX')./rangeX';
    Y_values = forwardSigmoidnet2v2(NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp, x', x_selectors_per_mid);
    Y_diff = Y_values - test_Y_flat;
    all_Y_values_no_met_data(:,pred_i) = Y_values;
    all_Y_diffs_no_met_data(:,pred_i) = Y_diff;
    all_residuals_no_met_data(pred_i) = Y_diff' * Y_diff;
    x_raw = [ (all_pred_fluxes_met_constraints(:,pred_i)'/rnullspace')'; all_pred_mets_met_constraints(:,pred_i) ];
    x = (x_raw-minX')./rangeX';
    Y_values = forwardSigmoidnet2v2(NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp, x', x_selectors_per_mid);
    Y_diff = Y_values - test_Y_flat;
    all_Y_values_met_constraints(:,pred_i) = Y_values;
    all_Y_diffs_met_constraints(:,pred_i) = Y_diff;
    all_residuals_met_constraints(pred_i) = Y_diff' * Y_diff;
    x_raw = [ (all_pred_fluxes_perfect_mets(:,pred_i)'/rnullspace')'; all_pred_mets_perfect_mets(:,pred_i) ];
    x = (x_raw-minX')./rangeX';
    Y_values = forwardSigmoidnet2v2(NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp, x', x_selectors_per_mid);
    Y_diff = Y_values - test_Y_flat;
    all_Y_values_perfect_mets(:,pred_i) = Y_values;
    all_Y_diffs_perfect_mets(:,pred_i) = Y_diff;
    all_residuals_perfect_mets(pred_i) = Y_diff' * Y_diff;
end

save("all_result_data.mat", "all_pred_fluxes_met_constraints", "all_pred_fluxes_no_met_data", "all_pred_fluxes_perfect_mets", ...
   "all_pred_mets_met_constraints", "all_pred_mets_no_met_data", "all_pred_mets_perfect_mets", ...
"all_true_fluxes", "all_true_mets", "minX", "mu_per_MID_and_tp", "NNs_per_MID_and_tp", "rangeX", "rnullspace", "sigma_per_MID_and_tp", ...
"test_per_MID_and_tp", "sample_decomp_ratios_rnd1", "sample_decomp_ratios_rnd2", "sample_met_pools_rnd", "total_samples_scaled_rnd", ...
"needed_mids_ix1", "needed_mids_ix2", "cfgNN", "cfgNN2", ...
"measured_decomp_mets1", "measured_decomp_mets2", "model", "relevant_met_is", "all_Y_values_met_constraints", "all_Y_values_no_met_data", "all_Y_values_perfect_mets", ...
"test_pred_mids", "x_selectors_per_mid");
