
function [error] = evalEstimatesForTestSamplesWithConstrainedMetConcs(config_script_NN, config_script_testsamples, rng_start, start_i, end_i)
% Base idea: learn NNs for every (relevant) timepoint of every measured
% metabolite/MID
% first calculate the relevant timepoints (PCA)
% then calculate the relevant input values (fluxes/met concs)
% learn a neural Network
% to estimate fluxes: unite all the NNs and their needed inputs to one
% function, apply lsqnonlin on it

error = 0;

if ~exist('config_script_testsamples', 'var')
    disp("call with appropriate configs eg. evalEstimatesForTestSamples('configAraCoreC','configEvalTestSamplesC')");
    error = -1;
    return;
end

if ~exist('rng_start', 'var')
    rng_start = mod(now*1e10,1e5);
end

if ~exist('start_i', 'var')
    start_i = 1;
end

if ~exist('end_i', 'var')
    end_i = 100;
end


load("evalClosestSamples_"+config_script_NN+"_"+config_script_testsamples+".mat","closest1000_samples");

cfgNN = eval(config_script_NN);
cfgTest = eval(config_script_testsamples);
load(cfgTest.sample_decomp_ratios_file,"sample_decomp_ratios");
sample_decomp_ratios_test = sample_decomp_ratios;
num_test_samples = size(sample_decomp_ratios_test,3);
load(cfgTest.met_pools_sample_file,"sample_met_pools");
sample_met_pools_test = sample_met_pools;
load(cfgTest.scaled_fd_sample_file,"total_samples_scaled");
total_samples_scaled_test = total_samples_scaled;

load(cfgNN.scaled_fd_sample_file,"total_samples_scaled");
load(cfgNN.met_pools_sample_file,"sample_met_pools");
load(cfgNN.model_file, "measured_emu_mids", "model", "needed_emus", "rnullspace");
measured_decomp_emu_mids = unique(extractBefore(measured_emu_mids,"[")+extractAfter(measured_emu_mids,"]"));

decomp_relevant_measured_mids_b = ismember(extractAfter(cfgNN.mid_name_input,"."), ...
                                            ["0", "1", "2", "3", "4", "5"]);
decomp_measured_mids = cfgNN.mid_name_input(decomp_relevant_measured_mids_b);
needed_mids_ix = zeros(size(decomp_measured_mids));
for  mid_i=1:length(decomp_measured_mids)
    needed_mids_ix(mid_i) = find(measured_decomp_emu_mids == decomp_measured_mids(mid_i));
end

[emu_mid_met_inx, ~, ~, ~] = create_emu_data(model, needed_emus);
relevant_met_is = unique(emu_mid_met_inx);
    
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
measured_decomp_mets = unique(extractBefore(measured_emu_mids,"["));
met_decomp_matrix = zeros(length(relevant_met_is), length(measured_decomp_mets));
met_decomp_min = zeros(length(measured_decomp_mets),1);
for  comp_met_i=1:length(relevant_met_is)
    measured_decomp_met_i = find(measured_decomp_mets == extractBefore(string(model.mets(relevant_met_is(comp_met_i))),"["));
    met_decomp_matrix(comp_met_i, measured_decomp_met_i) = rangeX(num_nspfactors+comp_met_i);
    met_decomp_min(measured_decomp_met_i) = met_decomp_min(measured_decomp_met_i) + minX(num_nspfactors+comp_met_i);
end
relative_constraint_factor = 1.5; % could also be an array with values for
% each measured decompartmentalized mjetabolite, baed on sd of the corresponding measurements
test_sds = 0.01*ones(length(needed_mids_ix),length(cfgNN.selected_timepoints));
for  mid_i=1:length(decomp_measured_mids)
    test_sds(mid_i, :) = cfgNN.met_mid_sds(cfgNN.met_mid_sds(:,1)==extractBefore(decomp_measured_mids(mid_i),":"), 2);
end

for  mid_i=1:length(decomp_measured_mids)
    for tp_i = 1:length(cfgNN.selected_timepoints)
        current_tp = cfgNN.selected_timepoints(tp_i);
        load(cfgNN.NN_output_dir + "/NN6_for_"+decomp_measured_mids(mid_i)+"_tp_+"+string(current_tp)+"_withMets_3L.mat",'compactMdl',...
            'filter_i','rand_indices','mu','sigma','test_result');

        NNs_per_MID_and_tp{mid_i,tp_i} = compactMdl;
        mu_per_MID_and_tp{mid_i,tp_i} = mu;
        sigma_per_MID_and_tp{mid_i,tp_i} = sigma;
        test_per_MID_and_tp{mid_i,tp_i} = test_result;
    end
end

rng(rng_start);
rand_samples = randperm(num_test_samples);

max_iterations = 10;
fitOptions = optimoptions(@lsqnonlin,'Display','iter', 'FunctionTolerance',1e-10, 'StepTolerance',1e-10, ...
                    'MaxIterations',max_iterations,'MaxFunctionEvaluations',1600000, 'UseParallel', true);

num_rand_samples = end_i - start_i + 1;
flux_diff_log2 = zeros(length(model.rxns),num_rand_samples);
pred_nsfs = zeros(num_nspfactors, num_rand_samples);
pred_fluxes2 = flux_diff_log2;
true_fluxes = total_samples_scaled_test(1:end-1,rand_samples(start_i:end_i));
met_diff_log2 = zeros(length(relevant_met_is),num_rand_samples);
pred_mets2 = met_diff_log2;
true_mets = sample_met_pools_test(relevant_met_is,rand_samples(start_i:end_i));
for rand_i = start_i:end_i
    local_i = rand_i - start_i + 1;
    test_i = rand_samples(rand_i);
    test_Y = squeeze(sample_decomp_ratios_test(needed_mids_ix,cfgNN.selected_timepoints, test_i));
    best_fit_samples = squeeze(closest1000_samples(test_i,:,1));
   
    measured_decomp_met_conc = zeros(length(measured_decomp_mets),1);
    for  comp_met_i=1:length(relevant_met_is)
        measured_decomp_met_i = find(measured_decomp_mets == extractBefore(string(model.mets(relevant_met_is(comp_met_i))),"["));
        measured_decomp_met_conc(measured_decomp_met_i) = measured_decomp_met_conc(measured_decomp_met_i) + ...
                                sample_met_pools_test(relevant_met_is(comp_met_i), test_i);
    end

    nsf_base = (rnsp_factor_samples(:,best_fit_samples(1:5))-minX(1:num_nspfactors)')./rangeX(1:num_nspfactors)';
    met_base = (sample_met_pools(relevant_met_is,best_fit_samples(1:5))-minX(num_nspfactors+1:end)')./rangeX(num_nspfactors+1:end)';
    
    [standardized_nsfs, standardized_met_concs, resnorm, residual] = fit_NN_to_NSFs_and_constrained_mets_Chi2( ...
        fitOptions, test_Y, test_sds, nsf_base, met_base, NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp, ...
        met_decomp_matrix, met_decomp_min, measured_decomp_met_conc, relative_constraint_factor);
    save_file_name = "evalResultsConstrainedMetsChi2_"+config_script_NN+"_"+config_script_testsamples+"_"+string(rng_start)+"_"+string(start_i)+ ...
        "_"+string(end_i)+".mat";
    predicted_nsfs = minX(1:num_nspfactors)' + standardized_nsfs'.*rangeX(1:num_nspfactors)';
    predicted_mets = minX(num_nspfactors+1:end)' + standardized_met_concs'.*rangeX(num_nspfactors+1:end)';
    predicted_fluxes = rnullspace*predicted_nsfs;
    pred_nsfs(:,local_i) = predicted_nsfs;
    pred_fluxes2(:,local_i) = predicted_fluxes;
    pred_mets2(:,local_i) = predicted_mets;
    flux_diff_log2(:,local_i) = log10(max(predicted_fluxes,1e-6))-log10(max(total_samples_scaled_test(1:end-1,test_i),1e-6));
    met_diff_log2(:,local_i) = log10(max(predicted_mets,1e-6))-log10(max(sample_met_pools_test(relevant_met_is,test_i),1e-6));
    disp(resnorm);
    save(save_file_name, "flux_diff_log2","met_diff_log2", "pred_fluxes2", "pred_nsfs", "pred_mets2", "true_fluxes", "true_mets");
end

%flux_diff_log2_good = abs(flux_diff_log2)<0.5;
%flux_diff_log_good = abs(flux_diff_log)<0.5;
%flux_diff_log2_stats = sum(flux_diff_log2_good,2);
%flux_diff_log_stats = sum(flux_diff_log_good,2);
%figure(); histogram(flux_diff_log2_stats,100);
%figure(); histogram(flux_diff_log_stats,100);
