function [error] = evalConfidenceIntervalsWithConstrainedMetConcs2(config_script_NN, config_script_NN2, rand_i, model_met_i)
% Base idea: learn NNs for every (relevant) timepoint of every measured
% metabolite/MID
% first calculate the relevant timepoints (PCA)
% then calculate the relevant input values (fluxes/met concs)
% learn a neural Networke
% to estimate fluxes: unite all the NNs and their needed inputs to one
% function, apply lsqnonlin on it
    %gpu1 = gpuDevice(1)

error = 0;

cfgNN = eval(config_script_NN);
cfgNN2 = eval(config_script_NN2);
load(cfgNN.model_file, "measured_emu_mids", "model", "needed_emus", "rnullspace");
measured_decomp_emu_mids = unique(extractBefore(measured_emu_mids,"[")+extractAfter(measured_emu_mids,"]"));
measured_decomp_mets1 = unique(extractBefore(measured_emu_mids,"["));
load(cfgNN2.model_file, "measured_emu_mids");
measured_decomp_emu_mids2 = unique(extractBefore(measured_emu_mids,"[")+extractAfter(measured_emu_mids,"]"));
measured_decomp_mets2 = unique(extractBefore(measured_emu_mids,"["));
measured_decomp_mets = unique([ measured_decomp_mets1; measured_decomp_mets2]);
load('all_result_dataCN.mat')

decomp_relevant_measured_mids_b = ismember(extractAfter(cfgNN.mid_name_input,"."), ...
                                            ["0", "1", "2", "3", "4", "5"]);
decomp_measured_mids = cfgNN.mid_name_input(decomp_relevant_measured_mids_b);
needed_mids_ix = zeros(size(decomp_measured_mids));
for  mid_i=1:length(decomp_measured_mids)
    needed_mids_ix(mid_i) = find(measured_decomp_emu_mids == decomp_measured_mids(mid_i));
end
decomp_relevant_measured_mids_b2 = ismember(extractAfter(cfgNN2.mid_name_input,"."), ...
                                            ["0", "1", "2", "3", "4", "5"]);
decomp_measured_mids2 = cfgNN2.mid_name_input(decomp_relevant_measured_mids_b2);
needed_mids_ix2 = zeros(size(decomp_measured_mids2));
for  mid_i=1:length(decomp_measured_mids2)
    needed_mids_ix2(mid_i) = find(measured_decomp_emu_mids2 == decomp_measured_mids2(mid_i));
end

test_sds = 0.01*ones(length(needed_mids_ix) + length(needed_mids_ix2),length(cfgNN.selected_timepoints));
for  mid_i=1:length(decomp_measured_mids)
    test_sds(mid_i, :) = cfgNN.met_mid_sds(cfgNN.met_mid_sds(:,1)==extractBefore(decomp_measured_mids(mid_i),":"), 2);
end
for  mid_i=1:length(decomp_measured_mids2)
    test_sds(length(decomp_measured_mids) + mid_i, :) = cfgNN2.met_mid_sds(cfgNN2.met_mid_sds(:,1)==extractBefore(decomp_measured_mids2(mid_i),":"), 2);
end

num_nspfactors = size(rnullspace,2);



relevant_decomp_met_names = extractBefore(string(model.mets(relevant_met_is)),"[");
met_decomp_matrix = zeros(length(relevant_met_is), length(measured_decomp_mets));
met_decomp_min = zeros(length(measured_decomp_mets),1);
for  comp_met_i=1:length(relevant_met_is)
    measured_decomp_met_i = find(measured_decomp_mets == extractBefore(string(model.mets(relevant_met_is(comp_met_i))),"["));
    met_decomp_matrix(comp_met_i, measured_decomp_met_i) = rangeX(num_nspfactors+comp_met_i);
    met_decomp_min(measured_decomp_met_i) = met_decomp_min(measured_decomp_met_i) + minX(num_nspfactors+comp_met_i);
end
relative_constraint_factor = 1.5; % could also be an array with values for
% each measured decompartmentalized mjetabolite, baed on sd of the corresponding measurements


max_iterations = 10;
fitOptions = optimoptions(@lsqnonlin,'Display','iter', 'FunctionTolerance',1e-10, 'StepTolerance',1e-10, ...
                    'MaxIterations',max_iterations,'MaxFunctionEvaluations',1600000, 'UseParallel', true);


   
measured_decomp_met_conc = zeros(length(measured_decomp_mets),1);
for  comp_met_i=1:length(relevant_met_is)
    measured_decomp_met_i = find(measured_decomp_mets == extractBefore(string(model.mets(relevant_met_is(comp_met_i))),"["));
    measured_decomp_met_conc(measured_decomp_met_i) = measured_decomp_met_conc(measured_decomp_met_i) + ...
                            sample_met_pools_rnd(relevant_met_is(comp_met_i), rand_i);
end


flux_base = all_pred_fluxes_met_constraints(:,rand_i);
nsfs_base_unnorm = (flux_base' / rnullspace')';
nsf_base = (nsfs_base_unnorm - minX(1:num_nspfactors)') ./ rangeX(1:num_nspfactors)';

met_base = (all_pred_mets_met_constraints(:,rand_i) - minX(num_nspfactors+1:end)') ./ rangeX(num_nspfactors+1:end)';

flux_corr_perfect_mets = diag(corr(all_pred_fluxes_perfect_mets',all_true_fluxes'));
[~, rxns_sorted_by_accuracy] = sort(flux_corr_perfect_mets);
nan_start = find(isnan(flux_corr_perfect_mets(rxns_sorted_by_accuracy)),1);

result_file_name = "results_met_confidence_intervals_metConstraints_"+string(model_met_i)+"_"+string(rand_i)+".mat";
if isfile(result_file_name)
    return
end

met_i = find(relevant_met_is == model_met_i);
Chi2_lim = 6.63;

[ci_lb, ci_ub, x_lb, x_ub] = calc_comp_met_confidence_interval_constrained_mets( ...
    fitOptions, all_Y_values_met_constraints(:,rand_i), test_sds, nsf_base, met_base, NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp, ...
    met_decomp_matrix, met_decomp_min, measured_decomp_met_conc, relative_constraint_factor, ...
    met_i, Chi2_lim, x_selectors_per_mid);

save(result_file_name, "ci_lb", "ci_ub", "x_ub", "x_lb");

    
