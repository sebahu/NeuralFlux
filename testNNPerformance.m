%function [] = testNNPerformance(mid_i)

addpath("application_core");
addpath("configs");

config_script_testsamples='configEvalTestSamplesC';
config_script_NN='configAraCoreC';

cfgTest = eval(config_script_testsamples);
load(cfgTest.met_pools_sample_file,"sample_met_pools");
sample_met_pools_test = sample_met_pools;
num_test_samples = size(sample_met_pools_test,2);
load(cfgTest.sample_decomp_ratios_file,"sample_decomp_ratios");
sample_decomp_ratios_test = sample_decomp_ratios;
load(cfgTest.scaled_fd_sample_file,"total_samples_scaled");
total_samples_scaled_test = total_samples_scaled;

cfgNN = eval(config_script_NN);
load(cfgNN.scaled_fd_sample_file,"total_samples_scaled");
load(cfgNN.met_pools_sample_file,"sample_met_pools");
load(cfgNN.model_file, "measured_emu_mids", "model", "needed_emus", "rnullspace");
measured_decomp_emu_mids = unique(extractBefore(measured_emu_mids,"[")+extractAfter(measured_emu_mids,"]"));
needed_mids_ix = zeros(size(cfgNN.mid_name_input));
for  mid_i=1:length(cfgNN.mid_name_input)
    needed_mids_ix(mid_i) = find(measured_decomp_emu_mids == cfgNN.mid_name_input(mid_i));
end

[emu_mid_met_inx, ~, ~, ~] = create_emu_data(model, needed_emus);
relevant_met_is = unique(emu_mid_met_inx);
    
rnsp_factor_samples = (total_samples_scaled(1:end-1,:)'/rnullspace')';
minX = min([rnsp_factor_samples',sample_met_pools(relevant_met_is,:)']);
rangeX = max([rnsp_factor_samples',sample_met_pools(relevant_met_is,:)']) - minX;
rangeX(rangeX<=0)=1;

rnsp_factor_samples_test = (total_samples_scaled_test(1:end-1,:)'/rnullspace')';
scaled_X = ([rnsp_factor_samples_test',sample_met_pools_test(relevant_met_is,:)']-minX)./rangeX;
NN_results_direct = zeros(length(cfgNN.mid_name_input), length(cfgNN.selected_timepoints), num_test_samples);

for mid_i = 1:length(cfgNN.mid_name_input)
    for tp_i = 1:length(cfgNN.selected_timepoints)
        current_tp = cfgNN.selected_timepoints(tp_i);
        load(fullfile(cfgNN.NN_output_dir,"NN6_for_"+cfgNN.mid_name_input(mid_i)+"_tp_+"+string(current_tp)+"_withMets_3L.mat"),'compactMdl',...
                'mu','sigma');
        NN_results_direct(mid_i,tp_i,:) = predict(compactMdl, scaled_X);
    end    
end
sim_results = squeeze(sample_decomp_ratios_test(needed_mids_ix,cfgNN.selected_timepoints,:));
resdiff = NN_results_direct - sim_results;
mid_stds = zeros(length(cfgNN.mid_name_input), length(cfgNN.selected_timepoints));
for mid_i = 1:length(cfgNN.mid_name_input)
    mid_stds(mid_i,:) = std(squeeze(resdiff(mid_i,:,:))');
end
save("NN_performance", "NN_results_direct", "mid_stds");


