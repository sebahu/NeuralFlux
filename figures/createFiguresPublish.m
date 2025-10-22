
if 1

load("cor_groups_n","cor_groups_n","cor_groups_subsystems");

load("all_result_dataCN.mat", "all_pred_fluxes_perfect_mets", "all_pred_fluxes_met_constraints", ...
   "all_true_fluxes", "all_true_mets", "minX", "rangeX", "rnullspace", ...
   "sample_decomp_ratios_rnd1", "sample_decomp_ratios_rnd2", "sample_met_pools_rnd", "total_samples_scaled_rnd", ...
"needed_mids_ix1", "needed_mids_ix2", "cfgNN", "cfgNN2", ...
"measured_decomp_mets1", "measured_decomp_mets2", "relevant_met_is", "all_Y_values_perfect_mets", ...
"test_pred_mids", "x_selectors_per_mid");

load("run_AraCoreC/preparedAraCore.mat", "model")
load("min_max_fluxes.mat");

high_corr_threshold = 0.85;
medium_corr_threshold = 0.3;

all_pred_fluxes_perfect_mets_CN = all_pred_fluxes_perfect_mets;
all_pred_fluxes_met_constraints_CN = all_pred_fluxes_met_constraints;
all_true_mets_CN = all_true_mets;
minX_CN = minX;
rangeX_CN = rangeX;
sample_met_pools_rnd_CN = sample_met_pools_rnd;
total_samples_scaled_rnd_CN = total_samples_scaled_rnd;
relevant_met_is_CN = relevant_met_is;
all_Y_values_perfect_mets_CN = all_Y_values_perfect_mets;
test_pred_mids_CN = test_pred_mids;
flux_corr_control = diag(corr(all_true_fluxes(:,[2:end,1])',all_true_fluxes'));
flux_corr_perfect_mets_CN = diag(corr(all_pred_fluxes_perfect_mets_CN',all_true_fluxes'));
flux_corr_met_constraints_CN = diag(corr(all_pred_fluxes_met_constraints_CN',all_true_fluxes'));
[~, rxns_sorted_by_accuracy_CN] = sort(flux_corr_perfect_mets_CN);

load("all_result_dataC.mat", "all_pred_fluxes_perfect_mets", "all_pred_fluxes_met_constraints", ...
   "all_true_mets", "minX", "rangeX", "rnullspace", ...
   "sample_decomp_ratios_rnd", "sample_met_pools_rnd", "total_samples_scaled_rnd", ...
"needed_mids_ix", "measured_decomp_mets", "relevant_met_is", "all_Y_values_perfect_mets", ...
"test_pred_mids");   
all_pred_fluxes_perfect_mets_C = all_pred_fluxes_perfect_mets;
all_pred_fluxes_met_constraints_C = all_pred_fluxes_met_constraints;
all_true_mets_C = all_true_mets;
minX_C = minX;
rangeX_C = rangeX;
sample_met_pools_rnd_C = sample_met_pools_rnd;
total_samples_scaled_rnd_C = total_samples_scaled_rnd;
relevant_met_is_C = relevant_met_is;
all_Y_values_perfect_mets_C = all_Y_values_perfect_mets;
test_pred_mids_C = test_pred_mids;
flux_corr_perfect_mets_C = diag(corr(all_pred_fluxes_perfect_mets_C',all_true_fluxes'));
[~, rxns_sorted_by_accuracy_C] = sort(flux_corr_perfect_mets_C);

load("all_result_dataN.mat", "all_pred_fluxes_perfect_mets", "all_pred_fluxes_met_constraints", ...
   "all_true_mets", "minX", "rangeX", "rnullspace", ...
   "sample_decomp_ratios_rnd", "sample_met_pools_rnd", "total_samples_scaled_rnd", ...
"needed_mids_ix", "measured_decomp_mets", "relevant_met_is", "all_Y_values_perfect_mets", ...
"test_pred_mids");   
all_pred_fluxes_perfect_mets_N = all_pred_fluxes_perfect_mets;
all_pred_fluxes_met_constraints_N = all_pred_fluxes_met_constraints;
all_true_mets_N = all_true_mets;
minX_N = minX;
rangeX_N = rangeX;
sample_met_pools_rnd_N = sample_met_pools_rnd;
total_samples_scaled_rnd_N = total_samples_scaled_rnd;
relevant_met_is_N = relevant_met_is;
all_Y_values_perfect_mets_N = all_Y_values_perfect_mets;
test_pred_mids_N = test_pred_mids;
flux_corr_perfect_mets_N = diag(corr(all_pred_fluxes_perfect_mets_N',all_true_fluxes'));
[~, rxns_sorted_by_accuracy_N] = sort(flux_corr_perfect_mets_N);

sum(sum(sample_met_pools_rnd_C ~= sample_met_pools_rnd_CN))
sum(sum(total_samples_scaled_rnd_C ~= total_samples_scaled_rnd_CN))

flux_corr_perfect_mets_diff = flux_corr_perfect_mets_CN - flux_corr_perfect_mets_C;

high_corr_subsystems = (unique(string(model.subSystems(flux_corr_perfect_mets_CN>high_corr_threshold))));

high_corr_corr_groups = unique(cor_groups_n(flux_corr_perfect_mets_CN>high_corr_threshold));

subsystems = unique(string(model.subSystems));

subsystem_blocked_n = zeros(length(subsystems),1);
subsystem_high_corr_CN_n = zeros(length(subsystems),1);
subsystem_medium_corr_CN_n = zeros(length(subsystems),1);
subsystem_high_corr_C_n = zeros(length(subsystems),1);
subsystem_medium_corr_C_n = zeros(length(subsystems),1);
subsystem_high_corr_N_n = zeros(length(subsystems),1);
subsystem_medium_corr_N_n = zeros(length(subsystems),1);
subsystem_total_n = zeros(length(subsystems),1);
for i = 1:length(subsystems)
    subsystem_blocked_n(i) = sum(isnan(flux_corr_perfect_mets_CN) & string(model.subSystems) == subsystems(i));
    subsystem_high_corr_CN_n(i) = sum(flux_corr_perfect_mets_CN>high_corr_threshold & string(model.subSystems) == subsystems(i));
    subsystem_medium_corr_CN_n(i) = sum(flux_corr_perfect_mets_CN>medium_corr_threshold & string(model.subSystems) == subsystems(i)) - subsystem_high_corr_CN_n(i);
    subsystem_high_corr_C_n(i) = sum(flux_corr_perfect_mets_C>high_corr_threshold & string(model.subSystems) == subsystems(i));
    subsystem_medium_corr_C_n(i) = sum(flux_corr_perfect_mets_C>medium_corr_threshold & string(model.subSystems) == subsystems(i)) - subsystem_high_corr_C_n(i);
    subsystem_high_corr_N_n(i) = sum(flux_corr_perfect_mets_N>high_corr_threshold & string(model.subSystems) == subsystems(i));
    subsystem_medium_corr_N_n(i) = sum(flux_corr_perfect_mets_N>medium_corr_threshold & string(model.subSystems) == subsystems(i)) - subsystem_high_corr_N_n(i);
    subsystem_total_n(i) = sum(string(model.subSystems) == subsystems(i));
end
subsystem_low_corr_CN_n = subsystem_total_n - subsystem_high_corr_CN_n - subsystem_medium_corr_CN_n - subsystem_blocked_n;
subsystem_low_corr_C_n = subsystem_total_n - subsystem_high_corr_C_n - subsystem_medium_corr_C_n - subsystem_blocked_n;
subsystem_low_corr_N_n = subsystem_total_n - subsystem_high_corr_N_n - subsystem_medium_corr_N_n - subsystem_blocked_n;

subsystem_high_ratio_CN = subsystem_high_corr_CN_n./subsystem_total_n;

load("subsystems","subsystems");
supersystems = unique(subsystems(:,2));
supersystems_total_n = zeros(size(supersystems));
supersystems_high_corr_CN_n = zeros(size(supersystems));
supersystems_medium_corr_CN_n = zeros(size(supersystems));
supersystems_high_corr_C_n = zeros(size(supersystems));
supersystems_medium_corr_C_n = zeros(size(supersystems));
supersystems_high_corr_N_n = zeros(size(supersystems));
supersystems_medium_corr_N_n = zeros(size(supersystems));
supersystems_blocked_n = zeros(size(supersystems));
for i=1:length(supersystems)
supersystems_total_n(i)=sum(subsystem_total_n(subsystems(:,2)==supersystems(i)));
supersystems_blocked_n(i)=sum(subsystem_blocked_n(subsystems(:,2)==supersystems(i)));
supersystems_high_corr_CN_n(i)=sum(subsystem_high_corr_CN_n(subsystems(:,2)==supersystems(i)));
supersystems_medium_corr_CN_n(i)=sum(subsystem_medium_corr_CN_n(subsystems(:,2)==supersystems(i)));
supersystems_high_corr_C_n(i)=sum(subsystem_high_corr_C_n(subsystems(:,2)==supersystems(i)));
supersystems_medium_corr_C_n(i)=sum(subsystem_medium_corr_C_n(subsystems(:,2)==supersystems(i)));
supersystems_high_corr_N_n(i)=sum(subsystem_high_corr_N_n(subsystems(:,2)==supersystems(i)));
supersystems_medium_corr_N_n(i)=sum(subsystem_medium_corr_N_n(subsystems(:,2)==supersystems(i)));
end
supersystems_low_corr_CN_n = supersystems_total_n - supersystems_high_corr_CN_n - supersystems_medium_corr_CN_n - supersystems_blocked_n;
supersystems_low_corr_C_n = supersystems_total_n - supersystems_high_corr_C_n - supersystems_medium_corr_C_n - supersystems_blocked_n;
supersystems_low_corr_N_n = supersystems_total_n - supersystems_high_corr_N_n - supersystems_medium_corr_N_n - supersystems_blocked_n;
supersystems_text = supersystems;
supersystems_text(supersystems == "") = "none assigned";

load('NN_performanceCN.mat', "mid_stdsCN");
figure(); set(gca,'fontsize',24);
set(gcf,"Color","w");
hold on
histogram(mid_stdsCN,0:0.001:0.035);
rectangle('Position', [ median(mid_stdsCN(:))-0.0001, 0, 0.0002, 250], ...
    'LineStyle', 'none', 'FaceColor', [1 0 0]);
ylabel("# of mass isotopomer/time point pairs");
xlabel("standard deviation predicted values - simulated value");
mean(mid_stdsCN(:))
median(mid_stdsCN(:))

hold off

data_NN = readtable("comp_effort_NNs.txt");
data_sim = readtable("comp_effort_sim.txt");

figure(); set(gcf,"Color","w");
ratios = data_sim.Var1./data_NN.Var1;
low_ratio = floor(min(ratios)/100)*100;
high_ratio = (floor((max(ratios)-1)/100)+1)*100;
h = histogram(ratios,low_ratio:25:high_ratio);
yticks([0:floor((max(h.BinCounts)-1)/4)+1:max(h.BinCounts), max(h.BinCounts)]);
rectangle('Position', [ median(data_sim.Var1./data_NN.Var1)-2, 0, 4, max(h.BinCounts)+3], ...
    'LineStyle', 'none', 'FaceColor', [1 0 0]);




flux_corr_perfect_mets_max = max([flux_corr_perfect_mets_CN, flux_corr_perfect_mets_C, flux_corr_perfect_mets_N],[],2);
[~, rxns_sorted_by_accuracy] = sort(flux_corr_perfect_mets_max);
high_accuracy_start = find(flux_corr_perfect_mets_CN(rxns_sorted_by_accuracy)>high_corr_threshold,1);
medium_accuracy_start = find(flux_corr_perfect_mets_CN(rxns_sorted_by_accuracy)>medium_corr_threshold,1);
nan_start = find(isnan(flux_corr_perfect_mets_CN(rxns_sorted_by_accuracy)),1);
figure(); set(gca,'fontsize',48);
set(gcf,"Color","w");
hold on
plot(NaN,NaN,"ksquare","MarkerSize",30);
plot(NaN,NaN,"kx","MarkerSize",30);

rectangle('Position', [ 0, -0.4, medium_accuracy_start-0.5, 1.4], ...
    'LineStyle', 'none', 'FaceColor', [0.9290 0.6940 0.1250]);
rectangle('Position', [ medium_accuracy_start-0.5, -0.4, high_accuracy_start-medium_accuracy_start, 1.4], ...
    'LineStyle', 'none', 'FaceColor', [0.8500 0.3250 0.0980]);
rectangle('Position', [ high_accuracy_start-0.5, -0.4, nan_start-high_accuracy_start, 1.4], ...
    'LineStyle', 'none', 'FaceColor', [0 0.4470 0.7410]);

rxns_for_cis = [108   133   184   205   207   209   217   219   230   257   260   266   267   239 ];
rxns_for_cis_in_sorted = find(ismember(rxns_sorted_by_accuracy(1:nan_start-1), rxns_for_cis))';
%scatter(rxns_for_cis_in_sorted,flux_corr_perfect_mets_max(rxns_sorted_by_accuracy(rxns_for_cis_in_sorted)),300,"square","g");
for r_x = rxns_for_cis_in_sorted
    rectangle('Position', [ r_x-0.5, -0.4, 1, 1.4], 'LineStyle', 'none', 'FaceColor', [1 1 1]);
end

b(1).FaceColor = [0 0.4470 0.7410]; % dark blue
b(2).FaceColor = [0.8500 0.3250 0.0980]; % dark orange
b(3).FaceColor = [0.9290 0.6940 0.1250];  % dark yellow
b(4).FaceColor = [0.4940 0.1840 0.5560];  % "dark purple "

plot([0, nan_start],[high_corr_threshold, high_corr_threshold], "k");
plot([0, nan_start],[medium_corr_threshold, medium_corr_threshold], "k");
plot([0, nan_start],[0, 0], "k");
xlim([0,nan_start]);

plot([high_accuracy_start-0.5, high_accuracy_start-0.5],[-0.4,1], "k");
plot([medium_accuracy_start-0.5, medium_accuracy_start-0.5],[-0.4,1], "k");
yticks([-0.4, 0, medium_corr_threshold, high_corr_threshold, 1]);
ylim([-0.4,1]);
xticks([medium_accuracy_start/2,  medium_accuracy_start/2 + high_accuracy_start/2, high_accuracy_start/2 + nan_start/2])
xticklabels([string(medium_accuracy_start-1), string(high_accuracy_start-medium_accuracy_start), ...
            string(nan_start-high_accuracy_start)])
axes(2).TickLength = [0,0];
ylabel(["Pearson correlation" "to true values"]);
xlabel("reactions");
scatter(1:nan_start-1,flux_corr_control(rxns_sorted_by_accuracy(1:nan_start-1)),400,"square","k", 'LineWidth', 1);
scatter(1:nan_start-1,flux_corr_perfect_mets_max(rxns_sorted_by_accuracy(1:nan_start-1)),400,'x',"k", 'LineWidth', 1);
legend(["random flux distribution",strcat("best estimated value,",string(newline),"using C and N data")], ...
"Position", [0.1,0.65,0.6,0.1], "Box", 0, "FontSize", 48 );

hold off


% num of high/medium/low acc rxns per supergroup
[~, supersystems_i_sorted] = sort(supersystems_high_corr_CN_n);
supersystems_data_CN = [supersystems_high_corr_CN_n, supersystems_medium_corr_CN_n, supersystems_low_corr_CN_n, supersystems_blocked_n];
supersystems_data_C = [supersystems_high_corr_C_n, supersystems_medium_corr_C_n, supersystems_low_corr_C_n, supersystems_blocked_n];
supersystems_data_N = [supersystems_high_corr_N_n, supersystems_medium_corr_N_n, supersystems_low_corr_N_n, supersystems_blocked_n];
supersystems_data3 = zeros(length(supersystems)*3, 4);
supersystems_data3(1:3:size(supersystems_data3,1),:)= supersystems_data_N(supersystems_i_sorted,:);
supersystems_data3(2:3:size(supersystems_data3,1),:)= supersystems_data_C(supersystems_i_sorted,:);
supersystems_data3(3:3:size(supersystems_data3,1),:)= supersystems_data_CN(supersystems_i_sorted,:);

figure(); hold on
set(gcf,"Color","w");
set(gca,'fontsize',30);
bar_indices0 = (1:3*length(supersystems))-1;
bar_indices = floor(bar_indices0/3)*4 + mod(bar_indices0,3) + 1;
bar_obj = barh(bar_indices, supersystems_data3,1,'stacked');
b(1).FaceColor = [0 0.4470 0.7410]; % dark blue
b(2).FaceColor = [0.8500 0.3250 0.0980]; % dark orange
b(3).FaceColor = [0.9290 0.6940 0.1250];  % dark yellow
b(4).FaceColor = [0.4940 0.1840 0.5560];  % "dark purple "

yticks(2:4:4*length(supersystems));
yticklabels("  "+supersystems_text(supersystems_i_sorted));
xlabel("Number of reactions");
%for line_i = 1:length(supersystems)-1
%    plot([0,200],[line_i*3 + 0.5, line_i*3 + 0.5],'k');
%end
lgd = legend(["High","Medium","Low","Blocked"],'Orientation','horizontal');
lgd.Box=0;
lgd.Position = [0.479861117576041 0.707575748464237 0.394270824749644 0.0406417101462256];
lgd.AutoUpdate = 0;
text([100,100,100], 39:3.5:46, ["N", "C", "C and N"],'VerticalAlignment','middle', 'FontSize',30);
plot([sum(supersystems_data3(34,:))+2,98], [45, 39], 'k');
plot([sum(supersystems_data3(35,:))+2,98], [46, 42.5], 'k');
plot([sum(supersystems_data3(36,:))+2,98], [47, 46], 'k');
% 34:36

hold off


cis_lb_CN = -1*ones(length(model.rxns),100);
cis_ub_CN = cis_lb_CN;
cis_lb_C = cis_lb_CN;
cis_ub_C = cis_lb_CN;
cis_lb_N = cis_lb_CN;
cis_ub_N = cis_lb_CN;
for rxn_i = 1:length(model.rxns)
    for rand_i = 1:100
        for label = ["C", "N", "CN" ]
            if isfile("results_ci"+label+"/results_confidence_intervals_metConstraints_"+string(rxn_i)+"_"+string(rand_i)+".mat")
                load("results_ci"+label+"/results_confidence_intervals_metConstraints_"+string(rxn_i)+"_"+string(rand_i), "ci_lb", "ci_ub");
                eval("cis_lb_"+label+"(rxn_i, rand_i) = ci_lb;");
                eval("cis_ub_"+label+"(rxn_i, rand_i) = ci_ub;");
            end
        end
    end
end

rxns_to_plot = [{[239   266  217  205 219 108   133   184      207   209 230   257   260  267]}];
rand_i_to_plot = [7 33 59 80 97];

plot_cis(rxns_to_plot, rand_i_to_plot, model, min_fluxes, cis_lb_CN, all_true_fluxes, cis_ub_CN, max_fluxes);

end
