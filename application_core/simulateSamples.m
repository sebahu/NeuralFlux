function [] = simulateSamples(model_file, split_dir, split_size, start_slot, end_slot, logsPerHour, simDurationHours, atom_label_fraction)
load(model_file, "model", "measured_emu_mids", "label_feed_rxn_ids", "label_feed_emus", ...
            "needed_emus", "needed_emu_map_mappings", "needed_emu_map_rxns");

% make this parameterizable
%pp = parpool(64);
%parfor slot_i = start_slot:end_slot
for slot_i = start_slot:end_slot
    disp("slot started: "+string(slot_i));
    slot_start_i = 1+(slot_i-1)*split_size;

    S = load(fullfile(string(split_dir),"samples_scaled_"+string(slot_start_i)+".mat"), 'samples_split_scaled');
    samples100_scaled = S.samples_split_scaled;
    S = load(fullfile(string(split_dir),"met_pools100_"+string(slot_start_i)+".mat"), 'met_pools100');
    met_pools100 = S.met_pools100;
    sample_ratios = single(zeros(length(measured_emu_mids),logsPerHour*simDurationHours+1,split_size)); %uint8
    sample_d_ratios = single(zeros(length(measured_emu_mids),logsPerHour*simDurationHours+1,split_size));
    for sample_i = 1:split_size %100
        simFluxesPerHour = samples100_scaled(1:end-1,sample_i); %scale it to ~ 6mgdw/gdw*h growth and scale to umol/gdw*h fluxes
        sim_fluxes = simFluxesPerHour/logsPerHour;
        met_pools = met_pools100(:,sample_i); %default: 20nmol/gdw, for umol/gdw*h fluxes

        [emu_mids, mapping_matrix, emu_mid_emu_is, emu_mid_met_inx, join_emu_mid_factor1_is, join_emu_mid_factor2_is] = ...
            create_d_ratio_from_ratio_matrix(model, needed_emus, needed_emu_map_mappings, needed_emu_map_rxns, sim_fluxes);

        [~,measured_mid_emu_mid_i] =ismember(measured_emu_mids, emu_mids);

        %parsave(strcat("100k_AraCore_samples_30_dir/mappingMatrix_",string(sample_i)),mapping_matrix);
        emu_pools = met_pools(emu_mid_met_inx);

        % find minimal step
        history_length = simDurationHours * logsPerHour;
        total_out_rate = sum(mapping_matrix(1:length(emu_pools),:),2)./emu_pools;
        max_flux_rate = max(total_out_rate);
        flux_rate_limit = 10000;
        if max_flux_rate > flux_rate_limit
            emu_pools(total_out_rate>flux_rate_limit) = emu_pools(total_out_rate>flux_rate_limit) .* total_out_rate(total_out_rate>flux_rate_limit)/flux_rate_limit;
            max_flux_rate = flux_rate_limit;
        end
        met_pools(emu_mid_met_inx) = emu_pools;

        steps_per_log = 20*ceil(max_flux_rate/10);

        normal_emus_num = length(emu_mid_met_inx);
        % this can either be parameterised or switched to a calculation
        % based on natural stable isotope occurrence
        emu_mid_start_values = calcStartMIDs(normal_emus_num, emu_mids, atom_label_fraction);

        [history_ratio, history_mids, mapping_matrix, history_d_ratio] = simMIDsEuler(sim_fluxes, ...
            history_length, mapping_matrix, label_feed_rxn_ids, label_feed_emus, met_pools, emu_mid_start_values, steps_per_log, ...
            emu_mid_met_inx, join_emu_mid_factor1_is, join_emu_mid_factor2_is);
        
        rounded_x = single(history_ratio(measured_mid_emu_mid_i,:)); %uint8(history_ratio(measured_mid_emu_mid_i,:)*100);
        sample_ratios(:,:,sample_i)=rounded_x;
        rounded_dx = single(history_d_ratio(measured_mid_emu_mid_i,:));
        sample_d_ratios(:,:,sample_i)=rounded_dx;
    end
    parsave(fullfile(string(split_dir),"mid_ratios_100_"+string(slot_start_i)+".mat"),sample_ratios);
    % at the moment, d_ratios are not used, thus also not stored
    % kept here for future use
    % parsave(fullfile(string(split_dir),"mid_d_ratios_100_"+string(slot_start_i)+".mat"),sample_d_ratios);    
end

end
