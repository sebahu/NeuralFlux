function [] = createMetConcSampling(model_file, total_samples, split_size, split_dir, all_samples_file)
    load(model_file, "model", "special_species_loc", "special_species_conc");
    mets_no_cmp = extractBefore(string(model.mets),"[");
    all_species_id = unique(mets_no_cmp);
    [~, met_loc] = ismember(mets_no_cmp, all_species_id);
    comp_mhpc_exist = zeros(length(all_species_id),4);
    comp_mhpc_exist(ismember(all_species_id+"[m]",string(model.mets)),1)=1;
    comp_mhpc_exist(ismember(all_species_id+"[h]",string(model.mets)),2)=1;
    comp_mhpc_exist(ismember(all_species_id+"[p]",string(model.mets)),3)=1;
    comp_mhpc_exist(ismember(all_species_id+"[c]",string(model.mets)),4)=1;
    link2 = zeros(length(mets_no_cmp),1);
    link2(endsWith(string(model.mets),"[m]"))=0;
    link2(endsWith(string(model.mets),"[h]"))=length(all_species_id);
    link2(endsWith(string(model.mets),"[p]"))=2*length(all_species_id);
    link2(endsWith(string(model.mets),"[c]"))=3*length(all_species_id);
    

    
    for slot_start_i = 1:split_size:total_samples
        m_part = 0.05+rand(split_size,1)/10;
        h_part = 0.15+rand(split_size,1)/4;
        p_part = 0.01+rand(split_size,1)/25;
        mhpc_part = [m_part, h_part, p_part, 1-m_part-h_part-p_part];
        rel_conc = permute((1+rand(length(all_species_id),4,split_size)) .* comp_mhpc_exist,[3,2,1]) .* mhpc_part;
        sum_rel_conc = squeeze(sum(rel_conc,2));
        rel_conc_2d = reshape(permute(rel_conc,[1,3,2]), split_size, []); % TODO change here, for-loop ?
        total_met_conc = 0.01 * (10 .^ rand(split_size,length(all_species_id)));
        total_met_conc(:,special_species_loc) = 20*total_met_conc(:,special_species_loc) .* special_species_conc';
        met_pools100 = (total_met_conc(:,met_loc) .* rel_conc_2d(:,met_loc+link2) ./ sum_rel_conc(:,met_loc))';
        
        save(string(split_dir)+"/met_pools100_"+string(slot_start_i)+".mat", 'met_pools100', '-v7.3');
        if slot_start_i == 1
            sample_met_pools = zeros(size(met_pools100,1),total_samples);
        end
        sample_met_pools(:,slot_start_i:(slot_start_i+split_size-1)) = met_pools100;
    end
    save(all_samples_file, "sample_met_pools", "-v7.3");
end
