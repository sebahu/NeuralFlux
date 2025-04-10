
function [] = decompartmentalizeSampleRatios(model_file, total_samples, split_size, split_dir, sample_ratios_file, sample_d_ratios_file, ...
                sample_decomp_ratios_file, sample_decomp_d_ratios_file)
    load(model_file, "measured_emu_mids", "measured_emus", "model", "atom_map_rxns", "atom_map_mapping", ...
         "corrected_atom_map_mapping", "needed_emus", "needed_emu_map_mappings", "needed_emu_map_rxns");
    %corrected_atom_map_mapping = create_corrected_atom_map_mapping(model, atom_map_rxns, atom_map_mapping);
    %[needed_emus, ~, ~] = create_needed_emus(atom_map_rxns, corrected_atom_map_mapping, measured_emus);

    [emu_mid_met_inx, ~, ~, emu_mids] = create_emu_data(model, needed_emus);
    load(string(split_dir)+"/met_pools100_1.mat", 'met_pools100');
    sample_met_pools = zeros(size(met_pools100,1),total_samples);
    for slot_start_i = 1:split_size:total_samples
        load(string(split_dir)+"/met_pools100_"+string(slot_start_i)+".mat", 'met_pools100');
        sample_met_pools(:,slot_start_i:(slot_start_i+split_size-1)) = met_pools100;
    end
    
    load(sample_ratios_file, "sample_ratios");    
    sample_decomp_ratios =  single(decompartmentalizeData(total_samples, sample_ratios, sample_met_pools, measured_emu_mids, emu_mids, emu_mid_met_inx));
    save(sample_decomp_ratios_file, "sample_decomp_ratios","-v7.3");
    clear sample_ratios;
    clear sample_decomp_ratios;
    %load(sample_d_ratios_file, "sample_d_ratios");
    %sample_decomp_d_ratios =  single(decompartmentalizeData(total_samples, sample_d_ratios, sample_met_pools, measured_emu_mids, emu_mids, emu_mid_met_inx));
    %save(sample_decomp_d_ratios_file, "sample_decomp_d_ratios","-v7.3");
end

function [sample_decomp_data, sample_decomp_met_pools] = decompartmentalizeData(total_samples, comp_data, sample_met_pools, measured_emu_mids, emu_mids, emu_mid_met_inx)
    measured_decomp_emu_mids = unique(extractBefore(measured_emu_mids,"[")+extractAfter(measured_emu_mids,"]"));
    %sample_decomp_ratios = zeros(length(measured_decomp_emu_mids),size(sample_ratios,2),total_samples,'uint8');
    sample_decomp_data = zeros(length(measured_decomp_emu_mids),size(comp_data,2),total_samples,'single');
    sample_decomp_met_pools = zeros(length(measured_decomp_emu_mids),total_samples);
    measured_emu_mid_met_inx = zeros(size(measured_emu_mids));
    for i = 1:length(measured_emu_mids)
        measured_emu_mid_met_inx(i) = emu_mid_met_inx(emu_mids == measured_emu_mids(i));
    end
    for decomp_emu_mid_i = 1 : length(measured_decomp_emu_mids)
        if startsWith(measured_decomp_emu_mids(decomp_emu_mid_i), "iCit:C#1,C#2,C#3,C#4,C#5,C#6.")
            continue
        end
        comp_measured_emu_mid_is = find(extractBefore(measured_emu_mids,"[")+extractAfter(measured_emu_mids,"]") == measured_decomp_emu_mids(decomp_emu_mid_i));
        comp_emu_mid_is = find(extractBefore(emu_mids,"[")+extractAfter(emu_mids,"]") == measured_decomp_emu_mids(decomp_emu_mid_i));
        comp_emu_mid_met_is = emu_mid_met_inx(comp_emu_mid_is);
        if length(comp_measured_emu_mid_is) == 1
            %sample_decomp_ratios(decomp_emu_mid_i,:,:) = sample_ratios(comp_measured_emu_mid_is,:,:);
            sample_decomp_data(decomp_emu_mid_i,:,:) = comp_data(comp_measured_emu_mid_is,:,:);
            sample_decomp_met_pools(decomp_emu_mid_i,:) = sample_met_pools(comp_emu_mid_met_is,:);
        else
            decomp_met_conc = sum(sample_met_pools(comp_emu_mid_met_is,:));
            sample_decomp_met_pools(decomp_emu_mid_i,:) = sum(sample_met_pools(comp_emu_mid_met_is,:));
            comp_met_fractions = sample_met_pools(comp_emu_mid_met_is,:)./decomp_met_conc;
            %tmp = squeeze(sum(double(permute(sample_ratios(comp_measured_emu_mid_is,:,:),[1,3,2])).*comp_met_fractions,1));
            %sample_decomp_ratios(decomp_emu_mid_i,:,:) = uint8(tmp');
            tmp = permute(comp_data(comp_measured_emu_mid_is,:,:),[1,3,2]);
            sample_decomp_data(decomp_emu_mid_i,:,:) = squeeze(sum(single(tmp).*single(comp_met_fractions)))';
        end
    end


end
