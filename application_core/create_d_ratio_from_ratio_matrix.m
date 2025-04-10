function [emu_mids, mapping_matrix, emu_mid_emu_is, emu_mid_met_inx, join_emu_mid_factor1_is, join_emu_mid_factor2_is] = ...
    create_d_ratio_from_ratio_matrix(model, needed_emus, needed_emu_map_mappings, needed_emu_map_rxns, sim_fluxes)

% this is the matrix used during simulation

% result: matrix: #emu_mids * #emu_mids+1
[emu_mid_met_inx, emu_met_inx, emu_mid_emu_is, emu_mids] = create_emu_data(model, needed_emus);

normal_emu_mids = emu_mids;

% now the mappings with special "+" join-EMUs
join_emu_mappings_is = find(contains(needed_emu_map_mappings,"+"));
join_emu_mid_factor1_is = [];
join_emu_mid_factor2_is = [];

for join_emu_mapping_ii = 1:length(join_emu_mappings_is)
    join_emu_mapping_i = join_emu_mappings_is(join_emu_mapping_ii);
    join_emu_mapping = needed_emu_map_mappings(join_emu_mapping_i);
    join_emu1 = extractBefore(join_emu_mapping,"+");
    join_emu2 = extractBefore(extractAfter(join_emu_mapping,"+"),"=");
    join_emu1_max_mid = sum(startsWith(normal_emu_mids,join_emu1+"."))-1;
    join_emu2_max_mid = sum(startsWith(normal_emu_mids,join_emu2+"."))-1;
    join_emu_max_mid = join_emu1_max_mid+join_emu2_max_mid; 
    for mid_sum = 0:join_emu_max_mid
        for mid1 = max(0,mid_sum-join_emu2_max_mid):min(mid_sum,join_emu1_max_mid)
            mid2 = mid_sum-mid1;
            new_emu_mid = join_emu1+"+"+join_emu2+"."+string(mid_sum)+","+string(mid1)+","+string(mid2);
            if ~ismember(new_emu_mid, emu_mids)
                emu_mids(end+1) = new_emu_mid;
                join_emu_mid_factor1_is(end+1) = find(normal_emu_mids==join_emu1+"."+string(mid1));
                join_emu_mid_factor2_is(end+1) = find(normal_emu_mids==join_emu2+"."+string(mid2));
            end
        end
    end
end

mapping_matrix = zeros(length(emu_mids),length(normal_emu_mids)+1);

for rxn_i = 1:length(model.rxns)
    rxn = model.rxns(rxn_i);
    current_map_is = find(needed_emu_map_rxns == rxn);
    emus_with_mappings = string([]);
    for map_ii = 1:length(current_map_is)
        map_i = current_map_is(map_ii);
        current_mapping = needed_emu_map_mappings(map_i);
        sourceEMU = extractBefore(current_mapping,"=");
        sourceMet = extractBefore(sourceEMU,":");
        source_i = string(model.mets) == sourceMet;
        if ~contains(sourceEMU,"+")
            emus_with_mappings(end+1) = sourceEMU;
        end
        destEMU = extractAfter(current_mapping,"=");
        destMet = extractBefore(destEMU,":");
        dest_i = string(model.mets) == destMet;
        emus_with_mappings(end+1)=destEMU;
        if(sum(source_i)+sum(dest_i) ==0)
            continue
        end
        S_value = model.S(source_i, rxn_i);
        % Workaround: we have two different cases, where S_value is different
        % from 1: a) two (or more) molecules of the same metabolite are involved. Then,
        % there should be as many mappings, thus we need to use 1 here b) a
        % fractional amount is given, as in biomass reaction - we take the
        % fraction - workaround solution: if it is a whole number, take 1
        if S_value == round(S_value)
            S_value=1;
        end
        % now enter the value for all related MIDs
        dest_MID_is = find(startsWith(emu_mids, destEMU+"."));
        for dest_ii = 1:length(dest_MID_is)
            dest_i = dest_MID_is(dest_ii);
            % this also handles the case with >1 emu_mid as source, eg.
            % A+B.1,0,1 and A+B.1,1,0 
            source_b = startsWith(emu_mids+",", sourceEMU + "." + extractAfter(emu_mids(dest_i),".")+",");
            mapping_matrix(source_b,dest_i) = mapping_matrix(source_b,dest_i)+S_value*sim_fluxes(rxn_i);
        end
    end
    % now handle all source EMUs (needed for reactions like GMPS_c and CTPS_c) without mapping as an export
    rxn_source_emus = needed_emus(model.S(emu_met_inx,rxn_i)<0);
    rxn_export_emus = setdiff(rxn_source_emus, emus_with_mappings);
    for emu_i = 1:length(rxn_export_emus)
        export_met_i = emu_met_inx(needed_emus == rxn_export_emus(emu_i));
        rxn_export_mids_b = startsWith( emu_mids, rxn_export_emus(emu_i)+"." );
        % it is export, S is negative, but we use positive values -> -
        mapping_matrix(rxn_export_mids_b,end) = mapping_matrix(rxn_export_mids_b,end) - model.S(export_met_i,rxn_i)*sim_fluxes(rxn_i);
    end
end


end
