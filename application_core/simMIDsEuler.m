function [history_ratio, history_mids, mapping_matrix, history_d_ratio] = simMIDsEuler(sim_fluxes, ...
    history_length, mapping_matrix, label_feed_rxn_ids, label_feed_emus, met_conc, emu_mid_start_ratios, flux_part, ...
    emu_mid_met_inx, join_emu_mid_factor1_is, join_emu_mid_factor2_is)
% Performs simulation of an experiment with stable isotope labeling
%
% USAGE:
%
%    [history_ratio, history_unlabeled, history_labeled, mapping_matrix, high_flux_atoms] = simlabeling(sim_fluxes, history_length, mapping_matrix, label_feed_rxn_ids, label_feed_atoms, total_atoms, flux_part)
%
% INPUT:
%    sim_fluxes:       fluxes of the reactions in model to be simulated, in
%                      model.rxns. unit is per (time interval of one
%                      history step)
%    history_length:   number of simulated history steps (* sim_fluxes is
%                      the total flux), for which the labeling values are
%                      stored and returned as result of the function
%    mapping_matrix: matrix of atom transitions, contains the atom
%                      transitions from mid pool "row" to mid pool
%                      "column", extra rows are for combined emus in join
%                      reactions. 
%    label_feed_rxn_ids: reactions which feed labeled atoms into the network
%    label_feed_emu_mids: emu mids input reactions, which are labeled
%    met_concs  :      met pool sizes
%    emu_mid_start_values: values for the MIDs at start
%    flux_part:        how many mini-steps shall be simulated per history step
%    join_emu_mid_factor1_is, join_emu_mid_factor2_is: the ids to the mids
%                       which make up the extra rows

emu_mid_ratios = emu_mid_start_ratios;
normal_emus_num = length(emu_mid_ratios);
emu_mid_met_totals = met_conc(emu_mid_met_inx);
emu_mid_pools = emu_mid_met_totals .* emu_mid_start_ratios;

history_mids = zeros(size(emu_mid_ratios,1),history_length+1);
history_d_ratio = history_mids;
history_ratio = history_mids;
history_mids(:,1) = emu_mid_pools;
history_ratio(:,1) = emu_mid_ratios;

%minimal_flux = 1E-10;
%sanitize fluxes
%sim_fluxes(abs(sim_fluxes)<minimal_flux)=0;
base_rxn_flux = sim_fluxes/flux_part;


%precalculating some values, which are reused every step
base_mapping_matrix_in = sparse((mapping_matrix(1:normal_emus_num,1:end-1)'/flux_part));
base_mapping_matrix_in_special = sparse((mapping_matrix(normal_emus_num+1:end,1:end-1)'/flux_part));
total_out_flux = sum(mapping_matrix(1:normal_emus_num,:)/flux_part,2);

emu_mid_import = zeros(size(emu_mid_ratios));
if ~isempty(label_feed_emus)
for feed_rxn_i = 1:length(label_feed_rxn_ids)
    dest=label_feed_emus(feed_rxn_i);
    % relying on selection of feed rxns to have pos. flux
    emu_mid_import(dest) = emu_mid_import(dest) + base_rxn_flux(label_feed_rxn_ids(feed_rxn_i));
end
end

% history_length steps with flux_part mini steps of 1/flux_part poolsize
for history_step = 1:history_length    
    for mini_step = 1:flux_part
        % first calculate all outgoing atoms
        emu_mid_out_flux = total_out_flux .* emu_mid_ratios;
        emu_mid_in_flux = base_mapping_matrix_in*emu_mid_ratios + ...
            base_mapping_matrix_in_special * (emu_mid_ratios(join_emu_mid_factor1_is).*emu_mid_ratios(join_emu_mid_factor2_is));

        emu_mid_change = emu_mid_in_flux + emu_mid_import - emu_mid_out_flux;       
        %emu_mid_change(abs(emu_mid_change)<1e-16)=0;
       
        emu_mid_pools = max(min(emu_mid_pools + emu_mid_change,emu_mid_met_totals),0);
        emu_mid_ratios = emu_mid_pools./emu_mid_met_totals;

        if mini_step == 1
            history_d_ratio(:,history_step) = flux_part*emu_mid_change./emu_mid_met_totals;
        end
    end
    history_mids(:,history_step+1) = emu_mid_pools;
    history_ratio(:,history_step+1) = emu_mid_ratios;
end

end