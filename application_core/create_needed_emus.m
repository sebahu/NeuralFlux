function [needed_emus, needed_emu_map_mappings, needed_emu_map_rxns] = create_needed_emus( ...
    atom_map_rxns, corrected_atom_map_mapping, corrected_atom_map_mapping_with_ids, measured_emus, model)
% this function calculates all needed emus, which are necessary to
% determine the state of the measured emus in an enrichment simulation
% input: atom_map_rxns: the rxn ids, corresponding to the mappings in
% corrected_atom_map_mapping
% corrected_atom_map_mapping: the atom mappings, one per atom pair, in the
% direction of the irreversible reactions
% (i.e.: substrate atom=product atom), format for atoms:
% <metabolite>:<element, e.g. C>#<atom number of this element in the
% metabolite according to InChI>
% corrected_atom_map_mapping_with_ids: the same, but also with molecule id,
% in case more than molecule of the same metabolite is involved in the
% reaction: <metabolite>$<ID>:<element ....
% measured_emus: the measured emus, containing several atoms of a
% metabolite molecule:  <metabolite>:[<element, e.g. C>#<atom number>,]*[<element, e.g. C>#<atom number>,]
%
% return values:
% needed_emus: the emus, which are necessary to
% determine the state of the measured emus in an enrichment simulation
% needed_emu_map_mappings: all mappings, <emu1>=<emu2>, which are necessary to
% determine the state of the measured emus in an enrichment simulation
% needed_emu_map_rxns: the rxn ids, corresponding to the mappings

corrected_atom_map_mapping_prod_atoms = extractAfter(corrected_atom_map_mapping,"=");
% starting point: the measured emus are obviously also needed
needed_emus = measured_emus;
needed_emu_map_mappings = string([]);
needed_emu_map_rxns = string([]);
last_handled_needed_emu_i = 0;

% add needed emus until no further are needed
while length(needed_emus) > last_handled_needed_emu_i
    last_handled_needed_emu_i = last_handled_needed_emu_i+1;
    % work on the next, not yet covered emu
    needed_emu = needed_emus(last_handled_needed_emu_i);
    prod_met = extractBefore(needed_emu,":");
    % create the individual needed atoms, as used in corrected_atom_map_mapping
    needed_atoms = prod_met + ":"+ split(extractAfter(needed_emu,":"),",");
    % all mappings, which have the metabolite of the needed emu also as
    % metabolite, are candidates to be necessary for the current emu
    candidate_mappings_b = startsWith(extractAfter(corrected_atom_map_mapping,"="),prod_met+":");
    % and the corresponding rxns
    needed_rxns = intersect(unique(atom_map_rxns(candidate_mappings_b)),string(model.rxns));
    % check each cadidate rxn
    for rxn_i = 1:length(needed_rxns)
        needed_rxn = needed_rxns(rxn_i);
        % only check mappings for the emu's atoms
        needed_mappings_b = atom_map_rxns == needed_rxn & ...
            ismember(corrected_atom_map_mapping_prod_atoms, needed_atoms);
        needed_mappings = corrected_atom_map_mapping(needed_mappings_b);
        % we could have >1 instance of the metabolite in this reaction -.
        % we need the molecule ids
        needed_mappings_with_id = corrected_atom_map_mapping_with_ids(needed_mappings_b);
        while ~ isempty(needed_mappings)
            current_prod_id = [];
            substr_atoms_with_id = string([]);
            handled_mappings_b = false(size(needed_mappings));
            for map_i = 1:length(needed_mappings)
                mapping_w_id = needed_mappings_with_id(map_i);
                % checking the molecule id
                prod_id = extractBefore(extractAfter(extractAfter(mapping_w_id+"$:","="),"$"),":");
                if isempty(current_prod_id)
                    % if we just start a new instance, remember its id
                    current_prod_id = prod_id;
                end
                % only work on mappings of the current metabolite instance
                if prod_id == current_prod_id
                    % add the current mappings's substrate atom to the
                    % substrate emu, mark the current mapping as done
                    substr_atoms_with_id(end+1) = extractBefore(mapping_w_id,"=");
                    handled_mappings_b(map_i) = true;
                end
                % current metabolite instance is also done, when enough
                % atoms are mapped
                if length(substr_atoms_with_id) == length(needed_atoms)
                    break;
                end
            end
            % remove handled mappings from list of candidates
            needed_mappings = needed_mappings(~handled_mappings_b);
            needed_mappings_with_id = needed_mappings_with_id(~handled_mappings_b);
            substr_mets_with_id = unique(extractBefore(substr_atoms_with_id, ":"));
            new_emu_mapping = "";
            % construct new needed emu and mappings, note: it could be a
            % join-emu, i.e. emu1+emu2, even from two different metabolites
            for substr_met_i = 1:length(substr_mets_with_id)
                substr_met_with_id = substr_mets_with_id(substr_met_i);
                substr_met_no_id = extractBefore(substr_met_with_id+"$","$");
                substr_atoms_only = extractAfter(substr_atoms_with_id(startsWith(substr_atoms_with_id,substr_met_with_id+":")),":");
                substr_atoms_numbers = double(extractAfter(substr_atoms_only,"#"));
                % atom numbers are ordered numerically in emus
                [~,atom_order] = sort(substr_atoms_numbers);
                substr_emu = substr_met_no_id + ":" + join(substr_atoms_only(atom_order),",");
                if ~ismember(substr_emu, needed_emus)
                    needed_emus(end+1) = substr_emu;
                end
                new_emu_mapping = new_emu_mapping + "+"+substr_emu;
            end
            new_emu_mapping = extractAfter(new_emu_mapping,1) + "=" + needed_emu;
            needed_emu_map_mappings(end+1,1) = new_emu_mapping;
            needed_emu_map_rxns(end+1,1) = needed_rxn;
        end
    end
end
needed_emus = sort(needed_emus);
%end

