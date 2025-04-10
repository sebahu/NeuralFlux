function [corrected_atom_map_mapping, corrected_atom_map_mapping_with_ids] = create_corrected_atom_map_mapping( ...
    model, atom_map_rxns, atom_map_mapping)

atom_name_prefix_length = 2;
% first, have the atom mappings in correct order, i.e.
% substrate-atom=product-atom (reminder: the model has reversable
% reactions split into forward and reversed reaction)
corrected_atom_map_mapping_with_ids = atom_map_mapping;
for map_i = 1:length(corrected_atom_map_mapping_with_ids)
    met1 = extractBefore(corrected_atom_map_mapping_with_ids(map_i),"$");
    rxn_i = find(string(model.rxns)==atom_map_rxns(map_i));
    met1_i = find("M_"+string(model.mets)==met1);    
    met1_s = model.S(met1_i,rxn_i);
    if met1_s>0
        corrected_atom_map_mapping_with_ids(map_i) = extractAfter(extractAfter(atom_map_mapping(map_i),"="), atom_name_prefix_length) ...
            + "=" + extractAfter(extractBefore(atom_map_mapping(map_i),"="),atom_name_prefix_length);
    else
        corrected_atom_map_mapping_with_ids(map_i) = extractAfter(extractBefore(atom_map_mapping(map_i),"="),atom_name_prefix_length) ...
            + "=" + extractAfter(extractAfter(atom_map_mapping(map_i),"="), atom_name_prefix_length);        
    end
end
corrected_atom_map_mapping = replace(corrected_atom_map_mapping_with_ids, {'$A','$B','$C','$D'},"");
end
