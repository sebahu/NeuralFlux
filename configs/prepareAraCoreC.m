function [] = prepareAraCoreC(output_model_file) % curve_data )


% biomass is around 570000 g/mol, opt biomass flux around 0.005
% 0.005 umol/gdw*h  * 570000 g/mol = 2850 ug/gdw*h

load(fullfile('models','AraCore_v2_1.mat'), 'model');


model_bio1 = 'Bio_opt';
model_bio1_rxn_index = find(string(model.rxns) == model_bio1);

% split reversible reactions (lb < 0)
rev_rxns = model.lb < 0;
rev_rxns_ids = string(model.rxns(rev_rxns));
rev_rxns_num = length(rev_rxns_ids);
num_rxns = length(model.rxns);

model.c(num_rxns+(1:rev_rxns_num)) = 0;
model.grRules(num_rxns+(1:rev_rxns_num)) = model.grRules(rev_rxns);
model.lb(num_rxns+(1:rev_rxns_num)) = 0;
model.ub(num_rxns+(1:rev_rxns_num)) = -model.lb(rev_rxns);
model.lb(rev_rxns) = 0;
model.rules(num_rxns+(1:rev_rxns_num)) = model.rules(rev_rxns);
model.rxnConfidenceScores(num_rxns+(1:rev_rxns_num)) = model.rxnConfidenceScores(rev_rxns);
model.rxnECNumbers(num_rxns+(1:rev_rxns_num)) = model.rxnECNumbers(rev_rxns);
model.rxnNames(num_rxns+(1:rev_rxns_num)) = cellstr(string(model.rxnNames(rev_rxns)) + " (reversed)");
model.rxnNotes(num_rxns+(1:rev_rxns_num)) = model.rxnNotes(rev_rxns);
model.subSystems(num_rxns+(1:rev_rxns_num)) = model.subSystems(rev_rxns);
rev_rxns_arr = regexprep(model.rxns(rev_rxns),"(_[a-z])$","_rev$1");
rev_rxns_arr(~endsWith(rev_rxns_arr, "_rev_"+lettersPattern(1))) = cellstr(string(rev_rxns_arr(~endsWith(rev_rxns_arr, "_rev_"+lettersPattern(1))))+"_rev");
model.rxns(num_rxns+(1:rev_rxns_num)) = rev_rxns_arr;

model.S(:,num_rxns+(1:rev_rxns_num)) = -model.S(:,rev_rxns);
model.A(:,num_rxns+(1:rev_rxns_num)) = -model.A(:,rev_rxns);

% limits for exports
model.ub(443:522)=0;
model.ub([443,447,451,455,459,463,471,475,479,483,487,491,495,499,503,507,511,515,519])=1;
model.ub(467) = 2; % Glu_c

%model.ub(443:522)=maxFlux(443:522)/100;


% mappings
atom_name_prefix_length = 2;
atom_C_id_table = readtable(fullfile('models','all_atoms.C.sorted.txt'), 'ReadVariableNames', false, 'Delimiter', ' ');
atom_names = extractAfter(string(atom_C_id_table.Var1),atom_name_prefix_length);

AtomTransitionRDT_table = readtable(fullfile("models","all_mapping.C.sorted.txt"),"Delimiter"," ", 'ReadVariableNames', false);
atom_map_rxns = string(AtomTransitionRDT_table.Var1);
atom_map_mapping = string(AtomTransitionRDT_table.Var2);

% add mappings for reverse reactions
atom_map_rxns_b = ismember(atom_map_rxns,rev_rxns_ids);
atom_map_rxns_rev = regexprep(atom_map_rxns(atom_map_rxns_b),"(_[a-z])$","_rev$1");
atom_map_rxns_rev(~endsWith(atom_map_rxns_rev, "_rev_"+lettersPattern(1))) = cellstr(string(atom_map_rxns_rev(~endsWith(atom_map_rxns_rev, "_rev_"+lettersPattern(1))))+"_rev");
atom_map_mapping = [ atom_map_mapping; atom_map_mapping(atom_map_rxns_b)];
atom_map_rxns = [ atom_map_rxns; atom_map_rxns_rev ];

% TODO

amino_acids = ["Gln"; "Asp"; "Glu"; "Asn"; "Ser"; "Cys"; "Thr"; "Gly"; "Met"; "Pro";
                    "Ala"; "Arg"; "Lys"; "His"; "Ile"; "Leu"; "Phe"; "Trp"; "Tyr"; "Val"];

amino_acid_emu_mets = reshape(amino_acids + ["[c]","[h]","[m]","[p]"],[],1);

additional_emu_mets = ["2PGA[c]"; "2PGA[h]"; "6PG[h]"; "Cit[c]"; "Cit[h]"; "Cit[m]"; "iCit[c]"; "iCit[h]"; "iCit[m]"; ...
    "CDP[c]"; "CTP[c]"; "GDP[c]"; "GDP[h]"; "GMP[c]"; "GTP[c]"; "GTP[h]"; "UDP[c]"; "UMP[c]"; "UTP[c]"; ...
    "F6P[c]"; "F6P[h]"; "FBP[c]"; "FBP[h]"; "Fum[c]"; "Fum[h]"; "Fum[m]"; "Glc[c]"; "Glc[h]"; "Ru5P[h]"; "RuBP[h]"; ...
    "THF[c]"; "THF[h]"; "THF[m]"];
% todo add "iCit[p]"
%"starch1[c]"; "starch1[h]"; "starch2[c]"; "starch2[h]"; "starch3[h]"; "starch5[h]"; "cellulose1[c]"; "cellulose2[c]"; "cellulose3[c]"];

emu_mets = [amino_acid_emu_mets; additional_emu_mets];

emu_atoms = zeros(size(emu_mets));
for emu_i = 1:length(emu_mets)
    emu_atoms(emu_i) = max(double(extractAfter(atom_names(startsWith(atom_names, emu_mets(emu_i)+":")),"#")));
end

measured_emus = emu_mets;
measured_emu_mids = string(zeros(sum(emu_atoms+1),1));

current_mid = 0;
for emu_i = 1:length(emu_mets)
    measured_emus(emu_i) = emu_mets(emu_i)+":"+strjoin("C#"+string(1:emu_atoms(emu_i)),",");
    measured_emu_mids(current_mid+(1:(emu_atoms(emu_i)+1))) = ...
        measured_emus(emu_i)+"."+string(0:emu_atoms(emu_i))';
    current_mid = current_mid + emu_atoms(emu_i) + 1;
end
%measured_emu_mids = measured_emu_mids(endsWith(measured_emu_mids,".0") | endsWith(measured_emu_mids,".1") | ...
%    endsWith(measured_emu_mids,".2") | endsWith(measured_emu_mids,".3"));

%nitrate and ammonium seem to be the source of choice of N in this setting ...
co2_feed_rxn_id = find(string(model.rxns) == 'Im_CO2');
label_feed_rxn_ids = [ co2_feed_rxn_id];

% END TODO

[corrected_atom_map_mapping, corrected_atom_map_mapping_with_ids] = create_corrected_atom_map_mapping(model, atom_map_rxns, atom_map_mapping);
[needed_emus, needed_emu_map_mappings, needed_emu_map_rxns] = create_needed_emus(atom_map_rxns, corrected_atom_map_mapping, ...
                                                                            corrected_atom_map_mapping_with_ids, measured_emus, model);
[~, ~, ~, emu_mids] = create_emu_data(model, needed_emus);

% TODO
co2_import_id = find(emu_mids == 'CO2[c]:C#1.1');
label_feed_emus = [ co2_import_id ];
% END TODO

rnullspace = null(model.S,'r');

mets_no_cmp = extractBefore(string(model.mets),"[");
all_species_id = unique(mets_no_cmp);
N_species_lit_conc_table = readtable(fullfile('models','N_relevant_species_no_cmp_literature_concentration.txt'), 'ReadVariableNames', false, 'Delimiter', '\t');
N_species_id = string(N_species_lit_conc_table.Var1);
[~,special_species_loc] = ismember(N_species_id,all_species_id);
special_species_conc = N_species_lit_conc_table.Var2 * 0.011; %lit values are nmol/gFW

save(output_model_file, "atom_map_rxns", "atom_map_mapping", "atom_names", "model", ...
    "rnullspace", "measured_emu_mids", "measured_emus", "label_feed_rxn_ids", "label_feed_emus", ...
    "special_species_loc", "special_species_conc", "needed_emus", "corrected_atom_map_mapping_with_ids", ...
    "corrected_atom_map_mapping", "needed_emu_map_mappings", "needed_emu_map_rxns", "emu_mids", "-v7.3");


end