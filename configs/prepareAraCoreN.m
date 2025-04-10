function [] = prepareAraCoreN(output_model_file) % curve_data )


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
atom_N_id_table = readtable(fullfile('models','all_atoms.N.sorted.txt'), 'ReadVariableNames', false, 'Delimiter', ' ');
atom_names = extractAfter(string(atom_N_id_table.Var1),atom_name_prefix_length);

AtomTransitionRDT_table = readtable(fullfile("models","all_mapping.N.sorted.txt"),"Delimiter"," ", 'ReadVariableNames', false);
atom_map_rxns = string(AtomTransitionRDT_table.Var1);
atom_map_mapping = string(AtomTransitionRDT_table.Var2);

% add mappings for reverse reactions
atom_map_rxns_b = ismember(atom_map_rxns,rev_rxns_ids);
atom_map_rxns_rev = regexprep(atom_map_rxns(atom_map_rxns_b),"(_[a-z])$","_rev$1");
atom_map_rxns_rev(~endsWith(atom_map_rxns_rev, "_rev_"+lettersPattern(1))) = cellstr(string(atom_map_rxns_rev(~endsWith(atom_map_rxns_rev, "_rev_"+lettersPattern(1))))+"_rev");
atom_map_mapping = [ atom_map_mapping; atom_map_mapping(atom_map_rxns_b)];
atom_map_rxns = [ atom_map_rxns; atom_map_rxns_rev ];

amino_acids = ["Gln"; "Asp"; "Glu"; "Asn"; "Ser"; "Cys"; "Thr"; "Gly"; "Met"; "Pro"; 
                    "Ala"; "Arg"; "Lys"; "His"; "Ile"; "Leu"; "Phe"; "Trp"; "Tyr"; "Val"];

amino_acid_atom_nums = [2; 1; 1; 2; 1; 1; 1; 1; 1; 1;
                        1; 4; 2; 3; 1; 1; 1; 2; 1; 1 ];

additional_emus = ["GABA[c]:N#1", "GABA[m]:N#1", "GTP[c]:N#1,N#2,N#3,N#4,N#5", "GTP[h]:N#1,N#2,N#3,N#4,N#5", "CTP[c]:N#1,N#2,N#3", "UTP[c]:N#1,N#2", "urea[m]:N#1,N#2"];
additional_emu_atom_nums = [1, 1, 5, 5, 3, 2, 2];

measured_emus = string(zeros(length(amino_acids)*4+length(additional_emus),1));
measured_emu_mids = string(zeros(sum(amino_acid_atom_nums+1)*4+sum(additional_emu_atom_nums+1),1));

current_mid = 0;
for aa_i = 1:length(amino_acids)
    measured_emus(aa_i*4+(-3:0)) = amino_acids(aa_i)+["[c]";"[h]";"[m]";"[p]"]+":"+strjoin("N#"+string(1:amino_acid_atom_nums(aa_i)),",");
    mids_added = 4*(amino_acid_atom_nums(aa_i)+1);
    measured_emu_mids(current_mid+(1:mids_added)) = ...
        reshape((amino_acids(aa_i)+["[c]";"[h]";"[m]";"[p]"]+":"+ ...
                strjoin("N#"+string(1:amino_acid_atom_nums(aa_i)),",")+"."+string(0:amino_acid_atom_nums(aa_i)))',1,mids_added)';
    current_mid = current_mid + mids_added;
end
for additional_emu_i = 1:length(additional_emus)
    measured_emus(length(amino_acids)*4+additional_emu_i) = additional_emus(additional_emu_i);
    mids_added = (additional_emu_atom_nums(additional_emu_i)+1);
    measured_emu_mids(current_mid+(1:mids_added)) = ...
        (additional_emus(additional_emu_i)+"."+string(0:additional_emu_atom_nums(additional_emu_i)))';
    current_mid = current_mid + mids_added;
end

%nitrate and ammonium seem to be the source of choice of N in this setting ...
nitrate_feed_rxn_id = find(string(model.rxns) == 'Im_NO3');
nh3_feed_rxn_id = find(string(model.rxns) == 'Im_NH4');
label_feed_rxn_ids = [ nitrate_feed_rxn_id, nh3_feed_rxn_id];


[corrected_atom_map_mapping, corrected_atom_map_mapping_with_ids] = create_corrected_atom_map_mapping(model, atom_map_rxns, atom_map_mapping);
[needed_emus, needed_emu_map_mappings, needed_emu_map_rxns] = create_needed_emus(atom_map_rxns, corrected_atom_map_mapping, ...
                                                                        corrected_atom_map_mapping_with_ids, measured_emus, model);
[~, ~, ~, emu_mids] = create_emu_data(model, needed_emus);

nitrate_import_id = find(emu_mids == 'NO3[c]:N#1.1');
nh3_import_id = find(emu_mids == 'NH4[c]:N#1.1');
label_feed_emus = [ nitrate_import_id, nh3_import_id ];

rnullspace = null(model.S,'r');

mets_no_cmp = extractBefore(string(model.mets),"[");
all_species_id = unique(mets_no_cmp);
N_species_lit_conc_table = readtable(fullfile('models','N_relevant_species_no_cmp_literature_concentration.txt'), 'ReadVariableNames', false, 'Delimiter', '\t');
N_species_id = string(N_species_lit_conc_table.Var1);
[~,special_species_loc] = ismember(N_species_id,all_species_id);
special_species_conc = N_species_lit_conc_table.Var2 * 0.011; %lit values are nmol/gFW

save(output_model_file, "atom_map_rxns", "atom_map_mapping", "atom_names", "model", ...
    "rnullspace", "measured_emu_mids", "measured_emus", "label_feed_rxn_ids", "label_feed_emus", ...
    "special_species_loc", "special_species_conc", "needed_emus", "corrected_atom_map_mapping_with_ids", "emu_mids", ...
    "corrected_atom_map_mapping", "needed_emu_map_mappings", "needed_emu_map_rxns", "-v7.3");


end