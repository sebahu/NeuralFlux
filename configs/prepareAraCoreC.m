function [] = prepareAraCoreC(output_model_file) % curve_data )
% a script to adapt the existing model to any needs of NeuralFlux
% most importantly: NeuralFlux handles reactions as irreversible
% reversible reactions need to to be split in forward and reverse reaction

% this specific script modifies the AraCore v2.1 model

% note: biomass is around 570000 g/mol, opt biomass flux around 0.005
% 0.005 umol/gdw*h  * 570000 g/mol = 2850 ug/gdw*h

load(fullfile('models','AraCore_v2_1.mat'), 'model');


% split reversible reactions (lb < 0)
% find reversible reactions and store some variables
rev_rxns = model.lb < 0;
rev_rxns_ids = string(model.rxns(rev_rxns));
rev_rxns_num = length(rev_rxns_ids);
num_rxns = length(model.rxns);

% the new reversed reactions are not part of the objective
model.c(num_rxns+(1:rev_rxns_num)) = 0;
% the new reversed reactions use the grRules of the reversible reactions
model.grRules(num_rxns+(1:rev_rxns_num)) = model.grRules(rev_rxns);
% the new reversed reactions are not reversible
model.lb(num_rxns+(1:rev_rxns_num)) = 0;
% the new reversed reactions have as upper boundary the reversed
% lower boundary of the reversible reactions
model.ub(num_rxns+(1:rev_rxns_num)) = -model.lb(rev_rxns);
% now, the original reversible reactions are made irreversible (forward)
model.lb(rev_rxns) = 0;
% the new reversed reactions use the rules of the reversible reactions
model.rules(num_rxns+(1:rev_rxns_num)) = model.rules(rev_rxns);
% the new reversed reactions use the rxnConfidenceScores of the reversible reactions
model.rxnConfidenceScores(num_rxns+(1:rev_rxns_num)) = model.rxnConfidenceScores(rev_rxns);
% the new reversed reactions get the rxnECNumbers of the reversible reactions
model.rxnECNumbers(num_rxns+(1:rev_rxns_num)) = model.rxnECNumbers(rev_rxns);
% the new reversed reactions get an appropriate name
model.rxnNames(num_rxns+(1:rev_rxns_num)) = cellstr(string(model.rxnNames(rev_rxns)) + " (reversed)");
% the new reversed reactions get the rxnNotes of the reversible reactions
model.rxnNotes(num_rxns+(1:rev_rxns_num)) = model.rxnNotes(rev_rxns);
% the new reversed reactions get the subSystems of the reversible reactions
model.subSystems(num_rxns+(1:rev_rxns_num)) = model.subSystems(rev_rxns);
% the new reversed reactions get the string "rev" inserted into their id (before the compartment)
rev_rxns_arr = regexprep(model.rxns(rev_rxns),"(_[a-z])$","_rev$1");
rev_rxns_arr(~endsWith(rev_rxns_arr, "_rev_"+lettersPattern(1))) = cellstr(string(rev_rxns_arr(~endsWith(rev_rxns_arr, "_rev_"+lettersPattern(1))))+"_rev");
model.rxns(num_rxns+(1:rev_rxns_num)) = rev_rxns_arr;

% the new reversed reactions get reversed entries in the Stoichimetric matrix
model.S(:,num_rxns+(1:rev_rxns_num)) = -model.S(:,rev_rxns);
model.A(:,num_rxns+(1:rev_rxns_num)) = -model.A(:,rev_rxns);

% limits for exports (indices are from knowledge of AraCore v2.1
% setting all complex exports to 0
model.ub(443:522)=0;
% allowing amino acid exports from the cell plams
model.ub([443,447,451,455,459,463,471,475,479,483,487,491,495,499,503,507,511,515,519])=1;
% allowing glucose export from the cell plasm
model.ub(467) = 2; % Glu_c

%model.ub(443:522)=maxFlux(443:522)/100;


% adding atom transition mapping information for the model
% this is not part of the AraCore v2.1 model, it is provided as additional text files
% first the atoms
atom_name_prefix_length = 2;
atom_C_id_table = readtable(fullfile('models','all_atoms.C.sorted.txt'), 'ReadVariableNames', false, 'Delimiter', ' ');
atom_names = extractAfter(string(atom_C_id_table.Var1),atom_name_prefix_length);

% now the mappings, format: reaction_id <substrate_atom_id>=<product_atom_id>
AtomTransitionRDT_table = readtable(fullfile("models","all_mapping.C.sorted.txt"),"Delimiter"," ", 'ReadVariableNames', false);
atom_map_rxns = string(AtomTransitionRDT_table.Var1);
atom_map_mapping = string(AtomTransitionRDT_table.Var2);

% add mappings for reverse reactions: duplicates of the mappings for the
% original reversible reaction, with adapted rxn id (inserted "rev")
atom_map_rxns_b = ismember(atom_map_rxns,rev_rxns_ids);
atom_map_rxns_rev = regexprep(atom_map_rxns(atom_map_rxns_b),"(_[a-z])$","_rev$1");
atom_map_rxns_rev(~endsWith(atom_map_rxns_rev, "_rev_"+lettersPattern(1))) = cellstr(string(atom_map_rxns_rev(~endsWith(atom_map_rxns_rev, "_rev_"+lettersPattern(1))))+"_rev");
atom_map_mapping = [ atom_map_mapping; atom_map_mapping(atom_map_rxns_b)];
atom_map_rxns = [ atom_map_rxns; atom_map_rxns_rev ];

% the measured EMU MIDs
% 1st, collect the metabolites
%  - all amino acids, in all compartents
amino_acids = ["Gln"; "Asp"; "Glu"; "Asn"; "Ser"; "Cys"; "Thr"; "Gly"; "Met"; "Pro";
                    "Ala"; "Arg"; "Lys"; "His"; "Ile"; "Leu"; "Phe"; "Trp"; "Tyr"; "Val"];

amino_acid_emu_mets = reshape(amino_acids + ["[c]","[h]","[m]","[p]"],[],1);

% some additional metablites from CBC and TCA cycle, as they are not present in all compartments
% they are specificially listed with compartment
additional_emu_mets = ["2PGA[c]"; "2PGA[h]"; "6PG[h]"; "Cit[c]"; "Cit[h]"; "Cit[m]"; "iCit[c]"; "iCit[h]"; "iCit[m]"; ...
    "CDP[c]"; "CTP[c]"; "GDP[c]"; "GDP[h]"; "GMP[c]"; "GTP[c]"; "GTP[h]"; "UDP[c]"; "UMP[c]"; "UTP[c]"; ...
    "F6P[c]"; "F6P[h]"; "FBP[c]"; "FBP[h]"; "Fum[c]"; "Fum[h]"; "Fum[m]"; "Glc[c]"; "Glc[h]"; "Ru5P[h]"; "RuBP[h]"; ...
    "THF[c]"; "THF[h]"; "THF[m]"];
% todo add "iCit[p]"
%"starch1[c]"; "starch1[h]"; "starch2[c]"; "starch2[h]"; "starch3[h]"; "starch5[h]"; "cellulose1[c]"; "cellulose2[c]"; "cellulose3[c]"];

emu_mets = [amino_acid_emu_mets; additional_emu_mets];

% converting into EMUs: take the longest EMU (containing n C atoms, where
% n is the max atom numer of the metabolite from the mappings
emu_atoms = zeros(size(emu_mets));
for emu_i = 1:length(emu_mets)
    emu_atoms(emu_i) = max(double(extractAfter(atom_names(startsWith(atom_names, emu_mets(emu_i)+":")),"#")));
end

measured_emus = emu_mets;
measured_emu_mids = string(zeros(sum(emu_atoms+1),1));

% as MIDs for the EMUs, take all, from M0 to Mn (n = max atom number of the EMU)
current_mid = 0;
for emu_i = 1:length(emu_mets)
    measured_emus(emu_i) = emu_mets(emu_i)+":"+strjoin("C#"+string(1:emu_atoms(emu_i)),",");
    measured_emu_mids(current_mid+(1:(emu_atoms(emu_i)+1))) = ...
        measured_emus(emu_i)+"."+string(0:emu_atoms(emu_i))';
    current_mid = current_mid + emu_atoms(emu_i) + 1;
end
%measured_emu_mids = measured_emu_mids(endsWith(measured_emu_mids,".0") | endsWith(measured_emu_mids,".1") | ...
%    endsWith(measured_emu_mids,".2") | endsWith(measured_emu_mids,".3"));

% For 13C labeling, CO2 is the only source ...
co2_feed_rxn_id = find(string(model.rxns) == 'Im_CO2');
label_feed_rxn_ids = [ co2_feed_rxn_id];


[corrected_atom_map_mapping, corrected_atom_map_mapping_with_ids] = create_corrected_atom_map_mapping(model, atom_map_rxns, atom_map_mapping);
[needed_emus, needed_emu_map_mappings, needed_emu_map_rxns] = create_needed_emus(atom_map_rxns, corrected_atom_map_mapping, ...
                                                                            corrected_atom_map_mapping_with_ids, measured_emus, model);
[~, ~, ~, emu_mids] = create_emu_data(model, needed_emus);

% TODO
co2_import_id = find(emu_mids == 'CO2[c]:C#1.1');
label_feed_emus = [ co2_import_id ];
% END TODO

% some precalcuations
% we need the "rational" nullspace quite often
rnullspace = null(model.S,'r');

% metabolite concentrations for some special species
mets_no_cmp = extractBefore(string(model.mets),"[");
all_species_id = unique(mets_no_cmp);
N_species_lit_conc_table = readtable(fullfile('models','N_relevant_species_no_cmp_literature_concentration.txt'), 'ReadVariableNames', false, 'Delimiter', '\t');
N_species_id = string(N_species_lit_conc_table.Var1);
[~,special_species_loc] = ismember(N_species_id,all_species_id);
special_species_conc = N_species_lit_conc_table.Var2 * 0.011; %lit values are nmol/gFW

% the values, which are actually stored and reused
save(output_model_file, "atom_map_rxns", "atom_map_mapping", "atom_names", "model", ...
    "rnullspace", "measured_emu_mids", "measured_emus", "label_feed_rxn_ids", "label_feed_emus", ...
    "special_species_loc", "special_species_conc", "needed_emus", "corrected_atom_map_mapping_with_ids", ...
    "corrected_atom_map_mapping", "needed_emu_map_mappings", "needed_emu_map_rxns", "emu_mids", "-v7.3");


end
