function [emu_mid_start_values] = calcStartMIDs(normal_emus_num, emu_mids, atom_label_fraction)
    % "6PG:C#1,C#2,C#3,C#4,C#5,C#6.4"
    emu_mid_start_values = zeros(normal_emus_num,1);
    emu_mid_start_values(endsWith(emu_mids(1:normal_emus_num),".0")) = 1;
    emus = unique(extractBefore(emu_mids,"."));
    for emu_i = 1:length(emus)
        zero_label = 1;
        emu = emus(emu_i);
        num_atoms = count(emu,"#");
        for labeled_atoms = 1:num_atoms
            label_fraction = (1-atom_label_fraction)^(num_atoms-labeled_atoms)*...
                (atom_label_fraction^labeled_atoms)*nchoosek(num_atoms,labeled_atoms);
            zero_label = zero_label - label_fraction;
            emu_mid_i = find(emu_mids == emu + "." + string(labeled_atoms));
            if ~isempty(emu_mid_i)
                emu_mid_start_values(emu_mid_i) = label_fraction;
            end            
        end
        emu_mid_i = find(emu_mids == emu + ".0");
        if ~isempty(emu_mid_i)
            emu_mid_start_values(emu_mid_i) = zero_label;
        end            
    end
end