function[emu_mid_met_inx, emu_met_inx, emu_mid_emu_is, emu_mids] = create_emu_data(model, needed_emus)
emu_mets = extractBefore(needed_emus, ":");
emu_met_inx = zeros(size(emu_mets));
for emu_i = 1:length(model.mets)
    met = string(model.mets(emu_i));
    if ismember(met, emu_mets)
        emu_met_inx(emu_mets == met) = emu_i;
    end
end

emu_mid_emu_is = [];
emu_mids = string([]);
for emu_i = 1:length(needed_emus)
    emu_atoms = split(extractAfter(needed_emus(emu_i),":"),",");
    for m_i = 0:length(emu_atoms)
        emu_mids(end+1)=needed_emus(emu_i)+"."+string(m_i);
        emu_mid_emu_is(end+1)=emu_i;
    end    
end
emu_mid_met_inx = emu_met_inx(emu_mid_emu_is);