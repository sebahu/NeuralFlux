function [] = unite_and_scale_fds(model_file, input_dir, prefixes, num_samples_per_prefix, outputname, max_input)
    load(model_file,"atom_map_mapping", "atom_map_rxns", "atom_names", "model");

    input_reactions = [438,439];
    total_samples = zeros(length(model.rxns)+1,length(prefixes)*(num_samples_per_prefix));
    for prefix = 1:length(prefixes)
        load(fullfile(input_dir,strcat('samples_totalfluxlimited_',prefixes(prefix),'.mat')), 'samples');
        total_samples(:,(prefix-1)*num_samples_per_prefix+(1:num_samples_per_prefix)) = samples;
    end
    if ~exist('max_input', 'var')
        max_input = max(sum(total_samples(input_reactions,:)));
    end

    input_sums = sum(total_samples(input_reactions, :),1);
    scale = (max_input./input_sums);
    total_samples_scaled = total_samples.*scale;
    save(outputname, 'total_samples_scaled','-v7.3');

end
