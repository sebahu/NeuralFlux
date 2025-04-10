function [] = learnNNsForOneMIDFromNspLong(mid_name_input, model_file, scaled_fd_samples_file, ...
    sample_decomp_ratios_file, met_sample_file, selected_timepoints, NN_dir, used_nsfs_is, used_relevant_met_is)
% Base idea: learn NNs for every (relevant) timepoint of every measured
% metabolite/MID: simple approach: predefined time points
% brute force: all fluxes and met conc are input

    disp(mid_name_input);
    if contains(mid_name_input,";")
        mid_name_array = split(string(mid_name_input),";");
    else
        mid_name_array = string(mid_name_input);        
    end

    load(model_file, "measured_emu_mids", "model", "needed_emus", "rnullspace");
    [emu_mid_met_inx, ~, ~, ~] = create_emu_data(model, needed_emus);
    measured_decomp_emu_mids = unique(extractBefore(measured_emu_mids,"[")+extractAfter(measured_emu_mids,"]"));

    needed_mids_ix = zeros(size(mid_name_array));
    for  mid_i=1:length(mid_name_array)
        needed_mids_ix(mid_i) = find(measured_decomp_emu_mids == mid_name_array(mid_i));
    end

    load(sample_decomp_ratios_file,"sample_decomp_ratios");
    needed_sample_decomp_ratios = sample_decomp_ratios(needed_mids_ix,selected_timepoints,:);
    clear sample_decomp_ratios;
    % can be helpfull for debugging
    %save(string(NN_dir)+"/data_needed_for_NN_"+mid_name_input+".mat");

    num_samples = size(needed_sample_decomp_ratios,3);
    used_samples = round(0.833*num_samples);
    
    load(scaled_fd_samples_file, 'total_samples_scaled');
    load(met_sample_file, 'sample_met_pools');
    
    
    relevant_met_is = unique(emu_mid_met_inx);
    rnsp_factor_samples = (total_samples_scaled(1:end-1,:)'/rnullspace')';
    if ~isempty(used_nsfs_is)
        rnsp_factor_samples = rnsp_factor_samples(used_nsfs_is, :);
    end
    if ~isempty(used_relevant_met_is)
        relevant_met_is = relevant_met_is(used_relevant_met_is);
    end
    
    minX = min([rnsp_factor_samples',sample_met_pools(relevant_met_is,:)']);
    rangeX = max([rnsp_factor_samples',sample_met_pools(relevant_met_is,:)']) - minX;
    rangeX(rangeX<=0)=1;
        
    layerSizes = min([216 36 6], length(minX));
        
    for  mid_i=1:length(mid_name_array)
        mid_idx = find(measured_decomp_emu_mids == mid_name_array(mid_i));
        disp("Wanted MID is #"+mid_idx);

        rand_indices=randperm(num_samples);

        X = single([rnsp_factor_samples(:,rand_indices(1:used_samples))',sample_met_pools(relevant_met_is,rand_indices(1:used_samples))']);

        scaledX = (X-minX)./rangeX;
        testX = (single([rnsp_factor_samples(:,rand_indices(used_samples+1:end))',sample_met_pools(relevant_met_is,rand_indices(used_samples+1:end))'])-minX)./rangeX;
        for tp_i = 1:length(selected_timepoints)
            current_tp = selected_timepoints(tp_i);
            %Y = single(squeeze(sample_decomp_ratios(mid_idx,current_tp,rand_indices(1:used_samples))));
            Y = single(squeeze(needed_sample_decomp_ratios(mid_i,tp_i,rand_indices(1:used_samples))));
            scaledY = Y; %/100;
            disp("size X:");
            disp(size(scaledX));
            disp(" size Y:");
            disp(size(scaledY));
            rng("default") % For reproducibility

            filter_i = 1:size(scaledX,2); %[I(1:10),filter_extra(B(11:20)>0.2)]
            filteredX = scaledX(:,filter_i);
            
            mu = single(mean(filteredX));
            sigma = sqrt(sum(((single(filteredX)-mu).^2)*(1/size(filteredX,1)))/(1-1/size(filteredX,1)));
            % TODO Layersize could be calculated or parameterized
            Mdl = fitrnet(filteredX, scaledY, 'LayerSizes', layerSizes, 'Activations','sigmoid', 'Verbose', 1, ...
                'Standardize', true, 'IterationLimit', 2000);
            compactMdl = compact(Mdl);
            %refY = single(squeeze(sample_decomp_ratios(mid_idx,current_tp,rand_indices(used_samples+1:end))));
            refY = single(squeeze(needed_sample_decomp_ratios(mid_i,tp_i,rand_indices(used_samples+1:end))));
            predY = predict(Mdl,testX(:,filter_i));
            test_result = sum(abs(refY-predY)<=0.02)/(num_samples-used_samples);
            save(string(NN_dir)+"/NN6_for_"+mid_name_array(mid_i)+"_tp_+"+string(current_tp)+"_withMets_3L.mat",'compactMdl',...
                'filter_i','rand_indices','mu','sigma','test_result','minX','rangeX','-v7.3');
            % 2-Layer NN seems to be inferior
%             Mdl = fitrnet(filteredX, scaledY, 'LayerSizes',[80 9], 'Activations','sigmoid', 'Verbose', 1, ...
%                 'Standardize', true, 'IterationLimit', 2000);
%             compactMdl = compact(Mdl);
%             predY = predict(Mdl,testX(:,filter_i));
%             test_result = sum(abs(refY-predY)<=0.02)/(num_samples-used_samples);
%             save(string(NN_dir)+"/NN6_for_"+mid_name_array(mid_i)+"_tp_+"+string(current_tp)+"_withMets_2L.mat",'compactMdl',...
%                 'filter_i','rand_indices','mu','sigma','test_result','minX','rangeX','-v7.3');
        end
    end
end
