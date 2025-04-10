function [] = sampleModel(model_file, model_biomass_rxn_id, minBioRate, number_of_samples, samples_dir, path_to_gurobi, path_to_cobra)
% this needs to be called for all planned minBioRate values and their
% corresponding number of samples
% the standard while developing was .3, .5, .7, .8, .9, .95, and 100000
% samples each. A finer graining might be good    addpath(fullfile(path_to_gurobi,'linux64', 'matlab'));

    addpath(path_to_cobra);
    addpath(path_to_gurobi);
    initCobraToolbox();
    changeCobraSolver('gurobi', 'all');

    load(model_file, 'model');
    model_bio1_rxn_index = find(string(model.rxns) == model_biomass_rxn_id);
    
    % add synchronization metabolite
    model.b(end+1)=0;
    model.csense(end+1) = 'E';
    model.metCharges(end +1)=0;
    model.metFormulas(end+1)={'X'};
    model.metisinchikeyID(end+1)={'X'};
    model.metNames(end+1)={'dummy met to sync fluxes'};
    model.mets(end+1)={'flux_sync'};
    model.metSEEDID(end+1)={'X'};

    % synchronizing fluxes, to be able to limit total sum of fluxes
    model.S(end+1,:)=1;
    model.A(end+1,:)=1;
    model.S(:,end+1)=0;
    model.A(:,end+1)=0;
    model.S(end,end)=-1;
    model.A(end,end)=-1;

    model.c(end+1) = 0;
    model.grRules(end+1) = cellstr("");
    model.lb(end+1) = 0;
    model.ub(end+1) = 1000000;
    model.rules(end+1) = cellstr("");
    model.rxnConfidenceScores(end+1) = 0;
    model.rxnECNumbers(end+1) = cellstr("");
    model.rxnNames(end+1) = cellstr("Sum of all fluxes");
    model.rxnNotes(end+1) = cellstr("");
    model.rxns(end+1) = cellstr("sum_fluxes");


    model.c=zeros(size(model.rxns));
    model.c(model_bio1_rxn_index)=1;
    [minFlux, maxFlux] = fluxVariability(model,99);
    totalMinFlux99 = minFlux(end);
    opt_bio1 = maxFlux(model_bio1_rxn_index);

    [minFlux, maxFlux] = fluxVariability(model,minBioRate);
    model.lb=minFlux;
    model.ub=maxFlux;
    model.ub(end)=2*totalMinFlux99;

    % limits for exports

    %model.ub(443:522)=maxFlux(443:522)/100;

    [P,model_post] = chrrParseModel(model);
    model_post.c(:)=1;
    [samples, roundedPolytope, minSampledFlux, maxSampledFlux] = chrrExpSampler(model_post, 731, number_of_samples);
    save(strcat(samples_dir,'/samples_totalfluxlimited_',string(minBioRate*100),".mat"), "samples", '-v7.3')
