function [standardized_nsfs, standardized_met_concs, resnorm, residual] = ...
    fit_NN_to_NSFs_and_known_mets_Chi2(fitOptions, sample_MIDs_to_estimate, sds_per_MID, nsf_base, met_base, ...
                                NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp)

    num_of_nsfs = size(nsf_base,1);
    num_of_metfs = size(met_base,1);
    
    test_Y_flat = reshape(sample_MIDs_to_estimate',[],1);
    sds_per_Y_flat = reshape(sds_per_MID',[],1);
    
    degrees_of_freedom = max(length(test_Y_flat)/2, length(test_Y_flat) - (num_of_metfs + num_of_nsfs));

    fun = @(x)(double(forwardSigmoidnet2(NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp, x)-test_Y_flat)./sds_per_Y_flat);
    num_vars = num_of_nsfs + num_of_metfs; % length(mu_per_MID_and_tp{1,1});
    lb = zeros(num_vars,1);
    lb(num_of_nsfs+1:end)=0.9*met_base;
    ub = ones(num_vars,1);
    ub(num_of_nsfs+1:end)=min(1.1*met_base+0.001,1);
    x0=double([mean(nsf_base,2); mean(met_base,2)])';

    maxTotalIterations = 10*fitOptions.MaxIterations;
    fitOptions.MaxIterations = 3;
    totalIterations = 0;
    
    while totalIterations < maxTotalIterations
        t = datetime; 
        [x,resnorm,residual,exitflag,output] = lsqnonlin(fun,x0,lb,ub,fitOptions);
        tdiff_ms = milliseconds(datetime-t); disp(tdiff_ms+" ms");
        if resnorm/degrees_of_freedom < 1
            break
        end
        x0 = x;
        totalIterations = totalIterations + 4;
    end
    standardized_nsfs = x(1:num_of_nsfs);
    standardized_met_concs = x(num_of_nsfs+1:end);     
end
