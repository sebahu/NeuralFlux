function [ci_lb, ci_ub, x_lb, x_ub] = ...
    calc_comp_met_confidence_interval_constrained_mets(fitOptions, sample_MIDs_to_estimate, sds_per_MID, nsf_base_norm, met_base_norm, ...
                                NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp, ...
                            met_decomp_matrix, met_decomp_min, measured_decomp_met_conc, ...
                            relative_constraint_factor, met_i, Chi2_lim, x_selectors_per_mid)

    num_of_nsfs = size(nsf_base_norm,1);
    num_of_metfs = size(met_base_norm,1);
        
    test_Y_flat = reshape(sample_MIDs_to_estimate',[],1);
    sds_per_Y_flat = reshape(sds_per_MID',[],1);
    
    degrees_of_freedom = length(test_Y_flat) - (num_of_metfs + num_of_nsfs);

    %barrierFun = @(x,xmin,xmax)(0.003*(log(0.25*((xmax-xmin).^2))-log(max(0,(x-xmin).*(xmax-x)))));
    barrierFun = @(x,min,max)(1E+3./(1+exp(1./(((x-((min+max)/2))./(max-min)).^2))));
    barrierMetUp = @(x,ub)(max(0,(x-ub)+1E+3./(1+exp(1E+6*(ub-x)))));
    barrierMetDown = @(x,lb)(max(0,(lb-x)+1E+3./(1+exp(1E+6*(x-lb)))));
    fun_up = @(x, met_lb)([double(forwardSigmoidnet2v2(NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp, x, x_selectors_per_mid)-test_Y_flat)./sds_per_Y_flat;
                barrierFun((x(num_of_nsfs+1:end)*met_decomp_matrix)'+met_decomp_min, ...
                        measured_decomp_met_conc./relative_constraint_factor, ...
                        relative_constraint_factor.*measured_decomp_met_conc);
                barrierMetDown(x(num_of_nsfs+met_i), met_lb)]);
    fun_down = @(x, met_ub)([double(forwardSigmoidnet2v2(NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp, x, x_selectors_per_mid)-test_Y_flat)./sds_per_Y_flat;
                barrierFun((x(num_of_nsfs+1:end)*met_decomp_matrix)'+met_decomp_min, ...
                        measured_decomp_met_conc./relative_constraint_factor, ...
                        relative_constraint_factor.*measured_decomp_met_conc);
                barrierMetUp(x(num_of_nsfs+met_i), met_ub)]);
    num_vars = num_of_nsfs + num_of_metfs; % length(mu_per_MID_and_tp{1,1});
    lb = zeros(num_vars,1);
    ub = ones(num_vars,1);
    
    maxTotalIterations = 48;
    fitOptions.MaxIterations = 7;    
    x0=double([nsf_base_norm; met_base_norm; ])';
    met_opt = x0(num_of_nsfs+met_i);
    %fitOptions.FunctionTolerance =  (1+Chi2_lim)*degrees_of_freedom;
    
    fitOptionsStart = fitOptions;
    fitOptionsStart.MaxIterations = 100;         
    
        
    % finding lower interval border
    ci_lb_test0 = 0;
    ci_lb = met_opt;
    ci_lb_test1 = ci_lb_test0;
    
    ci_lb_found = false;
    x0=double([nsf_base_norm; met_base_norm; ])';
    x_lb = x0;
    while ~ci_lb_found        
        disp([ci_lb, ci_lb_test1, ci_lb_test0]);
        totalIterations = 0;
        met_ub = ci_lb_test1;
        x0=double([nsf_base_norm; met_base_norm; ])';
        fun_start = @(x)([x'-x0'; barrierMetUp(x(num_of_nsfs+met_i), met_ub);
                            barrierMetDown(x(num_of_nsfs+met_i), 0.9*met_ub)]);
        [x,resnorm,residual,exitflag,output] = lsqnonlin(fun_start,x0,lb,ub,fitOptionsStart);
        x0=x;
        fun = @(x)(fun_down(x, met_ub));
        while totalIterations < maxTotalIterations
            t = datetime; 
            [x,resnorm,residual,exitflag,output,lambda, jacobian] = lsqnonlin(fun,x0,lb,ub,fitOptions);
            tdiff_ms = milliseconds(datetime-t); disp(tdiff_ms+" ms");
            if resnorm/degrees_of_freedom < 1+Chi2_lim
                ci_lb = ci_lb_test1;
                x_lb = x;
                ci_lb_test1 = (ci_lb_test1 + ci_lb_test0)/2;
                break
            end
            x0 = x;
            totalIterations = totalIterations + 8;
        end
        if totalIterations >= maxTotalIterations
            ci_lb_test0 = ci_lb_test1;
            ci_lb_test1 = (ci_lb + ci_lb_test0)/2;
        end
        if ci_lb/ci_lb_test0 < 1.1
            ci_lb_found = true;
            ci_lb = ci_lb_test0;
        end
    end
    
    
    %finding the upper border
    ci_ub_test0 = 1;
    ci_ub = met_opt;
    ci_ub_test1 = ci_ub_test0;
    
    ci_ub_found = false;
    x0=double([nsf_base_norm; met_base_norm; ])';
    x_ub = x0;
    while ~ci_ub_found
        disp([ci_ub, ci_ub_test1, ci_ub_test0]);
        totalIterations = 0;
        met_lb = ci_ub_test1;
        x0=double([nsf_base_norm; met_base_norm; ])';
        fun_start = @(x)([x'-x0'; barrierMetDown(x(num_of_nsfs+met_i), met_lb);
                            barrierMetUp(x(num_of_nsfs+met_i), 1.1*met_lb)]);
        [x,resnorm,residual,exitflag,output] = lsqnonlin(fun_start,x0,lb,ub,fitOptionsStart);
        x0=x;
        fun = @(x)(fun_up(x, met_lb));
        while totalIterations < maxTotalIterations
            t = datetime;
            [x,resnorm,residual,exitflag,output, lambda, jacobian] = lsqnonlin(fun,x0,lb,ub,fitOptions);
            tdiff_ms = milliseconds(datetime-t); disp(tdiff_ms+" ms");
            if resnorm/degrees_of_freedom < 1+Chi2_lim
                ci_ub = ci_ub_test1;
                x_ub = x;
                ci_ub_test1 = (ci_ub + ci_ub_test0)/2;
                break
            end
            x0 = x;
            totalIterations = totalIterations + 8;
        end
        if totalIterations >= maxTotalIterations
            ci_ub_test0 = ci_ub_test1;
            ci_ub_test1 = (ci_ub + ci_ub_test0)/2;
        end
        if ci_ub_test0/ci_ub < 1.1
            ci_ub_found = true;
            ci_ub = ci_ub_test0;
        end
    end
   
    

    
    
%end
