function output = forwardSigmoidnet2(NNs_per_MID_and_tp,mu_per_MID_and_tp,sigma_per_MID_and_tp,X_in)
num_mids = size(mu_per_MID_and_tp,1);
num_tps = size(mu_per_MID_and_tp,2);
output = zeros(num_mids*num_tps,1);
for mid_i = 1:num_mids
    for tp_i = 1:num_tps
        % 1. Initialize first layer input - no need to save intermediate results.
        net = NNs_per_MID_and_tp{mid_i,tp_i};
        mu = mu_per_MID_and_tp{mid_i,tp_i};
        sigma = sigma_per_MID_and_tp{mid_i,tp_i};
        W = net.LayerWeights;
        B = net.LayerBiases;
        n = numel(W);

        X = X_in - mu;
        X= X./sigma;
        X(mu==0 & sigma==0)=0;

        Xi = double(X');
        % 2. Do forward pass.
        for i = 1:n
            Wi = double(W{i});
            Bi = double(B{i});
            Xi = Wi*Xi + Bi;
            if i<n
                Xi = 1./(1+exp(-Xi));
            end
        end
        output((mid_i-1)*num_tps + tp_i) = Xi';
    end
end

end
