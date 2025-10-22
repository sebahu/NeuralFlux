% default for running the estimation (w/o CI, but with known met
% concentrations for C labeling
config_script_NN = 'configAraCoreC';
config_script_testsamples = 'configEvalTestSamplesC';
rng_start = 12345;
start_i = 1;
end_i = 100;

evalEstimatesForTestSamplesWithKnownCompMetConcs(config_script_NN, config_script_testsamples, rng_start, start_i, end_i)
