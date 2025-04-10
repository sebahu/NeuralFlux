function [] = workflowUniversal(config_script, gurobi_path, cobra_path)
    addpath("application_core");
    addpath("configs");
    if ~exist('cobra_path', 'var')
        disp("call with appropriate config and paths, eg. workflowUniversal('configMinTest','~/apps/gurobi903/linux64/matlab', '~/apps/cobratoolbox')");
        return;
    end

    cfg = eval(config_script);
    workflowSample(config_script, gurobi_path, cobra_path);
    workflowSimulate(config_script, 1, cfg.num_slots);
    workflowHandleSimulateResults(config_script);
    workflowLearnNNs(config_script);
end
