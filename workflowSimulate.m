function [error] = workflowSimulate(config_script, start_slot, end_slot)
    addpath("application_core");
    addpath("configs");
    error = 0;
    if ~exist('config_script', 'var')
        disp("call with appropriate config, start and end_slot eg. workflowSimulate('configMinTest',1,3)");
        error = -1;
        return;
    end

    cfg = eval(config_script);
    if start_slot < 1 || start_slot > end_slot || end_slot > cfg.num_slots
        disp("Invalid arguments for start/end_slot: "+string(start_slot)+"/"+string(end_slot));
        disp(" Slots range from 1 to "+string(cfg.num_slots)+ ", start_slot must be <= end_slot.");
        return;
    end

    simulateSamples(cfg.model_file, cfg.split_dir, cfg.slot_size, start_slot, end_slot, cfg.logsPerHour, cfg.simDurationHours, cfg.atomLabelFraction);
end
