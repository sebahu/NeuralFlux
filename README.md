# NeuralFlux
Neural network based approach to INST-MFA
NeuralFlux was developed with COBRA toolbox 3.0, Gurobi 9.03 and Matlab 2021a

#Overview of files
./ (this folder): matlab files to execute (parts of) the workflow to sample fluxdistributions, simulate labeling enrichment for them, and train neural networks with the resulting data.

workflowUniversal.m: Matlab script, executing the whole workflow. Parameters: config_script, gurobi_path, cobra_path
- config_script: name of a configuration script in configs, e.g. configMinTestC or configAraCoreC
- gurobi_path: path to matlab frontend of gurobi installation, e.g. '~/apps/gurobi903/linux64/matlab'
- cobra_path: path to cobra toolbox installation, e.g. '~/apps/cobratoolbox'
When workflowUniversal is executed, no paralleization takes place. This is only usable as a test run for the pipeline on a minimal config

workflowSample.m; workflowSimulate.m; workflowHandleSimulateResults.m; workflowLearnNNs.m: the individual steps of the workflow, using the same config_script as the overall workflow

Individual parameters:
workflowSample.m: 
- gurobi_path, cobra_path: as in workflowUniversal.m
- max_input: optional, controls the normalization of the samples. If no value is given, the global maximum of the flux value(s) which form the base of the normalization are used. Otherwise, the fluxes are normalized for the configured normalization base to have this value.
workflowHandleSimulateResults.m: no further parameters
workflowSimulate: start_slot, end_slot - select only a subset of sampels to simulate, for parallelization
workflowLearnNNs: start_mid_index, end_mid_index - select only a subset of NNs to learn, for parallelization

