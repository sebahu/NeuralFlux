# NeuralFlux
Neural network based approach to INST-MFA
NeuralFlux was developed with COBRA toolbox 3.0, Gurobi 9.03 and Matlab 2021a

# Overview of files
## ./ (this folder)
matlab files to execute (parts of) the workflow to sample fluxdistributions, simulate labeling enrichment for them, and train neural networks with the resulting data.

workflowUniversal.m: Matlab script, executing the whole workflow. Parameters: config_script, gurobi_path, cobra_path
- config_script: name of a configuration script in configs, e.g. configMinTestC or configAraCoreC
- gurobi_path: path to matlab frontend of gurobi installation, e.g. '~/apps/gurobi903/linux64/matlab'
- cobra_path: path to cobra toolbox installation, e.g. '~/apps/cobratoolbox'
When workflowUniversal is executed, no paralleization takes place. This is only usable as a test run for the pipeline on a minimal config

A single command, to test the pipeline up to the neural network training:
`matlab -batch "workflowUniversal('configMinTestC', '~/apps/gurobi903/linux64/matlab',  '~/apps/cobratoolbox')"`

workflowSample.m; workflowSimulate.m; workflowHandleSimulateResults.m; workflowLearnNNs.m: the individual steps of the workflow, using the same config_script as the overall workflow

Individual parameters:
workflowSample.m: 
- gurobi_path, cobra_path: as in workflowUniversal.m
- max_input: optional, controls the normalization of the samples. If no value is given, the global maximum of the flux value(s) which form the base of the normalization are used. Otherwise, the fluxes are normalized for the configured normalization base to have this value.
workflowHandleSimulateResults.m: no further parameters
workflowSimulate: start_slot, end_slot - select only a subset of sampels to simulate, for parallelization
workflowLearnNNs: start_mid_index, end_mid_index - select only a subset of NNs to learn, for parallelization

## application_core
The matlab implementations of the preprocessing. sampling, simulation and neural network training

The matlab implementations of the flux estimation and confidence interval calculations

The usage of the functionality of these implementations is used in the workflow implementations of this folder and in the evaluation implementation, which tests the
estimate and confidence interval calculation for a arge set of sampled test flux distributions

## configs
contains matlab scripts that contain the configuration for a NeuralFlux setup
and matlab scripts, that provide a metabolic model with all needed information

- configAraCoreC.m: configuration script, with setup information for 13C labeling
  with CO2 as label source, measured metabolites include all amino acids but also
  metabolites of CBC and TCA cycle. All settings are documented in this file.
  
- prepareAraCoreC.m: matlab script, that provides an enhanced instance
  of the AraCore v2.1 model. Main additions include the splitting of reversible
  reactions into forward and reversed reaction, addition of atom transition maps,
  and definition of measured EMU MIDs  
  
The other files are configuration files or provide a model for other scenarios

## evaluation



## figures




## jobs



## models
Data needed for the models used in NeuralFlux model preparation scripts.
So far only AraCore v2.1 and supporting information is included.

all_atoms.C.sorted.txt, all_atoms.N.sorted.txt: all ids of carbon (C) or nitrogen (N) atoms included atom transistion mappings for AraCore 
all_mapping.C.sorted.txt, all_mapping.N.sorted.txt: all atom transistion mappings of carbon (C) or nitrogen (N) atoms for AraCore 
all_mapping.sorted.txt: all atom transistion mappings for AraCore, may contain errors (which are fixed in C and N mapping files)
AraCore_v2_1.mat
N_relevant_species_no_cmp_literature_concentration.txt





## runtime_data
parent folder for all runtime directories as specified in the config scripts.
Only the parent folder is included in git, all content is ignored.



