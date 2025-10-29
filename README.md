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
This folder contains scripts used in the proof-of-concept evaluation of NeuralFlux

### The proof-of-concept
The proof-of-concepts first sees the creation of the neural netwoprk sets for
AraCore v2.1 for carbon and nitrogen mapping, using the same sampled flux distributions
and metabolite concentration distributions. As a result then exist:
- 600k sampled flux and compartmentalized metabolite concentration distributions
- simulated 13C and 15N labeling enrichment for the samples
- Neural networks to predict the 13C and 15N enrichments

As the next step, additional test flux and compartmentalized distributions are created, with slightly different settings
then for the training data, to guarantee total independence of the two sets of samples.
For those test samples, 13C and 15N enrichment is simulated.

In the flux estimation step, the starting values are derived from the 600k sampled flux distributions which have the closest enrichment
values to the "measured" enrichment values.
To facilitate an efficient evaluation, this calculation is done for all sampled test flux distributions and the indices of the 1000 closest
orignal samples are stored.

The flux estimation for the whole flux distribution is done for 100 pseudo-randomly selected
test flux distributions. The results are compared to calculate the correlation
between estimated and true value for all reactions.

From the reactions with high correlation, a representaive selection of reactions, which are not pairwise fully coupled
to any of the other selected reactions, is used for the confidence interval calculation.

From the 100 pseudo-randomly selected test samples, 5 are selected, that have low pairwise correlation with each other.
For these reaction/test sample pairs, the confidence intervals are calculated.



### The scripts

- evalEstimatesForTestSamplesWithConstrainedMetConcs2.m: calculates flux estimations for all reactions (i.e. a whole flux distribution)
for a set of test samples, selected pseudo-randomly from all test samples. It combines 13C and 15N labeling data.

    - collect all the neural networks
    - get starting values as a combination of the closest matches from the original 600k samples.
    - create function which calls the neural networks to get the enrichment values for the given parameters 
      (i.e. the current estimate for the flux values and compartmentalized metabolite concentrations) and implements barrier
      functions to keep the metabolite concentrations in the limits given by the measured absolute concentrations
    - call lsqnonlin with this function to perform the parameter optimization
    

- evalConfidenceIntervalsWithConstrainedMetConcs2.m: calculates the confidence intervals for one reaction of one test sample.
  the reaction is encoded by its rank in the correlation analysis of the previous step, and the test sample is identified
  by its position in the pseudo-random sorted test samples


collectData2.m
collectData.m
evalConfidenceIntervalsWithConstrainedMetConcs.m
evalEstimatesForTestSamplesWithConstrainedMetConcs2.m
evalEstimatesForTestSamplesWithConstrainedMetConcs.m
evalEstimatesForTestSamplesWithKnownCompMetConcs2.m
evalEstimatesForTestSamplesWithKnownCompMetConcsC.m
evalEstimatesForTestSamplesWithKnownCompMetConcs.m
evalEstimatesForTestSamplesWithKnownCompMetConcsN.m
evalMetConfidenceIntervalsWithConstrainedMetConcs2.m
evalMetConfidenceIntervalsWithConstrainedMetConcs.m
prepareEstimatesForTestSamples2.m
prepareEstimatesForTestSamples.m





## figures




## jobs
Slurm jobs that facilitate the parallization of the computational expensive
steps of the workflow on HPC nodes.

- start_workflow_jobs.sh, start_workflow_jobs_after_sampling.sh: shell scripts, that launch
  slurm jobs for the workflow to create the samples, simulate them and train neural networks for
  one experimental setup. The few parallelization parameters are stored in additional
  shell scripts (config...sh, ih the same folder), that are sourced

- configAraCoreC.sh, configAraCoreN.sh, configEvalTestSamplesC.sh, configMinTestC.sh, configMinTest.sh:
  config shell scripts for the accordingly named matlab configuration scripts from the config folder

- workflowSample.job, workflowSimulate.job, workflowHandleSimulateResults.job, workflowLearnNNs.job:
  slurm jobs to execute one of the steps of the NeuralFlux workflow (matlab scripts of the main folder)
  with the same parameters regarding the configuration and additional parallization paramters

- evalConfidenceIntervalsWithConstrainedMetConcsCNpart3.job:
  slurm job to execute a confidence interval calculation for a preselection of reactions, for one
  flux distribution of the test samples. Used for the proof-of-concept. The reactions were selected
  after a best flux estimation for 100 test flux distributions and then a correlation analysis

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



