#!/bin/sh
# parameter: $1 config file to source, $2 gurobi base path, $3 cobra path
if [ -z "$3" ]
then
echo 'Use with $1 config file to source, $2 gurobi base path, $3 cobra path'
echo "e.g. sh start_workflow_jobs.sh configMinTest.sh ~/apps/gurobi903 ~/apps/cobratoolbox"
exit 1
fi

source $1
res=$(sbatch workflowSample.job $config_matlab $2 $3)
job1_id=$(echo $res | sed 's/.*job //')

wait2_for=""
for slot_i in $(seq 1 $slots_per_job $num_slots)
do
    res=$(sbatch -d afterok:${job1_id} workflowSimulate.job $config_matlab $slot_i $(($slot_i + $slots_per_job - 1)) )
    job2_id=$(echo $res | sed 's/.*job //')
    wait2_for="${wait2_for}:$job2_id"
done
res=$(sbatch -d afterok${wait2_for} workflowHandleSimulateResults.job $config_matlab)
job3_id=$(echo $res | sed 's/.*job //')
for learn_i in $(seq 1 $num_mids_to_learn)
do
    res=$(sbatch -d afterok:${job3_id} workflowLearnNNs.job $config_matlab $learn_i $learn_i)
done
