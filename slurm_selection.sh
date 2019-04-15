#!/bin/sh
#SBATCH --job-name=selection        # Job name
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=evan.irving-pease@arch.ox.ac.uk # Where to send mail
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=1gb           # Memory per processor
#SBATCH --time=01:30:00             # Time limit hrs:min:sec (each cycle takes ~7e-5 seconds, so 5e7 = 3,500 sec = ~ 1hr)
#SBATCH --output=array_%A-%a.out    # Standard output and error log

# fixed params for running the models
POP_HISTORY='horse-DOM2-const.pop'
MCMC_CYCLES=50000000    # number of MCMC cycles to run
SAMPLE_FREQ=10000       # frequency of sampling from the posterior
OUTPUT_FREQ=1000        # frequency of printing output to the screen
ALLELE_FREQ=20          # fraction of the allele frequency to update during a trajectory update move
NUMBER_CPUS=16          # how many CPUs are available on each node
NUMBER_REPS=4           # number of independent replicates to run

pwd; hostname; date

# we need the GSL library to run selection
module load gsl

# get a list of all the input files
models=($(ls selection/*.input))

# set the number of runs that each SLURM task should do
per_task=$NUMBER_CPUS/$NUMBER_REPS  # arcus-b nodes have 16 CPUs, but we run 4 independent replicates of each model

# calculate the starting and ending values for this task based on the SLURM task and the number of runs per task
start_num=$(( ($SLURM_ARRAY_TASK_ID - 1) * $per_task + 1 ))
end_num=$(( $SLURM_ARRAY_TASK_ID * $per_task ))

# don't overshoot the list of models
if (( $end_num > ${#models[@]} )); then
    end_num=${#models[@]}
fi

# Print the task and run range
echo "This is SLURM_ARRAY_TASK $SLURM_ARRAY_TASK_ID, which will do models $start_num to $end_num with $NUMBER_REPS replicates each.\n"

# run the inputs for this task
for (( i=$start_num; i<=end_num; i++ )); do

    input=${models[$i-1]}
    if test "$input+isset}"; then

        # run the replicates for this input
        for (( n=1; n<=NUMBER_REPS; n++ )); do

            # compose the command to run
            cmd="sr -D $input -o data/${input%.*}-n${n} -P $POP_HISTORY -a -n $MCMC_CYCLES -s $SAMPLE_FREQ -f $OUTPUT_FREQ -F $ALLELE_FREQ -e $RANDOM"

            # run selection, and log the command
            echo $cmd | tee ${input%.*}.log
            eval $cmd >> ${input%.*}.log &
        done
    fi
done

# wait for all the background jobs to finish, or else the script ends and SLURM terminates the jobs
wait

date