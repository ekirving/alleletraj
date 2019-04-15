#!/bin/sh
#SBATCH --job-name=selection        # Job name
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=evan.irving-pease@arch.ox.ac.uk # Where to send mail
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=1gb           # Memory per processor
#SBATCH --time=01:30:00             # Time limit hrs:min:sec (each cycle takes ~7e-5 seconds, so 5e7 = 3,500 sec = ~ 1hr)
#SBATCH --output=array_%A-%a.out    # Standard output and error log

# This is an example script that combines array tasks with
# bash loops to process many short runs. Array jobs are convenient
# for running lots of tasks, but if each task is short, they
# quickly become inefficient, taking more time to schedule than
# they spend doing any work and bogging down the scheduler for
# all users.

pwd; hostname; date

# we need the GSL library to run selection
module load gsl

# get a list of all the input files
models=($(ls selection/*.input))

# set the number of runs that each SLURM task should do
PER_TASK=16  # arcus-b nodes have 16 CPUs

# calculate the starting and ending values for this task based on the SLURM task and the number of runs per task
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# don't overshoot the list of models
if (( $END_NUM > ${#models[@]} )); then
    END_NUM=${#models[@]}
fi

# Print the task and run range
echo "This is task $SLURM_ARRAY_TASK_ID, which will do models $START_NUM to $END_NUM."

POP_HISTORY='horse-DOM2-const.pop'
MCMC_CYCLES=50000000
SAMPLE_FREQ=10000
PRINTL_FREQ=1000
ALLELE_FREQ=20

# Run the loop of runs for this task.
for (( run=$START_NUM; run<=END_NUM; run++ )); do
    input=${models[$run-1]}
    if test "$input+isset}"; then

        cmd="sr -D $input -o ${input%.*} -P $POP_HISTORY -a -n $MCMC_CYCLES -s $SAMPLE_FREQ -f $PRINTL_FREQ -F $ALLELE_FREQ -e $RANDOM"

        # run selection, and log the command
        echo $cmd | tee ${input%.*}.log
        eval $cmd >> ${input%.*}.log &
    fi
done

# wait for all the background jobs to finish
wait

date