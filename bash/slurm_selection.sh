#!/bin/sh
#SBATCH --job-name=selection        # Job name
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=evan.irving-pease@arch.ox.ac.uk # Where to send mail
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=1gb           # Memory per processor
#SBATCH --time=10:00:00             # Time limit hrs:min:sec
#SBATCH --output=array_%A-%a.out    # Standard output and error log

# TODO remove when done testing
#SLURM_ARRAY_TASK_ID=1

# fixed params for running the models
NUMBER_CPUS=16          # how many CPUs are available on each node

pwd; hostname; date

# we need the GSL library to run selection
#module load gsl

# magic setting that lets us read the file properly
IFS=$'\r\n'

# get a list of all the input files
models=( $(<$1) )

# calculate the starting and ending values for this task based on the SLURM task and the number of runs per task
start_num=$(( ($SLURM_ARRAY_TASK_ID - 1) * $NUMBER_CPUS + 1 ))
end_num=$(( $SLURM_ARRAY_TASK_ID * $NUMBER_CPUS ))

# don't overshoot the list of models
if (( $end_num > ${#models[@]} )); then
    end_num=${#models[@]}
fi

# Print the task and run range
echo "This is SLURM_ARRAY_TASK $SLURM_ARRAY_TASK_ID, which will do models ${start_num}-${end_num}."

# run the inputs for this task
for (( i=$start_num-1; i<end_num; i++ )); do
    # get the command
    cmd=${models[$i]}

    # run it in the background
    ( eval ${cmd}; ) &
done

# wait for all the background jobs to finish, or else the script ends and SLURM terminates the jobs
wait

date