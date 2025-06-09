#!/bin/bash
#SBATCH --job-name=git_clone		                    # Job name
#SBATCH --partition=batch		                        # Partition (queue) name
#SBATCH --ntasks=1			                            # Single task job
#SBATCH --cpus-per-task=4		                        # Number of cores per task - match this to the num_threads used by BLAST
#SBATCH --mem=40gb			                            # Total memory for job
#SBATCH --time=2:00:00  		                        # Time limit hrs:min:sec
#SBATCH --output=/scratch/ac05869/%x_%j.out			# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --error=/scratch/ac05869/%x_%j.err			# Location of standard output and error log files (replace cbergman with your myid)
#SBATCH --mail-user=ac05869@uga.edu                    # Where to send mail (replace cbergman with your myid)
#SBATCH --mail-type=END,FAIL                            # Mail events (BEGIN, END, FAIL, ALL)

#clone git repository in home directory
WD='/scratch/ac05869'

module load GCC/13.3.0
git clone https://github.com/lh3/wgsim.git
cd wgsim/
gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
cd $WD
~/wgsim/wgsim #calls the program

