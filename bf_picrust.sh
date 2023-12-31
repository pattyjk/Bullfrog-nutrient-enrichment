#!/bin/bash

# Sample slurm submission script for the Chimera compute cluster
# Lines beginning with # are comments, and will be ignored by
# the interpreter.  Lines beginning with #SBATCH are directives
# to the scheduler.  These in turn can be commented out by
# adding a second # (e.g. ##SBATCH lines will not be processed
# by the scheduler).
#
#
# set name of job
#SBATCH --job-name=bf_picrust
#
# set the number of processors/tasks needed
#SBATCH -n 12

#set an account to use
#if not used then default will be used
# for scavenger users, use this format:
#BATCH --account=patrick.kearns
# for contributing users, use this format:
##SBATCH --account=

# set max wallclock time  DD-HH:MM:SS

# the default time will be 1 hour if not set
#SBATCH --time=00-08:00:00

# set a memory request
#SBATCH --mem=48gb

# Set filenames for stdout and stderr.  %j can be used for the jobid.
# see "filename patterns" section of the sbatch man page for
# additional options
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#

# set the partition where the job will run.  Multiple partitions can
# be specified as a comma separated list
# Use command "sinfo" to get the list of partitions
#SBATCH --partition=Intel6240
##SBATCH --partition=Intel6240,Intel6248,DGXA100

#When submitting to the GPU node, these following three lines are needed:
##SBATCH --gres=gpu:1
##SBATCH --export=NONE
#source /etc/profile
 

#Optional
# mail alert at start, end and/or failure of execution
# see the sbatch man page for other options
#SBATCH --mail-type=ALL
# send mail to this address
#SBATCH --mail-user=patrick.kearns@umb.edu

# Put your job commands here, including loading any needed
# modules or diagnostic echos.
source activate picrust2
cd /hpcstor6/scratch01/p/patrick.kearns/Bullfrog_nut_enrich
picrust2_pipeline.py -s rep_seqs/dna-sequences.fasta -o picrust_norm -i asv_table/feature-table.biom
picrust2_pipeline.py -s rep_seqs/dna-sequences.fasta -o picrust_no_norm -i asv_table/feature-table.biom --skip_norm


