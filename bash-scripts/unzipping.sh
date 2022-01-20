#!/bin/bash

#SBATCH --job-name=unzipping
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000MB
#SBATCH --nodes=1

gunzip $1
