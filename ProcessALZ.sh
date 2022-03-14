#!/bin/bash
#SBATCH --job-name=AlzheimersProteomics    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lucast3@gene.com     # Where to send mail	
#SBATCH --ntasks=1                # Run on a single CPU
#SBATCH --partition="defq"         #Run on Himem or defq partition
#SBATCH --mem=90gb                     # Job memory request
#SBATCH --time=40:05:00               # Time limit hrs:min:sec
#SBATCH --output=/gstore/scratch/u/lucast3/ALZ_out.json   # Standard output and error log
#SBATCH --cpus-per-task=7            # Number of CPU cores per task
#SBATCH --qos=long            # quality of service for long jobs
pwd; hostname; date

ml R/prd

Rscript /gstore/scratch/u/lucast3/Alzheimers_Proteomics/ProcessData.R
 
