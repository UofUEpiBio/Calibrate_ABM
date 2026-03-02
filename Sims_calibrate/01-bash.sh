#!/bin/sh
#SBATCH --job-name=Sims_calibrate
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1418987/sima/CopyOfSims_calibrate_paper/Sims_calibrate/02-output-%A-%a.out
#SBATCH --array=1-200
#SBATCH --partition=vegayon-np
#SBATCH --account=vegayon-np
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=Sims_calibrate
#SBATCH --ntasks=1
/uufs/chpc.utah.edu/sys/installdir/r8/R/4.4.0/lib64/R/bin/Rscript  /uufs/chpc.utah.edu/common/home/u1418987/sima/CopyOfSims_calibrate_paper/Sims_calibrate/00-rscript.r
