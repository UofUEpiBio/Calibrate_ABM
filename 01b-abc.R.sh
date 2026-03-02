#!/bin/sh
#SBATCH --job-name=Sims_calibrate
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-np
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
/uufs/chpc.utah.edu/sys/installdir/r8/R/4.4.0/lib64/R/bin/Rscript --vanilla /uufs/chpc.utah.edu/common/home/u1418987/sima/CopyOfSims_calibrate_paper/01b-abc.R
