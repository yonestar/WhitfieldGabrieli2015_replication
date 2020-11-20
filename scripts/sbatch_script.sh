#!/bin/bash

#SBATCH -p blanca-ics
#SBATCH -n 1
#SBATCH -t 400
#SBATCH --mem=30G

cd /projects/joas2631

matlab -nodisplay -nosplash -nodesktop -r "pathdef; cd /work/ics/data/projects/wagerlab/labdata/collab/SADTRP_reboot/scripts; data_analytic_rep_wagertools; quit"
