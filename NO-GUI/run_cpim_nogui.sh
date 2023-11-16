#!/bin/bash
#SBATCH -J CPIM_external_field
#SBATCH -n 1
#SBATCH -p slims
#SBATCH --output=output_%j.out
#SBATCH --error=errors_%j.err
#SBATCH --mail-user=ajlhomme@uc.cl
#SBATCH --mail-type=ALL

make
./CPIM_external_field "CPIM_external_field.txt"