#!/bin/bash

#SBATCH -N1
module load gcc/6.1.0
cd /home/mb144122/share_docker/preprocessorRamOptimised
./run_programm.sh

