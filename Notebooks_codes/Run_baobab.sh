#!/bin/bash                                                                                                              

#SBATCH -J name                                                                                                          
#SBATCH --get-user-env                                                                                                   
#SBATCH --ntasks=1                                                                                                     
#SBATCH --cpus-per-task=1                                                                                                
#SBATCH -p dpt-bigmem-EL7                                                                                                      
#SBATCH --output=slurm-test-%J.out                                                                                       
#SBATCH -t 7-00:00:00                                                                                                     
#SBATCH --mail-type=FAIL
#SBATCH --mem-per-cpu=100000                                                                                                

##cd ##ROOTDIR##                                                                                           

module load Anaconda3
PYTHONPATH=/home/hassani/scratch/Doppler_Project/pygadgetreader:$PYTHONPATH
export PYTHONPATH

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job ID is $SLURM_JOBID

python file.py -kind lightcone -loc ./../output/Doppler_RSD_kevolution/gevolution_boxsize_4032_ngrid_4608_w_0m9_cs2_1_05062020/output -name w0d9_cs1_gevolution_lightcone0_0 -ini 165 -fin 212 -output ./gevolution_boxsize_4032_ngrid_4608_w_0m9_cs2_1_05062020.txt -e 1000

