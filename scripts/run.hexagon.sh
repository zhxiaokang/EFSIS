#!/bin/bash                                                                                                                                                                                                                                   

#SBATCH --ntasks=32                                                                                                                                                                                                                           
module load R
#                                                                                                                                                                                                                                             
# fix for temporary filesystem                                                                                                                                                                                                                
export TMP=/work/users/$USER/
export TMPDIR=$TMP
#                                                                                                                                                                                                                                             
mkdir $TMP/data
mkdir $TMP/scripts

cp -r ../data/Leukemia $TMP/data/
cp ./efsis_diff_num_sel.R $TMP/scripts

cd $TMP

aprun -n1 -N1 -m30000M R --no-save < ./efsis_diff_num_sel.R
