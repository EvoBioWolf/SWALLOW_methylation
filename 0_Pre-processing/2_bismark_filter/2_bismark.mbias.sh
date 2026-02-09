#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4763mb
#SBATCH --time=24:00:00
#SBATCH -J bismark.mbias

for RUN in $(ls *M-bias.txt | sed 's/\..*//g' ); do 

sed -n '/CHG context (R1)/q;p' ${RUN}.*M-bias.txt | egrep -v 'methyl|=|context' \
| grep . | awk -v s=${RUN} '{OFS="\t"}{print s,$0,"R1","CpG"}' > ${RUN}.MBIAS1.txt

sed -n '/CpG context (R2)/,$p' ${RUN}.*M-bias.txt | sed -n '/CHG context (R2)/q;p'  | egrep -v 'methyl|=|context' \
| grep . | awk -v s=${RUN} '{OFS="\t"}{print s,$0,"R2","CpG"}' > ${RUN}.MBIAS2.txt

done
echo -e 'Sample\tPosition\tCountMethylated\tCountUnmethylated\tPercentMethylated\tCoverage\tRead\tContext' > MBias.header
cat MBias.header *MBIAS1.txt > Methylation_Bias-Input_R1.txt
cat MBias.header *MBIAS2.txt > Methylation_Bias-Input_R2.txt

#clean-up
rm MBias.header
rm *M-Bias.txt*
