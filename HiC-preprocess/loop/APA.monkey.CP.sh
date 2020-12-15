#!/bin/bash
#SBATCH -J apacp
#SBATCH -p gpu
#SBATCH -N 1 
#SBATCH -o apa_cp1loops_%j.out
#SBATCH -e app_cp1loops_%j.err
#SBATCH --no-requeue
#SBATCH -A lch3000_g1
#SBATCH --qos=lch3000gpu
#SBATCH -c 1
#time java -Xmx8g -jar /lustre1/lch3000_pkuhpc/liuyt/software/bin/juicer_tools.1.7.6_jcuda.0.8.jar  apa  -r 10000 -n 10 -w 5 -u  /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/loops/hicfiles/07456A-CP_allValidPairs.hic monkey.CP.merged_loops.bedpe /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/loops/results/APA/monkey-CP1-allloops-n10-w5
#time java -Xmx8g -jar /lustre1/lch3000_pkuhpc/liuyt/software/bin/juicer_tools.1.7.6_jcuda.0.8.jar  apa  -r 10000 -n 10 -w 5 -u  /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/loops/hicfiles/07456B-CP_allValidPairs.hic monkey.CP.merged_loops.bedpe /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/loops/results/APA/monkey-CP2-allloops-n10-w5
time java -Xmx8g -jar /lustre1/lch3000_pkuhpc/liuyt/software/bin/juicer_tools.1.7.6_jcuda.0.8.jar  apa  -r 10000 -n 10 -w 5 -u  /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/loops/hicfiles/11002B_CP_allValidPairs.hic monkey.CP.merged_loops.bedpe /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/loops/results/APA/monkey-CP3-allloops-n10-w5
#time java -Xmx8g -jar /lustre1/lch3000_pkuhpc/liuyt/software/bin/juicer_tools.1.7.6_jcuda.0.8.jar  apa  -r 10000 -n 10 -w 5 -u  /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/merged-hic/hicfiles/CP_cis_allValidPairs.hic monkey.CP.merged_loops.bedpe /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/loops/results/APA/monkey-CP-allloops-n10-w5
 
