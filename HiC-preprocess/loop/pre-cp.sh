#!/bin/bash

bash hicpro2juicebox_cw.sh -i /home/lch3000_pkuhpc/lustre1/liuyt/monkey-brain/merged-hic/validfiles/CP_cisValidPairs    -g genome.chrom.sizes -j /home/lch3000_pkuhpc/lustre1/liuyt/software/bin/juicer_tools.1.9.8_jcuda.0.8.jar -t /home/lch3000_pkuhpc/lustre1/liuyt/monkey-brain/merged-hic/results/CP-new  -o /home/lch3000_pkuhpc/lustre1/liuyt/monkey-brain/merged-hic/results/CP-new 
