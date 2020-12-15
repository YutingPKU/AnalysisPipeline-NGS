#!/bin/bash

time java -Xmx8g -jar /lustre1/lch3000_pkuhpc/liuyt/software/bin/juicer_tools.1.7.6_jcuda.0.8.jar  dump norm VC  /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/loops/hicfiles/10084A-MEF_allValidPairs.hic  X BP 10000 chrX_norm.txt  
