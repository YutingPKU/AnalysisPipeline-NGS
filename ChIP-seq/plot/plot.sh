#ref=/lustre/user/liclab/publicData/igenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/mm10.refseq
#ref=mm10.refseq.NM
ref=/lustre/user/liclab/zhangc/Taolab/guan/ChIP-seq/peaks/merge/repeatMasker/tRNA.plot.ref
for i in `ls ../mapping/*.rmdup.bam`
do
{
	out=`basename $i`
	python ChIPseqPlot_V1.3.py -f $ref -b $i -o ${out/.sort.rmdup.bam/}".tRNA_repeat" -u 5000


}&
done

# ngsplot
#ngs.plot.r -G mm10 -R tss -C 8cell_rep1.bam -O 8cell_rep1.tss -T ATAC -L 3000 -FL 300
#ngs.plot.r -G mm9 -R genebody -C atac.cfg -O ngsPlot_genebody -P 6 -L 3000

