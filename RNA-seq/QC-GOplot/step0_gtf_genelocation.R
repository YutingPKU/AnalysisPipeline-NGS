
gtffile = "/lustre/user/liclab/liuyt/monkey-brain/ref_genome/ensemble/Mmul_8.0.1/Macaca_mulatta.Mmul_8.0.1.89.chr.gtf"
setwd("/lustre/user/liclab/liuyt/monkey-brain/ref_genome/ensemble/Mmul_8.0.1/")
ens <- ensemblGenome()
read.gtf(ens, "Macaca_mulatta.Mmul_8.0.1.89.chr.gtf")
my_gene <- getGenePositions(ens)
my_gr <- with(my_gene, GRanges(seqid, IRanges(start, end), strand, id = gene_id, gn = gene_name))
save(my_gr, file = "/lustre/user/liclab/liuyt/monkey-brain/merged-hic/RNA-Seq/data/Macaca_mulatta.Mmul_8.0.1.89.chr.gene.GRanges.rds")

setwd("/lustre/user/liclab/publicData/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Archives/archive-2014-05-23-16-03-55/Genes/")
ens <- ensemblGenome()
read.gtf(ens, "genes.gtf")
my_gene <- getGenePositions(ens)
my_gr <- with(my_gene, GRanges(seqid, IRanges(start, end), strand, id = gene_id,gn = gene_name))
save(my_gr, file = "/lustre/user/liclab/liuyt/monkey-brain/merged-hic/RNA-Seq/data/Homo_sapiens.Ensembl.GRCh37.gene.GRanges.rds")

setwd("/lustre/user/liclab/publicData/igenomes/Mus_musculus/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/")
ens <- ensemblGenome()
read.gtf(ens, "genes.gtf")
my_gene <- getGenePositions(ens)
my_gr <- with(my_gene, GRanges(seqid, IRanges(start, end), strand, id = gene_id, gn = gene_name))
save(my_gr, file = "/lustre/user/liclab/liuyt/monkey-brain/merged-hic/RNA-Seq/data/Mus_musculus.UCSC.mm10.gene.GRanges.rds")
