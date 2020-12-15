library(polyester)
library(Biostrings)

setwd("/lustre/user/liclab/liuyt/monkey-brain/merged-hic/RNA-Seq-new/XSAnno/step3_simNGS")
fold_changes = matrix(c(rep(1,797832),rep(1,797832)), nrow=797832)
# FASTA annotation
fasta_file = "rheMac8.blatFilter.fa"
fasta = readDNAStringSet(fasta_file)

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(20 * width(fasta) / 100)

# simulation call:
simulate_experiment('rheMac8.blatFilter.fa', reads_per_transcript=readspertx, 
                    num_reps=c(1,1), fold_changes=fold_changes, outdir='simulated_reads_rheMac8_1') 
print("simulating is done!")