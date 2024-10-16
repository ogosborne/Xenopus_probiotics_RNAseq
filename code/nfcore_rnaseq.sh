# singularity v3.11.4
# nextflow v24.04.4.5917
# nf-core/rnaseq v3.14.0-gb89fac3
nextflow run nf-core/rnaseq \
-profile scw \
--fasta XENLA_10.1_genome.fa \
--gff Xlaevisv10.17.gene.gff \
--input sample_sheet.csv \
--outdir pipeline_results \
--pseudo_aligner salmon \
-resume
