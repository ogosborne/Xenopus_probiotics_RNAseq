# extract proteins with gffread v.0.12.7
gffread -g XENLA_10.1_genome.fa \
-y XENLA_10.17.prot.fa \
Xlaevisv10.17.gene.gff 
# annotate proteins with eggnog-mapper v.2.1.5
emapper.py \
-i XENLA_10.17.prot.fa \
--itype proteins \
--tax_scope 33208 \
--output_dir . \
--data_dir DBs/eggnogmapper \
--output emapper \
-m diamond \
--cpu 20 | tee emapper.log
