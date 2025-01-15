# wgcna
Run WGCNA
```r
Rscript 01_wgcna_step1_filter_exp_matrix.R --input_exp exp_file.txt --large_exp 1 --outdir .
Rscript 02_wgcna_step2_test_power.R --outdir . --cut_height 80000 --network_type unsigned
Rscript 03_wgcna_step3_construct_network.R --input_power 18 --network_type unsigned --cor_type pearson --maxBlockSize 100000 --interest_geneset fake_interesting_gene.tsv --outdir .
Rscript 04_wgcna_step4_treat_relationship.R --interest_geneset fake_interesting_gene.tsv --trait_data trait.txt --outdir .
Rscript 05_wgcna_step5_export_hubgene.R --outdir .
```
