QIIME basic pipeline

make Fasting-map.txt

split_libraries.py -b 8 -m Fasting_map.txt -f 011216BV799wF-full.fasta -q 011216BV799wF-full.qual -o split_out -p

nohup pick_open_reference_otus.py -i split_out/seqs.fna -o open_ref/ -m uclust -f -p parameters_open_reference.txt &

filter_otus_from_otu_table.py -i open_ref/otu_table_mc2_w_tax.biom -o otu_table_no_singletons.biom -n 5

filter_taxa_from_otu_table.py -i otu_table_no_singletons.biom -o filtered_otu_table.biom -n D_4__Mitochondria,D_1__Cyanobacteria,D_2__Chloroplast,Unassigned

biom summarize-table -i otu_table_mc2_w_tax.biom

biom summarize-table -i filtered_otu_table.biom

