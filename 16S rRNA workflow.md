QIIME basic pipeline

1- make Fasting-map.txt
2- split_libraries.py -b 8 -m Fasting_map.txt -f 011216BV799wF-full.fasta -q 011216BV799wF-full.qual -o split_out_new_ref -p
3- nohup pick_open_reference_otus.py -i split_out_new_ref/seqs.fna -o open_ref_marco/ -m uclust -f -p parameters_open_reference.txt &
or
nohup pick_open_reference_otus.py -i Only_rain_112018/split_out_w53/seqs.fna -o Only_rain_112018/open_ref_w53_112018/ -m uclust -f -p parameters_open_reference.txt &
nohup pick_open_reference_otus.py -i split_out_RSF2/seqs.fna -o openref_RSF2_041316/ -f -p parameters_open_reference.txt -m usearch61 &
filter_otus_from_otu_table.py -i Dianthus_All_Samples/Dianthus_Wolbachia/open_ref_DW/otu_table_mc2_w_tax.biom -o Dianthus_All_Samples/Dianthus_Wolbachia/open_ref_DW/No_singleton_5reads/otu_table_no_singletons.biom -n 5
4- filter_taxa_from_otu_table.py -i Dianthus_paper/open_ref_dianthus/No_singleton/otu_table_no_singletons.biom -o Dianthus_paper/open_ref_dianthus/filtered_otu_table.biom -n D_4__Mitochondria,D_1__Cyanobacteria,D_2__Chloroplast,D_5__Wolbachia,D_5__Massilia,D_5__Buchnera,Unassigned
5-  biom summarize-table -i otu_table_mc2_w_tax.biom
6
