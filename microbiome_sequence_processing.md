## Microbiome sequence processing
```
#load QIIME 2022.8
source activate qiime2-2023.5

#change to directory where sequences are in
cd  /hpcstor6/scratch01/p/patrick.kearns/Ribbitr_16s_DNA_RNA_Nut_Enrich

#load raw FASTQ reads into QIIME
qiime tools import --type EMPSingleEndSequences --input-path ./data --output-path nut_enrich_seqs.qza

#demultiplex reads
qiime demux emp-single \
  --i-seqs nut_enrich_seqs.qza \
 --m-barcodes-file nut_enrich_16S_map.txt \
 --m-barcodes-column BarcodeSequence \
  --o-per-sample-sequences nut_enrich_demux.qza \
  --o-error-correction-details  nut_enrich_demux-details.qza \
  --p-no-golay-error-correction \
  
#quality filer
qiime quality-filter q-score \
--i-demux  nut_enrich_demux.qza \
--o-filtered-sequences  nut_enrich_demux-filtered.qza \
--o-filter-stats  nut_enrich_demux-filter-stats.qza \
 --p-min-quality 30 

 
 #export filter stats
  qiime tools export --input-path  nut_enrich_demux-filter-stats.qza --output-path filt_stats
 
  #call ASVs with deblur
  qiime deblur denoise-16S \
  --i-demultiplexed-seqs  nut_enrich_demux-filtered.qza \
  --p-trim-length 120 \
  --p-jobs-to-start 24 \
  --o-representative-sequences  nut_enrich_rep-seqs-deblur.qza \
  --o-table  nut_enrich_table-deblur.qza \
   --o-stats  nut_enrich_deblur-stats.qza
 
 #make phylogenetic tree with fasttree
 qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences nut_enrich_rep-seqs-deblur.qza \
  --output-dir phylogeny-align-to-tree-mafft-fasttree \
  --p-n-threads 24
  
 #export deblur stats
 qiime tools export --input-path nut_enrich_deblur-stats.qza --output-path deblur_stats
  
 #export rep seqs
 qiime tools export --input-path nut_enrich_rep-seqs-deblur.qza --output-path rep_seqs
  
#export tree as NWK format
qiime tools export --input-path phylogeny-align-to-tree-mafft-fasttree/tree.qza --output-path tree
 
 #export OTU table to biom then to text file
 qiime tools export --input-path nut_enrich_table-deblur.qza --output-path asv_table
 biom convert -i asv_table/feature-table.biom --to-tsv -o asv_table.txt
 
 
 #assign taxonomy with sklearn and silva database
qiime feature-classifier classify-sklearn   --i-classifier silva-138-99-515-806-nb-classifier.qza   --i-reads  nut_enrich_rep-seqs-deblur.qza   --o-classification nut_enrich_taxonomy.qza 

 #export taxonomy file
 qiime tools export --input-path taxonomy.qzv --output-path taxonomy
 
#make stacked bar visualizations
qiime taxa barplot --o-visualization taxa_plot  --m-metadata-file nut_enrich_16S_map.txt  --i-taxonomy nut_enrich_taxonomy.qza --i-table nut_enrich_table-deblur.qza
 ```
