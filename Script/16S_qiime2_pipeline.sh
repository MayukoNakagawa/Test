## clear previous analysis files in Results/
rm Results/Fastqc_Report/*
rm Results/Phix-removed/*
rm Results/Illumina_Adapter_Removed/*
rm Results/Illumina_Adapter_Removed/stats/*
rm Results/Primer-Filtered/fastq/*
rm Results/Demux/*
rm Results/Denoised-paired/*
rm Results/Denoised-paired/view/*
rm Results/Taxonomic-analysis/*

### define directories ####
## Sample reads
raw_data=Sample_reads

## bbtools installed directory
ADAPTSEQ=Software/bbmap/resources

## Results
rm_phix=Results/Phix-removed
rm_Adapter=Results/Illumina_Adapter_Removed
rm_primer=Results/Primer-Filtered
demux=Results/Demux
denoised=Results/Denoised-paired
taxonomy=Results/Taxonomic-analysis
phylo_div=Results/Phylogenic_diversity

## reference dataset for taxonomy assignment
classifier=Reference_Data/Silva132-99/silva-132-99-515-806-nb-classifier.qza

###remove phix reads
echo "remove Phix reads"

cd $raw_data
rm list.txt
ls *.fastq.gz > list.txt

while read R1
do read R2

   echo "Woring $R1 $R2" 

   bbduk.sh -Xmx20g in1=$R1 in2=$R2 out1=../$rm_phix/$R1 out2=../$rm_phix/$R2 ref=../$ADAPTSEQ/phix174_ill.ref.fa.gz  k=31  hdist=1
   done< list.txt

###remove Illumina Adapter
echo "remove Illumina Adapter"

while read R1
do read R2

   echo "Woring $R1 $R2"

   bbduk.sh -Xmx20g in1=$R1 in2=$R2 out1=../$rm_Adapter/$R1 out2=../$rm_Adapter/$R2 ref=../$ADAPTSEQ/nextera.fa.gz ktrim=r k=23 mink=11 hdist=1 tpe tbo overwrite=true stats=../$rm_Adapter/stats/$R1.log
   done< list.txt

###remove primers

# Forward and Reverse Sequence Primers(515F-GTGCCAGCMGCCGCGGTAA, 806R-GGACTACHVGGGTWTCTAAT)
tag5='GTGCCAGCMGCCGCGGTAA'

tag5_R='GGACTACHVGGGTWTCTAAT'

while read R1
do read R2

   echo "Working on $R1 $R2"

   if [[$R1 == *.gz]]; then

       gzip -d $R1
       gzip -d $R2
   else
       A="$( echo $R1 | cut -d'_' -f1)"
       B="$( echo $R1 | cut -d'_' -f2)"
       out_name_R1=$A"_"$B"_L001_R1_001"
       out_name_R2=$A"_"$B"_L001_R2_001"
       echo $out_name_R1

       cutadapt -g $tag5 -G $tag5_R -o ../$rm_primer/$out_name_R1".fastq"  -p ../$rm_primer/$out_name_R2".fastq" $R1  $R2  -e 0.2 #–discard-untrimmed
       gzip ../$rm_primer/$out_name_R1".fastq"
       gzip ../$rm_primer/$out_name_R2".fastq"
       mv ../$rm_primer/*.gz ../$rm_primer/fastq/
   fi
done< list.txt

########IMPORT TO QIIME2##############
cd ../

# start to use qiime2
source activate qiime2-2019.1

#(the commands are changed from 2018 ver) 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $rm_primer/fastq \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $demux/demux-paired-end.qza

qiime demux summarize \
  --i-data $demux/demux-paired-end.qza \
  --o-visualization $demux/demux-paired-end-stat.qzv

###Denoising Data###
echo "Denoise Data"

#change the trunc-len depend on the length of reads
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs  $results/$input_demux/demux-paired-end.qza \
  --p-trunc-len-f 140 \
  --p-trunc-len-r 140 \
  --o-representative-sequences  $results/$denoised/rep-seqs-dada2.qza \
  --o-table  $results/$denoised/table-dada2.qza \
  --o-denoising-stats $results/$denoised/stats-dada2.qza\
  --p-n-threads 8 \
  --p-max-ee 2 \
  --p-trunc-q  2 \
  --verbose
  
echo "Visualize denoised paired reads"


qiime feature-table summarize \
  --i-table $denoised/table-dada2.qza \
  --o-visualization $denoised/view/table.qzv \
  --m-sample-metadata-file $metadata_dir/meta_data.tsv

qiime feature-table tabulate-seqs \
  --i-data $denoised/rep-seqs-dada2.qza \
  --o-visualization $denoised/view/rep-seqs.qzv

###Taxonomy classification with SILVA database #####
##classifier=Reference_Data/Silva132-99/silva-132-99-515-806-nb-classifier.qza

echo "taxonomy assignment with using SILVA132"

qiime feature-classifier classify-sklearn \
  --i-classifier $classifier \
  --i-reads $denoised/rep-seqs-dada2.qza \
  --o-classification $taxonomy/taxonomy.qza

echo "making taxonomy table"
qiime metadata tabulate \
  --m-input-file $taxonomy/taxonomy.qza \
  --o-visualization $taxonomy/taxonomy.qzv

echo "making taxonomy barplot"
qiime taxa barplot \
  --i-table $denoised/table-dada2.qza \
  --i-taxonomy $taxonomy/taxonomy.qza \
  --m-metadata-file metadata/meta_data.tsv \
  --o-visualization $taxonomy/taxa-bar-plots.qzv

###### Phylogenic tree ####
echo "making phylogenic tree"

cho "make aligned read data"
qiime alignment mafft \
  --i-sequences $denoised/rep-seqs-dada2.qza \
  --o-alignment $phylo_div/aligned-rep-seqs.qza

echo "do alignment"
qiime alignment mask \
  --i-alignment $phylo_div/aligned-rep-seqs.qza \
  --o-masked-alignment $phylo_div/masked-aligned-rep-seqs.qza

qiime phylogeny fasttree \
  --i-alignment $phylo_div/masked-aligned-rep-seqs.qza \
  --o-tree $phylo_div/unrooted-tree.qza

qiime phylogeny midpoint-root \
  --i-tree $phylo_div/unrooted-tree.qza \
  --o-rooted-tree $phylo_div/rooted-tree.qza
