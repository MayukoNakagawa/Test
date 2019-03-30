results=Results
fastqc=Fastqc_Report

# Input director contains *.fastq.gz
raw_data=Sample_reads

# Remove if the output is exist
#rm -rf $fastqc

# Create the outout folder
mkdir $results/$fastqc



cd $raw_data

ls *.fastq.gz > list.txt

while read R1

do read R2

echo "Working $R1 $R2" 

fastqc $R1 $R2 -o ../$results/$fastqc -t 4

done< list.txt




multiqc ../$results/$fastqc/*.zip -o $fastqc/
