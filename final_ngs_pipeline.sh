#!/bin/bash

mkdir ngs_assessment
cd ngs_assessment
mkdir ngs_pipeline
cd ngs_pipeline
mkdir data scripts results
cd results
mkdir fastqc_untrimmed_reads
cd ..
cd data
mkdir aligned_data trimmed_fastq untrimmed_fastq reference

cd untrimmed_fastq

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz

cd ..
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz 
mv hg19.fa.gz /home/ubuntu/ngs_assessment/ngs_pipeline/data/reference

cd untrimmed_fastq
zcat NGS0001.R1.fastq.qz > NGS0001.R1.fastq
zcat NGS0001.R2.fastq.qz > NGS0001.R2.fastq
# Expanding the files since the qz version was not being read by fastqc.

# Pre-trimming QC analysis
fastqc -t 4 *.fastq # generate fastqc analysis on all fastq files
mv *fastqc* ~/ngs_assessment/ngs_pipeline/results/fastqc_untrimmed_reads # moving the fastqc files to the fastqc_untrimmed_reads directory
# fastqc gives an analysis of the quality of the data. Performing fastqc on the raw data gives an idea of the quality of the sequences and the necessary trimming steps required.
cd ../../
cd results/fastqc_untrimmed_reads
for zip in *.zip; do unzip $zip; done

cd ../..
cd data/untrimmed_fastq

# Trimmomatic
trimmomatic PE -threads 4  -phred33 /home/ubuntu/ngs_assessment/ngs_pipeline/data/untrimmed_fastq/NGS0001.R1.fastq /home/ubuntu/ngs_assessment/ngs_pipeline/data/untrimmed_fastq/NGS0001.R2.fastq -baseout /home/ubuntu/ngs_assessment/ngs_pipeline/data/trimmed_fastq/NGS0001_trimmed_R ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:25 MINLEN:50	
#Trimmomatic is used to trim the raw reads to remove low quality reads and adapter sequences so that downstream analysis such as sequence alignment are proper. The above command uses the  PE input because the given data is of paired-end reads. Threads is 4 to run the job using all 4 processors available in the openstack virtual machine. Since the reads have the phred33 quality scoring, the phred33 argument is used. The two files NGS0001_R1.fastq and NGS0001_R2.fastq are the input files on which trimmomatic will be performed and the output will be generated in the trimmed_fastq directory and the file name will be NGS0001_trimmed_R. The ILLUMINACLIP parameter is used since the reads are from an Illumina run. This will remove the adapters which in this case are from the Nextera kit from the reads. TRAILING : 25 is used to cut bases that fall below the quality of 25 from the 3' end.  MINLEN:50 will drop any reads below length 50 so that short reads are removed.


# Post trimming QC
fastqc -t 4 /home/ubuntu/ngs_assessment/ngs_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_1P
fastqc -t 4 /home/ubuntu/ngs_assessment/ngs_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_2P


# Alignment
#indexing
bwa index /home/ubuntu/ngs_assessment/ngs_pipeline/data/reference/hg19.fa.gz
# The index command from the bwa tool is used to generate the index files which will be needed by the bwa aligner to perform the alignment.

# The read group info is required to identify the samples.
# Read group info
#ID: @11V6WR1:111:D1375ACXX:1
#SM : NGS0001
#PL : ILLUMINA
#LB : Nextera-NGS-blood
#PU : 11V6WR1
#DT : 2022-04-10

bwa mem -t 4 -v 1 -R '@RG\tID:@11V6WR1:111:D1375ACXX:1\tSM:NGS0001\tPL:ILLUMINA\tLB:Nextera-NGS-blood\tDT:2023-04-10\tPU:11V6WR1' -I 250,50 /home/ubuntu/ngs_assessment/ngs_pipeline/data/reference/hg19.fa.gz /home/ubuntu/ngs_assessment/ngs_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_1P /home/ubuntu/ngs_assessment/ngs_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_2P > /home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/NGS0001_R.sam
#The - R input is to add the read group information which is important for sample identification. -I input is for the mean and standard deviation of insert sizes which is required for paired reads. The alignment will occur between the hg19 reference file and the trimmed reads NGS0001_trimmed_R_1P and NGS0001_trimmed_R_2P. The output will be generated in the directory aligned_data as NGS0001_R.sam. 

samtools view -h -b /home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/NGS0001_R.sam > /home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/NGS0001_R.bam
# This command coverts the sam file to bam file (-b) while including the header in the output(-h).

samtools sort /home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/NGS0001_R.sam > /home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/NGS0001_sorted.bam
# This command sorts the aligments based on the coordinates and the output will be generated in the aligned_data directory as NGS0001_sorted.bam.

samtools index /home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/NGS0001_sorted.bam
# This command  will generate index for the NGS0001_sorted.bam file. This is required to extract alignments that overlap in one region and is also required to view files in genome viewers.

picard MarkDuplicates I=/home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/NGS0001_sorted.bam O=/home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/NGS0001_sorted_marked.bam M=/home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/marked_dup_metrics.txt
# The Markduplicates command in Picard is used to mark the duplicate reads present. This tool will go through the alignment data which is provided by using the command  "I=" in either a SAM or BAM format (in this case the bam file NGS0001_sorted.bam) and will locate the  duplicate reads. The "O=" command specifies the output file which is a bam file where the duplicate reads are marked. The "M=" command specifies the metrics file which gives a summary of the duplicate reads found.

samtools index /home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/NGS0001_sorted_marked.bam
# This command  will generate index for the NGS0001_sorted_marked.bam.

cd /home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
# The -F command is used to filter reads. Using command -q, 20 is set  as the minimum MAPQ score required. Filtering on bitwise flag is done on the following criteria : The read is unmapped ;  The alignment or this read is not primary ; The read fails platform/vendor quality checks ;  The read is a PCR or optical duplicate. The command 1796 indicates that the above  combination of bitwise flag filters are used. This command will generate the file NGS0001_sorted_filtered.bam

samtools index NGS0001_sorted_filtered.bam
#This command  will generate index for the NGS0001_sorted_filtered.bam


samtools idxstats NGS0001_sorted_filtered.bam > /home/ubuntu/ngs_assessment/ngs_pipeline/results/idxstats.txt
#Samtools idxstats command gives a summary of alignment statistics

samtools flagstat NGS0001_sorted_filtered.bam > /home/ubuntu/ngs_assessment/ngs_pipeline/results/flagstat_results.txt
#Samtools flagstat gives the count of the number of alignments for each FLAG type.

picard CollectInsertSizeMetrics I=NGS0001_sorted_filtered.bam O=/home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001_insert_size_metrics.txt
#Picard collectInsertSizeMetrics command gives the insert size statistics. (I= input file ; O= output file)

bedtools coverage -a /home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam -b /home/ubuntu/ngs_assessment/ngs_pipeline/data/annotation.bed > /home/ubuntu/ngs_assessment/ngs_pipeline/results/coverage_summary.txt
#bedtools coverage gives the depth and breadth of coverage of one file (NGS0001_sorted_filtered.bam) by comparing and looking for overlaps with another file aligned data against the annotated bed file.

# freebayes
zcat /home/ubuntu/ngs_assessment/ngs_pipeline/data/reference/hg19.fa.gz > /home/ubuntu/ngs_assessment/ngs_pipeline/data/reference/hg19.fa
# Using the zcat command to expand the reference file

samtools faidx /home/ubuntu/ngs_assessment/ngs_pipeline/data/reference/hg19.fa
#samtools faidx indexes reference sequences in fast format

freebayes -f /home/ubuntu/ngs_assessment/ngs_pipeline/data/reference/hg19.fa /home/ubuntu/ngs_assessment/ngs_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam > /home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001.vcf
#The freebayes -f command uses the reference genome file, the sorted and filtered bam file and gives a vcf output which contains the variant calls

bgzip /home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001.vcf
# bgzip command is used to compress the vcf file.

tabix -p vcf /home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001.vcf.gz
#The tabix tool is used to create an index

vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
	~/ngs_assessment/ngs_pipeline/results/NGS0001.vcf.gz > ~/ngs_assessment/ngs_pipeline/results/NGS0001_filtered.vcf
#vcffilter -f command is used to specify the filtering criteria. The filtering expression needs to be mentioned in double quotes ("") following the -f command followed by the vcf file and the output file.

bedtools intersect -header -wa -a ~/ngs_assessment/ngs_pipeline/results/NGS0001_filtered.vcf -b /home/ubuntu/ngs_assessment/ngs_pipeline/data/annotation.bed \
 	> ~/ngs_assessment/ngs_pipeline/results/NGS0001_filtered_new.vcf
#-header parameter gives the header for the first file (NGS0001_filtered.vcf). -wa is used so that each original feature in the first file (vcf) is written for each overlap. - a is used to compare each feature in the first file (vcf) to that in the second file (bed) in search of overlaps

bgzip /home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001_filtered_new.vcf
#bgzip command is used to compress the vcf file.

tabix -p vcf /home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001_filtered_new.vcf.gz
#The tabix tool is used to create an index

cd /home/ubuntu/ngs_assessment/ngs_pipeline/data/annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

chmod +x annotate_variation.pl

./convert2annovar.pl -format vcf4 /home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001_filtered_new.vcf.gz > /home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001_filtered_new.avinput
#The convert2annovar.pl -format vcf4 command is used to convert the vcf file into a file format required by annovar to generate annotations. It converts the genotype calling format into the annovar format.

./table_annovar.pl ~/ngs_assessment/ngs_pipeline/results/NGS0001_filtered_new.avinput humandb/ -buildver hg19 -out ~/ngs_assessment/ngs_pipeline/results/NGS0001_filtered_new.avinput -remove  -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout
# The table_annovar.pl command takes vcf files and annotates them using the hg 19 build version and gives the output. 

snpEff -Xmx8g -v GRCh37.75 /home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001_filtered_new.vcf.gz > /home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001_snpeffann_new.vcf -csvStats /home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001_snpeff_new_stats.csv -s /home/ubuntu/ngs_assessment/ngs_pipeline/results/NGS0001_snpeff_new_summary.html
#the snpEff command is used to annotate the filtered vcf file using the GRCh37.75 reference genome version. -Xmx8g is used to increase the heapspace to 8gb to enable the execution of the command. -csvStats is used to generate a csv output and -s generates a html output.
