#/bin/bash

## Trim input fastq data of Illumina adaptor sequence.

module load bowtie2/2.2.3
module load java/1.8.0_60

if [ -f ./ref/NR_001566.1.fasta ] ; then
	echo "Reference file exist. Proceed with the script."
else
	mkdir ./ref
	cd ./ref
	echo ">scaRNA13 3' ligation adaptor
AATCTGTAGTCTTGGAGCCGCACAGGGTTGGTGGTACCCTCGAGCACACC
AGACTTGCAGAAAAAGCATACTCCAGAGGAAGCTGAGGCATGCCTGCTCG
AGAGCCAGCTGTTCCATGTGCAATTTTCCTCTGATAGTTTCTGGTCACTG
TTGCCACGGTGATAATGACTGGGCTATGTCATTATCTATCCGCCAACAGT
AAGAGAAGCTTTGCAGTCGAGATATTGTTTAGCAGATGGAGTGTTTTCTG
TTGAACACTAAGTACTGCCACAAGT""CTGTAGGCACCATCAATCGTTACGTAG" > ref.fa
	bowtie2-build ref.fa ref
	cd ..
fi

if [ -d ./unused_trimmomatic_out ] ; then
		echo "Directories already exist."
else
	mkdir unused_trimmomatic_out
	mkdir unused_flash_out
	mkdir trimmomatic_pe_out
fi
LeadQualTrim=25
MinLen=20
TailQualTrim=$LeadQualTrim

for file in $1*R1_001.fastq.gz ; do
	filename=${file%*_R1_001.fastq*}
	input1=$file
	input2=${file/_R1_/_R2_}
	samplename=${file%_S*_*}
	peread1="trimts_""$LeadQualTrim""_$MinLen""_""$samplename""_R1_001.fastq"
	peread2="trimts_""$LeadQualTrim""_$MinLen""_""$samplename""_R2_001.fastq"
	seread1="trimts_""$LeadQualTrim""_$MinLen""_""$samplename""_se1_001.fastq"
	seread2="trimts_""$LeadQualTrim""_$MinLen""_""$samplename""_se2_001.fastq"

	echo "#!/bin/bash" > tempjob
	echo "#SBATCH -N 1" >> tempjob
	echo "#SBATCH -n 1" >> tempjob
	echo "#SBATCH -c 6" >> tempjob
	echo "#SBATCH -p batch" >> tempjob
	echo "#SBATCH --time=12:00:00" >> tempjob
	echo "#SBATCH --mem=6000" >> tempjob
	echo "#SBATCH -J trim.$samplename" >> tempjob
	echo "#SBATCH -o $samplename.trim_ts.o" >> tempjob
	echo "#SBATCH -e $samplename.trim_ts.e" >> tempjob
	echo "java -Xmx4g -jar $tucf/Tools/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 6 -phred33 $input1 $input2 $peread1 $seread1 $peread2 $seread2 ILLUMINACLIP:/cluster/home/a/t/atai01/tucf_genomics/Tools/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10:1:true LEADING:$LeadQualTrim TRAILING:$TailQualTrim SLIDINGWINDOW:4:15 MINLEN:$MinLen" >> tempjob
	echo "mv trim*$samplename*se*fastq ./unused_trimmomatic_out" >> tempjob
	echo "/cluster/tufts/tucf_genomics/Tools/FLASH-1.2.11/flash -f 175 -s 10 -m 50 -M 120 -t 6 --output-prefix=$samplename $peread1 $peread2" >> tempjob
	echo "mv *$samplename*not*.fastq *$samplename*.hist* ./unused_flash_out" >> tempjob
	echo "mv trimts_*$samplename*_R*_001.fastq ./trimmomatic_pe_out" >> tempjob
	echo "bowtie2 --no-unal -p 6 -x ./ref/ref -U $samplename.extendedFrags.fastq -S $samplename.sam" >> tempjob
	echo "./pipeline_v3_scaRNA13.pl $samplename.sam" >> tempjob
	cat tempjob
	sbatch tempjob
done
rm tempjob

#zip results.zip *.txt *.log




