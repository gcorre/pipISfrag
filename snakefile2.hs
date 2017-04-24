#	This code is a snakemake pipeline integrating different tools to process Illumina 
#	sequencing paired-end reads for vector Insertion site calling. It provides
#	file tracking and input/output management.
#
#	Copyright (C) 2017 GENETHON, Guillaume Corre
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# docker run --rm --hostname "fragIS" -w /home/ -it -v /d/:/home/ gc/frag:V2 bash

print("#~~~~~~~~Version 24/04/2017~~~~~~~~#");
print("");
print("#       #    ##    ######    #######");
print(" #     #     ##    ##        ##   ##");
print("  #   #      ##    ######    #######");
print("   # #       ##        ##    ##   ##");
print("    #        ##    ######    ##   ##");
print("");
print("#~~~~~~~~~~~~~GENETHON~~~~~~~~~~~~~#");
	
	
# Input files R1 and R2 must be in the current path.
# Be careful to sample name, use the pattern : 
#
#
#
#
#
#
#
#
#
#
#	
#	
#workdir: "/home/"

from snakemake.utils import min_version
min_version("3.11.2")


BOWTIE2_INDEX = "/home/references/genomes/homo_sapiens/hg19_GRCh37/index_bowtie2/GRCh37.75"
GENOME="hg19"

import os
HOME = os.path.expanduser('~')

from datetime import date
#TODAY=str(date.today())
TODAY="2017-04-03"

ANNOTATION="dataset/"+GENOME+"/"+TODAY+"/"



from glob import glob;
FILES = glob('*_R[0-9]_*.fastq.gz');
import re;
TRAILING=list(set([re.search('_([0-9]+).fastq.gz',FILES[0]).group(1) for w in FILES]))[0]
FILES=[re.sub('_R[1-2]_[0-9]+.fastq.gz', '',w) for w in FILES];
SAMPLE=list(set(FILES))[0];
TAG=re.search('TAG([0-9]+)', SAMPLE).group(1)


rule targets:
	input: SAMPLE+"_TAG"+TAG+"_report.pdf"
		

rule Prepare_folders:
		input: 
		output: annot=ANNOTATION
		run:
			import os;
			if not os.path.exists(output[0]):
				os.makedirs(output[0])

	
rule Demultiplex_TAGS:
	input: R1="{NAME}_R1_"+TRAILING+".fastq.gz", R2="{NAME}_R2_"+TRAILING+".fastq.gz"
	output: R1="01-Demultiplex/{NAME}_R1_TAG{TAG}.fastq.gz",R2="01-Demultiplex/{NAME}_R2_TAG{TAG}.fastq.gz"
	log: "log/demultiplex-merged.log"
	threads: 1
	message: "Demultiplexing according to TAGs list"
	shell: """
			 fastq-multx -B ../dataset/illuminaTAGs_eautils.txt -b -x -m 1 {input.R1} {input.R2} -o 01-Demultiplex/{wildcards.NAME}_R1_%.fastq.gz 01-Demultiplex/{wildcards.NAME}_R2_%.fastq.gz > {log};
			"""
			
			
rule Find_provirus:
	input: R1=rules.Demultiplex_TAGS.output.R1, R2=rules.Demultiplex_TAGS.output.R2
	output: R1Prov="02-provirus/{NAME}_R1_TAG{TAG}_Prov.fastq", R2Prov="02-provirus/{NAME}_R2_TAG{TAG}_Prov.fastq", R1noProv="02-provirus/{NAME}_R1_TAG{TAG}_noProv.fastq", R2noProv="02-provirus/{NAME}_R2_TAG{TAG}_noProv.fastq",
	log:"log/Trim-Provirus_in_R1.log"
	threads: 1
	message: "Look for the Provirus in R1"
	shell: """
			cutadapt -g aaaatctctagcagtggcgcccgaacag -O 28 -e 0.1 --no-trim --no-indels -o {output.R1Prov} -p {output.R2Prov} --untrimmed-paired-output {output.R2noProv} --untrimmed-output {output.R1noProv} {input.R1} {input.R2} > {log}
			"""
			
			
rule Find_provirus_skewer:
	input: R1=rules.Demultiplex_TAGS.output.R1, R2=rules.Demultiplex_TAGS.output.R2
	output: R1noProv="02-provirus/{NAME}_R1_TAG{TAG}-unassigned-pair1.fastq" , R2noProv="02-provirus/{NAME}_R1_TAG{TAG}-unassigned-pair2.fastq"
	threads: 1
	log:"log/Trim-Provirus_in_R1.log"
	shell:"""
			skewer -m head -r 0.1 -d 0 -l 20 -t 8 -k 28 -b -x aaaatctctagcagtggcgcccgaacag -y NNNNNNNNNNNNNNNNNNNNNNNNNNNN -o "02-provirus/{wildcards.NAME}_R1_TAG{wildcards.TAG}" {input.R1} {input.R2} 2> {log}
			"""
	
	
rule Find_LTR:
	input: R1=rules.Find_provirus_skewer.output.R1noProv, R2=rules.Find_provirus_skewer.output.R2noProv
	output: R1LTR="03-LTR/{NAME}_R1_TAG{TAG}_LTR.fastq",
			R2LTR="03-LTR/{NAME}_R2_TAG{TAG}_LTR.fastq", 
			R1noLTR="03-LTR/{NAME}_R1_TAG{TAG}_noLTR.fastq", 
			R2noLTR="03-LTR/{NAME}_R2_TAG{TAG}_noLTR.fastq"
	log:"log/Trim-LTR_in_R1.log"
	threads: 1
	message: "Look for the LTR in R1"
	shell: """
			cutadapt -g GTCTGTTGTGTGACTCTGGTAAC -O 23 -e 0.1 --no-indels -o {output.R1LTR} -p {output.R2LTR} --untrimmed-paired-output {output.R2noLTR} --untrimmed-output {output.R1noLTR} {input.R1} {input.R2} > {log}
			"""
			
rule Find_LTR_skewer:
	input: R1=rules.Find_provirus_skewer.output.R1noProv, R2=rules.Find_provirus_skewer.output.R2noProv
	output: R1LTR="03-LTR/{NAME}_R1_TAG{TAG}_LTR.fastq", 
			R2LTR="03-LTR/{NAME}_R2_TAG{TAG}_LTR.fastq", 
			R1noLTR="03-LTR/{NAME}_R1_TAG{TAG}_noLTR.fastq",
			R2noLTR="03-LTR/{NAME}_R2_TAG{TAG}_noLTR.fastq",
	log:"log/Trim-LTR_in_R1.log"
	threads: 1
	message: "Look for the LTR in R1"
	shell: """
			skewer -m head -r 0.1 -d 0 -l 20 -t 8 -k 23 -b -x GTCTGTTGTGTGACTCTGGTAAC -y NNNNNNNNNNNNNNNNNNNNNNN -o "03-LTR/{wildcards.NAME}_R1_TAG{wildcards.TAG}"
			"""

rule Find_Elong:
	input: R1=rules.Find_LTR.output.R1LTR, R2=rules.Find_LTR.output.R2LTR
	output: R1elong="04-ELONG/{NAME}_R1_TAG{TAG}_Elong.fastq", 
			R2elong="04-ELONG/{NAME}_R2_TAG{TAG}_Elong.fastq",
			R1noelong="04-ELONG/{NAME}_R1_TAG{TAG}_noElong.fastq",
			R2noelong="04-ELONG/{NAME}_R2_TAG{TAG}_noElong.fastq"
	log:"log/Trim-ELONG_in_R1.log"
	threads: 1
	message: "Look for the Elong in R1"
	shell: """
			cutadapt -g TAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTA -O 43 -e 0.16 --no-indels -o {output.R1elong} -p {output.R2elong} --untrimmed-paired-output {output.R2noelong} --untrimmed-output {output.R1noelong} {input.R1} {input.R2} > {log}
			"""

rule Start_with_GCA:
	input: R1=rules.Find_Elong.output.R1elong, R2=rules.Find_Elong.output.R2elong
	output: GCAR1="05-GCA/{NAME}_R1_TAG{TAG}_GCA.fastq",
			noGCAR1="05-GCA/{NAME}_R1_TAG{TAG}_noGCA.fastq",
			GCAR2="05-GCA/{NAME}_R2_TAG{TAG}_GCA.fastq", 
			noGCAR2="05-GCA/{NAME}_R2_TAG{TAG}_noGCA.fastq"
	message: "Filtering R1 reads that start with GCA"
	log: "log/Trim-GCA_starting_R1.log"
	threads: 1
	shell: """
			cutadapt -g ^GCA  -e 0 --no-indels -o {output.GCAR1} -p {output.GCAR2} --untrimmed-output {output.noGCAR1} --untrimmed-paired-output {output.noGCAR2} {input.R1} {input.R2} > {log};
			"""			

rule Find_Linker_R1:
	input: R1=rules.Start_with_GCA.output.GCAR1, R2=rules.Start_with_GCA.output.GCAR2
	output: R1Linker="06-LinkerR1/{NAME}_R1_TAG{TAG}_LinkerR1.fastq",
			R2Linker="06-LinkerR1/{NAME}_R2_TAG{TAG}_LinkerR1.fastq", 
			R1noLinker="06-LinkerR1/{NAME}_R1_TAG{TAG}_noLinkerR1.fastq",
			R2noLinker="06-LinkerR1/{NAME}_R2_TAG{TAG}_noLinkerR1.fastq"
	log:"log/Trim-Linker_in_R1.log"
	threads: 1
	message: "Look for the Linker in R1"
	shell: """
			cutadapt -a GTCCCTTAAGCGGAGCCCT -O 8 -e 0.2 --no-indels -o {output.R1Linker} -p {output.R2Linker} --untrimmed-paired-output {output.R2noLinker} --untrimmed-output {output.R1noLinker} {input.R1} {input.R2} > {log}
			"""	
 	
rule Find_Linker_R2:
	input: R1=rules.Find_Linker_R1.output.R1noLinker, R2=rules.Find_Linker_R1.output.R2noLinker
	output: R1Linker="06-LinkerR2/{NAME}_R1_TAG{TAG}_LinkerR2.fastq", 
			R2Linker="06-LinkerR2/{NAME}_R2_TAG{TAG}_LinkerR2.fastq",
			R1noLinker="06-LinkerR2/{NAME}_R1_TAG{TAG}_noLinkerR2.fastq", 
			R2noLinker="06-LinkerR2/{NAME}_R2_TAG{TAG}_noLinkerR2.fastq"
	log:"log/Trim-Linker_in_R2.log"
	threads: 1
	message: "Look for the Linker in R2"
	shell: """
			cutadapt -g AGGGCTCCGCTTAAGGGAC -O 19 -e 0.16 --no-indels -o {output.R2Linker} -p {output.R1Linker} --untrimmed-paired-output {output.R1noLinker} --untrimmed-output {output.R2noLinker} {input.R2} {input.R1} > {log}
			"""				
#######################################################################################
#######################################################################################			
## Process R1 reads that contain the Linker:
#######################################################################################
#######################################################################################
	
rule Collapse_identical_reads:
	input: rules.Find_Linker_R1.output.R1Linker
	output: file = "07-collapsed/{NAME}_R1_TAG{TAG}_LinkerR1_trimmed.fastq"
	message: "Collapsing R1 with linker identical reads before mapping"
	log:"log/collapse_R1wLinker.log"
	threads: 1
	shell: """
			#fastx_collapser -v -i {input} -o {output} > {log}
			 seqcluster collapse -f {input} -d -o "07-collapsed/" -m 1
			"""

rule FilterSize_R1collapsed:
	input: rules.Collapse_identical_reads.output.file
	output: filtered="07-collapsed/{NAME}_R1_TAG{TAG}_LinkerR1_trimmed_filtered.fastq",
			tooshort="07-collapsed/{NAME}_R1_TAG{TAG}_LinkerR1_trimmed_2short.fastq"
	message:"Filtering reads shorter than 20bp"
	log:"log/collapsed_R1_filtersize.log"
	threads:1
	shell: """
			cutadapt -m 20 -o {output.filtered} --too-short-output {output.tooshort} {input} > {log} 
			"""
			
			
####
# Map reads containing the linker in R1, remove reads with multiple hits (XS:i: TAG) and reads with mismatches in the first 3 positions after the LTR (taking into account the orientation)
# Keep the header in sam file
#### 	
		
rule Map_R1_with_Linker :
	input: rules.FilterSize_R1collapsed.output.filtered
	output: mappedsam="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_mapped.sam",mapped="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_mapped.fastq", unmapped="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_unmapped.fastq", exact="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_mapped_exact.sam"
	log: met="log/mappingR1wLinker.log",summary="log/mappingR1wLinker_numbers.log"
	threads: 8
	message: "Mapping R1 reads that contain the linker and filter out those with mismatches in the 3 first positions after the LTR"
	shell: """
			bowtie2 -N 0 -L 25 -i S,25,0 --score-min L,0,-0.16 --gbar 10 -p {threads} -x {BOWTIE2_INDEX} --no-unal --un {output.unmapped} --al {output.mapped} --met-file {log.met} {input} -S {output.mappedsam} 2> {log.summary}
			awk -F "\\t" '(/^@/) || ($2==0 && !/XS:i:/ && !/MD:Z:[012][A-Za-z].*\t/) || ($2==16 && !/XS:i:/ && !/MD:Z:.*[A-Za-z][012]\t/) {{print $0}}' {output.mappedsam} > {output.exact}
			"""
# awk 'function abs(v) {return v < 0 ? -v : v};$0~/XS:i:/ {AS=substr($12,6,length($12)); XS=substr($13,6,length($13));{if((AS==0 && XS<0) || (abs(AS)+0 < abs(XS))){print $0}};next} {print $0}' TM024-1-03-library-TAG3_S2_L001_R1_TAG3_mapped.sam 

			
####
# convert SAM to BAM and discard bad quality alignments.
# Sort BAM and make index for visualization.
# Convert BAM to BED.
# Add read length and number of collapsed reads to the BED file.
# Convert Start/end positions to IS position depending on the strand.
# Sort the Bed file, make clusters by x nt and resort.
# Group by cluster and read length
# Compute sum of reads for each size and number of "variants" (reads of the same size mapping on same position but with mismatches).
# Group by cluster
# Compute number of reads and number of different fragment length for the same IS.
####

rule Filter_Map_R1:
	input: rules.Map_R1_with_Linker.output.exact
	output: bam=temp("08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_mapped.bam"),
			sortbam="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_mapped_sorted.bam",
			idx="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_mapped_sorted.bai",
			bed="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_mapped_sorted.bed",
			ISsortednum="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_IS_sortedNum.bed",
			IS=("08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_IS.bed"),
			ISsorted="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_IS_sorted.bed",
			IScollapsedSize="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_IS_collapsedBySizeCluster.bed",
			IScluster="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_IS_cluster.bed",
			ISclustersorted="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_IS_cluster_sorted.bed",
			IScollapsed="08-mappingR1wLinker/{NAME}_R1_TAG{TAG}_IS_Collapsed.bed"
	params: qual="10", collapsing_window="-1"
	threads: 1
	message:"Filtering mapped reads, converting to BED,Collapsing by 3nt and counting read abundance by IS"
	shell: """ 
			samtools view -q {params.qual} -bS {input} > {output.bam};
			samtools sort {output.bam} > {output.sortbam};
			samtools index {output.sortbam} {output.idx};
			bedtools bamtobed -i {output.sortbam} > {output.bed};
			awk 'function abs(a){{return ((a < 0) ? -a : a)}} OFS="\\t" {{split($4,a,"x");print $0,abs($3-$2),a[2]}}' {output.bed} > {output.ISsortednum};
			awk 'OFS="\\t" {{if($6=="-") {{print $1,$3-1,$3,$4,$5,$6,$7,$8}} else {{print $1,$2,$2+1,$4,$5,$6,$7,$8}}}}' {output.ISsortednum} > {output.IS};
			bedtools sort -i {output.IS} > {output.ISsorted};
			bedtools cluster -s  -d {params.collapsing_window} -i {output.ISsorted} > {output.IScluster};
			sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -k 7,7n {output.IScluster} > {output.ISclustersorted};
			bedtools groupby -i {output.ISclustersorted}  -g 9,7 -c 1,2,3,4,5,6,7,8,8,9 -o distinct,mode,mode,collapse,median,distinct,distinct,sum,count,distinct  | cut -f 3- > {output.IScollapsedSize};
			bedtools groupby -i {output.IScollapsedSize}  -g 10 -c 1,2,3,4,5,6,8,7 -o distinct,mode,mode,collapse,median,distinct,sum,count_distinct | cut -f 2-  > {output.IScollapsed}
			"""
			
#######################################################################################
#######################################################################################
# Process paired reads with the linker in R2 but not in R1	
# Amplicon may be too long to include the complete linker in R1.
# R1 was already trimmed before for LTR and elongation. R2 was already trimmed for the linker before.
#######################################################################################
 #######################################################################################
 
# look for the reverse-complement Elongation in R2 3'end.
rule Find_RC_LTR_R2:
	input: R1=rules.Find_Linker_R2.output.R1Linker,R2=rules.Find_Linker_R2.output.R2Linker
	output: trimmedLTRRCR1="06-RCLTRR2/{NAME}_R1_TAG{TAG}_trimmedElongRCR2.fastq", 
			trimmedLTRRCR2="06-RCLTRR2/{NAME}_R2_TAG{TAG}_trimmedElongRCR2.fastq"
	log: "log/Trim-RC-LTR_in_R2.log"
	threads: 1
	message: "Trimming R2 that have LTR"
	shell: """
			cutadapt -a TGCTAGAGATTTTCCACACTGACTAAAAGGGTCTGAGGGATCTCTA -q 20 -O 8 -e 0.25 -o {output.trimmedLTRRCR2} -p {output.trimmedLTRRCR1} {input.R2} {input.R1} > {log};
			"""

rule FilterSize_R1R2:
	input: R1=rules.Find_RC_LTR_R2.output.trimmedLTRRCR1, R2=rules.Find_RC_LTR_R2.output.trimmedLTRRCR2
	output: R1filtered="06-RCLTRR2/{NAME}_R1_TAG{TAG}_filtered.fastq",
			R1tooshort="06-RCLTRR2/{NAME}_R1_TAG{TAG}_2short.fastq",
			R2filtered="06-RCLTRR2/{NAME}_R2_TAG{TAG}_filtered.fastq",
			R2tooshort="06-RCLTRR2/{NAME}_R2_TAG{TAG}_2short.fastq"
	message:"Filtering any R1R2 read pairs if any is shorter than 20bp"
	log:"log/R1R2_filtersize.log"
	threads:1
	shell: """
			cutadapt -m 20 -o {output.R1filtered} -p {output.R2filtered} --too-short-output {output.R1tooshort}  --too-short-paired-output {output.R2tooshort} --pair-filter any {input.R1} {input.R2} > {log} 
			"""			


####################################################
# Map paired end reads, trim 3' end to remove uncut linker parts, tolerate dovetail and max fragments length 800bp. Keep only concordant mapping (same Chr, <800bp).		
	
rule Map_R1_R2_pairs:
	input: R1=rules.FilterSize_R1R2.output.R1filtered, R2=rules.FilterSize_R1R2.output.R2filtered
	output: mapped="08-mappingR1R2/{NAME}_TAG{TAG}_mapped.sam", R1unal="08-mappingR1R2/{NAME}_R1_TAG{TAG}_unaligned.fastq"
	message: "Mapping R1 and R2 reads as pairs"
	params: fragLenght="800", Trim3="5"
	log: met="log/mapping_R1R2.log",summary="log/mapping_R1R2_numbers.log"
	threads: 8
	shell: """
			bowtie2 -p {threads}  -N 0 -L 25 -i S,25,0 --score-min L,0,-0.16 --gbar 10 -x {BOWTIE2_INDEX} -1 {input.R1} -2 {input.R2} -S {output.mapped} --trim3 {params.Trim3} --no-unal --al-conc "08-mappingR1R2/{wildcards.NAME}_R%_TAG{wildcards.TAG}_aligned.fastq" --un-conc "08-mappingR1R2/{wildcards.NAME}_R%_TAG{wildcards.TAG}_unaligned.fastq" --dovetail -X {params.fragLenght} --no-mixed --met-file {log.met} 2> {log.summary};
			"""

#########################################################################
# Filter out pairs of reads for which the mapping quality score is below 20, meaning p= 0.01 that the reported alignment may occured elsewhere.
# Then keep only R1 reads from filtered out reads and co
# 
			
rule Filter_Map_R1R2_unique:
	input: rules.Map_R1_R2_pairs.output.mapped
	output: badqual = "08-mappingR1R2/{NAME}_TAG{TAG}_bad_qual.sam",
			confident="08-mappingR1R2/{NAME}_TAG{TAG}_qual_gt20.sam",
			R1_badqual= "08-mappingR1R2/{NAME}_R1_TAG{TAG}_qual_lt20.sam" ,
			R1_badqualfq= "08-mappingR1R2/{NAME}_R1_TAG{TAG}_qual_lt20.fastq"
	params: qual="20"
	message: "filtering bad quality alignments that may be multiple hits, keep R1 reads and convert to fastq."
	log: "log/mapping_R1R2_filtering.log"
	threads :1
	shell: """
			samtools view -h -q {params.qual} {input} -o {output.confident} -U {output.badqual};
			samtools view -h -f 64 {output.badqual} > {output.R1_badqual};
			samtools flagstat {output.R1_badqual} > {log};
			samtools view {output.R1_badqual} -b | samtools fastq - > {output.R1_badqualfq};
			"""
			
#########################################################################
# Sort SAM alignement by read name and convert to SAM then BEDPE with mate1 first for reproductibility of filtering steps.  
# Sort sam by position and make an index for IGV .
# Extract IS position, fragment length from BEDPE.
# Could do the last steps together but it may be interesting to filter out length with low reads count ?.


rule Filter_Map_R1R2:
	input: sam = rules.Filter_Map_R1R2_unique.output.confident
	output: bam = "08-mappingR1R2/{NAME}_TAG{TAG}.bam",
			bed = "08-mappingR1R2/{NAME}_TAG{TAG}.bed",
			IS= "08-mappingR1R2/{NAME}_TAG{TAG}_IS.bed",
			idx= "08-mappingR1R2/{NAME}_TAG{TAG}.bai",
			ISsorted="08-mappingR1R2/{NAME}_TAG{TAG}_sorted.bed",
			IScluster="08-mappingR1R2/{NAME}_TAG{TAG}_clusters.bed",
			IScollapsedSize="08-mappingR1R2/{NAME}_TAG{TAG}_IS_collapsedBySizeCluster.bed",
			ISclustersorted="08-mappingR1R2/{NAME}_TAG{TAG}_IS_cluster_sorted.bed",
			IScollapsed="08-mappingR1R2/{NAME}_TAG{TAG}_IS_Collapsed.bed"
			
			
			
	threads: 1
	message: "Sort SAM by read name, convert to bed, find IS and calculate fragments size from R1 and R2 extremities"
	shell: """
			samtools sort {input.sam} -o {output.bam} -O BAM -n;
			bedtools bamtobed -bedpe -mate1 -i {output.bam} > {output.bed};
			samtools sort {input.sam} -o {output.bam} -O BAM;
			samtools index {output.bam} {output.idx};
			awk 'OFS="\\t" {{if($9=="+"){{print $1,$2,$2+1,$7,$8,$9, $6-$2}} else {{print $1,$3-1,$3,$7,$8,$9,$3-$5}}}}' {output.bed} > {output.IS};
			bedtools sort -i {output.IS} > {output.ISsorted};
			bedtools cluster -s -d -1 -i {output.ISsorted} > {output.IScluster};
			sort -k 8,8n -k 7,7n {output.IScluster} > {output.ISclustersorted};
			bedtools groupby -i {output.ISclustersorted} -g 8,7 -c 1,2,3,4,5,6,7,7,8 -o distinct,mode,mode,collapse,median,distinct,distinct,count,distinct  | cut -f 3- > {output.IScollapsedSize};
			bedtools groupby -i {output.IScollapsedSize} -g 9 -c 1,2,3,4,5,6,8,7 -o distinct,mode,mode,collapse,median,distinct,sum,count_distinct | cut -f 2-  > {output.IScollapsed}
			"""
			
			
#######################################################################################
#######################################################################################
# Process R1 reads that do not contain the linker sequence nor in their R2 counterpart
# These reads are only use for a qualitative anaylsis of IS, not for the quantification.	
#######################################################################################
#######################################################################################
rule Merge_R1alones:
	input: A=rules.Find_Linker_R2.output.R1noLinker, 
			B=rules.Filter_Map_R1R2_unique.output.R1_badqualfq,
			C = rules.Map_R1_R2_pairs.output.R1unal, 
			D=rules.Map_R1_with_Linker.output.unmapped
	output: "07-MergeR1only/{NAME}_R1_TAG{TAG}_merged.fastq"
	message: "Merging R1 reads without Linker and R1 from pairs that were mapped multiple times and R1 reads from pairs that do not map concordantly."
	threads:1
	log: 
	shell: """
			cat {input.A} {input.B} {input.C} {input.D} > {output} 
			"""
			
			
## This step will be removed when using the sonication method.
	
rule Cut_TTAA_R1_noLinker:
	input: rules.Merge_R1alones.output
	output: TTAA="07-TTAA_cutR1alone/{NAME}_R1_TAG{TAG}_TTAA-cut.fastq", 
			TTAA_cut="07-TTAA_cutR1alone/{NAME}_R1_TAG{TAG}_TTAA-cut20bp.fastq",
			tooshort="07-TTAA_cutR1alone/{NAME}_R1_TAG{TAG}_TTAA-cut20bp_2short.fastq"
	message: "Cutting TTAA sequences that are still present and filter reads shorter than 20bp"
	threads: 1
	log: "log/Cut_TTAA_Merge_R1.log"
	shell: """
			  awk '{{if(NR %4 ==2) {{x=index($0,"TTAA");if(x>0){{print substr($0,1,x)}}else{{print $0}}}} else {{if(NR % 4 ==0 && x>0) {{print substr($0,1,x)}} else{{print $0}}}}}}' {input} > {output.TTAA};
			  cutadapt -m 20 -o {output.TTAA_cut} --too-short-output {output.tooshort} {output.TTAA} > {log}
			"""	


			
rule Collapse_identical_R1_noLinkerReads:
	input: rules.Cut_TTAA_R1_noLinker.output.TTAA_cut
	output: file = "07-collapse_R1Merged/{NAME}_R1_TAG{TAG}_TTAA-cut20bp_trimmed.fastq"
	message: "Collapsing identical reads before mapping"
	log:"log/Collapse_R1_Merge.log"
	threads: 1
	shell: """
			#fastx_collapser -v -i {input} -o {output} > {log}
			 seqcluster collapse -f {input} -d -o "07-collapse_R1Merged/" -m 0
			"""	

			
rule Map_R1_noLinker:
		input: R1=rules.Collapse_identical_R1_noLinkerReads.output.file
		output: mappedsam="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_mapped.sam",mapped="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_mapped.fastq",unmapped="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_unmapped.fastq", exact="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_mapped_exact.sam"
		message:"Mapping R1 reads without Linker in R1 nor R2 and R1 that were unmapped in previous steps"
		threads: 8
		log: met="log/mapping_Merged_R1only.log",summary="log/mappingR1alone_numbers.log"
		shell: """
				bowtie2 -N 0 -L 25 -i S,25,0 --score-min L,0,-0.16 --gbar 10 -p {threads} -x {BOWTIE2_INDEX} --trim3 5 --no-unal --un {output.unmapped} --al {output.mapped} --met-file {log.met} {input.R1} -S {output.mappedsam} 2> {log.summary}
				awk -F "\\t" '(/^@/) || ($2==0 && !/XS:i:/ && !/MD:Z:[012][A-Za-z].*\t/) || ($2==16 && !/XS:i:/ && !/MD:Z:.*[A-Za-z][012]\t/) {{print $0}}' {output.mappedsam} > {output.exact}
				"""


			
rule Filter_Map_R1_noLinker:
	input: rules.Map_R1_noLinker.output.exact
	output: bam="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_mapped.bam",
			sortbam="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_mapped_sorted.bam",
			idx="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_mapped_sorted.bai",
			bed="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_mapped_sorted.bed",
			IS="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_IS.bed",
			ISsorted="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_IS_sorted.bed",
			ISsortednum="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_IS_sortedNum.bed",
			IScollapsedSize="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_IS_collapsedBySizeCluster.bed",
			IScluster="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_IS_cluster.bed",
			ISclustersorted="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_IS_cluster_sorted.bed",
			IScollapsed="08-mapping_MergedR1NoLinker/{NAME}_R1_TAG{TAG}_IS_Collapsed.bed"
	params: qual="10"
	threads: 1
	message:"Filtering mapped reads, converting to BED, Collapsing by 3nt and counting read abundance by IS"
	shell: """ 
			samtools view -q {params.qual} -bS {input} > {output.bam};
			samtools sort {output.bam} > {output.sortbam};
			samtools index {output.sortbam} {output.idx};
			bedtools bamtobed -i {output.sortbam} > {output.bed};
			awk 'function abs(a){{return ((a < 0) ? -a : a)}} OFS="\\t" {{split($4,a,"x");print $0,0,a[2]}}' {output.bed} > {output.ISsortednum};
			awk 'OFS="\\t" {{if($6=="-") {{print $1,$3-1,$3,$4,$5,$6,$7,$8}} else {{print $1,$2,$2+1,$4,$5,$6,$7,$8}}}}' {output.ISsortednum} > {output.IS};
			bedtools sort -i {output.IS} > {output.ISsorted};
			bedtools cluster -s  -i {output.ISsorted} > {output.IScluster};
			sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -k 7,7n {output.IScluster} > {output.ISclustersorted};
			bedtools groupby -i {output.ISclustersorted}  -g 9,7 -c 1,2,3,4,5,6,7,8,8,9 -o distinct,mode,mode,collapse,median,distinct,distinct,sum,count,distinct  | cut -f 3- > {output.IScollapsedSize};
			bedtools groupby -i {output.IScollapsedSize}  -g 10 -c 1,2,3,4,5,6,8,7 -o distinct,mode,mode,collapse,median,distinct,sum,count_distinct | cut -f 2-  > {output.IScollapsed}
			"""
		
		
		
# convert chromosom name from 'x' to 'chrx' for compatibility.
# sort merged bed file and collapse IS at the same position.
	
rule QualIS:
	input: R1full=rules.Filter_Map_R1.output.IScollapsed, 
		R1alone=rules.Filter_Map_R1_noLinker.output.IScollapsed, 
		R1R2=rules.Filter_Map_R1R2.output.IScollapsed
	output: merged=("09-qualitativeIS/{NAME}_TAG{TAG}_qualIS.bed"), 
			merged_corrected="09-qualitativeIS/{NAME}_TAG{TAG}_qualIS_merged.bed",
			merged_sorted = "09-qualitativeIS/{NAME}_TAG{TAG}_qualIS_merged_sorted.bed",
			merged_sorted_collapsed= "09-qualitativeIS/{NAME}_TAG{TAG}_qualIS_merged_sorted_collapsed.bed"
	message: "Merging and collapsing IS from each step"
	threads: 8
	log: "log/qualitativeIS.log"
	shell: """
			cat  {input.R1full} {input.R1alone} {input.R1R2} > {output.merged};
			awk '{{print "chr"$0}}' {output.merged} > {output.merged_corrected};
			sort -k1,1 -k2,2n {output.merged_corrected} > {output.merged_sorted};
			bedtools merge -s -d -1 -c 1,2,3,4,5,6,7 -o distinct,mode,mode,collapse,median,distinct,sum -i {output.merged_sorted} | cut -f 5- | bedtools sort -i - > {output.merged_sorted_collapsed};
			"""
	
			

#rule annotateIS_:
#	input: rules.QualIS.output.merged_corrected
#	output: slop="13-qualitativeIS/{NAME}_TAG{TAG}_qualIS_merged_slop400.bed",slopsorted="13-qualitativeIS/{NAME}_TAG{TAG}_qualIS_merged_slop400_sorted.bed",bigCoverage = "13-qualitativeIS/{NAME}_TAG{TAG}_qualIS_H3K27me3.stab",bedcoverage = "13-qualitativeIS/{NAME}_TAG{TAG}_qualIS_H3K27me3.bed"
#	log: "log/H3K27me3.log"
#	shell:
#		"""
#		bedtools slop -i {input} -g ../references/hg19.chrom.sizes -b 200 | cut -f 1-4 > {output.slop};
#		sort -k1,1 -k2,2n -k6,6 {output.slop} > {output.slopsorted};
#		#bigWigAverageOverBed ../dataset/epigenetics/GSM621664_BI.Mobilized_CD34_Primary_Cells.H3K27me3.UW_RO_01536.bw {output.slopsorted} {output.bigCoverage} -stats={log} -bedOut={output.bedcoverage} -minMax
#		"""
	
rule QuantIS:
	input:R1full=rules.Filter_Map_R1.output.IScollapsed,
		R1R2=rules.Filter_Map_R1R2.output.IScollapsed
	output: merged="09-quantitativeIS/{NAME}_TAG{TAG}_quantIS_merged.bed", 
			merged_corrected="09-quantitativeIS/{NAME}_TAG{TAG}_quantIS_merged_corrected.bed",
			collapsed = "09-quantitativeIS/{NAME}_TAG{TAG}_quantIS_collapsed.bed"
	threads: 1
	shell: """
			cat {input.R1full} {input.R1R2} > {output.merged};
			awk '{{print "chr"$0}}' {output.merged} > {output.merged_corrected};
			bedtools sort -i {output.merged_corrected} | bedtools merge -d -1 -s -c 1,2,3,4,5,6,7,8 -o distinct,mode,mode,collapse,median,distinct,sum,count_distinct -i - | cut -f 5- > {output.collapsed}
			"""
			
			
			
			
############################################################################
############################################################################
############################################################################
############################################################################
# Download annotation files from different sources	: UCSC known genes, Ensembl genes, Refseq Genes, GencodeV19 and the HGNC gene symbol database
#
############################################################################
	
rule GetAnnotations_UCSC:
	input: rules.Prepare_folders.output.annot
	output: folder = ANNOTATION+"UCSC/",ChrSize= ANNOTATION+"UCSC/hg19.chrom.sizes",genes=ANNOTATION+"UCSC/knownGene.txt"
	log: "log/UCSC_remote_wget.log"
	threads: 1
	message: "Downloading annotation from UCSC"
	shell : """
			xargs -i wget -nc -P {output.folder} '{{}}'  < ../Scripts/hg19_annotationUCSC_urls.txt 2> {log};
			gunzip {output.folder}/*.gz
			"""
			
rule GetAnnotations_HGNC:
	input:  rules.Prepare_folders.output.annot
	output: HGNC=ANNOTATION+"HGNC/HGNC.txt"
	threads: 1
	message: "Downloading annotation from HGNC"
	shell : """
			perl ../Scripts/retreive_HUGO.pl > {output.HGNC};
			"""
			
rule GetAnnotations_Ensembl:	
	input: rules.Prepare_folders.output.annot
	output: folder =ANNOTATION+"Ensembl/", ChrSize= ANNOTATION+"Ensembl/hg19.chrom.sizes", genes=ANNOTATION+"Ensembl/ensGene.txt"
	threads: 1
	log: "log/Ensembl_remote_wget.log"
	message: "Downloading annotation from Ensembl"
	shell:"""
			xargs -i wget -nc -P {output.folder} '{{}}'  < ../Scripts/hg19_annotationEnsembl_urls.txt 2> {log};
			gunzip {output.folder}/*.gz
			"""

rule GetAnnotations_RefSeq:
	input:  rules.Prepare_folders.output.annot
	output: folder = ANNOTATION+"RefSeq/", genes = ANNOTATION+"RefSeq/refFlat.txt", ChrSize= ANNOTATION+"RefSeq/hg19.chrom.sizes"
	threads: 1
	log: "log/RefSeq_remote_wget.log"
	message: "Downloading annotation from RefSeq"
	shell: """
			 xargs -i wget -nc -P {output.folder} '{{}}'  < ../Scripts/hg19_annotationRefSeq_urls.txt 2> {log};
			 gunzip {output.folder}/*.gz
			"""
			
			
rule GetAnnotations_GenCode:
	input:  rules.Prepare_folders.output.annot
	output: folder = ANNOTATION+"GencodeV19/"
	threads: 1
	log: "log/Gencode_remote_wget.log"
	message: "Downloading annotation from GencodeV19"
	shell: """
			xargs -i wget -nc -P {output.folder} '{{}}'  < ../Scripts/hg19_annotationGencodeV19_urls.txt 2> {log};
			gunzip {output.folder}/*.gz
			"""	
			
############################################################################
############################################################################
############################################################################
############################################################################
# Convert refFlat table to bed, convert exon start-end to bed with exon rank in transcript numbered according to strand !!!!!

						
rule FormatAnnotationRefSeq:
	input: rules.GetAnnotations_RefSeq.output.folder,
			refFlat=rules.GetAnnotations_RefSeq.output.genes, 
			gensize= rules.GetAnnotations_RefSeq.output.ChrSize
	output: reflatBed=ANNOTATION+"RefSeq/refFlat.bed",
			reflatTSS=ANNOTATION+"RefSeq/refFlatTSS.bed",
			reflatTSSSorted=ANNOTATION+"RefSeq/refFlatTSSsorted.bed",
			reflatCore=ANNOTATION+"RefSeq/refFlatCoreProm.bed",
			reflatProx=ANNOTATION+"RefSeq/refFlatProxProm.bed",
			reflatDist=ANNOTATION+"RefSeq/refFlatDistProm.bed",
			reflatExons=ANNOTATION+"RefSeq/refFlatExons.bed"
	threads: 1
	log:
	message: "From RefSeq refFlat table, convert to BED, extract TSS, core, proximal and distal promoter"
	shell: """
			awk 'OFS="\\t" {{print $3,$5,$6,$2,0,$4}}' {input.refFlat} | bedtools sort -i - > {output.reflatBed};
			awk 'OFS="\\t" {{if($6 == "+") {{print $1,$2,$2+1,$4,$5,$6}} else {{print $1,$3-1,$3,$4,$5,$6}}}}' {output.reflatBed} > {output.reflatTSS};
			bedtools sort -i {output.reflatTSS} > {output.reflatTSSSorted};
			bedtools flank -s -l 40 -r 0 -g {input.gensize} -i {output.reflatTSSSorted} > {output.reflatCore};
			bedtools flank -s -l 160 -r 0 -g {input.gensize} -i {output.reflatCore} > {output.reflatProx};
			bedtools flank -s -l 800 -r 0 -g {input.gensize} -i {output.reflatProx} > {output.reflatDist};
			awk 'OFS="\\t" {{split($10,a,",");split($11,b,","); for (i = 1; i <= $9; ++i) {{if($4=="-") {{k=$9-i+1}} else {{k=i}}; print $3,a[i],b[i],$2"_exon_"k,$1,$4}}}}' {input.refFlat} | bedtools sort -i - > {output.reflatExons}
			"""
	
rule FormatAnnotationUCSC:
	input: rules.GetAnnotations_UCSC.output.folder,
			known=rules.GetAnnotations_UCSC.output.genes,
			gensize= rules.GetAnnotations_UCSC.output.ChrSize
	output: UCSCBed=ANNOTATION+"UCSC/knownGenes.bed",
			UCSCTSS=ANNOTATION+"UCSC/knownGeneTSS.bed",
			UCSCTSSSorted=ANNOTATION+"UCSC/UCSCTSSsorted.bed",
			UCSCCore=ANNOTATION+"UCSC/UCSCCoreProm.bed",
			UCSCProx=ANNOTATION+"UCSC/UCSCProxProm.bed",
			UCSCDist=ANNOTATION+"UCSC/UCSCDistProm.bed",
			UCSCExons=ANNOTATION+"UCSC/UCSCExons.bed"
			
	threads: 1
	shell: """
			awk 'OFS="\\t" {{print $2,$4,$5,$1,0,$3}}' {input.known} | bedtools sort -i - > {output.UCSCBed};
			awk 'OFS="\\t" {{if($6 == "+") {{print $1,$2,$2+1,$4,$5,$6}} else {{print $1,$3-1,$3,$4,$5,$6}}}}' {output.UCSCBed} > {output.UCSCTSS};
			bedtools sort -i {output.UCSCTSS} > {output.UCSCTSSSorted};
			bedtools flank -s -l 40 -r 0 -g {input.gensize} -i {output.UCSCTSSSorted} > {output.UCSCCore};
			bedtools flank -s -l 160 -r 0 -g {input.gensize} -i {output.UCSCCore} > {output.UCSCProx};
			bedtools flank -s -l 800 -r 0 -g {input.gensize} -i {output.UCSCProx} > {output.UCSCDist};
			awk 'OFS="\\t" {{split($9,a,",");split($10,b,","); for (i = 1; i <= $8; ++i) {{if($3=="-") {{k=$8-i+1}} else {{k=i}}; print $2,a[i],b[i],$1"_exon_"k,$1,$3}}}}' {input.known} | bedtools sort -i - > {output.UCSCExons}
			"""
#join  -16 -21 -e NULL -t $'\t' <(sort -k 6,6 knownCanonical.txt) <(sort -k 1,1 kgXref.txt) > annotatedUCSC.txt -a1




rule FormatAnnotationHGNC:
	input: rules.GetAnnotations_HGNC.output
	output: touch(ANNOTATION+"HGNC/temp.bed")
	threads: 1
	
	
	
	
rule FormatAnnotationEnsembl:
	input: rules.GetAnnotations_Ensembl.output.folder,
			EnsGenes=rules.GetAnnotations_Ensembl.output.genes,
			gensize= rules.GetAnnotations_Ensembl.output.ChrSize
	output: EnsBed=ANNOTATION+"Ensembl/EnsGenes.bed",
			EnsTSS=ANNOTATION+"Ensembl/EnsGeneTSS.bed",
			EnsTSSSorted=ANNOTATION+"Ensembl/EnsTSSsorted.bed",
			EnsCore=ANNOTATION+"Ensembl/EnsCoreProm.bed",
			EnsProx=ANNOTATION+"Ensembl/EnsProxProm.bed",
			EnsDist=ANNOTATION+"Ensembl/EnsDistProm.bed",
			EnsExons=ANNOTATION+"Ensembl/EnsExons.bed"
			
	threads: 1
	shell: """
			awk 'OFS="\\t" {{print $3,$5,$6,$2,$13,$4}}' {input.EnsGenes} | bedtools sort -i - > {output.EnsBed};
			awk 'OFS="\\t" {{if($6 == "+") {{print $1,$2,$2+1,$4,$5,$6}} else {{print $1,$3-1,$3,$4,$5,$6}}}}' {output.EnsBed} > {output.EnsTSS};
			bedtools sort -i {output.EnsTSS} > {output.EnsTSSSorted};
			bedtools flank -s -l 40 -r 0 -g {input.gensize} -i {output.EnsTSSSorted} > {output.EnsCore};
			bedtools flank -s -l 160 -r 0 -g {input.gensize} -i {output.EnsCore} > {output.EnsProx};
			bedtools flank -s -l 800 -r 0 -g {input.gensize} -i {output.EnsProx} > {output.EnsDist};
			awk 'OFS="\\t" {{split($10,a,",");split($11,b,","); for (i = 1; i <= $9; ++i) {{if($4=="-") {{k=$9-i+1}} else {{k=i}}; print $3,a[i],b[i],$2"_exon_"k,$13,$4}}}}' {input.EnsGenes} | bedtools sort -i - > {output.EnsExons}
			
			"""	

rule FormatAnnotationGencode:
	input: rules.GetAnnotations_GenCode.output.folder
	output: touch(ANNOTATION+"GencodeV19/temp.bed")	
	threads: 1
				
		
############################################################################
############################################################################
############################################################################
############################################################################
# Annotate IS with different databases
	
rule AnnotateISRefSeq:
	input: IS=rules.QualIS.output.merged_sorted_collapsed,
			genes = rules.FormatAnnotationRefSeq.output.reflatBed,
			TSS = rules.FormatAnnotationRefSeq.output.reflatTSSSorted,
			CoreP=rules.FormatAnnotationRefSeq.output.reflatCore,
			ProxP=rules.FormatAnnotationRefSeq.output.reflatProx,
			DistP=rules.FormatAnnotationRefSeq.output.reflatDist,
			Exons=rules.FormatAnnotationRefSeq.output.reflatExons
			
	output: all="10-Annotation/RefSeq/{NAME}_TAG{TAG}_IntragenicIS.txt",  
			distTSS = "10-Annotation/RefSeq/{NAME}_TAG{TAG}_distTSS_IS.txt",
			IntronExon = "10-Annotation/RefSeq/{NAME}_TAG{TAG}_IntronExonIS.txt",
			IntraGeneIS = "10-Annotation/RefSeq/{NAME}_TAG{TAG}_IntragenicIS.bed"
	threads: 1
	shell: """
			bedtools intersect -a {input.IS} -b {input.genes} {input.CoreP} {input.ProxP} {input.DistP} -wa -wb -loj -filenames > {output.all};
			bedtools closest -a {input.IS} -b {input.TSS} -D b > {output.distTSS};
			awk 'BEGIN{{OFS="\\t"}} /refFlat.bed/ {{print $1,$2,$3,$4,$5,$6,$7}} ' {output.all} | sort -u > {output.IntraGeneIS};
			bedtools intersect -a {output.IntraGeneIS} -b {input.Exons} -wa -wb -loj -filenames > {output.IntronExon};
			"""

			
rule AnnotateISUCSC:
	input: IS=rules.QualIS.output.merged_sorted_collapsed,
			genes = rules.FormatAnnotationUCSC.output.UCSCBed,
			TSS = rules.FormatAnnotationUCSC.output.UCSCTSSSorted,
			CoreP=rules.FormatAnnotationUCSC.output.UCSCCore,
			ProxP=rules.FormatAnnotationUCSC.output.UCSCProx,
			DistP=rules.FormatAnnotationUCSC.output.UCSCDist,
			Exons=rules.FormatAnnotationUCSC.output.UCSCExons
			
	output: IntraGenes="10-Annotation/UCSC/{NAME}_TAG{TAG}_IntragenicIS.txt", 
			InterGenes="10-Annotation/UCSC/{NAME}_TAG{TAG}_IntergenicIS.txt", 
			distTSS = "10-Annotation/UCSC/{NAME}_TAG{TAG}_distTSS_IS.txt",
			IntronExon = "10-Annotation/UCSC/{NAME}_TAG{TAG}_IntronExonIS.txt",
			InterProm = "10-Annotation/UCSC/{NAME}_TAG{TAG}_PromoterInter.txt",
			IntraGeneIS = "10-Annotation/UCSC/{NAME}_TAG{TAG}_IntragenicIS.bed"	
	threads: 1
	shell: """
			bedtools intersect -a {input.IS} -b {input.genes} -wa -wb -filenames > {output.IntraGenes};
			bedtools intersect -a {input.IS} -b {input.genes} -wa -v -wb -filenames > {output.InterGenes};
			bedtools closest -a {input.IS} -b {input.TSS} -D b > {output.distTSS};
			cut -f 1-7 {output.IntraGenes} | sort -u > {output.IntraGeneIS};
			bedtools intersect -a {output.IntraGeneIS} -b {input.Exons} -wa -wb -loj -filenames > {output.IntronExon};
			bedtools intersect -a {output.InterGenes} -b {input.CoreP} {input.ProxP} {input.DistP} -wa -wb -loj -filenames > {output.InterProm};
			"""

rule AnnotateISEnsembl:
	input: IS=rules.QualIS.output.merged_sorted_collapsed,
			genes = rules.FormatAnnotationEnsembl.output.EnsBed,
			TSS = rules.FormatAnnotationEnsembl.output.EnsTSSSorted,
			CoreP=rules.FormatAnnotationEnsembl.output.EnsCore,
			ProxP=rules.FormatAnnotationEnsembl.output.EnsProx,
			DistP=rules.FormatAnnotationEnsembl.output.EnsDist,
			Exons=rules.FormatAnnotationEnsembl.output.EnsExons
			
	output: IntraGenes="10-Annotation/Ensembl/{NAME}_TAG{TAG}_IntragenicIS.txt", 
			InterGenes="10-Annotation/Ensembl/{NAME}_TAG{TAG}_IntergenicIS.txt", 
			distTSS = "10-Annotation/Ensembl/{NAME}_TAG{TAG}_distTSS_IS.txt",
			IntronExon = "10-Annotation/Ensembl/{NAME}_TAG{TAG}_IntronExonIS.txt",
			InterProm = "10-Annotation/Ensembl/{NAME}_TAG{TAG}_PromoterInter.txt",
			IntraGeneIS = "10-Annotation/Ensembl/{NAME}_TAG{TAG}_IntragenicIS.bed"	
	threads: 1
	shell: """
			bedtools intersect -a {input.IS} -b {input.genes} -wa -wb -filenames > {output.IntraGenes};
			bedtools intersect -a {input.IS} -b {input.genes} -wa -v -wb -filenames > {output.InterGenes};
			bedtools closest -a {input.IS} -b {input.TSS} -D b > {output.distTSS};
			cut -f 1-7 {output.IntraGenes} | sort -u > {output.IntraGeneIS};
			bedtools intersect -a {output.IntraGeneIS} -b {input.Exons} -wa -wb -loj -filenames > {output.IntronExon};
			bedtools intersect -a {output.InterGenes} -b {input.CoreP} {input.ProxP} {input.DistP} -wa -wb -loj -filenames > {output.InterProm};
			"""

rule AnnotateISGencode:
	input: IS=rules.QualIS.output.merged_sorted_collapsed, DB=rules.FormatAnnotationGencode.output
	output: touch("10-Annotation/Gencode/{NAME}_TAG{TAG}_reportGencode.pdf")
	threads: 1	
	
##########################################################################


rule QC_stat:
	input: rules.Demultiplex_TAGS.output.R1, 
			rules.Demultiplex_TAGS.output.R2,
			rules.Find_Linker_R1.output.R1Linker,
			rules.Find_Linker_R1.output.R2Linker,
			rules.Find_Elong.output.R1elong,
			rules.Find_Elong.output.R2elong,
			rules.Find_Linker_R2.output.R1Linker,
			rules.Find_Linker_R2.output.R1Linker,
			rules.Map_R1_R2_pairs.output.mapped,
			rules.Map_R1_noLinker.output.mapped,
			rules.Map_R1_with_Linker.output.mapped
	output: QC="QC_{NAME}_TAG{TAG}/"
	message: "Gathering statistics from sam and fastq files"
	shell: """
			mkdir {output};
			find -type f -regex ".+\.sam$" -exec  sh -c 'samtools flagstat {{}} > {{}}.samstat' \;
			find -type f -regex ".+\.fastq\(.gz\)?$" -exec sh -c 'fastq-stats -s 0 {{}} > {{}}.fqstat' \;
			"""
			
#find -type f -regex ".+\.fastq.gz$" -exec sh -c 'fastq-stats -s 0 {{}} > {{}}.fqstat' \;
#fastqc -t 8 {input} --adapters ../dataset/contaminants/adapter_list.txt -o {output.QC};
#multiqc . -f;


rule MakeReport:
	input:  rules.QuantIS.output.collapsed, rules.QualIS.output.merged_sorted_collapsed, rules.QC_stat.output
	output: touch("{NAME}_TAG{TAG}_report.pdf")
	threads: 1
	
onsuccess:
    print("Workflow finished, no error")