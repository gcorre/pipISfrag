#	This code is a snakemake pipeline integrating different tools to process Illumina 
#	sequencing paired-end reads for vector Insertion site calling. It provides
#	file tracking and input/output management.
#

#	Copyright (C) 2017 GENETHON, Guillaume Corre
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

print("#~~~~~~~~Version 03/15/2017~~~~~~~~#");
print("");
print("#       #    ##    ######    #######");
print(" #     #     ##    ##        ##   ##");
print("  #   #      ##    ######    #######");
print("   # #       ##        ##    ##   ##");
print("    #        ##    ######    ##   ##");
print("");
print("#~~~~~~~~~~~~~GENETHON~~~~~~~~~~~~~#");
	
#workdir: "/home/"

BOWTIE2_INDEX = "/home/references/genomes/homo_sapiens/hg19_GRCh37/index_bowtie2/GRCh37.75"

from datetime import date
TODAY=str(date.today())

from glob import glob;
FILES = glob('*_R[0-9]_*.fastq.gz');
import re;
TRAILING=list(set([re.search('_([0-9]+).fastq.gz',FILES[0]).group(1) for w in FILES]))[0]
FILES=[re.sub('_R[1-2]_[0-9]+.fastq.gz', '',w) for w in FILES];
SAMPLE=list(set(FILES))[0];
TAG=re.search('TAG([0-9]+)', SAMPLE).group(1)

rule targets:
	input: "13-Report/"+SAMPLE+"_TAG"+TAG+"_report.pdf","../dataset/UCSC/hg19/"+TODAY+"/"

rule Demultiplex_TAGS:
	input: R1="{NAME}_R1_"+TRAILING+".fastq.gz", R2="{NAME}_R2_"+TRAILING+".fastq.gz"
	output: R1="01-Demultiplex/{NAME}_R1_TAG{TAG}.fastq.gz",R2="01-Demultiplex/{NAME}_R2_TAG{TAG}.fastq.gz"
	log: "log/demultiplex-merged.log"
	threads: 8
	message: "Demultiplexing according to TAGs list"
	shell: """
			 fastq-multx -B ../dataset/illuminaTAGs_eautils.txt -b -x -m 1 {input.R1} {input.R2} -o 01-Demultiplex/{wildcards.NAME}_R1_%.fastq.gz 01-Demultiplex/{wildcards.NAME}_R2_%.fastq.gz > {log};
			"""
			
rule Find_provirus:
	input: R1=rules.Demultiplex_TAGS.output.R1, R2=rules.Demultiplex_TAGS.output.R2
	output: R1Prov="02-provirus/{NAME}_R1_TAG{TAG}_Prov.fastq", R2Prov="02-provirus/{NAME}_R2_TAG{TAG}_Prov.fastq", R1noProv="02-provirus/{NAME}_R1_TAG{TAG}_noProv.fastq", R2noProv="02-provirus/{NAME}_R2_TAG{TAG}_noProv.fastq",
	log:"log/Provirus.log"
	threads: 8
	message: "Look for the Provirus in R1"
	shell: """
			cutadapt -g aaaatctctagcagtggcgcccgaacag -O 28 -e 0.1 --no-trim --no-indels -o {output.R1Prov} -p {output.R2Prov} --untrimmed-paired-output {output.R2noProv} --untrimmed-output {output.R1noProv} {input.R1} {input.R2} > {log}
			"""
			
rule Find_LTR:
	input: R1=rules.Find_provirus.output.R1noProv, R2=rules.Find_provirus.output.R2noProv
	output: R1LTR="03-LTR/{NAME}_R1_TAG{TAG}_LTR.fastq", R2LTR="03-LTR/{NAME}_R2_TAG{TAG}_LTR.fastq", R1noLTR="03-LTR/{NAME}_R1_TAG{TAG}_noLTR.fastq", R2noLTR="03-LTR/{NAME}_R2_TAG{TAG}_noLTR.fastq",
	log:"log/LTR.log"
	threads: 8
	message: "Look for the LTR in R1"
	shell: """
			cutadapt -g GTCTGTTGTGTGACTCTGGTAAC -m 20 -O 23 -e 0.1 --no-indels -o {output.R1LTR} -p {output.R2LTR} --untrimmed-paired-output {output.R2noLTR} --untrimmed-output {output.R1noLTR} {input.R1} {input.R2} > {log}
			"""

rule Find_Elong:
	input: R1=rules.Find_LTR.output.R1LTR, R2=rules.Find_LTR.output.R2LTR
	output: R1elong="04-ELONG/{NAME}_R1_TAG{TAG}_Elong.fastq", R2elong="04-ELONG/{NAME}_R2_TAG{TAG}_Elong.fastq", R1noelong="04-ELONG/{NAME}_R1_TAG{TAG}_noElong.fastq", R2noelong="04-ELONG/{NAME}_R2_TAG{TAG}_noElong.fastq",
	log:"log/ELONG.log"
	threads: 8
	message: "Look for the Elong in R1"
	shell: """
			cutadapt -g TAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTA -m 20 -O 43 -e 0.16 --no-indels -o {output.R1elong} -p {output.R2elong} --untrimmed-paired-output {output.R2noelong} --untrimmed-output {output.R1noelong} {input.R1} {input.R2} > {log}
			"""

rule Start_with_GCA:
	input: R1=rules.Find_Elong.output.R1elong, R2=rules.Find_Elong.output.R2elong
	output: GCAR1="09-GCA/{NAME}_R1_TAG{TAG}_GCA.fastq", noGCAR1="09-GCA/{NAME}_R1_TAG{TAG}_noGCA.fastq",GCAR2="09-GCA/{NAME}_R2_TAG{TAG}_GCA.fastq", noGCAR2="09-GCA/{NAME}_R2_TAG{TAG}_noGCA.fastq"
	message: "Filtering R1 reads that start with GCA"
	log: "log/GCA-R1.log"
	shell: """
			cutadapt -g ^GCA  -m 20 -e 0 --no-indels -o {output.GCAR1} -p {output.GCAR2} --untrimmed-output {output.noGCAR1} --untrimmed-paired-output {output.noGCAR2} {input.R1} {input.R2} > {log};
			"""			

rule Find_Linker_R1:
	input: R1=rules.Start_with_GCA.output.GCAR1, R2=rules.Start_with_GCA.output.GCAR2
	output: R1Linker="05-Linker/{NAME}_R1_TAG{TAG}_Linker.fastq", R2Linker="05-Linker/{NAME}_R2_TAG{TAG}_Linker.fastq", R1noLinker="05-Linker/{NAME}_R1_TAG{TAG}_noLinker.fastq", R2noLinker="05-Linker/{NAME}_R2_TAG{TAG}_noLinker.fastq",
	log:"log/LinkerR1.log"
	threads: 8
	message: "Look for the Linker in R1"
	shell: """
			cutadapt -a GTCCCTTAAGCGGAGCCCT -m 20 -O 8 -e 0.2 --no-indels -o {output.R1Linker} -p {output.R2Linker} --untrimmed-paired-output {output.R2noLinker} --untrimmed-output {output.R1noLinker} {input.R1} {input.R2} > {log}
			"""	
 	
rule Find_Linker_R2:
	input: R1=rules.Find_Linker_R1.output.R1noLinker, R2=rules.Find_Linker_R1.output.R2noLinker
	output: R1Linker="06-LinkerR2/{NAME}_R1_TAG{TAG}_Linker.fastq", R2Linker="06-LinkerR2/{NAME}_R2_TAG{TAG}_Linker.fastq", R1noLinker="06-LinkerR2/{NAME}_R1_TAG{TAG}_noLinker.fastq", R2noLinker="06-LinkerR2/{NAME}_R2_TAG{TAG}_noLinker.fastq"
	log:"log/LinkerR2.log"
	threads: 8
	message: "Look for the Linker in R2"
	shell: """
			cutadapt -g AGGGCTCCGCTTAAGGGAC -O 19 -e 0.16 --no-indels -o {output.R2Linker} -p {output.R1Linker} --untrimmed-paired-output {output.R1noLinker} --untrimmed-output {output.R2noLinker} {input.R2} {input.R1} > {log}
			"""				

#######################################################################################			
## Process R1 reads that contain the Linker:
#######################################################################################

####
# This step will be removed when using the sonication method. Just here to compare with VISA-MseI
####
#rule Cut_TTAA:
#	input: rules.Find_Linker_R1.output.R1Linker
#	output: TTAA="08-TTAA_cut/{NAME}_R1_TAG{TAG}_TTAA-cut.fastq",TTAA_cut="08-TTAA_cut/{NAME}_R1_TAG{TAG}_TTAA-cut20bp.fastq"
#	message: "Cutting TTAA sequences that are still present and remove sequences that are too short"
#	threads: 8
#	log: "log/TTAA.log"
#	shell: """
#			  awk '{{if(NR %4 ==2) {{x=index($0,"TTAA");if(x>0){{print substr($0,1,x)}}else{{print $0}}}} else {{if(NR % 4 ==0 && x>0) {{print substr($0,1,x)}} else{{print $0}}}}}}' {input} > {output.TTAA}
#			  cutadapt -m 20 -o {output.TTAA_cut} {output.TTAA} > {log}
#			"""			

			
rule Collapse_identical_reads:
	input: rules.Find_Linker_R1.output.R1Linker
	output: file = "10-collapsed/{NAME}_R1_TAG{TAG}_Linker_trimmed.fastq"
	message: "Collapsing identical reads before mapping"
	log:"log/collapsingR1.log"
	shell: """
			#fastx_collapser -v -i {input} -o {output} > {log}
			 seqcluster collapse -f {input} -d -o "10-collapsed/" -m 1
			"""

####
# Map reads containing the linker in R1, remove reads with multiple hits (XS:i: TAG) and reads with mismatches in the first 3 positions after the LTR (taking into account the orientation)
# Keep the header in sam file
#### 	
		
rule Map_R1_with_Linker :
	input: rules.Collapse_identical_reads.output.file
	output: mapped="11-mapping/{NAME}_R1_TAG{TAG}_mapped.sam",unmapped="11-mapping/{NAME}_R1_TAG{TAG}_unmapped.fastq", exact="11-mapping/{NAME}_R1_TAG{TAG}_mapped_exact.sam"
	log: met="log/mappingR1.log",summary="log/mapping_numbers.log"
	threads: 8
	message: "Mapping R1 reads that contain the linker and filter out those with mismatches in the 3 first positions after the LTR"
	shell: """
			bowtie2 -N 1 -L 25 -i S,25,0 --score-min L,0,-0.15 --gbar 10 -p {threads} -x {BOWTIE2_INDEX} --no-unal --un {output.unmapped} --met-file {log.met} {input} -S {output.mapped} 2> {log.summary}
			awk -F "\\t" '(/^@/) || ($2==0 && !/XS:i:/ && !/MD:Z:[012][A-Za-z].*\t/) || ($2==16 && !/XS:i:/ && !/MD:Z:.*[A-Za-z][012]\t/) {{print $0}}' {output.mapped} > {output.exact}
			"""

			
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
	output: bam=temp("12-Annotation/{NAME}_R1_TAG{TAG}_mapped.bam"),
			sortbam="12-Annotation/{NAME}_R1_TAG{TAG}_mapped_sorted.bam",
			idx="12-Annotation/{NAME}_R1_TAG{TAG}_mapped_sorted.bai",
			bed="12-Annotation/{NAME}_R1_TAG{TAG}_mapped_sorted.bed",
			ISsortednum="12-Annotation/{NAME}_R1_TAG{TAG}_IS_sortedNum.bed",
			IS=temp("12-Annotation/{NAME}_R1_TAG{TAG}_IS.bed"),
			ISsorted="12-Annotation/{NAME}_R1_TAG{TAG}_IS_sorted.bed",
			IScollapsedSize="12-Annotation/{NAME}_R1_TAG{TAG}_IS_collapsedBySizeCluster.bed",
			IScluster="12-Annotation/{NAME}_R1_TAG{TAG}_IS_cluster.bed",
			ISclustersorted="12-Annotation/{NAME}_R1_TAG{TAG}_IS_cluster_sorted.bed",
			IScollapsed="12-Annotation/{NAME}_R1_TAG{TAG}_IS_Collapsed.bed"
	params: qual="10", collapsing_window="-1"
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
			

#################################################################################
# Process paired reads with the linker in R2 but not in R1	
# Amplicon may be too long to include the complete linker in R1.
# R1 was already trimmed before for LTR and elongation. R2 was already trimmed for the linker before.
####
 
# look for the reverse-complement Elongation in R2 3'end.
rule Find_RC_LTR_R2:
	input: R1=rules.Find_Linker_R2.output.R1Linker,R2=rules.Find_Linker_R2.output.R2Linker
	output: trimmedLTRRCR1="07-cleaning2/{NAME}_R1_TAG{TAG}_trimmedElongRC.fastq", 
			trimmedLTRRCR2="07-cleaning2/{NAME}_R2_TAG{TAG}_trimmedElongRC.fastq"
	log: "log/RCLTR_R2.log"
	threads: 8
	message: "Trimming R2 that have LTR"
	shell: """
			cutadapt -a TGCTAGAGATTTTCCACACTGACTAAAAGGGTCTGAGGGATCTCTA -q 20 -O 8 -m 20 -e 0.25 -o {output.trimmedLTRRCR2} -p {output.trimmedLTRRCR1} {input.R2} {input.R1} > {log};
			"""	
			
####
# This step will be removed when moving to the sonication method.
		
#rule Cut_TTAA_in_R1R2:
	#input: R1=rules.Find_RC_LTR_R2.output.trimmedLTRRCR1,R2=rules.Find_RC_LTR_R2.output.trimmedLTRRCR2
#	output: TTAAR1="10-TTAA_cut/{NAME}_R1_TAG{TAG}_TTAA-cut.fastq", 
	#		TTAAR2="10-TTAA_cut/{NAME}_R2_TAG{TAG}_TTAA-cut.fastq",
#			qualR1="10-TTAA_cut/{NAME}_R1_TAG{TAG}.fastq",
#			qualR2="10-TTAA_cut/{NAME}_R2_TAG{TAG}.fastq"
			
#	message: "Cutting TTAA sequences that are still present, trimming low quality 3'"
#	threads: 8
	#log: "log/TTAAR1R2_quality.log"
#	shell: """
#			  awk '{{if(NR %4 ==2) {{x=index($0,"TTAA");if(x>0){{print substr($0,1,x)}}else{{print $0}}}} else {{if(NR % 4 ==0 && x>0) {{print substr($0,1,x)}} else{{print $0}}}}}}' {input.R1} > {output.TTAAR1};
#			  awk '{{if(NR %4 ==2) {{x=index($0,"TTAA");if(x>0){{print substr($0,1,x)}}else{{print $0}}}} else {{if(NR % 4 ==0 && x>0) {{print substr($0,1,x)}} else{{print $0}}}}}}' {input.R2} > {output.TTAAR2};
#			  cutadapt -q 20 -m 20 -A XXXXXXX -o {output.qualR1} -p {output.qualR2} {output.TTAAR1} {output.TTAAR2} > {log};
#			"""	

####################################################
# Mapped paired end reads, trim 3' end to remove uncut linker parts, tolerate dovetail and max fragments length 800bp. Keep only concordant mapping (same Chr, <800bp).		
	
rule Map_R1_R2_pairs:
	input: R1=rules.Find_RC_LTR_R2.output.trimmedLTRRCR1, R2=rules.Find_RC_LTR_R2.output.trimmedLTRRCR2
	output: mapped="11-mappingR1R2/{NAME}_TAG{TAG}_mapped.sam", R1unal="11-mappingR1R2/{NAME}_R1_TAG{TAG}_unaligned.fastq"
	message: "Mapping R1 and R2 reads as pairs"
	params: fragLenght="800", Trim3="5"
	log: met="log/mappingR1R2.log",summary="log/mappingR1R2_numbers.log"
	threads: 8
	shell: """
			bowtie2 -p {threads} -x {BOWTIE2_INDEX} -1 {input.R1} -2 {input.R2} -S {output.mapped} --trim3 {params.Trim3} --no-unal --al-conc "11-mappingR1R2/{wildcards.NAME}_R%_TAG{wildcards.TAG}_aligned.fastq" --un-conc "11-mappingR1R2/{wildcards.NAME}_R%_TAG{wildcards.TAG}_unaligned.fastq" --dovetail -X {params.fragLenght} --no-mixed --met-file {log.met} 2> {log.summary};
			"""

#########################################################################
# Filter out pairs of reads for which the mapping quality score is below 20, meaning p= 0.01 that the reported alignment may occured elsewhere.
# Then keep only R1 reads from filtered out reads and co
# 
			
rule Filter_Map_R1R2_unique:
	input: rules.Map_R1_R2_pairs.output.mapped
	output: badqual = "11-mappingR1R2/{NAME}_TAG{TAG}_bad_qual.sam",confident="11-mappingR1R2/{NAME}_TAG{TAG}_qual_gt20.sam",R1_badqual= "11-mappingR1R2/{NAME}_R1_TAG{TAG}_qual_lt20.sam" ,R1_badqualfq= "11-mappingR1R2/{NAME}_R1_TAG{TAG}_qual_lt20.fastq"
	params: qual="20"
	message: "filtering bad quality alignments that may be multiple hits, keep R1 reads and convert to fastq."
	log: "log/R1R2_filtering.log"
	threads :8
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


rule Filter_Map_R1R2:
	input: sam = rules.Filter_Map_R1R2_unique.output.confident
	output: bam = "12-Annotation_R1R2/{NAME}_TAG{TAG}.bam",
			bed = "12-Annotation_R1R2/{NAME}_TAG{TAG}.bed",
			IS= "12-Annotation_R1R2/{NAME}_TAG{TAG}_IS.bed",
			idx= "12-Annotation_R1R2/{NAME}_TAG{TAG}.bai",
			ISsorted="12-Annotation_R1R2/{NAME}_TAG{TAG}_sorted.bed",
			IScluster="12-Annotation_R1R2/{NAME}_TAG{TAG}_clusters.bed",
			IScollapsedSize="12-Annotation_R1R2/{NAME}_TAG{TAG}_IS_collapsedBySizeCluster.bed",
			ISclustersorted="12-Annotation_R1R2/{NAME}_TAG{TAG}_IS_cluster_sorted.bed",
			IScollapsed="12-Annotation_R1R2/{NAME}_TAG{TAG}_IS_Collapsed.bed"
			
			
			
	threads: 8
	message: "Sort SAM by read name, convert to bed, find IS and calculate fragments size from R1 and R2 extremities"
	shell: """
			samtools sort {input.sam} -o {output.bam} -O BAM -n;
			bedtools bamtobed -bedpe -mate1 -i {output.bam} > {output.bed};
			samtools sort {input.sam} -o {output.bam} -O BAM;
			samtools index {output.bam} {output.idx};
			awk 'OFS="\\t" {{if($9=="+"){{print $1,$2,$2+1,$7,$8,$9, $6-$2}} else {{print $1,$3,$3+1,$7,$8,$9,$3-$5}}}}' {output.bed} > {output.IS};
			bedtools sort -i {output.IS} > {output.ISsorted};
			bedtools cluster -s -d -1 -i {output.ISsorted} > {output.IScluster};
			sort -k 8,8n -k 7,7n {output.IScluster} > {output.ISclustersorted};
			bedtools groupby -i {output.ISclustersorted} -g 8,7 -c 1,2,3,4,5,6,7,7,8 -o distinct,mode,mode,collapse,median,distinct,distinct,count,distinct  | cut -f 3- > {output.IScollapsedSize};
			bedtools groupby -i {output.IScollapsedSize} -g 9 -c 1,2,3,4,5,6,8,7 -o distinct,mode,mode,collapse,median,distinct,sum,count_distinct | cut -f 2-  > {output.IScollapsed}
			"""

###################################################################################################
# Process R1 reads that do not contain the linker sequence nor in their R2 counterpart
# These reads are only use for a qualitative anaylsis of IS, not for the quantification.	

rule Merge_R1alones:
	input: A=rules.Find_Linker_R2.output.R1noLinker, B=rules.Filter_Map_R1R2_unique.output.R1_badqualfq, C = rules.Map_R1_R2_pairs.output.R1unal, D=rules.Map_R1_with_Linker.output.unmapped
	output: "11-R1only/{NAME}_R1_TAG{TAG}_merged.fastq"
	message: "Merging R1 reads without Linker and R1 from pairs that were mapped multiple times and R1 reads from pairs that do not map concordantly."
	threads:8
	log: 
	shell: """
			cat {input.A} {input.B} {input.C} {input.D} > {output} 
			"""
			
			
## This step will be removed when using the sonication method.
	
rule Cut_TTAA_R1_noLinker:
	input: rules.Merge_R1alones.output
	output: TTAA="08-TTAA_cutR1alone/{NAME}_R1_TAG{TAG}_TTAA-cut.fastq",TTAA_cut="08-TTAA_cutR1alone/{NAME}_R1_TAG{TAG}_TTAA-cut20bp.fastq"
	message: "Cutting TTAA sequences that are still present and filter reads shorter than 20bp"
	threads: 8
	log: "log/TTAAR1alone.log"
	shell: """
			  awk '{{if(NR %4 ==2) {{x=index($0,"TTAA");if(x>0){{print substr($0,1,x)}}else{{print $0}}}} else {{if(NR % 4 ==0 && x>0) {{print substr($0,1,x)}} else{{print $0}}}}}}' {input} > {output.TTAA};
			  cutadapt -m 20 -o {output.TTAA_cut} {output.TTAA} > {log}
			"""	
			
rule Collapse_identical_R1_noLinkerReads:
	input: rules.Cut_TTAA_R1_noLinker.output.TTAA_cut
	output: file = "10-collapsedR1alone/{NAME}_R1_TAG{TAG}_TTAA-cut20bp_trimmed.fastq"
	message: "Collapsing identical reads before mapping"
	log:"log/collapsingR1alone.log"
	shell: """
			#fastx_collapser -v -i {input} -o {output} > {log}
			 seqcluster collapse -f {input} -d -o "10-collapsedR1alone/" -m 0
			"""	

			
rule Map_R1_noLinker:
		input: R1=rules.Collapse_identical_R1_noLinkerReads.output.file
		output: mapped="11-mappingR1alone/{NAME}_R1_TAG{TAG}_mapped.sam",unmapped="11-mappingR1alone/{NAME}_R1_TAG{TAG}_unmapped.fastq", exact="11-mappingR1alone/{NAME}_R1_TAG{TAG}_mapped_exact.sam"
		message:"Mapping R1 reads without Linker in R1 nor R2"
		threads: 8
		log: met="log/mappingR1alone.log",summary="log/mappingR1alone_numbers.log"
		shell: """
				bowtie2 -N 1 -L 25 -i S,25,0 --score-min L,0,-0.15 --gbar 10 -p {threads} -x {BOWTIE2_INDEX} --trim3 5 --no-unal --un {output.unmapped} --met-file {log.met} {input.R1} -S {output.mapped} 2> {log.summary}
				awk -F "\\t" '(/^@/) || ($2==0 && !/XS:i:/ && !/MD:Z:[012][A-Za-z].*\t/) || ($2==16 && !/XS:i:/ && !/MD:Z:.*[A-Za-z][012]\t/) {{print $0}}' {output.mapped} > {output.exact}
				"""
			

			
rule Filter_Map_R1_noLinker:
	input: rules.Map_R1_noLinker.output.exact
	output: bam="12-AnnotationR1alone/{NAME}_R1_TAG{TAG}_mapped.bam",
			sortbam="12-AnnotationR1alone/{NAME}_R1_TAG{TAG}_mapped_sorted.bam",
			idx="12-AnnotationR1alone/{NAME}_R1_TAG{TAG}_mapped_sorted.bai",
			bed="12-AnnotationR1alone/{NAME}_R1_TAG{TAG}_mapped_sorted.bed",
			IS="12-AnnotationR1alone/{NAME}_R1_TAG{TAG}_IS.bed",
			ISsorted="12-AnnotationR1alone/{NAME}_R1_TAG{TAG}_IS_sorted.bed",
			ISsortednum="12-AnnotationR1alone/{NAME}_R1_TAG{TAG}_IS_sortedNum.bed",
			IScollapsedSize="12-AnnotationR1alone/{NAME}_R1_TAG{TAG}_IS_collapsedBySizeCluster.bed",
			IScluster="12-AnnotationR1alone/{NAME}_R1_TAG{TAG}_IS_cluster.bed",
			ISclustersorted="12-AnnotationR1alone/{NAME}_R1_TAG{TAG}_IS_cluster_sorted.bed",
			IScollapsed="12-AnnotationR1alone/{NAME}_R1_TAG{TAG}_IS_Collapsed.bed"
	params: qual="10"
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
	output: merged=temp("13-qualitativeIS/{NAME}_TAG{TAG}_qualIS.bed"), 
			merged_corrected="13-qualitativeIS/{NAME}_TAG{TAG}_qualIS_merged.bed",
			merged_sorted = "13-qualitativeIS/{NAME}_TAG{TAG}_qualIS_merged_sorted.bed",
			merged_sorted_collapsed= "13-qualitativeIS/{NAME}_TAG{TAG}_qualIS_merged_sorted_collapsed.bed"
	message: "Merging and collapsing IS from each step"
	threads: 8
	log: "log/qualIS.log"
	shell: """
			cat  {input.R1full} {input.R1alone} {input.R1R2} > {output.merged};
			awk '{{print "chr"$0}}' {output.merged} > {output.merged_corrected};
			sort -k1,1 -k2,2n {output.merged_corrected} > {output.merged_sorted};
			bedtools merge -s -d -1 -c 1,2,3,4,5,6,7 -o distinct,mode,mode,collapse,median,distinct,sum -i {output.merged_sorted} | cut -f 5- > {output.merged_sorted_collapsed};
			"""
	
			

rule annotateIS_:
	input: rules.QualIS.output.merged_corrected
	output: slop="13-qualitativeIS/{NAME}_TAG{TAG}_qualIS_merged_slop400.bed",slopsorted="13-qualitativeIS/{NAME}_TAG{TAG}_qualIS_merged_slop400_sorted.bed",bigCoverage = "13-qualitativeIS/{NAME}_TAG{TAG}_qualIS_H3K27me3.stab",bedcoverage = "13-qualitativeIS/{NAME}_TAG{TAG}_qualIS_H3K27me3.bed"
	log: "log/H3K27me3.log"
	shell:
		"""
		bedtools slop -i {input} -g ../references/hg19.chrom.sizes -b 200 | cut -f 1-4 > {output.slop};
		sort -k1,1 -k2,2n -k6,6 {output.slop} > {output.slopsorted};
		#bigWigAverageOverBed ../dataset/epigenetics/GSM621664_BI.Mobilized_CD34_Primary_Cells.H3K27me3.UW_RO_01536.bw {output.slopsorted} {output.bigCoverage} -stats={log} -bedOut={output.bedcoverage} -minMax
		"""
	
rule QuantIS:
	input:R1full=rules.Filter_Map_R1.output.IScollapsed,
		R1R2=rules.Filter_Map_R1R2.output.IScollapsed
	output: touch("13-quantitativeIS/{NAME}_TAG{TAG}_qualIS_merged.bed")

	
rule MakeReport:
	input: rules.QuantIS.output, rules.QualIS.output.merged_sorted_collapsed
	output: touch("13-Report/{NAME}_TAG{TAG}_report.pdf")
	
	
rule GetAnnotations_UCSC_HGNC:
	input: 
	output: HGNC="../dataset/UCSC/hg19/"+TODAY+"/HGNC_.txt",UCSC="../dataset/UCSC/hg19/"+TODAY+"/"
	log: "../dataset/UCSC/hg19/"+TODAY+"/ucsc_remote_wget.log"
	shell : """
			perl ../Scripts/retreive_HUGO.pl > {output.HGNC};
			wget -N -i ../Scripts/hg19_annotation_urls.txt -P {output.UCSC} -o {log} -x;
			"""
			
	
	