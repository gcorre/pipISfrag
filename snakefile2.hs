workdir: "/home/"
#from snakemake.utils import report
#TAG="1"
	
BOWTIE2_INDEX = "/home/references/genomes/homo_sapiens/hg19_GRCh37/index_bowtie2/GRCh37.75"

### JJJJJ

from glob import glob;
FILES = glob('*R[0-9]_[0-9]*.fastq.gz');
import re;
TRAILING=list(set([re.search('_([0-9]+).fastq.gz',FILES[0]).group(1) for w in FILES]))[0]
FILES=[re.sub('_R[1-2]_[0-9]+.fastq.gz', '',w) for w in FILES];
SAMPLE=list(set(FILES))[0];
TAG=re.search('TAG([0-9]+)', SAMPLE).group(1)

rule targets:
	input: "12-Annotation_R1R2/"+SAMPLE+"_TAG"+TAG+"_IS.bed", "12-Annotation/"+SAMPLE+"_R1_TAG"+TAG+"_IS_sorted.bed"

	
rule Demultiplex:
	input: R1="{NAME}_R1_"+TRAILING+".fastq.gz", R2="{NAME}_R2_"+TRAILING+".fastq.gz"
	output: R1="01-Demultiplex/{NAME}_R1_TAG{TAG}.fastq.gz",R2="01-Demultiplex/{NAME}_R2_TAG{TAG}.fastq.gz"
	log: "log/demultiplex-merged.log"
	threads: 8
	message: "Demultiplexing according to TAGs list"
	shell: """
			 fastq-multx -B ../dataset/illuminaTAGs_eautils.txt -b -x -m 1 {input.R1} {input.R2} -o 01-Demultiplex/{wildcards.NAME}_R1_%.fastq.gz 01-Demultiplex/{wildcards.NAME}_R2_%.fastq.gz > {log};
			"""
			
			
rule provirus:
	input: R1=rules.Demultiplex.output.R1, R2=rules.Demultiplex.output.R2
	output: R1Prov="02-provirus/{NAME}_R1_TAG{TAG}_Prov.fastq", R2Prov="02-provirus/{NAME}_R2_TAG{TAG}_Prov.fastq", R1noProv="02-provirus/{NAME}_R1_TAG{TAG}_noProv.fastq", R2noProv="02-provirus/{NAME}_R2_TAG{TAG}_noProv.fastq",
	log:"log/Provirus.log"
	threads: 8
	message: "Look for the Provirus in R1"
	shell: """
			cutadapt -g aaaatctctagcagtggcgcccgaacag -O 28 -e 0.1 --no-trim --no-indels -o {output.R1Prov} -p {output.R2Prov} --untrimmed-paired-output {output.R2noProv} --untrimmed-output {output.R1noProv} {input.R1} {input.R2} > {log}
			"""
			
rule LTR:
	input: R1=rules.provirus.output.R1noProv, R2=rules.provirus.output.R2noProv
	output: R1LTR="03-LTR/{NAME}_R1_TAG{TAG}_LTR.fastq", R2LTR="03-LTR/{NAME}_R2_TAG{TAG}_LTR.fastq", R1noLTR="03-LTR/{NAME}_R1_TAG{TAG}_noLTR.fastq", R2noLTR="03-LTR/{NAME}_R2_TAG{TAG}_noLTR.fastq",
	log:"log/LTR.log"
	threads: 8
	message: "Look for the LTR in R1"
	shell: """
			cutadapt -g GTCTGTTGTGTGACTCTGGTAAC -O 23 -e 0.1 --no-indels -o {output.R1LTR} -p {output.R2LTR} --untrimmed-paired-output {output.R2noLTR} --untrimmed-output {output.R1noLTR} {input.R1} {input.R2} > {log}
			"""
	
rule Elong:
	input: R1=rules.LTR.output.R1LTR, R2=rules.LTR.output.R2LTR
	output: R1elong="04-ELONG/{NAME}_R1_TAG{TAG}_Elong.fastq", R2elong="04-ELONG/{NAME}_R2_TAG{TAG}_Elong.fastq", R1noelong="04-ELONG/{NAME}_R1_TAG{TAG}_noElong.fastq", R2noelong="04-ELONG/{NAME}_R2_TAG{TAG}_noElong.fastq",
	log:"log/ELONG.log"
	threads: 8
	message: "Look for the Elong in R1"
	shell: """
			cutadapt -g TAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTA -O 43 -e 0.16 --no-indels -o {output.R1elong} -p {output.R2elong} --untrimmed-paired-output {output.R2noelong} --untrimmed-output {output.R1noelong} {input.R1} {input.R2} > {log}
			"""

rule GCA:
	input: R1=rules.Elong.output.R1elong, R2=rules.Elong.output.R2elong
	output: GCAR1="09-GCA/{NAME}_R1_TAG{TAG}_GCA.fastq", noGCAR1="09-GCA/{NAME}_R1_TAG{TAG}_noGCA.fastq",GCAR2="09-GCA/{NAME}_R2_TAG{TAG}_GCA.fastq", noGCAR2="09-GCA/{NAME}_R2_TAG{TAG}_noGCA.fastq"
	message: "Filtering R1 reads that start with GCA"
	log: "log/GCA-R1.log"
	shell: """
			cutadapt -g ^GCA  -m 20 -e 0 --no-indels -o {output.GCAR1} -p {output.GCAR2} --untrimmed-output {output.noGCAR1} --untrimmed-paired-output {output.noGCAR2} {input.R1} {input.R2} > {log};
			"""			

rule LinkerR1:
	input: R1=rules.GCA.output.GCAR1, R2=rules.GCA.output.GCAR2
	output: R1Linker="05-Linker/{NAME}_R1_TAG{TAG}_Linker.fastq", R2Linker="05-Linker/{NAME}_R2_TAG{TAG}_Linker.fastq", R1noLinker="05-Linker/{NAME}_R1_TAG{TAG}_noLinker.fastq", R2noLinker="05-Linker/{NAME}_R2_TAG{TAG}_noLinker.fastq",
	log:"log/LinkerR1.log"
	threads: 8
	message: "Look for the Linker in R1"
	shell: """
			cutadapt -a GTCCCTTAAGCGGAGCCCT -O 8 -e 0.2 --no-indels -o {output.R1Linker} -p {output.R2Linker} --untrimmed-paired-output {output.R2noLinker} --untrimmed-output {output.R1noLinker} {input.R1} {input.R2} > {log}
			"""				
## If R1 doesn't contains the linker, look in R2, else keep only R1.


rule cut_TTAA:
	input: rules.LinkerR1.output.R1Linker
	output: TTAA="08-TTAA_cut/{NAME}_R1_TAG{TAG}_TTAA-cut.fastq"
	message: "Cutting TTAA sequences that are still present"
	threads: 8
	log: "log/TTAA.log"
	shell: """
			  awk '{{if(NR %4 ==2) {{x=index($0,"TTAA");if(x>0){{print substr($0,1,x)}}else{{print $0}}}} else {{if(NR % 4 ==0 && x>0) {{print substr($0,1,x)}} else{{print $0}}}}}}' {input} > {output.TTAA}
			"""			
			
rule collapse_reads:
	input: rules.cut_TTAA.output.TTAA
	output: file = "10-collapsed/{NAME}_R1_TAG{TAG}_TTAA-cut_trimmed.fastq"
	message: "Collapsing identical reads before mapping"
	log:"log/collapsingR1.log"
	shell: """
			#fastx_collapser -v -i {input} -o {output} > {log}
			 seqcluster collapse -f {input} -d -o "10-collapsed/" -m 0
			"""

			
rule mappingR1_full :
	input: rules.collapse_reads.output.file
	output: mapped="11-mapping/{NAME}_R1_TAG{TAG}_mapped.sam",unmapped="11-mapping/{NAME}_R1_TAG{TAG}_unmapped.fastq"
	log: met="log/mappingR1.log",summary="log/mapping_numbers.log"
	threads: 8
	message: "Mapping1 : Mapping with Bowtie2 of R1 reads starting with GCA"
	shell: """
			bowtie2 -N 1 -L 25 -i S,25,0 --score-min L,0,-0.15 --gbar 10 -p {threads} -x {BOWTIE2_INDEX} --no-unal --un {output.unmapped} --met-file {log.met} {input} -S {output.mapped} 2> {log.summary}
			"""

rule Filter_MapR1:
	input: rules.mappingR1_full.output.mapped
	output: bam="12-Annotation/{NAME}_R1_TAG{TAG}_mapped.bam",
			sortbam="12-Annotation/{NAME}_R1_TAG{TAG}_mapped_sorted.bam",
			idx="12-Annotation/{NAME}_R1_TAG{TAG}_mapped_sorted.bai",
			bed="12-Annotation/{NAME}_R1_TAG{TAG}_mapped_sorted.bed",
			IS="12-Annotation/{NAME}_R1_TAG{TAG}_IS.bed",
			ISsorted="12-Annotation/{NAME}_R1_TAG{TAG}_IS_sorted.bed",
			ISsortednum="12-Annotation/{NAME}_R1_TAG{TAG}_IS_sortedNum.bed",
			IScollapsedSize="12-Annotation/{NAME}_R1_TAG{TAG}_IS_collapsedBySizeCluster.bed",
			IScluster="12-Annotation/{NAME}_R1_TAG{TAG}_IS_cluster.bed",
			ISclustersorted="12-Annotation/{NAME}_R1_TAG{TAG}_IS_cluster_sorted.bed",
			IScollapsed="12-Annotation/{NAME}_R1_TAG{TAG}_IS_Collapsed.bed"
	params: qual="10"
	message:"Filtering mapped reads, converting to BED,Collapsing by 3nt and counting read abundance by IS"
	shell: """ 
			samtools view -q {params.qual} -bS {input} > {output.bam};
			samtools sort {output.bam} > {output.sortbam};
			samtools index {output.sortbam} {output.idx};
			bedtools bamtobed -i {output.sortbam} > {output.bed};
			awk 'function abs(a){{return ((a < 0) ? -a : a)}} OFS="\\t" {{split($4,a,"x");print $0,abs($3-$2),a[2]}}' {output.bed} > {output.ISsortednum};
			awk 'OFS="\\t" {{if($6=="-") {{print $1,$3-1,$3,$4,$5,$6,$7,$8}} else {{print $1,$2,$2+1,$4,$5,$6,$7,$8}}}}' {output.ISsortednum} > {output.IS};
			bedtools sort -i {output.IS} > {output.ISsorted};
			bedtools cluster -s  -i {output.ISsorted} > {output.IScluster};
			sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -k 7,7n {output.IScluster} > {output.ISclustersorted};
			bedtools groupby -i {output.ISclustersorted}  -g 9,7 -c 1,2,3,4,5,6,7,8,8,9 -o distinct,mode,mode,collapse,median,distinct,distinct,sum,count,distinct  | cut -f 3- > {output.IScollapsedSize};
			bedtools groupby -i {output.IScollapsedSize}  -g 10 -c 1,2,3,4,5,6,8,9,7 -o distinct,mode,mode,collapse,median,distinct,sum,sum,count_distinct | cut -f 2-  > {output.IScollapsed}
			"""
			
	
rule LinkerR2:
	input: R1=rules.LinkerR1.output.R1noLinker, R2=rules.LinkerR1.output.R2noLinker
	output: R1Linker="06-LinkerR2/{NAME}_R1_TAG{TAG}_Linker.fastq", R2Linker="06-LinkerR2/{NAME}_R2_TAG{TAG}_Linker.fastq", R1noLinker="06-LinkerR2/{NAME}_R1_TAG{TAG}_noLinker.fastq", R2noLinker="06-LinkerR2/{NAME}_R2_TAG{TAG}_noLinker.fastq"
	log:"log/LinkerR2.log"
	threads: 8
	message: "Look for the Linker in R2"
	shell: """
			cutadapt -g AGGGCTCCGCTTAAGGGAC -O 19 -e 0.16 --no-indels -o {output.R2Linker} -p {output.R1Linker} --untrimmed-paired-output {output.R1noLinker} --untrimmed-output {output.R2noLinker} {input.R2} {input.R1} > {log}
			"""	
	

rule RCLTR_R2:
	input: R1=rules.LinkerR2.output.R1Linker,R2=rules.LinkerR2.output.R2Linker
	output: trimmedLTRRCR1="07-cleaning2/{NAME}_R1_TAG{TAG}_trimmedElongRC.fastq", 
			trimmedLTRRCR2="07-cleaning2/{NAME}_R2_TAG{TAG}_trimmedElongRC.fastq"
	log: "log/RCLTR_R2.log"
	threads: 8
	message: "Trimming R2 that have LTR"
	shell: """
			cutadapt -a TGCTAGAGATTTTCCACACTGACTAAAAGGGTCTGAGGGATCTCTA -O 8 -m 20 -e 0.25 -o {output.trimmedLTRRCR2} -p {output.trimmedLTRRCR1} {input.R2} {input.R1} > {log};
			"""	
			
rule cut_TTAAR1R2:
	input: R1=rules.RCLTR_R2.output.trimmedLTRRCR1,R2=rules.RCLTR_R2.output.trimmedLTRRCR2
	output: TTAAR1="10-TTAA_cut/{NAME}_R1_TAG{TAG}_TTAA-cut.fastq", 
			TTAAR2="10-TTAA_cut/{NAME}_R2_TAG{TAG}_TTAA-cut.fastq",
			qualR1="10-TTAA_cut/{NAME}_R1_TAG{TAG}.fastq",
			qualR2="10-TTAA_cut/{NAME}_R2_TAG{TAG}.fastq"
			
	message: "Cutting TTAA sequences that are still present, trimming low quality 3'"
	threads: 8
	log: "log/TTAAR1R2_quality.log"
	shell: """
			  awk '{{if(NR %4 ==2) {{x=index($0,"TTAA");if(x>0){{print substr($0,1,x)}}else{{print $0}}}} else {{if(NR % 4 ==0 && x>0) {{print substr($0,1,x)}} else{{print $0}}}}}}' {input.R1} > {output.TTAAR1};
			  awk '{{if(NR %4 ==2) {{x=index($0,"TTAA");if(x>0){{print substr($0,1,x)}}else{{print $0}}}} else {{if(NR % 4 ==0 && x>0) {{print substr($0,1,x)}} else{{print $0}}}}}}' {input.R2} > {output.TTAAR2};
			  cutadapt -q 20 -m 20 -A XXXXXXX -o {output.qualR1} -p {output.qualR2} {output.TTAAR1} {output.TTAAR2} > {log};
			"""	

			
rule Mapping_R1R2:
	input: R1=rules.cut_TTAAR1R2.output.qualR1, R2=rules.cut_TTAAR1R2.output.qualR2
	output: mapped="11-mappingR1R2/{NAME}_TAG{TAG}_mapped.sam", bam="11-mappingR1R2/{NAME}_TAG{TAG}_mapped.bam",sortbam = "11-mappingR1R2/{NAME}_TAG{TAG}_mapped_sorted.sam",  idx="11-mappingR1R2/{NAME}_TAG{TAG}_mapped_sorted.bai"
	message: "Mapping R1 and R2 reads as pairs"
	params: qual="10"
	log: met="log/mappingR1R2.log",summary="log/mappingR1R2_numbers.log"
	threads: 8
	shell: """
			bowtie2 -p {threads} -x {BOWTIE2_INDEX} -1 {input.R1} -2 {input.R2} -S {output.mapped} --no-unal --al-conc "11-mappingR1R2/{wildcards.NAME}_R%_TAG{wildcards.TAG}_aligned.fastq" --un-conc "11-mappingR1R2/{wildcards.NAME}_R%_TAG{wildcards.TAG}_unaligned.fastq" --dovetail -X 800 --no-mixed --met-file {log.met} 2> {log.summary}
			"""
#cuicui			
rule Filter_MapR1R2:
	input: sam = rules.Mapping_R1R2.output.mapped
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
			bedtools bamtobed -i {output.bam} > {output.bed};
			awk '{{ ORS = (NR%2 ? "\\t" : RS) }} 1' {output.bed} | awk 'OFS="\\t" {{if($6=="+"){{print $1,$2,$2+1,substr($4, 1, length($4)-2),$9-$2,$6 }} else {{print $1,$3,$3+1,substr($4, 1, length($4)-2),$3-$8,$6}}}}' > {output.IS};
			samtools view -bS {input.sam} | samtools sort - -o {output.bam} -O bam;
			samtools index {output.bam} {output.idx};
			bedtools sort -i {output.IS} > {output.ISsorted};
			bedtools cluster -s  -i {output.ISsorted} > {output.IScluster};
			sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 -k 5,5n {output.IScluster} > {output.ISclustersorted};
			bedtools groupby -i {output.ISclustersorted}  -g 7,5 -c 1,2,3,4,6,7,5,5, -o distinct,mode,mode,collapse,distinct,distinct,distinct,count  | cut -f 3- > {output.IScollapsedSize};
			bedtools groupby -i {output.IScollapsedSize}  -g 6 -c 1,2,3,4,5,7,8 -o distinct,mode,mode,collapse,distinct,count_distinct,sum | cut -f 2-  > {output.IScollapsed}
			"""

			
				



			
			
			