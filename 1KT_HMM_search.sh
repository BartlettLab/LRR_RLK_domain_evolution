#!/bin/bash

#BSUB -W 1:00             	# How much time does your job need (HH:MM)
#BSUB -R rusage[mem=1500]	# How much memory
#BSUB -R span[hosts=1]		# Keep on one CPU cluster
#BSUB -n 1                	# Where X is in the set {1..X}
#BSUB -J search         	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file
#BSUB -q short            	# Which queue to use {short, long, parallel, GPU, interactive}

module load MAFFT/7.313
module load hmmer/3.1b2
module load blast/2.2.22
module load noisy/1.5.12 
module load iqtree/1.6.3
# exit when any command fails
set -e

run="CLV1"

<<COMMENT
##this section will take a list of genes, make an HMM, search all 1kt and n best hit collect sequences, then scan them for domain content



cd /home/jm33a/domain_evolution/1KT_searches/$run/

echo "constructing database for input sequences, collecting and aligning inputs..."
    cat $(ls ~/supertree/databases/primary_transcript_databases/*.oneline.fa) > combined_onelines.fa
    grep -w -A 1 -f *geneIDs.txt --no-group-separator combined_onelines.fa | awk '{print $1}' > $run.input_seqs.fa #collect input clade seqs
    rm combined_onelines.fa
    linsi $run.input_seqs.fa > $run.input_align.fasta
    
    

echo "make hmm from input alignment..."
	hmmbuild -o hmmout.txt $run.model.hmm $run.input_align.fasta
	rm hmmout.txt

echo "run hmm searches..."
	#search_databases=~/supertree/databases/primary_transcript_databases/*.oneline.fa #genomes from LRR_RLK paper
	
	search_databases=/project/uma_madelaine_bartlett/JarrettMan/sequence_databases/1KP_seqs/seqs/*/*.prots.out #1 thousand transcriptomes
	for current_transcriptome in $search_databases; do
	   #current_transcriptome=/project/uma_madelaine_bartlett/JarrettMan/sequence_databases/1KP_seqs/seqs/AQXA-Oenothera_laciniata-MTJ_169/AQXA-SOAPdenovo-Trans-assembly.prots.out
		file_name=${current_transcriptome##*/} #store file_name as the variable after the last '/', which returns the file name only
		species_ID=$(echo $file_name | awk '{print substr($0,0,4)}') #only the first 4 characters of the file name, which is the species ID
		hmmsearch -o hmmout.txt --noali --tblout $run.HMMtable.output.txt $run.model.hmm $current_transcriptome #scan all sequences from a species with the HMM file
        rm hmmout.txt

        for current_hit in {4..5} #top 5 hits are on lines 4 to 8, always start with line 4
        do
            #current_hit=4
            gene_to_collect=$(sed -n "${current_hit}p" < $run.HMMtable.output.txt | awk '{print $1}') #numerical ID of best match
		    species_and_record_ID=">${species_ID}_${gene_to_collect}" #concatenate the species ID and the geneID with an underscore between, using fasta format >
	        #print the species/gene ID and sequence to the top hits list in fasta format
	    	echo $species_and_record_ID >> $run.top_hits.seqs.fa #add gene species/gene ID to running list
	    	grep -w -A 1 $gene_to_collect $current_transcriptome | tail -n 1 >> $run.top_hits.seqs.fa #add gene's sequence to running list
    	done
	done

echo "Building Pfam domain table from genes" 	##in case need to regenerate pfam searchable database: $ hmmpress Pfam-A.hmm top_hits.seqs.fa
	hmmscan --noali -o pfamout.temp --cut_tc --tblout $run.pfamout.tsv ~/pfam/hmmfiles/Pfam-A.hmm $run.top_hits.seqs.fa #this requires specialized rhmmer read_tblout in R. See if can get around this
	hmmscan --noali --cut_tc --tblout $run.pfamout.tsv ~/pfam/hmmfiles/Pfam-A.hmm $run.top_hits.seqs.fa
	rm pfamout.temp
echo "done"




exit


COMMENT

#####next step, use R script pfam_scan.R to process this list to only genes with both domains, and only the gene IDs. output is $run.both_domains_IDs.txt
Rscript /home/jm33a/domain_evolution/sort_full_LRR_RLKs_olnly.R $run





###next use the filtered list of gene IDs to collect sequence, make an alignment and tree. 
fgrep -w --no-group-separator -A 1 -f $run.hits_with_both_domains_IDs.txt $run.top_hits.seqs.fa > $run.hits_with_both_domains_seqs.fa #pull only seqs from R screening, dumpt them in round 2 folder
mkdir tree_from_hits
cp $run.hits_with_both_domains_seqs.fa tree_from_hits
cd tree_from_hits
cat /home/jm33a/domain_evolution/scaffold_seqs.fa $run.hits_with_both_domains_seqs.fa >> $run.hits_and_scaffold_seqs.fa




mafft --thread 16 $run.hits_and_scaffold_seqs.fa > $run.hits_and_scaffold_align.fasta #quick version

echo "running noisy to clean up alignment..."
	noisy -s --noconstant $run.hits_and_scaffold_align.fasta
	#remove noisy log files
	rm *typ.eps
	rm *sta.gr
	rm *idx.txt



echo "Running IQtree with auto find model of evo and inference of site rate evolution.."
	iqtree -s *out.fas -bb 1000 -fast -nt 16 -m JTT+F+R9 #quick version
	#iqtree -s *out.fas -bb 1000 -wsr -nt 16 -m JTT+F+R9 #full version, infer site rates


echo "cleaning up.."
		#remove iqtree log files
	rm *splits.nex
	rm *ckp.gz
	rm *bionj
	rm *mldist
	rm *.log
	rm *.iqtree
	rm *model.gz
	rm *uniqueseq.phy

exit







