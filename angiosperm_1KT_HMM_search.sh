#!/bin/bash

#BSUB -W 72:00             	# How much time does your job need (HH:MM)
#BSUB -R rusage[mem=500]	# How much memory
#BSUB -R span[hosts=1]		# Keep on one CPU cluster
#BSUB -n 17                	# Where X is in the set {1..X}
#BSUB -J search     	# Job Name
#BSUB -o out.%J           	# Append to output log file
#BSUB -e err.%J           	# Append to error log file
#BSUB -q long            	# Which queue to use {short, long, parallel, GPU, interactive}

module load MAFFT/7.313
module load hmmer/3.1b2
module load blast/2.2.22
module load noisy/1.5.12
module load iqtree/1.6.3
module load R/3.6.1_packages/tidyverse/1.3.0 gcc/8.1.0

run="CLV1"



#<<COMMENT


##this section will take a list of genes, make an HMM, search all 1kt and n best hit collect sequences, then scan them for domain content
cd /home/jm33a/domain_evolution/1KT_searches/$run/

echo "constructing database for input sequences, collecting and aligning inputs..."
    cat $(ls ~/supertree/databases/primary_transcript_databases/*.oneline.fa) > combined_onelines.fa
    grep -w -A 1 -f *geneIDs.txt --no-group-separator combined_onelines.fa | awk '{print $1}' > $run.input_seqs.fa #collect input clade seqs
    rm combined_onelines.fa
    linsi $run.input_seqs.fa > $run.input_align.fasta

echo "make $run hmm from input alignment..."
	hmmbuild -o hmmout.txt $run.model.hmm $run.input_align.fasta
	rm hmmout.txt


echo "collect the top 2 hits for each of the species listed in flowering_plants_species_IDs"
while read species_ID
do
    echo "Processing species $species_ID"
    directory=$(find /project/uma_madelaine_bartlett/JarrettMan/sequence_databases/1KP_seqs/seqs -name "*$species_ID*" -type d) #find directory of the species
    protein_database_path=$directory/*.prots.out #this is now the path to species' protein database
   
    test -f $protein_database_path #test if file exists. Write 0 to $? if yes, 1 if no
    if [ $? -eq 0 ] #only perform the search when a database exists, otherwise report the error
    then
        hmmsearch -o hmmout.txt --noali --tblout table.output.txt $run.model.hmm $protein_database_path #scan all sequences from a species with the HMM model
    
        #next collec the sequence of the top n hits found in the search
        for hit_number in {4..5} #top hit is on line 4, specify how many lines after that
        do
            hit=$(sed -n "${hit_number}p" < table.output.txt | awk '{print $1}') #gene ID from hit number line
            hit_full_ID=">${species_ID}_${hit}" #concatenate the species ID and the geneID with an underscore between, using fasta format >
            echo $hit_full_ID >> $run.top_hits.seqs.fa #add gene full gene ID to fasta sequences
            grep -w -A 1 $hit $protein_database_path | tail -n 1 >> $run.top_hits.seqs.fa #add gene's sequence to running list
        done
    else
        echo "$species_ID does not have a peptide database"
    fi
    
done < /home/jm33a/domain_evolution/flowering_plants_species_IDs

#clean up a bit before proceeding
rm hmmout.txt
rm table.output.txt

echo "removing hits with sequences that were not found"
    grep -v "#" $run.top_hits.seqs.fa > top_hits_clean.fa
    rm $run.top_hits.seqs.fa
    cat top_hits_clean.fa > $run.top_hits.seqs.fa
    rm top_hits_clean.fa


#this section to remove any genes for which both an LRR and RLK domain are not found
echo "Building Pfam domain table from genes" 	##in case need to regenerate pfam searchable database: $ hmmpress Pfam-A.hmm top_hits.seqs.fa
	hmmscan --noali -o pfamout.temp --cut_tc --tblout $run.pfamout.tsv ~/pfam/hmmfiles/Pfam-A.hmm $run.top_hits.seqs.fa 
	rm pfamout.temp
echo "use R script to process this list to only genes with both domains, and only the gene IDs. output is $run.both_domains_IDs.txt"
    Rscript /home/jm33a/domain_evolution/sort_full_LRR_RLKs_olnly.R $run



echo "use the filtered list of gene IDs to collect sequence, make an alignment and tree"
fgrep -w --no-group-separator -A 1 -f $run.hits_with_both_domains_IDs.txt $run.top_hits.seqs.fa > $run.hits_with_both_domains_seqs.fa #pull only seqs from R screening
mkdir tree_from_hits
cp $run.hits_with_both_domains_seqs.fa tree_from_hits
cd tree_from_hits
cat /home/jm33a/domain_evolution/scaffold_seqs.fa $run.hits_with_both_domains_seqs.fa >> $run.hits_and_scaffold_seqs.fa #add the backbone scaffold sequences to tree

mafft --thread 16 $run.hits_and_scaffold_seqs.fa > $run.hits_and_scaffold_align.fasta #quick version

echo "running noisy to clean up alignment..."
	noisy -s --noconstant $run.hits_and_scaffold_align.fasta
	#remove noisy log files
	rm *typ.eps
	rm *sta.gr
	rm *idx.txt



echo "Running IQtree with auto find model of evo and inference of site rate evolution.."
	iqtree -s *out.fas -fast -nt 16 -m LG+F+R7 #quick version
	#iqtree -s *out.fas -bb 1000 -wsr -nt 16 -m LG+F+R7 #full version, infer site rates


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








