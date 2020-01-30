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

<<COMMENT
echo "building database to find priors..."
    cat $(ls ~/supertree/databases/primary_transcript_databases/*.oneline.fa) > combined_onelines.fa
    database=combined_onelines.fa
echo "collecting prior sequences from database..."
    grep -w -A 1 -f clade_input_IDs.txt --no-group-separator $database | awk '{print $1}' > input_seqs.fa #collect sequences of inputs
echo "Aligning prior sequences..."
    mafft input_seqs.fa > input_alignment.fasta #align sequences
echo "building an HMM from the alignment..."
    hmmbuild -o hmmout.txt model.hmm input_alignment.fasta #build HMM of clade


echo "collect the top 2 hits for each of the species listed in flowering_plants_species_IDs"
while read species_ID
do
    echo "Processing species $species_ID"
    directory=$(find /project/uma_madelaine_bartlett/JarrettMan/sequence_databases/1KP_seqs/seqs -name "*$species_ID*" -type d) #find directory of the species
    protein_database_path=$directory/*.prots.out #this is now the path to species' protein database
   
    test -f $protein_database_path #test if file exists. Write 0 to $? if yes, 1 if no
    if [ $? -eq 0 ] #only perform the search when a database exists, otherwise report the error
    then
        hmmsearch -o hmmout.txt --noali --tblout table.output.txt model.hmm $protein_database_path #scan all sequences from a species with the HMM model
    
        #next collec the sequence of the top n hits found in the search
        for hit_number in {4..5} #top hit is on line 4, specify how many lines after that
        do
            hit=$(sed -n "${hit_number}p" < table.output.txt | awk '{print $1}') #gene ID from hit number line
            hit_full_ID=">${species_ID}_${hit}" #concatenate the species ID and the geneID with an underscore between, using fasta format >
            echo $hit_full_ID >> top_hits.seqs.fa #add gene full gene ID to fasta sequences
            grep -w -A 1 $hit $protein_database_path | tail -n 1 >> top_hits.seqs.fa #add gene's sequence to running list
        done
    else
        echo "$species_ID does not have a peptide database"
    fi
    
done < flowering_plants_species_IDs

#clean up a bit before proceeding
rm $database
rm hmmout.txt
rm table.output.txt

#species_ID=WAIL #nothing
#species_ID=WKSU #something


#next use all hits to infer a tree
echo "removing hits with sequences that were not found"
    grep -v "#" top_hits.seqs.fa > top_hits_clean.fa
    rm top_hits.seqs.fa
    cat top_hits_clean.fa > top_hits.seqs.fa
    rm top_hits_clean.fa


echo "adding original clade sequences to hits from 1KT, storing as tree_seqs.fa"
    cat top_hits.seqs.fa input_seqs.fa > tree_seqs.fa
echo "Aligning..."
    mafft tree_seqs.fa > tree_alignment.fasta
echo "Filtering with noisy..."
    noisy -s tree_alignment.fasta 

COMMENT

echo "inferring tree..."
#    iqtree -s *out.fas -bb 1000 -nt 16 -m LG+F+R7
    iqtree -s *out.fas -fast -nt 16 -m LG+F+R7
    


#remove iqtree log files
rm *splits.nex
rm *ckp.gz
rm *bionj
#rm *treefile
rm *mldist
rm *.log
rm *.iqtree
rm *model.gz
rm *uniqueseq.phy
rm *idx.txt
rm *sta.gr
rm *typ.eps



