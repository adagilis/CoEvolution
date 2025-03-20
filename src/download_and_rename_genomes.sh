#!/bin/bash

mkdir seqs

while getopts t: option
do
    case "${option}"
	in
        t) taxon=${OPTARG};;
    esac
done

datasets download genome taxon $taxon --assembly-level chromosome --include protein,gff3 --filename $taxon.zip 


unzip $taxon.zip

cd ncbi_dataset/data
for i in GCF*;
do
    mv $i/protein.faa $i.faa
    mv $i/genomic.gff $i.gff
done

#Rename files to species names. This may overwrite some files if multiple assemblies per species were downloaded, use at own discretion.

for i in *.faa;
    do esearch -db assembly -query ${i/.faa/} | efetch -format docsum | xtract -pattern DocumentSummary -element Organism | awk -v i=$i '{print "mv " i " " substr($1,1,1)"_"$2".faa"}';
done >> cmds.txt

#Separated out in case you want to comment the below out to check before executing manually. Good when you want specific name changes, etc.
cat cmds.txt | sh

mv *.faa ../../seqs
mv *.gff ../../seqs

cd ../..

rm -rf ncbi_datasets/

