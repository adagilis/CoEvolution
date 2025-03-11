#!/bin/bash

while getopts t:o: option
do
    case "${option}"
        t) taxon=${OPTARG};;
        o) output=${OPTARG};;
    esac
done

datasets download genome taxon $taxon --assembly-level chromosome --include protein,gtf,gff --filename $taxon.zip 

#Rename files to species names. This may overwrite some files if multiple assemblies per species were downloaded, use at own discretion.

for i in *.faa;
    do esearch -db assembly -query ${i/.faa/} | efetch -format docsum | xtract -pattern DocumentSummary -element Organism | awk -v i=$i '{print "mv " i " " substr($1,1,1)"_"$2".faa"}';
done >> cmds.txt

#Separated out in case you want to comment the below out to check before executing manually. Good when you want specific name changes, etc.
cat cmds.txt | sh