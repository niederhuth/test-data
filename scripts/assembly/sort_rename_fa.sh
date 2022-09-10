#!/bin/bash

#Ugly script to sort and rename fasta files
#Use bash sort_rename_fa.sh <input_fasta> <cutoff_length>

cat $1 |\
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
awk -v a=$2 '$1 >= a' |\
sort -k1,1rn |\
cut -f 2- |\
tr "\t" "\n" |\
sed s/^\>/old_name\=/ |\
awk -v i=9000000 '/old_name\=/{print ">contig_" ++i, $0; next}{print}' |\
sed s/contig_9/contig_0/

