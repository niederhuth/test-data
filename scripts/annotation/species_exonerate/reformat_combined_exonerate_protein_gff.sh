#Set path for scripts
path1=$(pwd | sed s/data.*/scripts/)

#name for output
name="protein_alignments"

#Combine the outputs from exonerate
for i in *
do
	cat ${i}/${i} >> protein_alignments
	echo ${i}
done

#Reformat combiend exonerate gff output
perl ${path1}/annotation/reformat_exonerate_protein_gff.pl --input_gff ${name} --output_gff tmp.gff

#Sort
gff3_sort -g tmp.gff -og ${name}.gff
rm tmp.gff

#Modify for maker
cat ${name}.gff | sed 's/mRNA/protein_match/g' | \
sed 's/exon/match_part/g' | sed s/protein2genome/protein_gff\:protein2genome/ > ${name}_maker_input.gff
