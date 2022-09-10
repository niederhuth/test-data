#Combine and reformat exonerate data into gff

path1=$(pwd | sed s/data.*/scripts/)
a=$(pwd | sed s/.*\\///)

if [ -f target_chunk_1_query_chunk_1 ]
then
	mkdir exonerate_output
	mv target_chunk_*_query_chunk_* exonerate_output
fi

mkdir tmp
cd exonerate_output

for i in target_chunk_*_query_chunk_*
do
	sed '1,2d' ${i} | grep -v "\-\-\ completed\ exonerate\ analysis" > ../tmp/${i}.tmp
done
cd ..
cat tmp/*tmp > ${a}

#Reformat exonerate output
perl ${path1}/annotation/pl/reformat_exonerate_protein_gff.pl --input_gff ${a} --output_gff tmp.gff

#Sort the gff file
gff3_sort -g tmp.gff -og ${a}.gff
rm tmp.gff

#Rename for maker
cat ${a}.gff | sed 's/mRNA/protein_match/g' | \
sed 's/exon/match_part/g' | sed s/protein2genome/protein_gff\:protein2genome/ > ${name}_maker_input.gff

rm -R tmp

echo "Done"
