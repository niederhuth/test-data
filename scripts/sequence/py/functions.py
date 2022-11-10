from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils import GC123
from Bio.Seq import Seq


#Calculate the GC, GC1, GC2, and GC3 for each sequence in a CDS fasta file
#GC1, GC2, & GC3 are codon specific metrics and the input is expected to be CDS sequences
#fasta is the input fasta file
#output, will write results to specified file, otherwise will return to screen
def gc123(fasta,output=()):
	#Create an array with a header line
	GCtable = [['Transcript','GC','GC1','GC2','GC3']]
	#parse and read over fasta file
	for i in SeqIO.parse(open(fasta, "r"), "fasta"):
		#use Biopython GC123 to calculate GC content and codon GC content
		gc = GC123(i.seq)
		#for each sequence, add these values to the array
		GCtable = GCtable + [[i.id,str(gc[0]),str(gc[1]),str(gc[2]),str(gc[3])]]
	#output results
	if output:
		with open(output, 'w') as out_file:
			out_file.writelines('\t'.join(i2) + '\n' for i2 in GCtable)
	else:
		return(GCtable)

