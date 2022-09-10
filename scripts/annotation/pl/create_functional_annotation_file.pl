#! /usr/bin/perl -w 

# create_functional_annotation_file_2020.pl

# This is an update to the previous functional annotation script.
# Functional descriptions will be based on blast hits or alternatively
# Pfam domain descriptions.

# ./create_functional_annotation_file_2020.pl --protein_fasta /data/run/plantajo/finalized_asm/scoparia/5_renamed_maker_fastas_gffs/fastas_gff_GreaterThan50aa/scoparia_v1.0_proteins_len_50aa.final.fasta  --model_annot /data/run/kchilds/rice_protein_functional_descriptions.txt  --model_blast /data/run/plantajo/finalized_asm/scoparia/5_renamed_maker_fastas_gffs/fastas_gff_GreaterThan50aa/add_functional_annot/scoparia_proteins_vs_rice.outfmt_default.blastp  --pfam_results_file  /data/run/plantajo/finalized_asm/scoparia/5_renamed_maker_fastas_gffs/fastas_gff_GreaterThan50aa/add_functional_annot/scoparia_proteins.final.fasta.iprscan  --max_hits 5  --output /data/run/kchilds/test.txt

# 15 June 2020
# Kevin Childs

# This script parses blast output from searches against UniRef
# and generates a table of annotation text based on the hits
# by doing some regex cleanup of the UniRef headers.

# 27 March 2009
# Kevin Childs

# The script has been modified to read a file of iprscan results instead of
# the rpsblast results vs pfam domains.  Only pfam domains will be extracted
# from the iprscan data.

# 22 September 2010
# Kevin Childs

# Change the script to take functional annotations primarily from arabidopsis functional descriptions, followed by
# grabbing functional descriptions from uniref, followed by using pfam domains as functional descriptions.

# 16 April 2014
# Kevin Childs

#module purge
#module load BioPerl/1.6.924
#protein=/data/run/ejennings/renaming_pd_genes/pd_purged_pilon_Greater35k.all.maker_MSGS_wOutTEs_deFusion_renamed.proteins.fasta
#annot=/data/run/ejennings/TAIR10_short_functional_descriptions.txt
#blast=/data/run/ejennings/pd_purged_pilon_Greater35k.all.maker_MSGS_wOutTEs_deFusion_renamed.proteins_TAIR10_blast.out
#pfam=/data/run/ejennings/pd_maker_renamed_pfam_prot_domains.out

#perl /data/run/ejennings/scripts/create_functional_annotation_file_2020.pl --protein_fasta $protein --model_annot $annot --model_blast $blast --pfam_results_file $pfam --max_hits 5 --output test.txt

use strict;
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;
use Cwd qw{ abs_path };
use Carp;

my $model_blast_file;
my $protein_fasta_file;
my $model_annot_file;
my $pfam_results_file;
my $output_file;
my $max_hits;
my $help;

my $result = GetOptions(
			'model_annot|a=s'    => \$model_annot_file,
			'model_blast|b=s'    => \$model_blast_file,
			'pfam_results_file|p=s'     => \$pfam_results_file,
			'max_hits|m=i'        => \$max_hits,
			'protein_fasta|f=s'        => \$protein_fasta_file,
			'output|o=s'          => \$output_file,
			'help|h'              => \$help,
);

if (defined($help)) {
	die "\n$0  --model_annot model_genome_annotation_file  --model_blast model_genome_blast_file  --pfam_results_file pfam_results_file  --max_hits max_number_uniref_hits_to_use  --output  path_to_output_file\n\n";
}
unless (defined($protein_fasta_file) && -e $protein_fasta_file) {
	die "\nMust provide the protein fasta file with --protein_fasta flag\n\n";
}
unless (defined($model_blast_file) && -e $model_blast_file) {
	die "\nMust provide raw blast results file with --model_blast flag\n\n";
}
unless (defined($model_annot_file) && -e $model_annot_file) {
	die "\nMust provide model genome annotation file with --model_annot flag\n\n";
}
#unless (defined($pfam_results_file) && -e $pfam_results_file) {
#    die "\nMust provide raw pfam results file with --pfam_results_file flag\n\n";
#}
unless (defined($max_hits)) {
	die "\nMust provide the max number of blast hits to consider with --max_hits flag\n\n";
}
unless (defined($output_file) && !(-e $output_file)) {
	die "\nMust provide output file name with --output flag\n\n";
}

$protein_fasta_file = abs_path($protein_fasta_file);
$model_blast_file = abs_path($model_blast_file);
$model_annot_file = abs_path($model_annot_file);
$pfam_results_file = abs_path($pfam_results_file);
$output_file = abs_path($output_file);

# Determine which gene predictions had transcript evidence.
my %has_transcript_support;
my $in  = Bio::SeqIO->new(-file => "$protein_fasta_file" ,
			  -format => 'Fasta');
while ( my $seq = $in->next_seq() ) {
	my $seq_id = $seq->display_id();
	my $desc = $seq->desc();
	# >Csco26g01220.1 protein AED:0.00 eAED:0.00 QI:0|-1|0|1|-1|1|1|0|69
	# QI:0|-1|0|1|-1|1|1|0|69
	my @elems = split " ", $desc;
	my @qi = split /\|/, $elems[3];
	#print "$seq_id\n";
	#print "$desc\n";
	#print "$elems[0]\n";
	#print "$elems[3]\n";
	#print "$qi[2]\n";
	if ($qi[2] > 0) {
		$has_transcript_support{$seq_id} = 1;
	}
	else {
		$has_transcript_support{$seq_id} = 0;
	}
}



# Read in the model genome functional descriptions.
# This file should consist of protein IDs and functional descriptions separated by tabs.
# The protein IDs should match the IDs that are found in the model genome blast results.
my %model_annots;
open IN, $model_annot_file or die "\nUnable to read the $model_annot_file file.\n\n";
while (my $line = <IN>) {
	chomp $line;
	my ($gene_id, $annot) = split "\t", $line;
	my $short_gene_id;
	if ($gene_id =~ /(\w+)\.\d+/) {
		$short_gene_id = $1;
	}
	$annot =~ s/\, putative\, expressed//;
	$annot =~ s/\, expressed//;
	if ($annot =~ /^Expressed gene of unknown function/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	if ($annot =~ /^Hypothetical gene of unknown function/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	if ($annot =~ /^expressed gene/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	if ($annot =~ /^expressed protein/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	if ($annot =~ /unknown function domain containing protein/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	if ($annot =~ /unknown function containing protein/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	if ($annot =~ /protein of unknown function containing protein/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	if ($annot =~ /domain of unknown function DUF\d+ domain containing protein/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	if ($annot =~ /domain of unknown function, DUF\d+ domain containing protein/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	if ($annot =~ /protein of unknown function, DUF\d+ domain containing protein/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	if ($annot =~ /protein of unknown function DUF\d+ domain containing protein/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	if ($annot =~ /domain of unknown function family protein/) {
		$annot = "Homology to rice gene $short_gene_id";
	}
	$model_annots{$gene_id} = $annot;
}
close IN;

# Functional annotations from the model genome blast results will be prefered over pfam annotations.
# Start with the model genome blast results.
my ($query_name, $subject_name, %annotations, %blast_genes);
my $searchio = Bio::SearchIO->new( -format => 'blast',
					-file   => $model_blast_file,               
    );
while ( my $result = $searchio->next_result() ) {
	my $query_name = $result->query_name;
	my $annotation;
	my $blast_hit;
	my $hit_counter = 0;
	while( my $hit = $result->next_hit ) {
		my $subject_name = $hit->name;
		$hit_counter++;
		$annotation = $model_annots{$subject_name};
		$blast_hit = $subject_name;
		if (defined($annotation) || $hit_counter >= $max_hits) {
			if (defined($annotation)) {
				$annotations{$query_name} = "Arabidopsis blast: $annotation";
			}
			last;
		}        
    }
	if (!defined($annotation) && defined($has_transcript_support{$query_name}) && $has_transcript_support{$query_name} == 1) {
		# If there was no model genome homology but there was transcript support,
		# use this generic description.
		# This may be replaced by a Pfam description.
		$annotation = "Expressed gene of unknown function";
		if (!defined($blast_hit)) {
			$blast_hit = "NA";
		}
		$annotations{$query_name} = $annotation;
	} elsif (!defined($annotation)) {
		# If there was no model genome homology but there was no transcript support,
		# use this generic description.
		# This may be replaced by a Pfam description.
		$annotation = "Hypothetical gene of unknown function";
		if (!defined($blast_hit)) {
			$blast_hit = "NA";
		}
		$annotations{$query_name} = $annotation;
	}
	$blast_genes{$query_name} = $blast_hit;
}

# Functional annotations from the best pfam domain match
# will be used if the model genome annotation was poor or absent.
my %pfam_matched_genes;
open PFAM, $pfam_results_file or die "\nUnable to open pfam results file: $pfam_results_file\n\n";
while ( my $line = <PFAM> ) {
	chomp $line;
	if ($line =~ /^#/) {
		next;
	}
	my @elems = split "\t", $line;
	my $annotation_short = $elems[5];
	my $annotation_long = join(" ", @elems[12,]);
	my $query_name = $elems[0];
	my $score = $elems[8];
	if (exists($pfam_matched_genes{$query_name})) {
		# We have already provided this gene with a pfam domain annotation.
		# We don't need to reassign a weaker pfam domain.
		next;
    }
	if ($annotation_long =~ /^(.+?)\s+Molecular Function:/) {
		$annotation_long = $1;
	}
	if ($annotation_long =~ /rotein of unknown function,* \(*(DUF\d+)\)*/) {
		$annotation_long = $1;
	}
	if ($annotation_long =~ /omain of unknown function,* \(*(DUF\d+)\)*/) {
		$annotation_long = $1;
	}
	$annotation_long =~ s/ domain$//;
	if ($score > 1e-5) {
		# Let's only work with scores better than 1e-5.
		next;
	}
	if ($query_name =~ /__124_/) {
		# The pipes in the rice names have to be re-added.
		$query_name =~ s/__124_/\|/g;
	}
	if (!exists($annotations{$query_name}) || $annotations{$query_name} =~ /unknown function$/) {
		$annotations{$query_name} = "PFAM: $annotation_short domain containing protein";
		$pfam_matched_genes{$query_name} = 1;
	}
}

# Print everything out to a file.
open OUT, ">$output_file" || die "\nUnable to open $output_file for writing.\n\n";
foreach my $id (keys(%annotations)) {
	print OUT "$id\t$blast_genes{$id}\t$annotations{$id}\n";
}
close OUT;

exit;


