#! /usr/bin/perl
# reformat_exonerate_transcript_gff.pl
# This will reformat the exonerate transcript alignments so that they will display correctly in IGV.
# 19 March 2021
# Kevin Childs

use Getopt::Long;
use strict;
use warnings;

my $usage = "\n$0\n    --input_gff  input_gff_from_exonerate\n" .
                  "    --output_gff   reformated_output_gff\n" .
                  "    [--help]\n\n";

my ($input_gff, $output_gff);
my $help;

Getopt::Long::GetOptions( "input_gff=s" => \$input_gff,
                          "output_gff=s" => \$output_gff,
                          "help" => \$help) || die;

if (defined($help)) {
    print $usage;
    exit;
}

if (!defined($input_gff) ||  !(-e $input_gff) ||
    !defined($output_gff) ||  (-e $output_gff)) {
    die $usage;
}

open OLD_GFF, $input_gff or die "\nUnable to open $input_gff for writing.\n\n";
open NEW_GFF, ">$output_gff" or die "\nUnable to open $output_gff for reading.\n\n";

my $id = 0;

while (my $line = <OLD_GFF>) {
    chomp $line;
    #
    if ($line =~ /^#/ || $line =~ /^$/ || $line =~ /^>/) {
	   next;
    }
    #
    my @elems = split "\t", $line;
    #
    if ($elems[2] eq 'similarity' || $elems[2] eq 'splice3' || $elems[2] eq 'cds' ||
	$elems[2] eq 'splice5' || $elems[2] eq 'intron') {
	   next;
    }
    #
    if ($elems[2] eq 'gene' || $elems[2] eq 'transcript') { $elems[2] = "mRNA"; }
    #
    if ($elems[2] eq 'exon') { $elems[2] = "exon"; }
    #
    $elems[1] = "est2genome";
    #
    if ($elems[2] eq 'mRNA') {
        ++$id;
        my $name = "unknown";
        if ($elems[8] =~ /sequence ([\w\.]+) \;/) {
            $name = $1;
        } elsif ($elems[8] =~ /transcript_id \"([\w\.]+)\"\;/) {
            $name = $1;
        }
        $elems[8] = "ID=$id; Name=$name";
    } else {
        $elems[8] = "Parent=$id";
    }
    #
    my $new_line = join("\t", @elems);
    #
    print NEW_GFF "$new_line\n";
}

exit;

