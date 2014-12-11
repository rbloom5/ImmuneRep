#!/usr/bin/perl -w

# Author:  Jacob Glanville 
# Contact: jacob <dot> glanville <at> pfizer <dot> com
#
# Glanville J, Zhai W, Berka J et al. Precise determination of the diversity 
# of a combinatorial antibody library gives insight into the human immunoglobulin 
# repertoire. Proc Natl Acad Sci USA. 2009;106:20216–20221

use Getopt::Long;
use vars qw($rootdir);
BEGIN{
  use File::Basename;
  use Cwd 'abs_path';
  $rootdir = dirname(dirname(abs_path($0)));
};
use lib "$rootdir/lib";
use VDJFasta;
use strict;

################################ Arguments ################################

my ($seqfile,$outfile)=GatherOptions();

################################ Inputs    ###########################

my $fasta=VDJFasta->new();
   $fasta->loadSeqs($seqfile);

################################ Analysis  ################################

   $fasta->restackRange(94,102);
   $fasta->restackRange(25,35);
   $fasta->restackRange(46,56);
   $fasta->writeSeqs($outfile);

############################### Subroutines ###############################

sub GatherOptions {
  my $file   =   ""; # sequence file
  my $outfile=   "";
  GetOptions(
    "--file=s"  => \$file,
    "--outfile=s" => \$outfile
  );
  unless(-f $file ){
    printUsage();
  }
  #if($outfile=~m/^ *$/){
  #  printUsage();
  #}
  return ($file,$outfile);
}

sub printUsage {
    print "\nUsage: $0\n";
    print " --file     amino acid VDJ hmm alignment\n";
    print " --outfile  tabbed results file\n";
    exit;
}
