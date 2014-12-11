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

my ($seqfile)=GatherOptions();

################################ Inputs    ################################

my $igfasta=VDJFasta->new();
   $igfasta->loadSeqs("$seqfile");

################################ Analysis  ################################

   $igfasta->aascFvc2m();

############################### Subroutines ###############################

sub GatherOptions {
  my $file   =   ""; # sequence file
  GetOptions(
    "--file=s"  => \$file
  );
  unless(-f $file ){
    print "\nUsage: $0\n";
    print " --file     amino acid sequence file in fasta format\n";
    exit;
  }
  return ($file);
}

