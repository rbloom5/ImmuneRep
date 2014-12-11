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

my ($sims,$vseg,$dseg,$jseg,$cseg,$vshms,$dshms,$jshms,$vtrim,$ctrim)=GatherOptions();

################################ Inputs    ################################

my $fasta_v=VDJFasta->new();
   $fasta_v->loadSeqs($vseg);
my $fasta_j=VDJFasta->new();
   $fasta_j->loadSeqs($jseg);
my $fasta_c=VDJFasta->new();
   $fasta_c->loadSeqs($cseg);

# only load d-seg if a file has been offered
my $fasta_d=VDJFasta->new();
if(-f $dseg){
   $fasta_d->loadSeqs($dseg);
}

my $vgene_count=$fasta_v->getSeqCount();
my $dgene_count=$fasta_d->getSeqCount();
my $jgene_count=$fasta_j->getSeqCount();
my $cgene_count=$fasta_c->getSeqCount();

################################ Analysis  ################################

for(my $v=0;$v<$vgene_count;$v++){
  for(my $j=0;$j<$jgene_count;$j++){
    for(my $c=0;$c<$cgene_count;$c++){
      for(my $s=0;$s<$sims;$s++){
        #my $rand_j=int(rand($fasta_j->getSeqCount()));
        my $rand_d=int(rand($fasta_d->getSeqCount()));
        #my $rand_c=int(rand($fasta_c->getSeqCount()));

        # determine the boundary effects
        my $v3prime_trim=0; #int(rand(3));
        my $d5prime_trim=0; #int(rand(3));
        my $d3prime_trim=0; #int(rand(3));
        my $d5prime_np=int(rand(6));
        my $d3prime_np=int(rand(6));
        my $j5prime_trim=0;

        # construct the header
        my $catseq_header  = $fasta_v->getHeader($v);
          $catseq_header  =~ s/ .*//;
        if(-f $dseg){
          my $catseq_d       = $fasta_d->getHeader($rand_d);
             $catseq_d       =~ s/ .*//;
             $catseq_header .= $catseq_d;
        }
        my $catseq_j       = $fasta_j->getHeader($j);
          $catseq_j       =~ s/ .*//;
          $catseq_header .= $catseq_j;
        my $catseq_c       = $fasta_c->getHeader($c);
          $catseq_c       =~ s/ .*//;
          $catseq_header .= $catseq_c;

        # construct the sequence
        my $catseq_sequence  = $fasta_v->getMutatedSeq($v,$vshms,$vtrim,$v3prime_trim);
          $catseq_sequence .= np_addition($d5prime_np);

          if(-f $dseg){
            $catseq_sequence .= $fasta_d->getMutatedSeq($rand_d,$dshms,$d5prime_trim,$d3prime_trim);
            $catseq_sequence .= np_addition($d3prime_np);
          }
          $catseq_sequence .= $fasta_j->getMutatedSeq($j,$jshms,$j5prime_trim,          0);
          $catseq_sequence .= $fasta_c->getMutatedSeq($c,          0,          0,$ctrim);

        # print the sequence
        print $catseq_header . "\n";
        print $catseq_sequence . "\n";
      }
    }
  }
}

############################### Subroutines ###############################

sub np_addition {
  my($additions)=@_;
  my $out="";
  my @nuc=("C","G","T","A");
  my $add=int(rand($additions));
  for(my $x=0;$x<$add;$x++){
    my $n=int(rand(4));
    $out.=$nuc[$n];
  }
  return $out;
}




sub GatherOptions {
  my $sims   =    0; # simulations per V-segment
  my $vseg   =   "";
  my $dseg   =   "";
  my $jseg   =   "";
  my $cseg   =   "";
  my $vshms  =    0; # shms in V-segment
  my $dshms  =    0; # shms in D-segment
  my $jshms  =    0; # shms in J-segment
  my $vtrim  =    0; # nucleotides cut off the 5-prime end of the V-segment
  my $ctrim  =    0; # nucleotides cut off 3-prime of C-segment
  my $help   =    0;

  GetOptions(
    "--sims=s" => \$sims,
    "--vseg=s" => \$vseg,
    "--dseg=s" => \$dseg,
    "--jseg=s" => \$jseg,
    "--cseg=s" => \$cseg,
    "--vshms=s" => \$vshms,
    "--dshms=s" => \$dshms,
    "--jshms=s" => \$jshms,
    "--vtrim=s" => \$vtrim,
    "--ctrim=s" => \$ctrim,
    "--help=s"  => \$help
  );

  if($help || ($sims == 0) ){
    print "\nUsage: $0\n";
    print " --sims     simulations per V-segment\n";
    print " --vseg     vsegment database\n";
    print " --dseg     vsegment database\n";
    print " --jseg     vsegment database\n";
    print " --cseg     vsegment database\n";
    print " --vshms    shms in V-segment\n";
    print " --dshms    shms in D-segment\n";
    print " --jshms    shms in J-segment\n";
    print " --vtrim    nucleotides cut off 5'prime end of V-segment\n";
    print " --ctrim    nucleotides cut off 3'prime end of C-segment\n";
    exit;
  }
  return ($sims,$vseg,$dseg,$jseg,$cseg,$vshms,$dshms,$jshms,$vtrim,$ctrim);
}

