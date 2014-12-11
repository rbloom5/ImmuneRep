#!/usr/bin/perl -w
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

# Author:  Jacob Glanville 
# Contact: jacob <dot> glanville <at> pfizer <dot> com
# 
# Glanville J, Zhai W, Berka J et al. Precise determination of the diversity 
# of a combinatorial antibody library gives insight into the human immunoglobulin 
# repertoire. Proc Natl Acad Sci USA. 2009;106:20216–20221

############################### Arguments ###############################

my ($seqfile,$annotate_dna,$annotate_aa,$annotate_c2m,$verbose,$vdb,$ddb,$jdb,$cdb)=GatherOptions();

############################### Inputs    ###############################

my $fasta=VDJFasta->new();
   $fasta->loadSeqs($seqfile);
my $corename=stripFastaSuffix($seqfile);
   
############################### Analysis  ###############################

   # germline classification
   if($verbose){ print "Running constant domain classification...\n"; }
   $fasta->getGermlines("C","$corename.C.germdata.txt","dna",$cdb);
   if($verbose){ print "Running J-segment classification...\n"; }
   $fasta->getGermlines("J","$corename.J.germdata.txt","dna",$jdb);
   if($verbose){ print "Running V-segment classification...\n"; }
   $fasta->getGermlines("V","$corename.V.germdata.txt","dna",$vdb);
   if($verbose){ print "Running D-segment classification...\n"; }
   $fasta->getGermlines("D","$corename.D.germdata.txt","dna",$ddb); 
   $fasta->writeDsegCoords("$corename.VDJ.coords.txt");

   if($verbose){ print "Translating all frames...\n"; }
   $fasta->translateAllFramesToFile("$corename.aa.fa");

my $aafasta=VDJFasta->new();
   $aafasta->loadSeqs("$corename.aa.fa");
   if($verbose){ print "Scoring aa sequences for immunoglobulin content...\n"; }
my @scoredHits=$aafasta->igScore();
   if($verbose){ print "Extracting ig-bearing subset...\n"; }
   $aafasta->printSeqSubsetbyList(\@scoredHits,"$corename.wIgs.fa");
   undef($aafasta);

my $igfasta=VDJFasta->new();
   $igfasta->loadSeqs("$corename.wIgs.fa");
   if($verbose){ print "Aligning all sequences to VDJ HMM...\n"; }
my $c2m=$igfasta->aascFvc2m();
my $c2mfasta=VDJFasta->new();
   $c2mfasta->loadSeqs($c2m);
   if($verbose){ print "Identifying translated CDR-H3 with VDJ HMM-boundary QC motifs...\n"; }
   $c2mfasta->qcH3("$corename.H3.acc.txt");
   if($verbose){ print "Identifying translated CDR-L3 with VDJ HMM-boundary QC motifs...\n"; }
   $c2mfasta->qcL3("$corename.L3.acc.txt");
   if($verbose){ print "Identifying translated CDR-H1 with VDJ HMM-boundary QC motifs...\n"; }
   $c2mfasta->qcH1("$corename.H1.acc.txt");
   if($verbose){ print "Identifying translated CDR-H2 with VDJ HMM-boundary QC motifs...\n"; }
   $c2mfasta->qcH2("$corename.H2.acc.txt");
   if($verbose){ print "Identifying translated CDR-L1 with VDJ HMM-boundary QC motifs...\n"; }
   $c2mfasta->qcL1("$corename.L1.acc.txt");
   if($verbose){ print "Identifying translated CDR-L2 with VDJ HMM-boundary QC motifs...\n"; }
   $c2mfasta->qcL2("$corename.L2.acc.txt");

   if($annotate_dna){
     # add header annotations to dna
     if($verbose){ print "Add all annotations to $corename.VDJ.H3.L3.CH1.fa...\n"; }
     $fasta->addHeaderAnnotation("$corename.V.germdata.txt");
     $fasta->addHeaderAnnotation("$corename.D.germdata.txt");
     $fasta->addHeaderAnnotation("$corename.J.germdata.txt");
     $fasta->addHeaderAnnotation("$corename.H3.acc.txt");
     $fasta->addHeaderAnnotation("$corename.L3.acc.txt");
     $fasta->addHeaderAnnotation("$corename.C.germdata.txt");
     $fasta->addHeaderAnnotation("$corename.VDJ.coords.txt"); 
     $fasta->addHeaderAnnotation("$corename.H1.acc.txt");
     $fasta->addHeaderAnnotation("$corename.H2.acc.txt");
     $fasta->addHeaderAnnotation("$corename.L1.acc.txt");
     $fasta->addHeaderAnnotation("$corename.L2.acc.txt");
     $fasta->writeSeqs("$corename.VDJ.H3.L3.CH1.fa");
   } 

   if($annotate_aa){
     if($verbose){ print "Add all annotations to $corename.aa.VDJ.H3.L3.CH1.fa...\n"; }
     $igfasta->addHeaderAnnotation("$corename.V.germdata.txt");
     $igfasta->addHeaderAnnotation("$corename.D.germdata.txt");
     $igfasta->addHeaderAnnotation("$corename.J.germdata.txt");
     $igfasta->addHeaderAnnotation("$corename.H3.acc.txt");
     $igfasta->addHeaderAnnotation("$corename.L3.acc.txt");
     $igfasta->addHeaderAnnotation("$corename.C.germdata.txt");
     $igfasta->addHeaderAnnotation("$corename.VDJ.coords.txt");
     $igfasta->addHeaderAnnotation("$corename.H1.acc.txt");
     $igfasta->addHeaderAnnotation("$corename.H2.acc.txt");
     $igfasta->addHeaderAnnotation("$corename.L1.acc.txt");
     $igfasta->addHeaderAnnotation("$corename.L2.acc.txt");
     $igfasta->writeSeqs("$corename.aa.VDJ.H3.L3.CH1.fa");
   }

   if($annotate_c2m){
     if($verbose){ print "Add all annotations to $corename.aa.VDJ.H3.L3.CH1.c2m...\n"; }
     # add header annotations to translated alignment
     $c2mfasta->addHeaderAnnotation("$corename.V.germdata.txt");
     $c2mfasta->addHeaderAnnotation("$corename.D.germdata.txt");
     $c2mfasta->addHeaderAnnotation("$corename.J.germdata.txt");
     $c2mfasta->addHeaderAnnotation("$corename.H3.acc.txt");
     $c2mfasta->addHeaderAnnotation("$corename.L3.acc.txt");
     $c2mfasta->addHeaderAnnotation("$corename.C.germdata.txt");
     $c2mfasta->addHeaderAnnotation("$corename.VDJ.coords.txt");
     $c2mfasta->addHeaderAnnotation("$corename.H1.acc.txt");
     $c2mfasta->addHeaderAnnotation("$corename.H2.acc.txt");
     $c2mfasta->addHeaderAnnotation("$corename.L1.acc.txt");
     $c2mfasta->addHeaderAnnotation("$corename.L2.acc.txt");
     $c2mfasta->restackRange(94,102);
     $c2mfasta->restackRange(25,35);
     $c2mfasta->restackRange(46,56);
     $c2mfasta->writeSeqs("$corename.aa.VDJ.H3.L3.CH1.c2m");
   }

############################### subroutines ###############################

sub GatherOptions {
  my $file          =""; 
  my $annotate_dna  = 1;
  my $annotate_aa   = 0;
  my $annotate_c2m  = 0;
  my $verbose       = 0;
  my $vdb           ="";
  my $ddb           ="";
  my $jdb           ="";
  my $cdb           ="";

  GetOptions(
     "--file=s" => \$file,
     "--dna=s"  => \$annotate_dna,
     "--aa=s"   => \$annotate_aa,
     "--c2m=s"  => \$annotate_c2m,
     "--vdb=s"  => \$vdb,
     "--ddb=s"  => \$ddb,
     "--jdb=s"  => \$jdb,
     "--cdb=s"  => \$cdb,
     "--verbose=s" => \$verbose
  );

  unless(-f $file){
    print "\nUsage: $0\n";
    print "  --file=seqs.fa     input dna file in fasta\n";
    print "  --dna=1            output a header-annotated dna fasta file\n";
    print "  --aa=1             output a header-annotated aa fasta file\n";
    print "  --c2m=1            output a header-annotated c2m fasta file\n";
    print "  --verbose=1        produce verbose messages while running\n";
    print "  --vdb=\"alt.fa\"   alternative variable segment blast database\n";
    print "  --ddb=\"alt.fa\"   alternative D segment blast database\n";
    print "  --jdb=\"alt.fa\"   alternative J segment blast database\n";
    print "  --cdb=\"alt.fa\"   alternative constant segment blast database\n";
    exit;
  }
  return ($file,$annotate_dna,$annotate_aa,$annotate_c2m,$verbose,$vdb,$ddb,$jdb,$cdb);
}

sub stripFastaSuffix {
  my($corename)=@_;
   $corename=~s/.fa$//;
   $corename=~s/.fasta$//;
   $corename=~s/.fna$//;
   return $corename;
}

