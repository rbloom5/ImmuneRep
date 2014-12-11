package VDJFasta;
use strict;
use vars qw($rootdir);
use Cwd;

# Author:  Jacob Glanville 
# Contributing authors: Yun He & Marina Sirota 
# Contact: jacob <dot> glanville <at> pfizer <dot> com
#
# Glanville J, Zhai W, Berka J et al. Precise determination of the diversity 
# of a combinatorial antibody library gives insight into the human immunoglobulin 
# repertoire. Proc Natl Acad Sci USA. 2009;106:20216–20221

BEGIN {
  use File::Basename; 
  use Cwd 'abs_path';
  $rootdir = dirname(dirname(abs_path($0)));
};


# Constructor ----------------------------------------------
sub new {
  my ($class) = @_;
  my $self = {};

  my $raman = dirname(dirname(dirname(abs_path($0))));

  $self->{filename}     = "";
  $self->{headers}      = [];
  $self->{sequence}     = [];
  $self->{germline}     = [];
  $self->{nseqs}        =  0;
  $self->{mids}         = {};

  $self->{accVsegQstart}  = {}; # example: 124
  $self->{accVsegQend}    = {}; # example: 417
  $self->{accJsegQstart}  = {};
  $self->{accJsegQend}    = {};
  $self->{accDsegQstart}  = {};
  $self->{accDsegQend}    = {};
  $self->{accCsegQstart}  = {};
  $self->{accCsegQend}    = {};

  $self->{ranVseg}        =  0; # boolean
  $self->{ranDseg}        =  0;
  $self->{ranJseg}        =  0;
  $self->{ranCseg}        =  0;

  $self->{Vseg}           = {}; # example: IGHV1-46 293 3
  $self->{Dseg}           = {};
  $self->{Jseg}           = {};
  $self->{Cseg}           = {};

  $self->{V}              = []; # example: IGHV1-46 293 3
  $self->{D}              = [];
  $self->{J}              = [];
  $self->{C}              = [];
  $self->{H3}             = [];
  $self->{L3}             = [];
  $self->{coords}         = [];

  $self->{dnaVsegdb}    = $raman . "/vdjfasta/db/imgt.V.dna.fa"; 
  $self->{dnaJsegdb}    = $raman . "/vdjfasta/db/imgt.J.dna.fa";
  $self->{dnaDsegdb}    = $raman . "/vdjfasta/db/imgt.HD.dna.nr.fa";
  $self->{dnaCsegdb}    = $raman . "/vdjfasta/db/imgt.CH1.dna.nr.fa";

  #$self->{dnaVsegdb}    = "/tools/vdjfasta/db/imgt.V.dna.fa"; 
  #$self->{dnaJsegdb}    = "/tools/vdjfasta/db/imgt.J.dna.fa";
  #$self->{dnaDsegdb}    = "/tools/vdjfasta/db/imgt.HD.dna.nr.fa";
  #$self->{dnaCsegdb}    = "/tools/vdjfasta/db/imgt.CH1.dna.nr.fa";

  $self->{VhVkHMM}      = $raman . "/vdjfasta/db/Vh-linker-Vk.hmm";
  #$self->{VhVkHMM}      = "/tools/vdjfasta/db/Vh-linker-Vk.hmm";
 
  $self->{blast}        = $raman . "/ncbi-blast/bin/blastn";
  $self->{hmmsearch}    = $raman . "/hmmer/binaries/hmmsearch";
  $self->{hmmalign}     = $raman . "/hmmer/binaries/hmmalign";

  $self->{verbose}      = 0;

  bless $self,'VDJFasta';
  return $self;
}

# Methods --------------------------------------------------

sub parseHeaderAnnotations {
  my($self)=@_;
  print "stub\n";

  my $annotated=0;

  for(my $x=0;$x<scalar(@{$self->{headers}});$x++){
    my @fields = split(/;/,${$self->{headers}}[$x]);
    if(scalar(@fields)<6){
      print "Not annotated!\n";
    }else{
      ${$self->{V}}[$x]     = $fields[1];  
      ${$self->{D}}[$x]     = $fields[2];
      ${$self->{J}}[$x]     = $fields[3];
      ${$self->{H3}}[$x]    = $fields[4];
      ${$self->{L3}}[$x]    = $fields[5];
      ${$self->{C}}[$x]     = $fields[6];
      ${$self->{coords}}[$x]= $fields[7];
      ${$self->{H1}}[$x]    = $fields[8];
      ${$self->{H2}}[$x]    = $fields[9];
      ${$self->{L1}}[$x]    = $fields[10];
      ${$self->{L2}}[$x]    = $fields[11];
      $annotated=1;
    }
  }
  return $annotated;
}

sub reportGermFreqs {
  my($self,$domain,$minlength,$min_shm,$max_shm,$isRedundant)=@_;
  print "stub\n";
}

sub getGermlines {
  my($self,$domain,$outfile,$sequence_type,$custom_reference_db)=@_;
  my $blast_type="blastn";
  if($sequence_type ne "dna"){
    # $blast_type="tblastn";
    print "This functionality is not yet active. Go away.\n";
    exit; 
  }
  my $reference_db="";
  my $evalue="";
  if($domain eq "V"){
    $reference_db=$self->{dnaVsegdb};
    $evalue=10e-10;
  }elsif($domain eq "J"){
    $reference_db=$self->{dnaJsegdb};
    $evalue=0.001;
  }elsif($domain eq "D"){
    $reference_db=$self->{dnaDsegdb};
    $evalue=0.01;
  }elsif($domain eq "C"){
    $reference_db=$self->{dnaCsegdb};
    $evalue=0.001;
  }else{
    print " (V J D C)\n";
    exit;
  }
  if(defined($custom_reference_db)){
    if(-f $custom_reference_db){
      $reference_db=$custom_reference_db;
    }
  }

  if($domain eq "D"){
    $self->batchBlastParseDsegClassifier($reference_db,$evalue,$outfile,$blast_type,$domain);
  }else{  
    $self->batchBlastParseClassifier($reference_db,$evalue,$outfile,$blast_type,$domain);
  }
}

sub batchBlastParseClassifier {
  my($self,$db,$e,$outfile,$blast_type,$domain)=@_;
  my @lines=$self->getBlastLines($self->{filename},$db,$e,$blast_type);

  my %accessions=();
  my %acc_tophit=();

  for(my $x=0;$x<scalar(@lines);$x++){
    my($acc,$match,$pid,$alnlength,$mismatch,$gaps,$qstart,$qend,$mstart,$mend,$eval,$bit)=split(/\t/,$lines[$x]);
    my $hitname=$self->parse_hitname($match);

    if(!(defined($acc_tophit{$acc}))){    
      $acc_tophit{$acc} = $hitname;
      my $prob = 1;
      $accessions{$acc} = $hitname . " " . $alnlength . " " . $mismatch . " ";
      if($domain eq "V"){
        ${$self->{accVsegQstart}}{$acc}  = $qstart;
        ${$self->{accVsegQend}}{$acc}    = $qend;
      }elsif($domain eq "J"){
        ${$self->{accJsegQstart}}{$acc}  = $qstart;
        ${$self->{accJsegQend}}{$acc}    = $qend;
      }elsif($domain eq "C"){
        ${$self->{accCsegQstart}}{$acc}  = $qstart;
        ${$self->{accCsegQend}}{$acc}    = $qend;
      }
    }elsif($acc_tophit{$acc} ne $hitname){
      my $prob = 1;
      if(defined($accessions{$acc})){
        my($best_match,$best_alnlength,$best_mismatch)=split(/ /,$accessions{$acc});
        $prob = $self->miscall_probability($best_alnlength,$alnlength,$best_mismatch,$mismatch);
        $prob = (int(1000 * $prob))/1000;
      }
      if($prob > 0.01){
        $accessions{$acc} .= $hitname . " " . $alnlength . " " . $mismatch . " ";
      }
    }
  }

  my @keys=keys %accessions;
  open(OUT,">$outfile");
  for(my $k=0;$k<scalar(@keys);$k++){
    print OUT $keys[$k] . "\t" . $accessions{$keys[$k]} . "\n";
    ${$self->{$domain . "seg"}}{$keys[$k]}=$accessions{$keys[$k]};
  }
  print OUT "WRITEDONEFLAG" . "\t" . "\n";
  close(OUT);
  $self->{"ran" . $domain . "seg"}=1;
}

sub writeDsegCoords {
  my($self,$outfile)=@_;
  my @accs= keys %{$self->{accDsegQstart}};

  open(OUTFILE,">$outfile");
  for(my $a=0;$a<scalar(@accs);$a++){
    print OUTFILE  $accs[$a] . "\t"
          . ${$self->{accVsegQstart}}{$accs[$a]} . " "
          . ${$self->{accVsegQend}}{$accs[$a]} . " "
          . ${$self->{accDsegQstart}}{$accs[$a]} . " "
          . ${$self->{accDsegQend}}{$accs[$a]} . " "
          . ${$self->{accJsegQstart}}{$accs[$a]} . " "
          . ${$self->{accJsegQend}}{$accs[$a]} . "\n";
  }
  close(OUTFILE);
}

sub isDsegSandwiched {
  my($self,$acc,$dseg_qstart,$dseg_qend)=@_;

  my $answer=1;

  # if not defined, then exit
  unless(defined(${$self->{accVsegQstart}}{$acc}) && defined(${$self->{accVsegQend}}{$acc})
         && ${$self->{accJsegQstart}}{$acc} && ${$self->{accJsegQend}}{$acc} ){
    return 0;
  }

  # 1. get middle of each segment
  my $v_average= (${$self->{accVsegQstart}}{$acc} + ${$self->{accVsegQend}}{$acc}) /2;
  my $j_average= (${$self->{accJsegQstart}}{$acc} + ${$self->{accJsegQend}}{$acc}) /2;
  my $d_average= ($dseg_qstart + $dseg_qend) / 2;
  
  # 2. determine orientation
  my $orientation="forward";
  my $vseg_3prime=${$self->{accVsegQend}}{$acc};
  my $jseg_5prime=${$self->{accJsegQstart}}{$acc};
  if( $v_average > $j_average ){
    $orientation="reverse";
    $vseg_3prime=${$self->{accVsegQstart}}{$acc};
    $jseg_5prime=${$self->{accJsegQend}}{$acc};
  }

  # get lower
  my $lower_d=$dseg_qstart;
  my $higher_d=$dseg_qend;
  if($dseg_qend<$lower_d){
    $lower_d=$dseg_qend;
    $higher_d=$dseg_qend;
  }

  # make sure D-seg is sandwiched
  if( ($v_average < $d_average) & ($d_average < $j_average) ) {
    # state is forward
    my $v_overlap= ${$self->{accVsegQend}}{$acc} - $lower_d;
    my $j_overlap= $higher_d - ${$self->{accJsegQstart}}{$acc};
    if($v_overlap>6){
      $answer=0;
    }
    if($j_overlap>6){
      $answer=0;
    }
  }elsif( ($v_average > $d_average) & ($d_average > $j_average)){
    # state is reversed
    my $v_overlap= $higher_d - ${$self->{accVsegQstart}}{$acc};
    my $j_overlap= ${$self->{accJsegQend}}{$acc} - $lower_d;
    if($v_overlap>6){
      $answer=0;
    }
    if($j_overlap>6){
      $answer=0;
    }
  }else{
    $answer=0;
  }

  return $answer;
}

sub batchBlastParseDsegClassifier {
  my($self,$db,$e,$outfile,$blast_type,$domain)=@_;
  unless($self->{ranVseg} && $self->{ranJseg}){
    print "V-seg and J-seg need to run first. Exiting...\n";
    exit;
  }
  my @lines=$self->getBlastLines($self->{filename},$db,$e,$blast_type);
 
  my %accessions=();
  my %acc_tophit=();

  for(my $x=0;$x<scalar(@lines);$x++){
    my($acc,$match,$pid,$alnlength,$mismatch,$gaps,$qstart,$qend,$mstart,$mend,$eval,$bit)=split(/\t/,$lines[$x]);
    my $hitname=$self->parse_hitname($match);

    if($self->isDsegSandwiched($acc,$qstart,$qend)){
      if(!(defined($acc_tophit{$acc}))){
        $acc_tophit{$acc} = $hitname;
        my $prob = 1;
        $accessions{$acc} = $hitname . " " . $alnlength . " " . $mismatch . " ";
        ${$self->{accDsegQstart}}{$acc}  = $qstart;
        ${$self->{accDsegQend}}{$acc}    = $qend;
      }elsif($acc_tophit{$acc} ne $hitname){
        my $prob = 1;
        if(defined($accessions{$acc})){
          my($best_match,$best_alnlength,$best_mismatch)=split(/ /,$accessions{$acc});
          $prob = $self->miscall_probability($best_alnlength,$alnlength,$best_mismatch,$mismatch);
          $prob = (int(1000 * $prob))/1000;
        }
        if($prob > 0.01){
          $accessions{$acc} .= $hitname . " " . $alnlength . " " . $mismatch . " ";
        }
      }
    }
  }

  my @keys=keys %accessions;
  open(OUT,">$outfile");
  for(my $k=0;$k<scalar(@keys);$k++){
    print OUT $keys[$k] . "\t" . $accessions{$keys[$k]} . "\n";
    ${$self->{$domain . "seg"}}{$keys[$k]}=$accessions{$keys[$k]};
  }
  print OUT "WRITEDONEFLAG" . "\t" . "\n";
  close(OUT);
  $self->{"ran" . $domain . "seg"}=1;
}

sub parse_hitname {
  my($self,$hitname)=@_;
  $hitname=~s/\*.*//;
  if($hitname=~m/IGH[GMEA]/){
    $hitname=~s/[0-9]*$//;
    $hitname=~s/P$//;
  }
  return $hitname;
}

sub miscall_probability{
  my($self,$sequence_length,$sequence_length2,$mutations1,$mutations2)=@_;
  my $length_diff=0;
  if($sequence_length2<$sequence_length){
    $length_diff=$sequence_length - $sequence_length2;
  }
  my $specific_mutations=$mutations2 - $mutations1 + $length_diff;
  my $total_mutations=$mutations2 + $specific_mutations;
  my $result=$self->nCk( ($sequence_length - $specific_mutations),($total_mutations - $specific_mutations) );
  my $result2=$self->nCk($sequence_length,$total_mutations);
  my $probability=($result/$result2);
  return $probability;
}

sub nCk {
  my($self,$n,$k)=@_;
  my $result=0; 
  #   n!
  # -----
  # k!(n-k)!
  my $nminusk=$n-$k;
  my @n_factorial=$self->factorial($n);
  my @nmink_factorial=$self->factorial($nminusk);

  my $numerator=1;
  my $denominator=1;
  for(my $x=1;$x<scalar(@n_factorial);$x++){
    if( $n_factorial[$x] && $nmink_factorial[$x] ){
      $n_factorial[$x]=0;
    }else{
      $numerator*=$x;
    }
  }
  for(my $x=$k;$x>1;$x--){
    $denominator *= $x;
  }
  $result=($numerator/$denominator);
  return $result;
}

sub factorial {
  my($self,$n)=@_;
  my @factorial=();
  while($n>0){
    $factorial[$n]=1;
    $n--;
  }
  return @factorial;
}

sub aascFvc2m {
  my($self)=@_;
  my $name = $self->{filename};
     $name =~ s/.[Ff][Nn]*[Aa][Ss]*[Tt]*[Aa]*$//;
  my $unique_file = $name . ".unique.fa";
  my $stockfile   = $name . ".stock";
  my $a2m_file    = $name . ".Vh-gs-Vk.a2m";
  my $c2m_file    = $name . ".Vh-gs-Vk.c2m";
  $self->makeAccessionsUnique();
  $self->writeSeqs($unique_file);
  $self->scFvAlign($unique_file,$stockfile);
  $self->stockholmToFasta($stockfile,$a2m_file);
  $self->cropToModel($a2m_file,$c2m_file);
  return $c2m_file;
}

sub cropToModel {
  my($self,$a2m,$c2m)=@_;

  open(OUT,">$c2m");
  my $seqname="";
  my $seqseq="";

  my $x=-1;
  open MSA, $a2m;
  while (<MSA>) {
   chomp;
   if (/^>/) {
     my $tempname = $_;

     if($seqname ne ""){
       print OUT $seqname . "\n";
       $seqseq=~s/^[\.a-z]*//;
       $seqseq=~s/[\.a-z]*$//;
       print OUT $seqseq  . "\n";
       $seqseq="";
     }

     $seqname=$tempname;
     $x++;
   } else{
     my $tempseq = $_;
     $seqseq .= $tempseq;
   }
  }
  close(MSA);
  $seqseq=~s/^[\.a-z]*//;
  $seqseq=~s/[\.a-z]*$//;
  print OUT $seqname . "\n";
  print OUT $seqseq  . "\n";
  close(OUT);
}

sub stockholmToFasta {
  my($self,$stockfile,$outfile)=@_;
  my %seqs = ();
  my $start=0;

  open(MSA,$stockfile);
  while(<MSA>){
    my $line=$_;
    chop($line);
    if($line =~ m/#/){
      $start=1;
    }
    if($start){
      if( !(($line =~ m/^#/) || ($line =~ m/^ *$/) || ($line =~ m/^\/\//)) ){
        my @fields = split(/  */,$line);
        if(exists($seqs{$fields[0]})){
          $seqs{$fields[0]} .= $fields[1];
        } else {
          $seqs{$fields[0]} = $fields[1];
        }
      }
    }
  }
  close(MSA);

  open(OUT,">$outfile");
  while ( my($key,$value) = each(%seqs)){
    print OUT ">$key\n";
    print OUT "$value\n";
  }
  close(OUT);
}

sub makeAccessionsUnique {
  my($self)=@_;
  my %accessions_already_seen=();
  
  for(my $s=0;$s<$self->{nseqs};$s++){
    ${$self->{headers}}[$s]=~s/ /_/g;
    if(defined($accessions_already_seen{${$self->{headers}}[$s]})){
      while(defined($accessions_already_seen{${$self->{headers}}[$s]})){
        ${$self->{headers}}[$s].= "." . int(rand(100000));
      }
    }
    $accessions_already_seen{${$self->{headers}}[$s]}=1;
  }
}

sub scFvAlign {
  my($self,$seqfile,$outfile)=@_;  
  my $cmd = $self->{hmmalign} 
          . "--amino"  #--allcol
          . $self->{VhVkHMM} 
          . " $seqfile > $outfile";
  `$cmd`;
}

sub igScore {
  my($self)=@_;

  # score all frames with hmm
  my $cmd = $self->{hmmsearch} 
          . " -E 1e-10 -A 0 " 
          . $self->{VhVkHMM} 
          . " " . $self->{filename} 
          . " > " . $self->{filename} 
          . ".1e-10.score.txt";
  `$cmd`;
 
  # generate hmmscore hits
  my @scored_hits=();
  open(FILE,$self->{filename} . ".1e-10.score.txt");
  my @lines=<FILE>;
  close(FILE);
  chomp(@lines);
  my $recording=0;
  for(my $x=0;$x<scalar(@lines);$x++){
    if($lines[$x]=~m/^>>/){
      $lines[$x]=~s/^>>* *//;
      $lines[$x]=~s/  *$//;
      push @scored_hits,$lines[$x];
    }
  }
  return @scored_hits;
}

sub addHeaderAnnotation {
  my($self,$annotation_file)=@_;

  # load file (accession	annotation)
  my %acc2cdr=();
  open(FILE,$annotation_file);
  my @lines=<FILE>;
  close(FILE);
  chomp(@lines);
  for(my $i=0;$i<scalar(@lines);$i++){
    my @fields=split(/\t/,$lines[$i]);
    if(defined($fields[1])){
      $fields[1]=~s/ *$//;
    }else{
      $fields[1]="";
    }
    $fields[0]=~s/^> *//;
    $fields[0]=~s/ *//;
    $acc2cdr{$fields[0]}=$fields[1];
  }

  # apply annotations
  for(my $n=0;$n<$self->{nseqs};$n++){
    ${$self->{headers}}[$n] .= ";"; 
    if(defined($acc2cdr{$self->getAccession($n)})){
      ${$self->{headers}}[$n] .= $acc2cdr{$self->getAccession($n)};
    }
  }
}

sub getBlastLines {
  my($self,$seqfile,$blastdb,$eval,$blast_type)=@_;
  my $command=$self->{blast} . " -task blastn-short -query $seqfile -db $blastdb -outfmt 6 -evalue $eval";
  my @lines=`$command`;
  return @lines;
}

sub getHeader {
  my($self,$c)=@_;
  return ${$self->{headers}}[$c];
}

sub getAccession {
  my($self,$c)=@_;
  my $header=$self->getHeader($c);
  $header=~s/ .*//;
  $header=~s/;.*//;
  $header=~s/^>//;
  $header=~s/_frame_[0123-]* *$//;
  return $header;
}

sub loadSeqs {
  my ($self,$file)=@_;
  if(-e $file){
    $self->{filename}=$file;
  } else {
    print "File does not exist.\n";
    return 0;
  }

  # load seqs into filename, headers and sequence
  my $x=-1;
  open FASTA, $self->{filename};
  while (<FASTA>) {
    chomp;
    if (/^>/) {
      push @{$self->{headers}},$_;
      $x++;
    } else {
      ${$self->{sequence}}[$x] .= $_;
    }
  }
  $x++;
  $self->{nseqs}=$x;
}

sub addSeq {
  my($self,$header,$sequence)=@_;
  push @{$self->{headers}},$header;
  push @{$self->{sequence}},$sequence;
}

sub getSeqLen {
  my($self,$c)=@_;
  my $seq=${$self->{sequence}}[$c];
  my $count=length($seq);
  return $count;
}

sub printSeqs {
  my($self)=@_;
  for(my $s=0;$s<$self->{nseqs};$s++){
    $self->printSeq($s);
  }
}

sub printSeqRange {
  my($self,$start,$stop)=@_;
  for(my $s=$start;$s<$stop;$s++){
    $self->printSeq($s);
  }
}

sub printSeq {
  my($self,$s)=@_;
  print ${$self->{headers}}[$s] . "\n";
  print ${$self->{sequence}}[$s] . "\n";
}

sub printUniqueSeqs {
  my($self)=@_;
  my $id=1;
  for(my $s=0;$s<$self->{nseqs};$s++){
    my $header=">$id-";
    $id++;
    my $header_part2 = ${$self->{headers}}[$s];
       $header_part2 =~ s/^>//;
       $header .= $header_part2;
    print $header . "\n";
    print ${$self->{sequence}}[$s] . "\n";
  }
}

sub getCDRHashFromHeaders {
  my($self,$cdr3)=@_;
  # assumes the headers have been encoded with H3 in field 5
  my $cdr3_field=4;
  if($cdr3 eq "L3"){
    $cdr3_field=5;
  }
  print "Fields is $cdr3_field\n";
  my %cdr_hash=();

  for(my $s=0;$s<$self->{nseqs};$s++){
    my @fields=split(/;/,${$self->{headers}}[$s]);
    if(defined($fields[$cdr3_field])){
      if(length($fields[$cdr3_field]) > 0){
        if(defined($cdr_hash{$fields[$cdr3_field]})){
          $cdr_hash{$fields[$cdr3_field]}.="," . $s;
        }else{
          $cdr_hash{$fields[$cdr3_field]} = $s;
        }
      }
    }
  }
  return %cdr_hash;
}

sub printSeqSubsetbyList {
  my($self,$accession_list,$outfile)=@_;

  my %selected_accessions=();
  my @lines=@$accession_list;

  open(OUTFILE,">$outfile");

  for(my $x=0;$x<scalar(@lines);$x++){
    my @data=split(/\t/,$lines[$x]);
    #print "Setting up |" . $data[0] . "|\n";
    $selected_accessions{$data[0]}=1;
  } 
 
  for(my $s=0;$s<$self->{nseqs};$s++){
    my $acc=${$self->{headers}}[$s];
    $acc=~s/ .*//;
    $acc=~s/^>//;
    #print "Checking on |$acc|\n";
    if($selected_accessions{$acc}){
      print OUTFILE ${$self->{headers}}[$s] . "\n";
      print OUTFILE ${$self->{sequence}}[$s] . "\n";
    }
  }
  close(OUTFILE);
}

sub getSeq{
  my($self,$seq)=@_;
  return (${$self->{headers}}[$seq],${$self->{sequence}}[$seq]);
}

sub getSeqCount{
  my($self)=@_;
  return $self->{nseqs};
}

sub writeSeq {
  my($self,$seq,$file)=@_;
  open(FILE,">$file");
  print FILE ${$self->{headers}}[$seq]  . "\n";
  print FILE ${$self->{sequence}}[$seq] . "\n";
  close(FILE);
  return $file;
}

sub writeSeqs {
  my($self,$file)=@_;
  open(FILE,">$file");
  for(my $s=0;$s<$self->{nseqs};$s++){
    print FILE ${$self->{headers}}[$s]  . "\n";
    print FILE ${$self->{sequence}}[$s] . "\n";
  }
  close(FILE);
}

sub writeSeqRange {
  my($self,$file,$start,$stop)=@_;
  open(FILE,">$file");
  if($stop>$self->{nseqs}){
    $stop=$self->{nseqs};
  }
  for(my $s=$start;$s<$stop;$s++){
    print FILE ${$self->{headers}}[$s]  . "\n";
    print FILE ${$self->{sequence}}[$s] . "\n";
  }
  close(FILE);
}

sub getMutatedSeq {
  my($self,$seq,$nmutations,$trim5prime,$trim3prime)=@_;
  my $mutseq=${$self->{sequence}}[$seq];

  # trime 5 prime and 3 prime
  if($trim5prime){
    $mutseq=substr($mutseq,$trim5prime);
  }
  if($trim3prime){
    $mutseq=substr($mutseq,0,(0 - $trim3prime));
  }

  # perform mutations
  my @chars=('C','A','T','G');
  my @res=split(/ */,$mutseq);
  my $mut_down=$nmutations;
  while($mut_down>0){
    my $position=int(rand(scalar(@res)));
    my $new=$chars[int(rand(4))];
    while($new eq uc($res[$position])){
      $new=$chars[int(rand(4))];
    }
    $res[$position]=$new;
    $mut_down--;
  }
  $mutseq=join('',@res);

  return $mutseq;
}

sub writeMutatedSeq {
  my($self,$seq,$file,$nmutations,$trim5prime,$trim3prime)=@_;

  my $mutseq=$self->getMutatedSeq($seq,$nmutations,$trim5prime,$trim3prime);
  
  open(FILE,">$file");
  print FILE ${$self->{headers}}[$seq]  . " with $nmutations mutations\n";
  print FILE $mutseq . "\n";
  close(FILE);
  return $file;
}



sub hmmSearch {
  my($self,$seqfile,$hmm,$eval)=@_;
  my $hmmsearch=$self->{hmmsearch};
  my @results=`$hmmsearch -E $eval $hmm $seqfile`;
  my @frames=();
  my $print=1;
  for(my $x=0;$x<scalar(@results);$x++){
    if($results[$x]=~m/^frame/){
      if($results[$x]=~m/domain/){
        my $line=$results[$x];
        $line=~s/^frame//;
        $line=~s/: .*//;
        chomp($line);
        push @frames,$line;
      }
    }
  }
  return @frames;
}

sub getPercentGc {
  my($self,$c)=@_;
  my ($header,$oligo)=$self->getSeq($c);
  my @chars=split(/ */,$oligo);
  my $gc=0;
  for(my $g=0;$g<scalar(@chars);$g++){
    if($chars[$g]=~m/[gcGC]/){
      $gc++;
    }
  }
  my $percent_gc=0;
  if(scalar(@chars)>0){
    $percent_gc=((  int(1000 * ($gc/scalar(@chars))) )/10);
  }
  return $percent_gc;
}

sub rounded {
  my ($self,$value)=@_;
  my $rounded=((  int(10 * $value) )/10);
  return $rounded;
}

sub getReverseStrand {
  my($self,$c)=@_;

  my $seq=uc($self->getSequence($c)); 
  my @nucs=split(/ */,$seq);
  my $len=scalar(@nucs);
  my $reverse="";
  for(my $x=($len-1);$x>=0;$x--){
    if($nucs[$x] eq "T"){
      $reverse.="A";
    }elsif($nucs[$x] eq "A"){
      $reverse.="T";
    }elsif($nucs[$x] eq "G"){
      $reverse.="C";
    }elsif($nucs[$x] eq "C"){
      $reverse.="G";
    }else{
      $reverse.="X";
    }
  }
  return $reverse;
}

sub getSequence {
  my($self,$seq)=@_;
  return ${$self->{sequence}}[$seq];
}

sub setSequence {
  my($self,$c,$sequence)=@_;
  return ${$self->{sequence}}[$c]=$sequence;
}

sub codon2aa {
  my($codon) = @_;

  if ( $codon =~ /^TC/i ) { return 'S' } # Serine
  elsif ( $codon =~ /TT[CT]/i ) { return 'F' } # Phenylalanine
  elsif ( $codon =~ /TT[AG]/i ) { return 'L' } # Leucine
  elsif ( $codon =~ /TA[CT]/i ) { return 'Y' } # Tyrosine
  elsif ( $codon =~ /TA[AG]/i ) { return 'X' } # Stop
  elsif ( $codon =~ /TG[CT]/i ) { return 'C' } # Cysteine
  elsif ( $codon =~ /TGA/i ) { return 'X' } # Stop
  elsif ( $codon =~ /TGG/i ) { return 'W' } # Tryptophan
  elsif ( $codon =~ /^CT/i ) { return 'L' } # Leucine
  elsif ( $codon =~ /^CC/i ) { return 'P' } # Proline
  elsif ( $codon =~ /CA[CT]/i ) { return 'H' } # Histidine
  elsif ( $codon =~ /CA[AG]/i ) { return 'Q' } # Glutamine
  elsif ( $codon =~ /^CG/i ) { return 'R' } # Arginine
  elsif ( $codon =~ /AT[ACT]/i ) { return 'I' } # Isoleucine
  elsif ( $codon =~ /ATG/i ) { return 'M' } # Methionine
  elsif ( $codon =~ /^AC/i ) { return 'T' } # Threonine
  elsif ( $codon =~ /AA[CT]/i ) { return 'N' } # Asparagine
  elsif ( $codon =~ /AA[AG]/i ) { return 'K' } # Lysine
  elsif ( $codon =~ /AG[CT]/i ) { return 'S' } # Serine
  elsif ( $codon =~ /AG[AG]/i ) { return 'R' } # Arginine
  elsif ( $codon =~ /^GT/i ) { return 'V' } # Valine
  elsif ( $codon =~ /^GC/i ) { return 'A' } # Alanine
  elsif ( $codon =~ /GA[CT]/i ) { return 'D' } # Aspartic Acid
  elsif ( $codon =~ /GA[AG]/i ) { return 'E' } # Glutamic Acid
  elsif ( $codon =~ /^GG/i ) { return 'G' } # Glycine
  elsif ( $codon =~ m/N/ ) { return 'X' }   # N - unknown character
  else { return "X"; } # bad codon
}

sub dna2aa {
  my($self,$seq,$frame)=@_;

  # header
  my $header=$self->getHeader($seq);
     $header=~s/ .*//;
     $header.="_frame_$frame";

  # frame options: 
  my @nucs=();
  my $startframe=$frame;
  if($frame>=0){
    @nucs=split(/ */,$self->getSequence($seq));
  }else{
    @nucs=split(/ */,$self->getReverseStrand($seq));
    $startframe=-1 - $frame;
  }

  # get translation
  my $translation="";
  for(my $n=$startframe;$n<scalar(@nucs);$n+=3){
    unless(defined($nucs[($n+1)])){
      $nucs[($n+1)]="N";
    }
    unless(defined($nucs[($n+2)])){
      $nucs[($n+2)]="N";
    }
    my $codon=$nucs[$n] . $nucs[($n+1)] . $nucs[($n+2)];
    $translation.=codon2aa($codon);
  }
  return($header,$translation);
}

sub translateAllFramesToFile {
  my($self,$outfile)=@_;
  open(FILE,">$outfile");
  for(my $s=0;$s<$self->{nseqs};$s++){ 
    for(my $f=-3;$f<3;$f++){
      my($header,$seq)=$self->dna2aa($s,$f);
      print FILE $header .  "\n";
      print FILE $seq . "\n";  
    }
  }
  close(FILE);  
}

sub restackRange {
  my($self,$range_start,$range_stop)=@_;
  # restack the residues within this range
  # n_right_stack_residues: the number of residues to stack to the end of the range
  # all others are stacked from left to right

  # first determine the longest segment in the range
  my $longest=0;
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $segment=$self->getHMMCOLRange($n,$range_start,$range_stop);
       $segment=~s/\.*//g;
       $segment=~s/^[A-Z-]//;
       $segment=~s/[A-Z-]$//;
    if(length($segment)>$longest){
      $longest=length($segment);
    }
  }
  
  # next proceed through all sequences, restacking the range and filling in the extra
  # space with enough dots to fit the $longest sequence
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $segment_before=$self->getHMMCOLRange($n,0,$range_start);
    my $segment_after =$self->getHMMCOLRange($n,$range_stop,length($self->getSeq($n)));
       $segment_before=~s/[a-z\.]*$//;
       $segment_after=~s/^[a-z\.]*//;

    # obtain segment
    my $segment=$self->getHMMCOLRange($n,$range_start,$range_stop);
       $segment=~s/\.*//g;
       $segment=~s/^[A-Z-]//; # remove boundary position, as it is represented in segment_before
       $segment=~s/[A-Z-]$//; # remove boundary position, as it is represented in segment_after
    # get hmmcol count
    my $hmmcols=$segment;
       $hmmcols=~s/[a-z]//g;
    my $hmmcol_count=length($hmmcols);
    # split in half
    my $halfway = int ( (length($segment) + 1)/2 ); #half way, round up
    my $segment_first_half=substr($segment,0,$halfway);
    my $segment_second_half=substr($segment,$halfway); # +1 I think
    # create middle gap
    my $middle_gap="";
    my $length_difference=$longest - length($segment);
    for(my $x=0;$x<$length_difference;$x++){
      $middle_gap.=".";
    }
    # create the new string
    my $structurally_correct_loop=lc($segment_first_half . $middle_gap . $segment_second_half);
    # reassign the HMM columns
       $structurally_correct_loop=~s/\-/./g;
    my @chars=split(/ */,$structurally_correct_loop);
    my $hmmcol_halfway=int ($hmmcol_count/2);
    for(my $p=(scalar(@chars)-1);(($hmmcol_halfway>0)&&($p>=0));$p--){
      $chars[$p]=uc($chars[$p]);
      $chars[$p]=~s/\./-/;
      $hmmcol_halfway--;
      $hmmcol_count--;
    }
    for(my $p=0;(($hmmcol_count>0)&&($p<scalar(@chars)));$p++){
      $chars[$p]=uc($chars[$p]);
      $chars[$p]=~s/\./-/;
      $hmmcol_count--;
    }
    $structurally_correct_loop="";
    for(my $x=0;$x<scalar(@chars);$x++){
      $structurally_correct_loop.=$chars[$x];
    }
    #print $segment_before . "\t" .  $structurally_correct_loop . "\t" . $segment_after . "\n";
    $self->setSequence($n,($segment_before . $structurally_correct_loop . $segment_after)); 
  }
}

sub getHMMCOLRange{
  my($self,$s,$startcol,$stopcol)=@_;
  my $sequence=$self->getSeq($s);
  my $subseq="";
  my @residues=split(/ */,$sequence);
  my $current_hmmcol=0;
  for(my $x=0;( ($x<scalar(@residues)) && ($current_hmmcol<=$stopcol));$x++){
    if( ($current_hmmcol >= $startcol) && ($current_hmmcol <= $stopcol) ){
      $subseq.=$residues[$x];
    }
    if($residues[$x]=~m/[A-Z-]/){
      $current_hmmcol++;
    }
  }
  return $subseq;
}

sub scoreMotif {
  my($self,$seq,$motif)=@_;
  # a sequence "TAVVYC" and a motif array "[ATY]","[CV]","[MTH]"
  $seq=~s/[a-z]//g;
  # next score states
  my @residues=split(/ */,$seq);

  my $score=0;
  for(my $m=0;$m<scalar(@$motif);$m++){
    if(defined($residues[$m])){
      my $position_motif=$$motif[$m];
      if($residues[$m]=~m/$position_motif/){
        $score++;
      }
    }
  }
  return $score;
}

sub qcH3 {
  my($self,$outfile)=@_;
  my $range1_start=85;
  my $range1_stop=94;
  my $range2_start=102;
  my $range2_stop=112;
  my $cdr_start=94;
  my $cdr_stop=102;
  my @range1_motif=("[IGRKTDQ]","[ASTVP]","[EDAS]","D","[TSA]","[AG]","[VTILF]","Y","[YF]","C");
  my @range2_motif=("W","[GS]","[QPR]","G","[TASM]","[LTM]","[VLIA]","[TAINS]","[VAI]","[SAP]","[SAP]");

  $self->qcCDR($outfile,$range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
}

sub qcL3 {
  my($self,$outfile)=@_;
  my $range1_start=207;
  my $range1_stop=215;
  my $range2_start=226;
  my $range2_stop=237;
  my $cdr_start=215;
  my $cdr_stop=226;
  my @range1_motif=("[AEPST]","[DE]","[D]","[AEFILSTV]","[AGV]","[DILMSTV]","[Y]","[YFH]","[C]");
  my @range2_motif=("[F]","[G]","[GQPT]","[G]","[T]","[KQR]","[LV]","[DET]","[IV]","[KL]");

  $self->qcCDR($outfile,$range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
}

sub qcH1 {
  my($self,$outfile)=@_;

  my $range1_start=13;
  my $range1_stop=23;
  my $range2_start=35;
  my $range2_stop=45;
  my $cdr_start=26;
  my $cdr_stop=34;
  my @range1_motif=("[P]","[GST]","[ADEGQRST]","[ST]", "[MLV]",    "[QFKRST]","[ILMV]","[ST]", "[C]",    "[AEKSTV]");
  my @range2_motif=("[W]","[AFIV]","[CHKQR]",   "[HKQ]","[AKPRSTVN]","[HPST]", "[EGS]",  "[KNQ]","[AGKRS]","[FLP]");

  $self->qcCDR($outfile,$range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
}

sub qcH2 {
  my($self,$outfile)=@_;

  my $range1_start=35;
  my $range1_stop=45;
  my $range2_start=65;
  my $range2_stop=75;
  my $cdr_start=48;
  my $cdr_stop=59;
  my @range1_motif=("[W]","[AFIV]","[CHKQR]",   "[HKQ]","[AKPRSTVN]","[HPST]", "[EGS]",  "[KNQ]","[AGKRS]","[FLP]");
  my @range2_motif=("[KQR]","[AFILTV]","[AISTKV]","[FILMV]", "[STFDN]", "[AKLRV]","[DE]","[DKNT]","[APS]","[IKST]");

  $self->qcCDR($outfile,$range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
}

sub qcL1 {
  my($self,$outfile)=@_;

  my $range1_start=141;
  my $range1_stop=151;
  my $range2_start=162;
  my $range2_stop=172;
  my $cdr_start=151;
  my $cdr_stop=161;
  my @range1_motif=("[LAGST]","[AILPSTV]","[GRADEQ]","[DEKPQRST]", "[PKTQSR]", "[AIV]",    "[KRST]","[ILM]","[ST]", "[C]");
  my @range2_motif=("[W]","[FHLYV]","[FLEQ]","[HEQS]","[HKQR]","[ANPQTS]", "[RDHEG]", "[PGHKQSTE]","[NALPSTV]","[FIPV]");

  $self->qcCDR($outfile,$range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
}

sub qcL2 {
  my($self,$outfile)=@_;

  my $range1_start=162;
  my $range1_stop=172;
  my $range2_start=184;
  my $range2_stop=193;
  my $cdr_start=177;
  my $cdr_stop=183;
  my @range1_motif=("[W]","[FHLYV]","[FLEQ]","[HEQS]","[HKQR]","[ANPQTS]", "[RDHEG]", "[EPGHKQST]","[NALPSTV]","[FIPV]");
  my @range2_motif=("[GDE]","[ITV]","[LPS]","[ADESV]","[RPM]","[F]","[ST]","[AGSV]","[S]");

  $self->qcCDR($outfile,$range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,\@range1_motif,\@range2_motif);
}


sub qcCDR {
  my($self,$outfile,$range1_start,$range1_stop,$range2_start,$range2_stop,$cdr_start,$cdr_stop,$range1_motif,$range2_motif)=@_;
  # confirms boundaries of CDRs are passing motif requirements
  open(OUT,">$outfile");
  for(my $n=0;$n<$self->{nseqs};$n++){
    my $range1_subseq=$self->getHMMCOLRange($n,$range1_start,$range1_stop);
       $range1_subseq=~s/[a-z\.]//g;
       $range1_subseq=uc($range1_subseq);
    my $range2_subseq=$self->getHMMCOLRange($n,$range2_start,$range2_stop);
       $range2_subseq=~s/[a-z\.]//g;
       $range2_subseq=uc($range2_subseq);
    my $cdr_subseq=$self->getHMMCOLRange($n,$cdr_start,$cdr_stop);
       $cdr_subseq=~s/[\.\-]//g;
       $cdr_subseq=uc($cdr_subseq);
    my $range1_score = $self->scoreMotif($range1_subseq,\@$range1_motif);
    my $range2_score = $self->scoreMotif($range2_subseq,\@$range2_motif);
    #print "testing $outfile\t" . $range1_subseq . "\t" . $range2_subseq . "\t" . $cdr_subseq . "\t" . $range1_score . "\t" . $range2_score . "\n";
    if($range1_score>6){
      if($range2_score>6){
         print OUT $self->getAccession($n) . "\t" . $cdr_subseq . "\n";
       }
     }
   }
  close(OUT);
}

sub printDNAseqs {
  my($self)=@_;
  for(my $n=0;$n<$self->{nseqs};$n++){
    if($self->isDNA($n)){
      $self->printSeq($n);
    }
  }
}

sub printAAseqs {
  my($self)=@_;
  for(my $n=0;$n<$self->{nseqs};$n++){
    unless($self->isDNA($n)){
      $self->printSeq($n);
    }
  }
}

sub isDNA {
  my($self,$c)=@_;
  my $seq=uc($self->getSequence($c));
     $seq=~s/\.*//g;
     $seq=~s/\-*//g;
     $seq=~s/X//g;
  my $chars=length($seq);
     $seq=~s/[ACGTN]*//g;
  my $non_nuc=length($seq);
  my $non_nuc_fraction=$non_nuc/$chars;
  if($non_nuc_fraction<0.10){
    return 1;
  }else{
    return 0;
  }
}

sub deMultiplexReads {
  my($self,$outfile,$midlength)=@_;

  open(OUT,">$outfile");
  for(my $s=0;$s<$self->{nseqs};$s++){
    print OUT $self->getAccession($s) . "\t" . $self->getMid($s,$midlength) . "\n";
  }
  close(OUT);
}

sub getMid {
  my($self,$seq,$midlength)=@_;
  my($header,$sequence)=$self->getSeq($seq);
  my $midrange=substr($sequence,0,$midlength);

  my %mids = $self->getMIDhash($midlength);
  if(defined(${$self->{mids}}{$midrange})){
    return ${$self->{mids}}{$midrange};
  }else{
    return "$midrange";
  }
}

sub loadMIDs {
  my($self,$midfile)=@_;

  unless(-f $midfile){
    print "Error! No midfile $midfile was detected!\n";
    exit;
  }else{
    open(FILE,$midfile);
    my @lines=<FILE>;
    chomp(@lines);
    close(FILE);
    
    for(my $x=0;$x<scalar(@lines);$x++){
      my @fields=split(/\t/,$lines[$x]);
      unless(scalar(@fields) == 2){
        print "Error! Wrong number of fields on line $x in midfile $midfile: " . $lines[$x] . "\n";
        print "Format should be: midname<tab>midseq\n";
        exit;
      }else{
        my $midname=$fields[0];
        my $midseq =$fields[1];
        ${$self->{mids}}{$midseq}=$midname;
      }
    }
  }
}

sub splitSeqHeader {
  my($self,$seq)=@_;
  my @header_fields=split(/;/,${$self->{headers}}[$seq]);
  return @header_fields;
}

sub printInSilicoSort {
  my($self,$minlength,$minshm,$maxshm)=@_;
  for(my $s=0;$s<$self->{nseqs};$s++){
    my @header_fields=$self->splitSeqHeader($s);
    # >GOIECX001C9I8F...;IGHV4-30-4/31 44 0 1;IGHD3-10 18 1 1;IGHJ4 41 0 1;CARTLSYASGSYDYW;;IGHM 88 0 1 
    if(defined($header_fields[1])){
      my @vseg_data=split(/ /,$header_fields[1]);
      if(defined($vseg_data[2])){
        if($vseg_data[1]>$minlength){
          if($vseg_data[2]>=$minshm){
            if($vseg_data[2]<=$maxshm){
              if(scalar(@vseg_data)<5){
                print ${$self->{headers}}[$s] . "\n";
                print ${$self->{seqs}}[$s] . "\n";
              }
            }
          }
        }
      }
    }
  }
}

sub printUniqueClones {
  my($self)=@_;
  # take the subset of unique clones
  my %clones_already_seen=();
  for(my $s=0;$s<$self->{nseqs};$s++){
    my @header_fields=$self->splitSeqHeader($s);
    if(defined($header_fields[4])){
      my @vseg_data=split(/ /,$header_fields[1]);
      my @jseg_data=split(/ /,$header_fields[3]);
      my $h3=$header_fields[4];
      my $clone=$vseg_data[0] . $h3 . $jseg_data[0];
      unless(defined($clones_already_seen{$clone})){
        print ${$self->{headers}}[$s] . "\n";
        print ${$self->{seqs}}[$s] . "\n";
      }
      $self->addCloneVariants(\%clones_already_seen,$vseg_data[0],$h3,$jseg_data[0]);
    }
  }
}

sub addCloneVariants {
  my($self,$clones_hash,$vseg,$h3,$jseg)=@_;
  # registers all distance 1 variants in a hash of CDR3's that have been previously seen
  my($hash,$germline,$sequence)=@_;
  my @aa=("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");
  my @chars=split(/ */,$h3);
  #print $sequence . "\n";
  for(my $position=1;$position<scalar(@chars);$position++){
    my $left_string=substr($h3,0,$position);
    my $right_string=substr($h3,($position+1),(scalar(@chars)-$position));
    for(my $a=0;$a<scalar(@aa);$a++){
      # altering position $position to residue $a
      $$clones_hash{$vseg . $left_string . $aa[$a]  . $right_string . $jseg}=1;
    }
  }
}

1;
