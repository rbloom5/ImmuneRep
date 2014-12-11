use Cwd;
use File::Basename;

#my $pwd = cwd();
#my $raman = dirname(dirname(cwd())) . "/hello";
#print dirname(dirname(cwd()));
#print "\n";
#print $raman;
#print "\n"

my $raman = dirname(dirname(cwd()));
print $raman . "/ncbi-blast/bin/blastn";
print "\n";